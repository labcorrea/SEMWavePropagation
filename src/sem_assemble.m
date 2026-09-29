function model = sem_assemble(mesh, material, mode)
%SEM_ASSEMBLE Lumped GLL mass and isotropic elastic stiffness.
% mode: 'longitudinal' or 'transverse' (1D), 'SH' or 'PSV' (2D), 'full' (3D).
dim = mesh.dim;
mode = lower(char(mode));
allowed = {{'longitudinal','transverse'},{'sh','psv'},{'full'}};
assert(any(strcmp(mode,allowed{dim})), 'Mode is incompatible with mesh dimension');
E = material.E; nu = material.nu; rho = material.rho;
validateattributes(E, {'numeric'}, {'scalar','positive','finite'});
validateattributes(rho, {'numeric'}, {'scalar','positive','finite'});
validateattributes(nu, {'numeric'}, {'scalar','finite','>',-1,'<',0.5});
mu = E/(2*(1+nu));
lambda = E*nu/((1+nu)*(1-2*nu));
if strcmp(mode,'psv') || strcmp(mode,'full')
    ndof = dim;
else
    ndof = 1;
end
nnode = size(mesh.nodes,1);
nlocal = size(mesh.conn,2);
ne = size(mesh.conn,1);
nd = nnode*ndof;
qweight = ones(nlocal,1)*prod(mesh.h/2);
for a = 1:dim
    qweight = qweight.*mesh.weights(mesh.localSubs{a});
end
% At each collocated GLL point, all other tensor-product factors are Kronecker deltas.
G = zeros(nlocal,nlocal,dim);
for a = 1:dim
    for q = 1:nlocal
        same = true(1,nlocal);
        for b = 1:dim
            if b ~= a
                same = same & (mesh.localSubs{b}(:)' == mesh.localSubs{b}(q));
            end
        end
        G(q,same,a) = mesh.D(mesh.localSubs{a}(q),mesh.localSubs{a}(same))/ (mesh.h(a)/2);
    end
end
Klocal = zeros(nlocal*ndof);
for q = 1:nlocal
    grad = squeeze(G(q,:,:));
    if dim == 1
        grad = reshape(grad,nlocal,1);
    end
    for a = 1:ndof
        ia = a:ndof:nlocal*ndof;
        for b = 1:ndof
            ib = b:ndof:nlocal*ndof;
            if ndof == 1
                if strcmp(mode,'longitudinal')
                    stiffness = E;
                else
                    stiffness = mu;
                end
                block = stiffness*(grad*grad');
            else
                block = lambda*(grad(:,a)*grad(:,b)') + ...
                    mu*(grad(:,b)*grad(:,a)');
                if a == b
                    block = block + mu*(grad*grad');
                end
            end
            Klocal(ia,ib) = Klocal(ia,ib) + qweight(q)*block;
        end
    end
end
rows = zeros(ne*numel(Klocal),1);
cols = rows; vals = rows;
mass = zeros(nd,1);
cursor = 0;
for e = 1:ne
    ids = reshape((mesh.conn(e,:)-1)*ndof+(1:ndof)',[],1);
    [rr,cc] = ndgrid(ids,ids);
    count = numel(Klocal);
    rows(cursor+(1:count)) = rr(:);
    cols(cursor+(1:count)) = cc(:);
    vals(cursor+(1:count)) = Klocal(:);
    cursor = cursor+count;
    for a = 1:ndof
        mass(ids(a:ndof:end)) = mass(ids(a:ndof:end)) + rho*qweight;
    end
end
K = sparse(rows,cols,vals,nd,nd);
M = spdiags(mass,0,nd,nd);
model = struct('mesh',mesh,'material',material,'mode',mode, ...
    'ndof',ndof,'M',M,'K',(K+K')/2, ...
    'cp',sqrt((lambda+2*mu)/rho),'cs',sqrt(mu/rho));
[model.dtCritical,model.lambdaMax] = sem_critical_dt(model.K,model.M);
if strcmp(mode,'longitudinal')
    model.waveSpeed = sqrt(E/rho);
elseif strcmp(mode,'transverse') || strcmp(mode,'sh')
    model.waveSpeed = model.cs;
else
    model.waveSpeed = model.cp;
end
end
