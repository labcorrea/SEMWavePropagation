function result = sem_solve_pml(model, pml, source, tEnd, options)
%SEM_SOLVE_PML Split-field velocity-stress PML with explicit RK4 time stepping.
% The PML acts on first spatial derivatives; sigma=0 recovers elasticity.
% A zero PML width leaves the corresponding face traction-free.
if nargin < 5, options = struct(); end
if ~isfield(options,'dt') || isempty(options.dt)
    error('sem_solve_pml:MissingTimeStep', ...
        'Choose options.dt explicitly. Eigenvalue limit: %.6g s.',model.dtCritical);
end
validateattributes(options.dt, {'numeric'}, {'scalar','positive','finite'});
validateattributes(tEnd, {'numeric'}, {'scalar','positive','finite'});
if options.dt > model.dtCritical
    error('sem_solve_pml:UnstableTimeStep', ...
        'Chosen dt = %.6g s exceeds eigenvalue limit %.6g s.', ...
        options.dt,model.dtCritical);
end
sigma = pml.sigma;
mesh = model.mesh;
nnode = size(mesh.nodes,1);
dim = mesh.dim;
ndof = model.ndof;
assert(isequal(size(sigma),[nnode dim]) && all(isfinite(sigma(:))) && ...
    all(sigma(:) >= 0), 'Invalid PML damping profile');
if ~isfield(options,'snapshotEvery') || isempty(options.snapshotEvery)
    options.snapshotEvery = max(1,ceil(tEnd/options.dt/40));
end
validateattributes(options.snapshotEvery, {'numeric'}, {'scalar','integer','positive'});
nSteps = ceil(tEnd/options.dt);
dt = tEnd/nSteps;
% RK4's negative-real-axis limit is about 2.785; stay below it for damping.
if max(sigma(:))*dt > 2.5
    error('sem_solve_pml:UnstableDamping', ...
        'Chosen dt is too large for PML damping; require max(sigma)*dt <= 2.5.');
end
position = source.position(:)';
assert(numel(position) == dim, 'Source position dimension mismatch');
assert(source.direction >= 1 && source.direction <= ndof && ...
    source.direction == floor(source.direction), 'Invalid source.direction');
[~,sourceNode] = min(sum((mesh.nodes-position).^2,2));
assert(all(mesh.nodes(sourceNode,:) >= pml.physicalBounds(:,1)'-1e-12) && ...
    all(mesh.nodes(sourceNode,:) <= pml.physicalBounds(:,2)'+1e-12), ...
    'Source must lie in the physical domain');

operators = sem_pml_operators(mesh);
B = operators.B;
qweight = operators.weights;
sigmaQ = sigma(operators.globalNodes,:);
nqp = numel(qweight);
mass = full(diag(model.M));
nodalMass = mass(1:ndof:end);
mu = model.material.E/(2*(1+model.material.nu));
lambda = model.material.E*model.material.nu/ ...
    ((1+model.material.nu)*(1-2*model.material.nu));
if ndof == 1
    nstress = dim;
    if strcmp(model.mode,'longitudinal')
        scalarStiffness = model.material.E;
    else
        scalarStiffness = mu;
    end
else
    nstress = ndof^2;
    scalarStiffness = 0;
end
u = zeros(nnode,ndof);
vSplit = zeros(nnode,ndof,dim);
stressSplit = zeros(nqp,nstress,dim);
saveSteps = unique([0:options.snapshotEvery:nSteps,nSteps]);
U = zeros(nnode*ndof,numel(saveSteps));
forceSamples = sem_source_samples(source,(0:2*nSteps)*dt/2);
saveIndex = 1;
for step = 0:nSteps
    if step == saveSteps(saveIndex)
        U(:,saveIndex) = reshape(u',[],1);
        saveIndex = saveIndex+1;
    end
    if step == nSteps, break; end
    [k1u,k1v,k1s] = rates(vSplit,stressSplit,forceSamples(2*step+1));
    [k2u,k2v,k2s] = rates(vSplit+dt*k1v/2,stressSplit+dt*k1s/2,forceSamples(2*step+2));
    [k3u,k3v,k3s] = rates(vSplit+dt*k2v/2,stressSplit+dt*k2s/2,forceSamples(2*step+2));
    [k4u,k4v,k4s] = rates(vSplit+dt*k3v,stressSplit+dt*k3s,forceSamples(2*step+3));
    u = u+dt*(k1u+2*k2u+2*k3u+k4u)/6;
    vSplit = vSplit+dt*(k1v+2*k2v+2*k3v+k4v)/6;
    stressSplit = stressSplit+dt*(k1s+2*k2s+2*k3s+k4s)/6;
    if any(~isfinite(u(:)))
        error('Non-finite PML wavefield: reduce options.dt');
    end
end
result = struct('time',saveSteps*dt,'U',U,'dt',dt, ...
    'dtCritical',model.dtCritical,'sourceNode',sourceNode, ...
    'sourceDof',(sourceNode-1)*ndof+source.direction);

    function [du,dv,ds] = rates(vParts,sParts,force)
        v = sum(vParts,3);
        stress = sum(sParts,3);
        du = v;
        dv = zeros(size(vParts));
        ds = zeros(size(sParts));
        gradient = zeros(nqp,ndof,dim);
        for j = 1:dim
            gradient(:,:,j) = B{j}*v;
        end
        for j = 1:dim
            if ndof == 1
                dv(:,1,j) = -(B{j}'*(qweight.*stress(:,j)))./nodalMass ...
                    -sigma(:,j).*vParts(:,1,j);
                ds(:,j,j) = scalarStiffness*gradient(:,1,j) ...
                    -sigmaQ(:,j).*sParts(:,j,j);
            else
                for a = 1:ndof
                    dv(:,a,j) = -(B{j}'*(qweight.*stress(:,a+(j-1)*ndof)))./nodalMass ...
                        -sigma(:,j).*vParts(:,a,j);
                end
                for a = 1:ndof
                    for b = 1:ndof
                        idx = a+(b-1)*ndof;
                        rate = zeros(nqp,1);
                        if a == b, rate = rate+lambda*gradient(:,j,j); end
                        if a == j, rate = rate+mu*gradient(:,b,j); end
                        if b == j, rate = rate+mu*gradient(:,a,j); end
                        ds(:,idx,j) = rate-sigmaQ(:,j).*sParts(:,idx,j);
                    end
                end
            end
        end
        dv(sourceNode,source.direction,1) = ...
            dv(sourceNode,source.direction,1)+force/mass((sourceNode-1)*ndof+source.direction);
    end
end
