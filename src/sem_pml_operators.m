function operators = sem_pml_operators(mesh)
%SEM_PML_OPERATORS Elementwise GLL gradients and quadrature weights.
% Stress remains discontinuous at shared nodes; this is required to recover
% the same SEM stiffness as sem_assemble when PML damping is zero.
dim = mesh.dim;
nnode = size(mesh.nodes,1);
nlocal = size(mesh.conn,2);
ne = size(mesh.conn,1);
nqp = ne*nlocal;
qweight = ones(nlocal,1)*prod(mesh.h/2);
for a = 1:dim
    qweight = qweight.*mesh.weights(mesh.localSubs{a});
end
operators.weights = repmat(qweight,ne,1);
operators.globalNodes = reshape(mesh.conn',[],1);
operators.B = cell(1,dim);
for a = 1:dim
    local = zeros(nlocal);
    for q = 1:nlocal
        same = true(1,nlocal);
        for b = 1:dim
            if b ~= a
                same = same & (mesh.localSubs{b}(:)' == mesh.localSubs{b}(q));
            end
        end
        local(q,same) = mesh.D(mesh.localSubs{a}(q), ...
            mesh.localSubs{a}(same))/(mesh.h(a)/2);
    end
    rows = zeros(ne*nlocal^2,1);
    cols = rows; vals = rows;
    for e = 1:ne
        ids = mesh.conn(e,:);
        [rr,cc] = ndgrid((e-1)*nlocal+(1:nlocal),ids);
        slot = (e-1)*nlocal^2+(1:nlocal^2);
        rows(slot) = rr(:);
        cols(slot) = cc(:);
        vals(slot) = local(:);
    end
    operators.B{a} = sparse(rows,cols,vals,nqp,nnode);
end
end
