function [dtCritical, lambdaMax] = sem_critical_dt(K, M)
%SEM_CRITICAL_DT Central-difference stability limit from K*v=lambda*M*v.
% K and M must be sparse; M is the positive diagonal GLL mass matrix.
assert(issparse(K) && issparse(M) && isequal(size(K),size(M)), ...
    'K and M must be sparse matrices of the same size');
m = full(diag(M));
assert(all(isfinite(m) & m > 0), 'M must have a positive diagonal');
assert(nnz(M) == numel(m), 'M must be diagonal');
invSqrtM = spdiags(1./sqrt(m),0,numel(m),numel(m));
A = invSqrtM*K*invSqrtM;
A = (A+A')/2;
if size(A,1) <= 200
    lambdaMax = max(real(eig(full(A))));
else
    opts = struct('tol',1e-10,'maxit',1000);
    lambdaMax = real(eigs(A,1,'lm',opts));
end
assert(isfinite(lambdaMax) && lambdaMax > 0, 'Invalid largest eigenvalue');
dtCritical = 2/sqrt(lambdaMax);
end
