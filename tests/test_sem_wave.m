function test_sem_wave()
% Small physical invariants and a free-free bar eigenfrequency check.
addpath(fullfile(fileparts(mfilename('fullpath')),'..','src'));
material = struct('E',70e9,'nu',0.33,'rho',2700);
modes = {'longitudinal','transverse','SH','PSV','full'};
bounds = {[0 1],[0 1],[0 1;0 0.8],[0 1;0 0.8], ...
    [0 1;0 0.8;0 0.6]};
elements = {5,5,[2 2],[2 2],[1 1 1]};
for k = 1:numel(modes)
    mesh = sem_mesh(bounds{k},elements{k},3);
    model = sem_assemble(mesh,material,modes{k});
    expectedMass = material.rho*prod(bounds{k}(:,2)-bounds{k}(:,1))*model.ndof;
    assert(issparse(model.K) && issparse(model.M));
    assert(nnz(model.M) == size(model.M,1));
    massDiagonal = full(diag(model.M));
    assert(abs(sum(massDiagonal)/expectedMass-1) < 1e-12);
    assert(all(massDiagonal > 0));
    assert(abs(model.dtCritical-2/sqrt(model.lambdaMax)) < 1e-12*model.dtCritical);
    assert(norm(model.K-model.K','fro')/norm(model.K,'fro') < 1e-12);
    for component = 1:model.ndof
        rigid = zeros(size(model.M,1),1);
        rigid(component:model.ndof:end) = 1;
        assert(norm(model.K*rigid)/norm(model.K,1) < 1e-11);
    end
    assert(size(mesh.cells,2) == 2^mesh.dim);
end
mesh = sem_mesh([0 1],5,3);
model = sem_assemble(mesh,material,'longitudinal');
frequencies2 = sort(real(eig(full(model.K),full(model.M))));
omega = sqrt(frequencies2(2));
assert(abs(omega/(pi*model.waveSpeed)-1) < 1e-3);
assert(abs(model.lambdaMax-max(frequencies2))/model.lambdaMax < 1e-10);
waveform = sem_hann_sine(100,2,1);
source = struct('position',0.5,'direction',1,'waveform',waveform);
assert(waveform.values(1) == 0 && waveform.values(end) == 0);
assert(all(sem_source_samples(source,[-1,waveform.time(end)+1]) == 0));
custom = struct('waveform',struct('time',[0 1 2],'values',[0 2 0]));
assert(isequal(sem_source_samples(custom,[-1 0.5 2.5]),[0 1 0]));
didError = false;
try
    sem_solve(model,source,0.01);
catch ME
    didError = strcmp(ME.identifier,'sem_solve:MissingTimeStep');
end
assert(didError, 'The user must choose a time step');
didError = false;
try
    sem_solve(model,source,0.01,struct('dt',1.01*model.dtCritical));
catch ME
    didError = strcmp(ME.identifier,'sem_solve:UnstableTimeStep');
end
assert(didError, 'A time step above the eigenvalue limit must fail');
result = sem_solve(model,source,0.01,struct('dt',0.8*model.dtCritical));
assert(all(isfinite(result.U(:))));
fprintf('SEM wave checks passed.\n');
end
