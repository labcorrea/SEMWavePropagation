function test_sem_pml()
% Check zero damping, half-space topology, and reflected wave reduction.
addpath(fullfile(fileparts(mfilename('fullpath')),'..','src'));
material = struct('E',1,'nu',0.25,'rho',1);
physicalBounds = [0 1];
thickness = [0.2 0.2];
mesh = sem_mesh([-0.2 1.2],70,2);
model = sem_assemble(mesh,material,'longitudinal');
pml = sem_pml_profile(model,physicalBounds,thickness);
waveform = sem_hann_sine(8,3,1);
source = struct('position',0.5,'direction',1,'waveform',waveform);
didError = false;
try
    sem_solve_pml(model,pml,source,0.1,struct('dt',1.01*model.dtCritical));
catch ME
    didError = strcmp(ME.identifier,'sem_solve_pml:UnstableTimeStep');
end
assert(didError, 'PML solver must reject steps above the eigenvalue limit');
options = struct('dt',0.4*model.dtCritical,'snapshotEvery',1);
absorbing = sem_solve_pml(model,pml,source,2.1,options);
zeroPml = pml;
zeroPml.sigma(:) = 0;
reflecting = sem_solve_pml(model,zeroPml,source,2.1,options);
[~,receiver] = min(abs(mesh.nodes(:,1)-0.65));
time = absorbing.time;
window = time > 1.25 & time < 1.85;
va = diff(absorbing.U(receiver,:))/absorbing.dt;
vr = diff(reflecting.U(receiver,:))/reflecting.dt;
ratio = max(abs(va(window(2:end))))/max(abs(vr(window(2:end))));
fprintf('1D PML reflected velocity ratio: %.3f\n',ratio);
assert(ratio < 0.5, 'PML did not sufficiently reduce reflected waves');

% The half-space top has no damping, while the side and bottom do.
halfBounds = [0 1;0 1];
halfThickness = [0.2 0.2;0 0.2];
halfMesh = sem_mesh([-0.2 1.2;0 1.2],[3 3],2);
halfModel = sem_assemble(halfMesh,material,'PSV');
halfPml = sem_pml_profile(halfModel,halfBounds,halfThickness);
top = abs(halfMesh.nodes(:,2)) < 1e-12;
assert(all(halfPml.sigma(top,2) == 0));
bottom = halfMesh.nodes(:,2) > 1;
assert(any(halfPml.sigma(bottom,2) > 0));

% With zero damping, the velocity-stress PML equations reproduce the
% undamped P-SV wavefield up to the time integrators' truncation error.
smallMesh = sem_mesh([0 1;0 0.8],[4 4],2);
smallModel = sem_assemble(smallMesh,material,'PSV');
noLayer = sem_pml_profile(smallModel,smallMesh.bounds,zeros(2,2));
vectorWaveform = sem_hann_sine(5,3,1);
vectorSource = struct('position',[0.5 0.4],'direction',2, ...
    'waveform',vectorWaveform);
smallOptions = struct('dt',0.08*smallModel.dtCritical,'snapshotEvery',1);
direct = sem_solve(smallModel,vectorSource,0.4,smallOptions);
split = sem_solve_pml(smallModel,noLayer,vectorSource,0.4,smallOptions);
relativeDifference = max(abs(direct.U(:)-split.U(:)))/max(abs(direct.U(:)));
fprintf('Zero-damping P-SV relative difference: %.3g\n',relativeDifference);
assert(relativeDifference < 0.01, 'Zero-damping PML disagrees with the SEM model');
fprintf('SEM PML checks passed.\n');
end
