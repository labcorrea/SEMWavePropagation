close all
clear variables
clc

% Full 3D elasticity with coupled P and S waves.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src'));

% Assemble the elastic cube with free boundaries.
material = struct('E',70e9,'nu',0.33,'rho',2700);
mesh = sem_mesh([0 0.04;0 0.04;0 0.04],[4 4 4],2);
model = sem_assemble(mesh,material,'full');

% Excite the x component at the center of the cube.
frequency = 1e5;
waveform = sem_hann_sine(frequency,2,1);
assert(sem_points_per_wavelength(mesh,material,model.mode,frequency) >= 5, ...
    'Mesh must have at least five GLL points per wavelength');
source = struct('position',[0.02 0.02 0.02],'direction',1,'waveform',waveform);

% Choose a time step below the stability limit and solve the response.
options.dt = 0.8*model.dtCritical;
result = sem_solve(model,source,3e-5,options);

% Export displacement snapshots for ParaView.
file = sem_export_paraview(model,result,fullfile(here,'output','wave_3d_full'));
fprintf('ParaView: %s\n',file);
sem_animate_wave(model,result);
