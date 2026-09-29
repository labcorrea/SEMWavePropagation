close all
clear variables
clc

% Plane-strain P-SV waves with coupled x and y displacements.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src'));

% Assemble a 2D mesh with free boundaries.
material = struct('E',70e9,'nu',0.33,'rho',2700);
mesh = sem_mesh([0 0.06;0 0.06],[9 9],2);
model = sem_assemble(mesh,material,'PSV');

% Excite the x component inside the domain.
frequency = 1e5;
waveform = sem_hann_sine(frequency,2,1);
assert(sem_points_per_wavelength(mesh,material,model.mode,frequency) >= 5, ...
    'Mesh must have at least five GLL points per wavelength');
source = struct('position',[0.03 0.03],'direction',1,'waveform',waveform);

% Choose a time step below the stability limit and solve the response.
options.dt = 0.8*model.dtCritical;
result = sem_solve(model,source,3e-5,options);

% Export displacement snapshots for ParaView.
file = sem_export_paraview(model,result,fullfile(here,'output','wave_2d_psv'));
fprintf('ParaView: %s\n',file);
sem_animate_wave(model,result);
