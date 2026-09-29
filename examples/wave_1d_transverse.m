close all
clear variables
clc

% Scalar transverse wave in a bar: displacement perpendicular to x.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src'));

% Assemble the bar in transverse shear mode.
material = struct('E',70e9,'nu',0.33,'rho',2700);
mesh = sem_mesh([0 0.10],15,3);
model = sem_assemble(mesh,material,'transverse');

% Apply a localized transverse force with a Hann-windowed pulse.
frequency = 1e5;
waveform = sem_hann_sine(frequency,2,1);
assert(sem_points_per_wavelength(mesh,material,model.mode,frequency) >= 5, ...
    'Mesh must have at least five GLL points per wavelength');
source = struct('position',0.03,'direction',1,'waveform',waveform);

% Choose a time step below the stability limit and solve the response.
options.dt = 0.8*model.dtCritical;
result = sem_solve(model,source,4e-5,options);

% Export displacement snapshots for ParaView.
file = sem_export_paraview(model,result,fullfile(here,'output','wave_1d_transverse'));
fprintf('ParaView: %s\n',file);
sem_animate_wave(model,result);
