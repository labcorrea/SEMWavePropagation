close all
clear variables
clc

% Longitudinal wave in a bar with PML at both ends.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src'));

% Extend the physical domain to include the absorbing layers.
material = struct('E',70e9,'nu',0.33,'rho',2700);
physicalBounds = [0 0.1];
thickness = [0.02 0.02];

mesh = sem_mesh([physicalBounds(:,1)-thickness(:,1), ...
                 physicalBounds(:,2)+thickness(:,2)],21,5);

model = sem_assemble(mesh,material,'longitudinal');
pml = sem_pml_profile(model,physicalBounds,thickness);

% Place the axial source inside the physical region.
frequency = 2e5;
waveform = sem_hann_sine(frequency,3,1);
assert(sem_points_per_wavelength(mesh,material,model.mode,frequency) >= 5, ...
    'Mesh must have at least five GLL points per wavelength');
source = struct('position',0.03,'direction',1,'waveform',waveform);

% Use a stable time step and solve with absorbing boundaries.
options.dt = 0.5*model.dtCritical;
result = sem_solve_pml(model,pml,source,3.5e-5,options);

% Export the time series for ParaView.
file = sem_export_paraview(model,result,fullfile(here,'output','wave_1d_longitudinal_pml'));
fprintf('ParaView: %s\n',file);
sem_animate_wave(model,result);
