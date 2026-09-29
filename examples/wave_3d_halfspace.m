close all
clear variables
clc

% 3D elastic waves in a half-space with a free surface at z=0.
% Positive z points into the solid.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src'));

% Zero PML thickness on the top face preserves the free surface.
material = struct('E',70e9,'nu',0.33,'rho',2700);
physicalBounds = [0 0.04;0 0.04;0 0.04];
thickness = [0.02 0.02;0.02 0.02;0 0.02];
mesh = sem_mesh([physicalBounds(:,1)-thickness(:,1), ...
                 physicalBounds(:,2)+thickness(:,2)],[8 8 6],2);
model = sem_assemble(mesh,material,'full');
pml = sem_pml_profile(model,physicalBounds,thickness);

% Apply a vertical force near the free surface.
frequency = 1e5;
waveform = sem_hann_sine(frequency,2,1);
assert(sem_points_per_wavelength(mesh,material,model.mode,frequency) >= 5, ...
    'Mesh must have at least five GLL points per wavelength');
source = struct('position',[0.02 0.02 0.01],'direction',3,'waveform',waveform);

% Solve with absorption on the sides and bottom.
options.dt = 0.5*model.dtCritical;
result = sem_solve_pml(model,pml,source,4e-5,options);

% Export the time series for ParaView.
file = sem_export_paraview(model,result,fullfile(here,'output','wave_3d_halfspace'));
fprintf('ParaView: %s\n',file);
sem_animate_wave(model,result);
