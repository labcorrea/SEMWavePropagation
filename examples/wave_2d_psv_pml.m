close all
clear variables
clc

% Plane-strain P-SV waves with PML on all sides.
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here,'..','src'));

% Include the physical region and four absorbing layers in the mesh.
material = struct('E',70e9,'nu',0.33,'rho',2700);
physicalBounds = [0 0.06;0 0.06];
thickness = [0.02 0.02;0.02 0.02];
mesh = sem_mesh([physicalBounds(:,1)-thickness(:,1), ...
                 physicalBounds(:,2)+thickness(:,2)],[13 13],2);
model = sem_assemble(mesh,material,'PSV');
pml = sem_pml_profile(model,physicalBounds,thickness);

% Excite the x component at the center of the physical region.
frequency = 1.5e5;
waveform = sem_hann_sine(frequency,3,1);
assert(sem_points_per_wavelength(mesh,material,model.mode,frequency) >= 5, ...
    'Mesh must have at least five GLL points per wavelength');
source = struct('position',[0.03 0.03],'direction',1,'waveform',waveform);

% Use a stable time step and solve with PML.
options.dt = 0.5*model.dtCritical;
result = sem_solve_pml(model,pml,source,3.5e-5,options);

% Export the time series for ParaView.
file = sem_export_paraview(model,result,fullfile(here,'output','wave_2d_psv_pml'));
fprintf('ParaView: %s\n',file);
sem_animate_wave(model,result);
