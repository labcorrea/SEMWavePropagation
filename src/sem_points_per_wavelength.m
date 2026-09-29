function points = sem_points_per_wavelength(mesh, material, mode, frequency)
%SEM_POINTS_PER_WAVELENGTH Conservative GLL resolution at the source frequency.
% Uses the largest adjacent-node gap and the slowest wave speed in the mode.
validateattributes(frequency,{'numeric'},{'scalar','positive','finite'});
mode = lower(char(mode));
if strcmp(mode,'longitudinal')
    speed = sqrt(material.E/material.rho);
else
    speed = sqrt(material.E/(2*(1+material.nu)*material.rho));
end
points = inf;
for a = 1:mesh.dim
    axisNodes = unique(mesh.nodes(:,a),'stable');
    points = min(points,(speed/frequency)/max(diff(axisNodes)));
end
end
