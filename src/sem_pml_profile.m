function pml = sem_pml_profile(model, physicalBounds, thickness, targetReflection, order)
%SEM_PML_PROFILE Polynomial split-field PML outside a rectangular interior.
% thickness is dim x 2: [lower upper] layer widths for each coordinate.
% A zero width leaves that face traction-free, as needed for a half-space.
if nargin < 4 || isempty(targetReflection), targetReflection = 1e-3; end
if nargin < 5 || isempty(order), order = 2; end
dim = model.mesh.dim;
assert(isequal(size(physicalBounds),[dim 2]), 'Invalid physicalBounds size');
assert(isequal(size(thickness),[dim 2]) && all(thickness(:) >= 0), ...
    'thickness must be a nonnegative dim x 2 array');
validateattributes(targetReflection, {'numeric'}, {'scalar','>',0,'<',1});
validateattributes(order, {'numeric'}, {'scalar','integer','positive'});
expectedBounds = [physicalBounds(:,1)-thickness(:,1), ...
                  physicalBounds(:,2)+thickness(:,2)];
assert(max(abs(expectedBounds(:)-model.mesh.bounds(:))) < ...
    1e-10*max(1,max(abs(expectedBounds(:)))), ...
    'Mesh bounds must include the specified PML layers');
sigma = zeros(size(model.mesh.nodes,1),dim);
for a = 1:dim
    for side = 1:2
        width = thickness(a,side);
        if width == 0, continue; end
        sigmaMax = -(order+1)*model.waveSpeed*log(targetReflection)/(2*width);
        if side == 1
            depth = max(0,(physicalBounds(a,1)-model.mesh.nodes(:,a))/width);
        else
            depth = max(0,(model.mesh.nodes(:,a)-physicalBounds(a,2))/width);
        end
        sigma(:,a) = sigma(:,a)+sigmaMax*min(depth,1).^order;
    end
end
pml = struct('sigma',sigma,'physicalBounds',physicalBounds, ...
    'thickness',thickness,'targetReflection',targetReflection,'order',order);
end
