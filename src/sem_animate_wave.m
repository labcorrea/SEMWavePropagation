function fig = sem_animate_wave(model, result, options)
%SEM_ANIMATE_WAVE Animate saved displacement snapshots in MATLAB.
% 1D: line plot. 2D: surface. 3D: three orthogonal slices.
% Scalar fields show signed displacement; vector fields show magnitude.
% options.framePause (seconds, default 0.05) controls playback speed.
% options.visible ('on' or 'off', default 'on') controls the figure.
if nargin < 3, options = struct(); end
if ~isfield(options,'framePause'), options.framePause = 0.05; end
if ~isfield(options,'visible'), options.visible = 'on'; end
validateattributes(options.framePause,{'numeric'},{'scalar','nonnegative','finite'});
assert(any(strcmp(options.visible,{'on','off'})), 'options.visible must be on or off');

mesh = model.mesh;
dim = mesh.dim;
nnode = size(mesh.nodes,1);
assert(size(result.U,1) == nnode*model.ndof && ...
    size(result.U,2) == numel(result.time), 'Result size does not match model');
if model.ndof == 1
    values = result.U;
    fieldName = 'Displacement';
else
    values = zeros(nnode,numel(result.time));
    for frame = 1:numel(result.time)
        displacement = reshape(result.U(:,frame),model.ndof,[])';
        values(:,frame) = sqrt(sum(displacement.^2,2));
    end
    fieldName = 'Displacement magnitude';
end
peak = max(abs(values(:)));
if peak == 0, peak = 1; end
limits = [0 peak];
if model.ndof == 1, limits = [-peak peak]; end

axesNodes = cell(1,dim);
shape = zeros(1,dim);
for a = 1:dim
    axesNodes{a} = unique(mesh.nodes(:,a),'stable');
    shape(a) = numel(axesNodes{a});
end
assert(prod(shape) == nnode, 'Expected a Cartesian mesh');
fig = figure('Name','SEM wave propagation','Visible',options.visible);
ax = axes('Parent',fig);
if dim == 1
    lineHandle = plot(ax,axesNodes{1},values(:,1),'LineWidth',1.5);
    xlabel(ax,'x (m)'); ylabel(ax,fieldName);
    xlim(ax,mesh.bounds(1,:)); ylim(ax,limits);
    grid(ax,'on');
elseif dim == 2
    x = axesNodes{1}; y = axesNodes{2};
    [X,Y] = meshgrid(x,y);
    height = 0.2*min(diff(mesh.bounds,1,2));
    frameValues = reshape(values(:,1),shape)';
    surfaceHandle = surf(ax,X,Y,height*frameValues/peak,frameValues, ...
        'EdgeColor','none');
    xlabel(ax,'x (m)'); ylabel(ax,'y (m)'); zlabel(ax,'Scaled displacement');
    xlim(ax,mesh.bounds(1,:)); ylim(ax,mesh.bounds(2,:));
    if model.ndof == 1, zlim(ax,[-height height]);
    else, zlim(ax,[0 height]); end
    view(ax,3); colorbar(ax); caxis(ax,limits);
else
    x = axesNodes{1}; y = axesNodes{2}; z = axesNodes{3};
    [X,Y,Z] = meshgrid(x,y,z);
    cuts = [x(round((numel(x)+1)/2)), ...
            y(round((numel(y)+1)/2)), ...
            z(round((numel(z)+1)/2))];
    xlabel(ax,'x (m)'); ylabel(ax,'y (m)'); zlabel(ax,'z (m)');
    xlim(ax,mesh.bounds(1,:)); ylim(ax,mesh.bounds(2,:)); zlim(ax,mesh.bounds(3,:));
    view(ax,3); colorbar(ax); caxis(ax,limits);
end
for frame = 1:numel(result.time)
    if ~isgraphics(fig), break; end
    if dim == 1
        set(lineHandle,'YData',values(:,frame));
    elseif dim == 2
        frameValues = reshape(values(:,frame),shape)';
        set(surfaceHandle,'ZData',height*frameValues/peak,'CData',frameValues);
    else
        frameValues = permute(reshape(values(:,frame),shape),[2 1 3]);
        delete(findobj(ax,'Type','Surface'));
        sliceHandles = slice(ax,X,Y,Z,frameValues,cuts(1),cuts(2),cuts(3));
        set(sliceHandles,'EdgeColor','none');
    end
    title(ax,sprintf('%s at t = %.3g s',fieldName,result.time(frame)));
    drawnow;
    if options.framePause > 0, pause(options.framePause); end
end
end
