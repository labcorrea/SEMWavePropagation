function result = sem_solve(model, source, tEnd, options)
%SEM_SOLVE Central-difference integration of M*u_tt + K*u = f(t).
% source.position is a point; source.direction is a displacement component.
% source.waveform contains sampled time and values arrays.
% options.dt is required and must not exceed model.dtCritical.
% options.snapshotEvery is optional. Boundaries are traction-free.
if nargin < 4, options = struct(); end
mesh = model.mesh;
validateattributes(tEnd, {'numeric'}, {'scalar','positive','finite'});
if ~isfield(options,'dt') || isempty(options.dt)
    error('sem_solve:MissingTimeStep', ...
        'Choose options.dt explicitly. Eigenvalue stability limit: %.6g s.',model.dtCritical);
end
validateattributes(options.dt, {'numeric'}, {'scalar','positive','finite'});
if options.dt > model.dtCritical
    error('sem_solve:UnstableTimeStep', ...
        'Chosen dt = %.6g s exceeds eigenvalue stability limit %.6g s.', ...
        options.dt,model.dtCritical);
end
if ~isfield(options,'snapshotEvery') || isempty(options.snapshotEvery)
    options.snapshotEvery = max(1,ceil(tEnd/options.dt/40));
end
validateattributes(options.snapshotEvery, {'numeric'}, {'scalar','integer','positive'});
nSteps = ceil(tEnd/options.dt);
dt = tEnd/nSteps;
position = source.position(:)';
assert(numel(position) == mesh.dim, 'Source position dimension mismatch');
assert(source.direction >= 1 && source.direction <= model.ndof && ...
    source.direction == floor(source.direction), 'Invalid source.direction');
[~,sourceNode] = min(sum((mesh.nodes-position).^2,2));
sourceDof = (sourceNode-1)*model.ndof+source.direction;
times = (0:nSteps)*dt;
forceSamples = sem_source_samples(source,times);
saveSteps = unique([0:options.snapshotEvery:nSteps,nSteps]);
U = zeros(size(model.M,1),numel(saveSteps));
u = zeros(size(model.M,1),1);
uPrev = u;
massDiagonal = full(diag(model.M));
saveIndex = 1;
for step = 0:nSteps
    if step == saveSteps(saveIndex)
        U(:,saveIndex) = u;
        saveIndex = saveIndex+1;
    end
    if step == nSteps, break; end
    acceleration = -model.K*u;
    acceleration(sourceDof) = acceleration(sourceDof)+forceSamples(step+1);
    acceleration = acceleration./massDiagonal;
    if step == 0
        uNext = u+0.5*dt^2*acceleration;
    else
        uNext = 2*u-uPrev+dt^2*acceleration;
    end
    uPrev = u;
    u = uNext;
    if any(~isfinite(u))
        error('Non-finite wavefield: reduce options.dt');
    end
end
result = struct('time',times(saveSteps+1),'U',U,'dt',dt, ...
    'dtCritical',model.dtCritical,'sourceNode',sourceNode,'sourceDof',sourceDof);
end
