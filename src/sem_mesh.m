function mesh = sem_mesh(bounds, nelem, p)
%SEM_MESH Cartesian tensor-product SEM grid with shared GLL nodes.
% bounds: dim x 2; nelem: number of elements per axis; p: polynomial order.
dim = size(bounds,1);
assert(dim >= 1 && dim <= 3 && size(bounds,2) == 2, 'bounds must be dim x 2, dim=1..3');
nelem = nelem(:)';
assert(numel(nelem) == dim && all(nelem >= 1) && all(nelem == floor(nelem)), 'Invalid nelem');
assert(all(bounds(:,2) > bounds(:,1)), 'Each bound must have positive length');
[xi,w,D] = sem_gll(p);
naxis = nelem*p+1;
axisNodes = cell(1,dim);
h = (bounds(:,2)-bounds(:,1))./nelem(:);
for a = 1:dim
    axisNodes{a} = zeros(naxis(a),1);
    for e = 1:nelem(a)
        ids = (e-1)*p+(1:p+1);
        axisNodes{a}(ids) = bounds(a,1)+h(a)*(e-1+(xi+1)/2);
    end
end
coord = cell(1,dim);
[coord{:}] = ndgrid(axisNodes{:});
nodes = zeros(prod(naxis),dim);
for a = 1:dim
    nodes(:,a) = coord{a}(:);
end

localSubs = cell(1,dim);
[localSubs{:}] = ndgrid(1:p+1);
for a = 1:dim
    localSubs{a} = localSubs{a}(:);
end
elemSubs = cell(1,dim);
elemAxes = arrayfun(@(n) 1:n, nelem, 'UniformOutput', false);
[elemSubs{:}] = ndgrid(elemAxes{:});
for a = 1:dim
    elemSubs{a} = elemSubs{a}(:);
end
conn = zeros(prod(nelem),(p+1)^dim);
for e = 1:prod(nelem)
    subs = cell(1,dim);
    for a = 1:dim
        subs{a} = (elemSubs{a}(e)-1)*p+localSubs{a};
    end
    if dim == 1
        conn(e,:) = subs{1};
    else
        conn(e,:) = sub2ind(naxis,subs{:});
    end
end

% Linear VTK subcells between GLL nodes preserve high-order geometry and data.
cellAxes = arrayfun(@(n) 1:n-1, naxis, 'UniformOutput', false);
cellSubs = cell(1,dim);
[cellSubs{:}] = ndgrid(cellAxes{:});
for a = 1:dim
    cellSubs{a} = cellSubs{a}(:);
end
switch dim
    case 1
        offsets = [0;1]; vtkType = 3;
    case 2
        offsets = [0 0;1 0;1 1;0 1]; vtkType = 9;
    case 3
        offsets = [0 0 0;1 0 0;1 1 0;0 1 0; ...
                   0 0 1;1 0 1;1 1 1;0 1 1]; vtkType = 12;
end
cells = zeros(numel(cellSubs{1}),size(offsets,1));
for c = 1:size(cells,1)
    for corner = 1:size(cells,2)
        subs = cell(1,dim);
        for a = 1:dim
            subs{a} = cellSubs{a}(c)+offsets(corner,a);
        end
        if dim == 1
            cells(c,corner) = subs{1};
        else
            cells(c,corner) = sub2ind(naxis,subs{:});
        end
    end
end
mesh = struct('dim',dim,'bounds',bounds,'nelem',nelem,'order',p, ...
    'xi',xi,'weights',w,'D',D,'h',h,'nodes',nodes,'conn',conn, ...
    'localSubs',{localSubs},'cells',cells,'vtkType',vtkType);
end
