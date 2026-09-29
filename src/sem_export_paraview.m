function pvdFile = sem_export_paraview(model, result, filenameBase)
%SEM_EXPORT_PARAVIEW Export snapshots to linear-subcell VTU files and a PVD index.
% filenameBase may contain a directory. Open the returned .pvd in ParaView.
mesh = model.mesh;
[outputDir,stem] = fileparts(filenameBase);
if isempty(outputDir), outputDir = pwd; end
if ~exist(outputDir,'dir'), mkdir(outputDir); end
nnode = size(mesh.nodes,1);
ncell = size(mesh.cells,1);
coordinates = [mesh.nodes,zeros(nnode,3-mesh.dim)];
pvdFile = fullfile(outputDir,[stem '.pvd']);
pvd = fopen(pvdFile,'w');
assert(pvd >= 0, 'Cannot create PVD file');
cleanup = onCleanup(@() fclose(pvd));
fprintf(pvd,'<?xml version="1.0"?>\n<VTKFile type="Collection" version="0.1" byte_order="LittleEndian">\n<Collection>\n');
for frame = 1:numel(result.time)
    name = sprintf('%s_%04d.vtu',stem,frame);
    file = fullfile(outputDir,name);
    fid = fopen(file,'w');
    assert(fid >= 0, 'Cannot create VTU file: %s',file);
    fprintf(fid,'<?xml version="1.0"?>\n<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n');
    fprintf(fid,'<UnstructuredGrid><Piece NumberOfPoints="%d" NumberOfCells="%d">\n',nnode,ncell);
    fprintf(fid,'<Points><DataArray type="Float64" NumberOfComponents="3" format="ascii">\n');
    fprintf(fid,'%.16g %.16g %.16g\n',coordinates');
    fprintf(fid,'</DataArray></Points>\n<Cells>\n');
    fprintf(fid,'<DataArray type="Int32" Name="connectivity" format="ascii">\n');
    fprintf(fid,'%d ',reshape((mesh.cells-1)',[],1));
    fprintf(fid,'\n</DataArray>\n<DataArray type="Int32" Name="offsets" format="ascii">\n');
    fprintf(fid,'%d ',(1:ncell)*size(mesh.cells,2));
    fprintf(fid,'\n</DataArray>\n<DataArray type="UInt8" Name="types" format="ascii">\n');
    fprintf(fid,'%d ',repmat(mesh.vtkType,1,ncell));
    fprintf(fid,'\n</DataArray></Cells>\n<PointData>\n');
    field = reshape(result.U(:,frame),model.ndof,[])';
    if model.ndof == 1
        fprintf(fid,'<DataArray type="Float64" Name="displacement" NumberOfComponents="1" format="ascii">\n');
        fprintf(fid,'%.16g\n',field');
        fprintf(fid,'</DataArray>\n');
    else
        vector = zeros(nnode,3);
        vector(:,1:model.ndof) = field;
        fprintf(fid,'<DataArray type="Float64" Name="displacement" NumberOfComponents="3" format="ascii">\n');
        fprintf(fid,'%.16g %.16g %.16g\n',vector');
        fprintf(fid,'</DataArray>\n');
    end
    fprintf(fid,'</PointData></Piece></UnstructuredGrid></VTKFile>\n');
    fclose(fid);
    fprintf(pvd,'<DataSet timestep="%.16g" group="" part="0" file="%s"/>\n',result.time(frame),name);
end
fprintf(pvd,'</Collection></VTKFile>\n');
end
