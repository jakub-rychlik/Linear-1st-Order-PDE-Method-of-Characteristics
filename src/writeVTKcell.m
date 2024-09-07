% Creates .vtk file which can be used in Paraview for rendering
% the graph of approximated solution
function writeVTKcell(filename,t,p,u)
    [np,dim]=size(p);
    [nt]=size(t,1);

    FID = fopen(strcat(filename,'.vtk'),'w+');
    fprintf(FID,'# vtk DataFile Version 2.0\nUnstructured Grid Example\nASCII\n');
    fprintf(FID,'DATASET UNSTRUCTURED_GRID\n');

    fprintf(FID,'POINTS %d float\n',np);
    s='%f %f %f \n';
    P=[p zeros(np,3-dim)];
    fprintf(FID,s,P');

    fprintf(FID,'\nCELLS %d %d\n',nt,nt*(dim+1));
    s='%d ';
    for k=1:dim
        s=horzcat(s,{' %d'});
    end
    s=cell2mat(horzcat(s,{' \n'}));
    fprintf(FID,s,[(dim)*ones(nt,1) t]');

    fprintf(FID,'\n\nCELL_TYPES %d\n',nt);
    s='%d ';
    fprintf(FID,s,5*ones(1,nt));

    fprintf(FID,'\n\nPOINT_DATA %s\nSCALARS u float \nLOOKUP_TABLE default\n',num2str(np));
    s='%f\n';
    fprintf(FID,s,u);

    fclose(FID);