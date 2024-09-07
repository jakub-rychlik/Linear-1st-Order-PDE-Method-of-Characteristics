% Adjusts vtk arrays and call function, which will create .vtk file for Paraview
function createVTK(chars,vtkPoints,vtkCells)
    vtkScalars = [];
    for k=1:length(chars)  % creates array of approximated solution values in vertices of triangles
        line = chars{k};
        for i=1:length(line)
            vtkScalars = [vtkScalars,line(i).sol];
        end
    end

    L=length(vtkScalars);  % adjusting dimensions of arrays
    vtkPoints = [vtkPoints(1:end-1,:)];
    i=1;
    while true
        if vtkCells(i,1) == L || vtkCells(i,2) == L || vtkCells(i,3) == L
            vtkCells = [vtkCells(1:i-1,:);vtkCells(i+1:end,:)];
        else
            i = i+1;
        end
        if i == length(vtkCells(:,1))
            if vtkCells(i,1) == L || vtkCells(i,2) == L || vtkCells(i,3) == L
                vtkCells = [vtkCells(1:end-1,:)];
            end
            break
        end
    end
    writeVTKcell('Paraview_data',vtkCells,vtkPoints,vtkScalars);
end