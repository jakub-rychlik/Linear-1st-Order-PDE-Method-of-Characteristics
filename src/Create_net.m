% Creates triangle net with vertices lying on the characteristics
% It also creates two array 'vtkPoints' and 'vtkCells' for generating .vtk file for Paraview 
function [new_chars_w_vert,vtkPoints,vtkCells] = Create_net(char_approx,Total_points)
    format long
    char_w_vertices = cell(1,length(char_approx));
    vtkPoints = [zeros(Total_points-1,3)];
    vtkCells = [];
    begin = 2;

    % first line
    L2 = length(char_approx{1}(1,:))-1;
    L3 = length(char_approx{2}(1,:));
    vtkPoints(1:length(char_approx{1}(1,:)),1:2) = (char_approx{1})';
    if length(char_approx{1}(1,:)) == 1  % there is only one point on first line
        char_w_vertices{1} = [struct('coor', char_approx{1}(:,1), 'prev_vert',[],'act_vert',[],'next_vert',[char_approx{2}(:,1),char_approx{2}(:,2)],'sol',0)];
        vtkCells = [1,2,3];
        
        L2 = length(char_approx{2}(1,:))-1;
        L3 = length(char_approx{3}(1,:));
        char_w_vertices{2} = [struct('coor', char_approx{2}(:,1), 'prev_vert',[char_approx{1}(:,1)],'act_vert',[char_approx{2}(:,2)],'next_vert',[char_approx{3}(:,1)],'sol',0)];
        vtkCells = [vtkCells;2,3,1 + L2+2];
        for p=2:min(L2,L3)
            char_w_vertices{2} = [char_w_vertices{2}, struct('coor', char_approx{2}(:,p),'prev_vert',[char_approx{1}(:,1),char_approx{1}(:,1)],'act_vert',[char_approx{2}(:,p-1),char_approx{2}(:,p+1)],'next_vert',[char_approx{3}(:,p-1),char_approx{3}(:,p)],'sol',0)];
            vtkCells = [vtkCells;1 + p,1 + L2+1+p-1,1 + L2+1+p];
            vtkCells = [vtkCells;1 + p,1 + p+1,1 + L2+1+p];
            vtkCells = [vtkCells;1,1 + p,1 + p+1];
        end

        if L2 == min(L2,L3) && L3 > min(L2,L3)
            char_w_vertices{2} = [char_w_vertices{2}, struct('coor', char_approx{2}(:,L2+1),'prev_vert',[char_approx{1}(:,1)],'act_vert',[char_approx{2}(:,L2)],'next_vert',[char_approx{3}(:,L2),char_approx{3}(:,L2+1)],'sol',0)];
            vtkCells = [vtkCells;1 + L2+1,1 + L2+1+L2,1 + L2+1+L2+1];

        elseif L2 > min(L2,L3) && L3 == min(L2,L3)
            for p=min(L2,L3)+1:L2
                char_w_vertices{2} = [char_w_vertices{2}, struct('coor', char_approx{2}(:,p),'prev_vert',[char_approx{1}(:,1),char_approx{1}(:,1)],'act_vert',[char_approx{2}(:,p-1),char_approx{2}(:,p+1)],'next_vert',[char_approx{3}(:,L3),char_approx{3}(:,L3)],'sol',0)];
                vtkCells = [vtkCells;1 + p,1 + p+1,1 + L2+1+L3];
            end
            char_w_vertices{2} = [char_w_vertices{2}, struct('coor', char_approx{2}(:,L2+1),'prev_vert',[char_approx{1}(:,1)],'act_vert',[char_approx{2}(:,L2)],'next_vert',[char_approx{3}(:,L3)],'sol',0)];
        end
        begin = 3;

    else  %more points on first line
        char_w_vertices{1} = [struct('coor', char_approx{1}(:,1), 'prev_vert',[],'act_vert',[char_approx{1}(:,2)],'next_vert',[char_approx{2}(:,1)],'sol',0)];
        vtkCells = [1,2,L2+2];

        if min(L2,L3) > 1
            for p=2:min(L2,L3)
                char_w_vertices{1} = [char_w_vertices{1}, struct('coor', char_approx{1}(:,p), 'prev_vert',[],'act_vert',[char_approx{1}(:,p-1),char_approx{1}(:,p+1)],'next_vert',[char_approx{2}(:,p-1),char_approx{2}(:,p)],'sol',0)];
                vtkCells = [vtkCells;p,L2+1+p-1,L2+1+p];
                vtkCells = [vtkCells;p,p+1,L2+1+p];
            end
        end

        if L2 == min(L2,L3) && L3 > min(L2,L3)
            char_w_vertices{1} = [char_w_vertices{1}, struct('coor', char_approx{1}(:,L2+1),'prev_vert',[],'act_vert',[char_approx{1}(:,L2)],'next_vert',[char_approx{2}(:,L2),char_approx{2}(:,L2+1)],'sol',0)];
            vtkCells = [vtkCells;L2+1,L2+1+L2,L2+1+L2+1];

        elseif L2 > min(L2,L3) && L3 == min(L2,L3)
            for p=min(L2,L3)+1:L2
                char_w_vertices{1} = [char_w_vertices{1}, struct('coor', char_approx{1}(:,p),'prev_vert',[],'act_vert',[char_approx{1}(:,p-1),char_approx{1}(:,p+1)],'next_vert',[char_approx{2}(:,L3),char_approx{2}(:,L3)],'sol',0)];
                vtkCells = [vtkCells;p,p+1,L2+1+L3];
            end
            char_w_vertices{1} = [char_w_vertices{1}, struct('coor', char_approx{1}(:,L2+1),'prev_vert',[],'act_vert',[char_approx{1}(:,L2)],'next_vert',[char_approx{2}(:,L3)],'sol',0)];
        end
    end

    prev_char = L2 + 1;
    % lines between
    for K=begin:length(char_approx)-1
        vtkPoints((prev_char + 1):(prev_char + length(char_approx{K}(1,:))),1:2) = (char_approx{K})';
        L1 = length(char_approx{K-1}(1,:))-1;  % length of previous char (-1)
        L2 = length(char_approx{K}(1,:))-1;  % length of actual one (-1)
        L3 = length(char_approx{K+1}(1,:));  % length of next one
        char_w_vertices{K} = [struct('coor', char_approx{K}(:,1), 'prev_vert',[char_approx{K-1}(:,1),char_approx{K-1}(:,2)],'act_vert',[char_approx{K}(:,2)],'next_vert',[char_approx{K+1}(:,1)],'sol',0)];
        vtkCells = [vtkCells;prev_char + 1,prev_char + 2,prev_char + L2+2];
        MIN = min([L1,L2,L3]);

        if MIN > 1  % creates regular triangles between chars
            for p=2:MIN
                char_w_vertices{K} = [char_w_vertices{K}, struct('coor', char_approx{K}(:,p), 'prev_vert',[char_approx{K-1}(:,p),char_approx{K-1}(:,p+1)],'act_vert',[char_approx{K}(:,p-1),char_approx{K}(:,p+1)],'next_vert',[char_approx{K+1}(:,p-1),char_approx{K+1}(:,p)],'sol',0)];
                vtkCells = [vtkCells;prev_char + p,prev_char + L2+1+p-1,prev_char + L2+1+p];
                vtkCells = [vtkCells;prev_char + p,prev_char + p+1,prev_char + L2+1+p];
            end
        end
        
        % special cases when one of the three chars ends earlier (and then also when one the of the remaining two ends earlier)
        if L1 == MIN && L2 == MIN && L3 == MIN
            char_w_vertices{K} = [char_w_vertices{K}, struct('coor', char_approx{K}(:,MIN+1),'prev_vert',[char_approx{K-1}(:,MIN+1)],'act_vert',[char_approx{K}(:,MIN)],'next_vert',[char_approx{K+1}(:,MIN)],'sol',0)];
        
        elseif L1 == MIN && L2 == MIN && L3 > MIN
            char_w_vertices{K} = [char_w_vertices{K}, struct('coor', char_approx{K}(:,MIN+1),'prev_vert',[char_approx{K-1}(:,MIN+1)],'act_vert',[char_approx{K}(:,MIN)],'next_vert',[char_approx{K+1}(:,MIN),char_approx{K+1}(:,MIN+1)],'sol',0)];
            vtkCells = [vtkCells;prev_char + MIN+1,prev_char + L2+1+MIN,prev_char + L2+1+MIN+1];
        
        elseif L1 > MIN && L2 == MIN && L3 == MIN
            char_w_vertices{K} = [char_w_vertices{K}, struct('coor', char_approx{K}(:,MIN+1),'prev_vert',[char_approx{K-1}(:,MIN+1),char_approx{K-1}(:,MIN+2)],'act_vert',[char_approx{K}(:,MIN)],'next_vert',[char_approx{K+1}(:,MIN)],'sol',0)];
        
        elseif L1 > MIN && L2 == MIN && L3 > MIN
            char_w_vertices{K} = [char_w_vertices{K}, struct('coor', char_approx{K}(:,MIN+1),'prev_vert',[char_approx{K-1}(:,MIN+1),char_approx{K-1}(:,MIN+2)],'act_vert',[char_approx{K}(:,MIN)],'next_vert',[char_approx{K+1}(:,MIN),char_approx{K+1}(:,MIN+1)],'sol',0)];
            vtkCells = [vtkCells;prev_char + MIN+1,prev_char + L2+1+MIN,prev_char + L2+1+MIN+1];
        
        elseif L1 == MIN && L2 > MIN && L3 == MIN
            for p=MIN+1:L2
                char_w_vertices{K} = [char_w_vertices{K}, struct('coor', char_approx{K}(:,p),'prev_vert',[char_approx{K-1}(:,MIN+1),char_approx{K-1}(:,MIN+1)],'act_vert',[char_approx{K}(:,p-1),char_approx{K}(:,p+1)],'next_vert',[char_approx{K+1}(:,MIN),char_approx{K+1}(:,MIN)],'sol',0)];
                vtkCells = [vtkCells;prev_char + p,prev_char + p+1,prev_char + L2+1+MIN];
                vtkCells = [vtkCells;prev_char,prev_char + p,prev_char + p+1];
            end
            char_w_vertices{K} = [char_w_vertices{K}, struct('coor', char_approx{K}(:,L2+1),'prev_vert',[char_approx{K-1}(:,MIN+1)],'act_vert',[char_approx{K}(:,L2)],'next_vert',[char_approx{K+1}(:,MIN)],'sol',0)];
        
        elseif L1 == MIN && L2 > MIN && L3 > MIN
            for p=MIN+1:min(L2,L3)
                char_w_vertices{K} = [char_w_vertices{K}, struct('coor', char_approx{K}(:,p),'prev_vert',[char_approx{K-1}(:,MIN+1),char_approx{K-1}(:,MIN+1)],'act_vert',[char_approx{K}(:,p-1),char_approx{K}(:,p+1)],'next_vert',[char_approx{K+1}(:,p-1),char_approx{K+1}(:,p)],'sol',0)];
                vtkCells = [vtkCells;prev_char + p,prev_char + L2+1+p-1,prev_char + L2+1+p];
                vtkCells = [vtkCells;prev_char + p,prev_char + p+1,prev_char + L2+1+p];
                vtkCells = [vtkCells;MIN+1,prev_char + p,prev_char + p+1];
            end

            if L2 == min(L2,L3) && L3 > min(L2,L3)
                char_w_vertices{K} = [char_w_vertices{K}, struct('coor', char_approx{K}(:,L2+1),'prev_vert',[char_approx{K-1}(:,MIN+1)],'act_vert',[char_approx{K}(:,L2)],'next_vert',[char_approx{K+1}(:,L2),char_approx{K+1}(:,L2+1)],'sol',0)];
                vtkCells = [vtkCells;prev_char + L2+1,prev_char + L2+1+L2,prev_char + L2+1+L2+1];
            
            elseif L2 > min(L2,L3) && L3 == min(L2,L3)
                for p=min(L2,L3)+1:L2
                    char_w_vertices{K} = [char_w_vertices{K}, struct('coor', char_approx{K}(:,p),'prev_vert',[char_approx{K-1}(:,MIN+1),char_approx{K-1}(:,MIN+1)],'act_vert',[char_approx{K}(:,p-1),char_approx{K}(:,p+1)],'next_vert',[char_approx{K+1}(:,L3),char_approx{K+1}(:,L3)],'sol',0)];
                    vtkCells = [vtkCells;prev_char + p,prev_char + p+1,prev_char + L2+1+L3];
                end
                char_w_vertices{K} = [char_w_vertices{K}, struct('coor', char_approx{K}(:,L2+1),'prev_vert',[char_approx{K-1}(:,MIN+1)],'act_vert',[char_approx{K}(:,L2)],'next_vert',[char_approx{K+1}(:,L3)],'sol',0)];
            end

        elseif L1 > MIN && L2 > MIN && L3 == MIN
            for p=MIN+1:min(L1,L2)
                char_w_vertices{K} = [char_w_vertices{K}, struct('coor', char_approx{K}(:,p),'prev_vert',[char_approx{K-1}(:,p),char_approx{K-1}(:,p+1)],'act_vert',[char_approx{K}(:,p-1),char_approx{K}(:,p+1)],'next_vert',[char_approx{K+1}(:,L3),char_approx{K+1}(:,L3)],'sol',0)];
                vtkCells = [vtkCells;prev_char + p,prev_char + p+1,prev_char + L2+1+L3];
            end

            if L1 == min(L1,L2) && L2 > min(L1,L2)
                for p=min(L1,L2)+1:L2
                    char_w_vertices{K} = [char_w_vertices{K}, struct('coor', char_approx{K}(:,p),'prev_vert',[char_approx{K-1}(:,L1+1),char_approx{K-1}(:,L1+1)],'act_vert',[char_approx{K}(:,p-1),char_approx{K}(:,p+1)],'next_vert',[char_approx{K+1}(:,L3),char_approx{K+1}(:,L3)],'sol',0)];
                    vtkCells = [vtkCells;prev_char + p,prev_char + p+1,prev_char + L2+1+L3];
                    vtkCells = [vtkCells;prev_char,prev_char + p,prev_char + p+1];
                end
                char_w_vertices{K} = [char_w_vertices{K}, struct('coor', char_approx{K}(:,L2+1),'prev_vert',[char_approx{K-1}(:,L1+1)],'act_vert',[char_approx{K}(:,L2)],'next_vert',[char_approx{K+1}(:,L3)],'sol',0)];
            
            elseif L1 > min(L1,L2) && L2 == min(L1,L2)
                char_w_vertices{K} = [char_w_vertices{K}, struct('coor', char_approx{K}(:,L2+1),'prev_vert',[char_approx{K-1}(:,L2+1),char_approx{K-1}(:,L2+2)],'act_vert',[char_approx{K}(:,L2)],'next_vert',[char_approx{K+1}(:,L3)],'sol',0)];
            end
        end
        prev_char = prev_char + L2+1;
    end

    
    line = char_approx{end};
    if length(vtkPoints((prev_char + 1):end,1)) == length(line(1,:))  % adding last values to vtk arrays
        vtkPoints((prev_char + 1):end,1:2) = (line)';
    else
        vtkPoints((prev_char + 1):end,1:2) = (line(:,1:end-1))';
        vtkPoints = [vtkPoints;(line(:,end))',0];
    end

    % last line
    if length(char_approx{end}(1,:)) == 1
        char_w_vertices{end} = [struct('coor', char_approx{end}(:,1), 'prev_vert',[char_approx{end-1}(:,1),char_approx{end-1}(:,2)],'act_vert',[],'next_vert',[],'sol',0)];
        
    else
        L1 = length(char_approx{end-1}(1,:))-1;
        L2 = length(char_approx{end}(1,:))-1;
        char_w_vertices{end} = [struct('coor', char_approx{end}(:,1), 'prev_vert',[char_approx{end-1}(:,1),char_approx{end-1}(:,2)],'act_vert',[char_approx{end}(:,2)],'next_vert',[],'sol',0)];
        
        if min(L1,L2) > 1
            for p=2:min(L1,L2)
                char_w_vertices{end} = [char_w_vertices{end}, struct('coor', char_approx{end}(:,p), 'prev_vert',[char_approx{end-1}(:,p),char_approx{end-1}(:,p+1)],'act_vert',[char_approx{end}(:,p-1),char_approx{end}(:,p+1)],'next_vert',[],'sol',0)];
            end
        end
        
        if L1 == min(L1,L2) && L2 > min(L1,L2)
            for p=min(L1,L2)+1:L2
                char_w_vertices{end} = [char_w_vertices{end}, struct('coor', char_approx{end}(:,p),'prev_vert',[char_approx{end-1}(:,L1+1),char_approx{end-1}(:,L1+1)],'act_vert',[char_approx{end}(:,p-1),char_approx{end}(:,p+1)],'next_vert',[],'sol',0)];
            end
            char_w_vertices{end} = [char_w_vertices{end}, struct('coor', char_approx{end}(:,L2+1),'prev_vert',[char_approx{end-1}(:,L1+1)],'act_vert',[char_approx{end}(:,L2)],'next_vert',[],'sol',0)];
        
        elseif L1 > min(L1,L2) && L2 == min(L1,L2)
            char_w_vertices{end} = [char_w_vertices{end}, struct('coor', char_approx{end}(:,L2+1),'prev_vert',[char_approx{end-1}(:,L2+1),char_approx{end-1}(:,L2+2)],'act_vert',[char_approx{end}(:,L2)],'next_vert',[],'sol',0)];
        end
    end

    vtkCells = vtkCells - 1;
    new_chars_w_vert = char_w_vertices;
end