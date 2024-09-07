% Finds closes vertex to the 'point' along all characteristics
function [close_vertex, index] = Find_closest_vertex(Chars,point)
    min_dist = 1000;
    close_vertex = [0;0];
    close_check = 0;
    index = [-1;-1];

    for k=1:length(Chars)
        if close_check==10  %too many iterations without any finding, point is probaly not in the Omega
            break
        end
        
        line = Chars{k};
        close_check = close_check + 1;

        for p=1:length(line)
            dist = norm(point-line(p).coor);
            if dist < min_dist  %new closer point found --> save it
                min_dist = dist;
                close_vertex = line(p);
                index = [k;p];  % number of char and the point on it for further usage
                close_check = 0;
            end
        end
    end
end