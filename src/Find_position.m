% Finds position of point V on i-th characteristics
function position = Find_position(V,i,chars)
    line = chars{i};
    for k=1:length(line)
        if ismembertol(line(k).coor,V)
            position = k;
            return
        end
    end
end