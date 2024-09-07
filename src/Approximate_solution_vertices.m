%Finds approximate solution in vertices of the trinagles from triangle net.
function new_char_w_vertices = Approximate_solution_vertices(char_w_vertices,parametrization_point_dist,ivc)
    global f c;
    new_char_w_vertices = char_w_vertices;
    if isequal(f,0) && isequal(c,0)  %if c=f=0 --> solution is constant along chars, we get solution easily from ivc
        for k=1:length(new_char_w_vertices)
            line = new_char_w_vertices{k};
            for i=1:length(line)
                line(i).sol = ivc(k);
            end
            new_char_w_vertices{k} = line;
        end

    elseif isequal(f,0) && ~isequal(c,0)  %c is non-zero --> approximate integrals of c*u along chars
        for k=1:length(new_char_w_vertices)
            line = new_char_w_vertices{k};
            point_dist = parametrization_point_dist{k};
            const = ivc(k);
            line(1).sol = const;

            if length(line) == 1  %case when char has only one point
                new_char_w_vertices{k} = line;
                continue
            elseif length(line) == 2  %case when char has two points
                c_integral = (point_dist(1)/2)*c(line(1).coor(1),line(1).coor(2))*line(1).sol;
                line(2).sol = (const - c_integral)/(1+c(line(2).coor(1),line(2).coor(2))*point_dist(1)/2);
            else
                c_integral = (point_dist(1)/2)*c(line(1).coor(1),line(1).coor(2))*line(1).sol;
                for i=2:length(line)-1
                    line(i).sol = (const - c_integral)/(1+c(line(i).coor(1),line(i).coor(2))*point_dist(i-1)/2);
                    c_integral = c_integral + c(line(i).coor(1),line(i).coor(2))*point_dist(i-1)/2 + c(line(i).coor(1),line(i).coor(2))*point_dist(i)/2;
                end
                line(end).sol = (const-c_integral)/(1+c(line(end).coor(1),line(end).coor(2))*point_dist(end)/2);
            end
            new_char_w_vertices{k} = line;
        end

    elseif ~isequal(f,0) && isequal(c,0)  %f is non-zero --> approximate integrals of f along chars
        for k=1:length(new_char_w_vertices)
            line = new_char_w_vertices{k};
            point_dist = parametrization_point_dist{k};
            const = ivc(k);
            line(1).sol = const;

            if length(line) == 1
                new_char_w_vertices{k} = line;
                continue
            elseif length(line) == 2
                f_integral = (point_dist(1)/2)*(f(line(2).coor(1),line(2).coor(2)) + f(line(1).coor(1),line(1).coor(2)));
                line(2).sol = const + f_integral;
            else
                f_integral = (point_dist(1)/2)*(f(line(2).coor(1),line(2).coor(2)) + f(line(1).coor(1),line(1).coor(2)));
                for i=2:length(line)-1
                    line(i).sol = const + f_integral;
                    f_integral = f_integral + (point_dist(i)/2)*(f(line(i+1).coor(1),line(i+1).coor(2)) + f(line(i).coor(1),line(i).coor(2)));
                end
                line(end).sol = const + f_integral;
            end
            new_char_w_vertices{k} = line;
        end

    else  %c and f are non-zero --> approximate integrals of c*u and f along chars
        for k=1:length(new_char_w_vertices)
            line = new_char_w_vertices{k};
            point_dist = parametrization_point_dist{k};
            const = ivc(k);
            line(1).sol = const;

            if length(line) == 1
                new_char_w_vertices{k} = line;
                continue
            elseif length(line) == 2
                f_integral = (point_dist(1)/2)*(f(line(2).coor(1),line(2).coor(2)) + f(line(1).coor(1),line(1).coor(2)));
                c_integral = (point_dist(1)/2)*c(line(1).coor(1),line(1).coor(2))*line(1).sol;
                line(2).sol = (const + f_integral - c_integral)/(1+c(line(2).coor(1),line(2).coor(2))*point_dist(1)/2);
            else
                f_integral = (point_dist(1)/2)*(f(line(2).coor(1),line(2).coor(2)) + f(line(1).coor(1),line(1).coor(2)));
                c_integral = (point_dist(1)/2)*c(line(1).coor(1),line(1).coor(2))*line(1).sol;
                for i=2:length(line)-1
                    line(i).sol = (const + f_integral - c_integral)/(1+c(line(i).coor(1),line(i).coor(2))*point_dist(i-1)/2);
                    c_integral = c_integral + c(line(i).coor(1),line(i).coor(2))*point_dist(i-1)/2 + c(line(i).coor(1),line(i).coor(2))*point_dist(i)/2;
                    f_integral = f_integral + (point_dist(i)/2)*(f(line(i+1).coor(1),line(i+1).coor(2)) + f(line(i).coor(1),line(i).coor(2)));
                end
                line(end).sol = (const + f_integral - c_integral)/(1 + c(line(end).coor(1),line(end).coor(2))*point_dist(end)/2);
            end
            new_char_w_vertices{k} = line;
        end
    end
end
