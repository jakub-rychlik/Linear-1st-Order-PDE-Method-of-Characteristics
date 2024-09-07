% Function approximates solution inside triangle using barycentric interpolation
function approx_solution = Approximate_solution_inside(Chars,P)
    approx_solution = NaN;
    max_iterations = 0;
    fixed_chars = Chars;
    
    if ~Set_def(P)  %if the point P is outside of Omega --> return Nan
        return
    end

    while isnan(approx_solution)  %loop until we find the triangle in which P lies
        [Vertex, index] = Find_closest_vertex(Chars,P);  %find the 1st/2nd/3rd/... closest point to P
        if ismembertol(Vertex.coor,P) %vertex of the triangle = P --> we know solution
            approx_solution = Vertex.sol;
            return
        end

        %Testing if P lies in one of the triangles which have 'Vertex' in common
        %If yes --> we approximate the solution in P
        if ~(length(Vertex.next_vert) == 0)
            if length(Vertex.next_vert(1,:)) == 2 && ~isequal(Vertex.next_vert(:,1),Vertex.next_vert(:,2))
                Check = Is_in_triangle(Vertex.coor,Vertex.next_vert(:,1),Vertex.next_vert(:,2),P);
                if Check
                    V_1 = fixed_chars{index(1)};
                    V_2 = fixed_chars{index(1)+1};
                    V_3 = fixed_chars{index(1)+1};
                    i_2 = Find_position(Vertex.next_vert(:,1),index(1)+1,fixed_chars);
                    i_3 = Find_position(Vertex.next_vert(:,2),index(1)+1,fixed_chars);
                    approx_solution = Barycentric_interp([V_1(index(2)).sol,V_2(i_2).sol,V_3(i_3).sol],Vertex.coor,Vertex.next_vert(:,1),Vertex.next_vert(:,2),P);
                    return
                end
            end

            if length(Vertex.act_vert(1,:)) == 2 && length(Vertex.next_vert(1,:)) == 2
                Check = Is_in_triangle(Vertex.coor,Vertex.act_vert(:,1),Vertex.next_vert(:,1),P);
                if Check
                    V_1 = fixed_chars{index(1)};
                    V_2 = fixed_chars{index(1)};
                    V_3 = fixed_chars{index(1)+1};
                    i_3 = Find_position(Vertex.next_vert(:,1),index(1)+1,fixed_chars);
                    approx_solution = Barycentric_interp([V_1(index(2)).sol,V_2(index(2)-1).sol,V_3(i_3).sol],Vertex.coor,Vertex.act_vert(:,1),Vertex.next_vert(:,1),P);
                    return
                end

                Check = Is_in_triangle(Vertex.coor,Vertex.act_vert(:,2),Vertex.next_vert(:,2),P);
                if Check
                    V_1 = fixed_chars{index(1)};
                    V_2 = fixed_chars{index(1)};
                    V_3 = fixed_chars{index(1)+1};
                    i_3 = Find_position(Vertex.next_vert(:,2),index(1)+1,fixed_chars);
                    approx_solution = Barycentric_interp([V_1(index(2)).sol,V_2(index(2)+1).sol,V_3(i_3).sol],Vertex.coor,Vertex.act_vert(:,2),Vertex.next_vert(:,2),P);
                    return
                end
            end

            if length(Vertex.next_vert(1,:)) == 1 && length(Vertex.act_vert(1,:)) == 1
                Check = Is_in_triangle(Vertex.coor,Vertex.act_vert(:,1),Vertex.next_vert(:,1),P);
                if Check
                    V_1 = fixed_chars{index(1)};
                    V_2 = fixed_chars{index(1)};
                    V_3 = fixed_chars{index(1)+1};
                    i_3 = Find_position(Vertex.next_vert(:,1),index(1)+1,fixed_chars);
                    approx_solution = Barycentric_interp([V_1(index(2)).sol,V_2(index(2)-1).sol,V_3(i_3).sol],Vertex.coor,Vertex.act_vert(:,1),Vertex.next_vert(:,1),P);
                    return
                end
            end

            if length(Vertex.next_vert(1,:)) == 2 && length(Vertex.act_vert(1,:)) == 1
                Check = Is_in_triangle(Vertex.coor,Vertex.act_vert(:,1),Vertex.next_vert(:,1),P);
                if Check
                    V_1 = fixed_chars{index(1)};
                    V_2 = fixed_chars{index(1)};
                    V_3 = fixed_chars{index(1)+1};
                    i_3 = Find_position(Vertex.next_vert(:,1),index(1)+1,fixed_chars);
                    approx_solution = Barycentric_interp([V_1(index(2)).sol,V_2(index(2)-1).sol,V_3(i_3).sol],Vertex.coor,Vertex.act_vert(:,1),Vertex.next_vert(:,1),P);
                    return
                end
            end
        end

        if ~(length(Vertex.prev_vert) == 0)
            if length(Vertex.prev_vert(1,:)) == 2 && ~isequal(Vertex.prev_vert(:,1),Vertex.prev_vert(:,2))
                Check = Is_in_triangle(Vertex.coor,Vertex.prev_vert(:,1),Vertex.prev_vert(:,2),P);
                if Check
                    V_1 = fixed_chars{index(1)};
                    V_2 = fixed_chars{index(1)-1};
                    V_3 = fixed_chars{index(1)-1};
                    i_2 = Find_position(Vertex.prev_vert(:,1),index(1)-1,fixed_chars);
                    i_3 = Find_position(Vertex.prev_vert(:,2),index(1)-1,fixed_chars);
                    approx_solution = Barycentric_interp([V_1(index(2)).sol,V_2(i_2).sol,V_3(i_3).sol],Vertex.coor,Vertex.prev_vert(:,1),Vertex.prev_vert(:,2),P);
                    return
                end
            end
    
            if length(Vertex.act_vert(1,:)) == 2 && length(Vertex.prev_vert(1,:)) == 2
                Check = Is_in_triangle(Vertex.coor,Vertex.act_vert(:,1),Vertex.prev_vert(:,1),P);
                if Check
                    V_1 = fixed_chars{index(1)};
                    V_2 = fixed_chars{index(1)};
                    V_3 = fixed_chars{index(1)-1};
                    i_3 = Find_position(Vertex.prev_vert(:,1),index(1)-1,fixed_chars);
                    approx_solution = Barycentric_interp([V_1(index(2)).sol,V_2(index(2)-1).sol,V_3(i_3).sol],Vertex.coor,Vertex.act_vert(:,1),Vertex.prev_vert(:,1),P);
                    return
                end

                Check = Is_in_triangle(Vertex.coor,Vertex.act_vert(:,2),Vertex.prev_vert(:,2),P);
                if Check
                    V_1 = fixed_chars{index(1)};
                    V_2 = fixed_chars{index(1)};
                    V_3 = fixed_chars{index(1)-1};
                    i_3 = Find_position(Vertex.prev_vert(:,2),index(1)-1,fixed_chars);
                    approx_solution = Barycentric_interp([V_1(index(2)).sol,V_2(index(2)+1).sol,V_3(i_3).sol],Vertex.coor,Vertex.act_vert(:,2),Vertex.prev_vert(:,2),P);
                    return
                end
            end

            if length(Vertex.prev_vert(1,:)) == 1 && length(Vertex.act_vert(1,:)) == 1
                Check = Is_in_triangle(Vertex.coor,Vertex.act_vert(:,1),Vertex.prev_vert(:,1),P);
                if Check
                    V_1 = fixed_chars{index(1)};
                    V_2 = fixed_chars{index(1)};
                    V_3 = fixed_chars{index(1)-1};
                    i_3 = Find_position(Vertex.prev_vert(:,1),index(1)-1,fixed_chars);
                    approx_solution = Barycentric_interp([V_1(index(2)).sol,V_2(index(2)-1).sol,V_3(i_3).sol],Vertex.coor,Vertex.act_vert(:,1),Vertex.prev_vert(:,1),P);
                    return
                end
            end

            if length(Vertex.prev_vert(1,:)) == 2 && length(Vertex.act_vert(1,:)) == 1
                if index(2) == 1
                    Check = Is_in_triangle(Vertex.coor,Vertex.act_vert(:,1),Vertex.prev_vert(:,2),P);
                    if Check
                        V_1 = fixed_chars{index(1)};
                        V_2 = fixed_chars{index(1)};
                        V_3 = fixed_chars{index(1)-1};
                        i_3 = Find_position(Vertex.prev_vert(:,2),index(1)-1,fixed_chars);
                        approx_solution = Barycentric_interp([V_1(index(2)).sol,V_2(index(2)-1).sol,V_3(i_3).sol],Vertex.coor,Vertex.act_vert(:,1),Vertex.prev_vert(:,2),P);
                        return
                    end
                else
                    Check = Is_in_triangle(Vertex.coor,Vertex.act_vert(:,1),Vertex.prev_vert(:,1),P);
                    if Check
                        V_1 = fixed_chars{index(1)};
                        V_2 = fixed_chars{index(1)};
                        V_3 = fixed_chars{index(1)-1};
                        i_3 = Find_position(Vertex.prev_vert(:,1),index(1)-1,fixed_chars);
                        approx_solution = Barycentric_interp([V_1(index(2)).sol,V_2(index(2)-1).sol,V_3(i_3).sol],Vertex.coor,Vertex.act_vert(:,1),Vertex.prev_vert(:,1),P);
                        return
                    end
                end
            end
        end

        if isnan(approx_solution)  %P did not lie in any of the triangles --> deleting Vertex from Chars and finding new one next iteration
            max_iterations = max_iterations + 1;
            if index(2) == length(Chars{index(1)})
                Chars{index(1)} = [Chars{index(1)}(1:index(2)-1)];
            else
                Chars{index(1)} = [Chars{index(1)}(1:index(2)-1),Chars{index(1)}(index(2)+1:end)];
            end
        end
        
        if max_iterations == 10
            break
        end
    end
end