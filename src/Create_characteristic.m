%Creates characteristics using Runge-Kutta method of 4th order with starting iteration step 'h'
function [new_line, point_dist] = Create_characteristic(init_cond,h,tol)
    format long
    global set_bdd;
    x = [init_cond];  %starting point of char on entry boundary
    i = 1;
    point_dist = [];

    while Set_def(x(:,i)) || Set_bdd(x(:,i))  %repeat while we are still inside Omega or on the boundary                     
        k_1 = B(x(1,i), x(2,i));
        k_2 = B(x(1,i)+0.5*h*k_1(1,1), x(2,i)+0.5*h*k_1(2,1));
        k_3 = B((x(1,i)+0.5*h*k_2(1,1)), (x(2,i)+0.5*h*k_2(2,1)));
        k_4 = B((x(1,i)+k_3(1,1)*h), (x(2,i)+k_3(2,1)*h));
        step = (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;  %step which we will add to the latest approximated point on char
        
        if norm(step) > tol  %the next point would be too far from the last one --> try to find closer one on char
            counter = 0;
            while norm(step) > tol
                j = 1;
                h = h - (j*(h/150));

                while norm(step) > tol
                    k_1 = B(x(1,i), x(2,i));
                    k_2 = B(x(1,i)+0.5*h*k_1(1,1), x(2,i)+0.5*h*k_1(2,1));
                    k_3 = B((x(1,i)+0.5*h*k_2(1,1)), (x(2,i)+0.5*h*k_2(2,1)));
                    k_4 = B((x(1,i)+k_3(1,1)*h), (x(2,i)+k_3(2,1)*h));
                    step = (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;
                    j = j+1;

                    if j == 150  %partition is too sparse, the closest point found will be added to char
                        counter = counter+1;
                        break
                    end
                end

                if counter==20  %too many tries between two points --> break inside cycle and move on
                    break
                end

                x = [x, x(:,end) + step];
                point_dist = [point_dist, h];
                i = i+1;
                h = 0.05;
                k_1 = B(x(1,i), x(2,i));
                k_2 = B(x(1,i)+0.5*h*k_1(1,1), x(2,i)+0.5*h*k_1(2,1));
                k_3 = B((x(1,i)+0.5*h*k_2(1,1)), (x(2,i)+0.5*h*k_2(2,1)));
                k_4 = B((x(1,i)+k_3(1,1)*h), (x(2,i)+k_3(2,1)*h));
                step = (1/6)*(k_1+2*k_2+2*k_3+k_4)*h;

                if ~Set_def(x(:,i))  %out of bounds --> break cycle
                    break
                end

                if ~(norm(step) > tol)  %satisfying value of 'step' found
                    x = [x, x(:,end) + step];
                    point_dist = [point_dist, h];
                    i = i+1;
                    break
                end
            end

        elseif Set_def(x(:,i)) || Set_bdd(x(:,i))  %'step' was good from the beginning, create new point on char and move on
            if (i == 1) && (norm(step) == 0) 
                break
            end
            x = [x, x(:,end) + step];
            point_dist = [point_dist, h];
            i = i+1;
            
        else  %out of bounds --> finish creating this char
            if i == 1 
                x = [x, x(:,end) + step];
                point_dist = [point_dist, h];
            end
            break
        end
    end
    new_line = x;
end