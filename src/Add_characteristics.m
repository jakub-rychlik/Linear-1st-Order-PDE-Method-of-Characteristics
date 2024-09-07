% Add_characteristics creates more characteristics if needed, so the estimates for distances between them are satisfied.
function [new_chars,new_en_bdd,new_point_dist] = Add_characteristics(chars,bdd,tol,bdd_steps,I,eps,parametrization_point_dist)
    added_lines = 0;
    h=0.05;  %iteration step for Runge-Kutta
    char_approx = chars;
    disc_en_bdd = bdd;
    new_point_dist = parametrization_point_dist;

    for K=2:length(disc_en_bdd(1,:))
        Lchar = char_approx{K-1+added_lines};
        Mchar = char_approx{K+added_lines};

        if norm(Lchar(:,end) - Mchar(:,end)) > tol  %checking whether the estime for distance of two lines is not satisfied
            L=bdd_steps(1,K-1);
            l=bdd_steps(2,K-1);
            M=bdd_steps(1,K);  
            m=bdd_steps(2,K);
            if L == M
                n = m-l;
            else
                n = eps-l+m;
            end

            lowest_quantity_of_new_char = ceil(norm(Lchar(:,end) - Mchar(:,end)))/tol;
            steps_ratio = (n)/(lowest_quantity_of_new_char);
            if floor(steps_ratio) > 0
                cons = eps;
            else
                cons = eps + floor(1/steps_ratio)*eps;
            end

            i = 1;
            interval_leap = floor((n*(1+floor(1/steps_ratio)))/lowest_quantity_of_new_char);  %new coefficient for creating new points on entry boundary
            while norm(Lchar(:,end) - Mchar(:,end)) > tol
                j = i*interval_leap;

                while norm(Lchar(:,end) - Mchar(:,end)) > tol
                    new_en_bdd_point = Entry_bdd(I(L)+l*(tol/eps)+(j*(tol/cons)));
                    [Mchar, point_dist] = Create_characteristic(new_en_bdd_point,h,tol);  %creating new char

                    if norm(Lchar(:,end) - Mchar(:,end)) > tol  %the estimate for new curve still not satisfied -->find some closer one
                        j = j-1;
                    else  %it is satisfied --> saving the entry boundary point and new char to our memory
                        disc_en_bdd = [disc_en_bdd(:,1:bdd_steps(3,K-1)+added_lines),new_en_bdd_point,disc_en_bdd(:,bdd_steps(3,K)+added_lines:end)];
                        char_approx = {char_approx{1:K-1+added_lines},Mchar,char_approx{K+added_lines:end}};
                        new_point_dist = {new_point_dist{1:K-1+added_lines},point_dist,new_point_dist{K+added_lines:end}};
                        added_lines = added_lines+1;
                        Lchar = Mchar;
                        Mchar = char_approx{K+added_lines};
                        i = i+1;
                        j = i*interval_leap;
                    end

                    if j == (i-1)*interval_leap | (I(L)+l*(tol/eps)+(j*(tol/cons)) > I(M)+m*(tol/eps))  %can't do any more steps, saving the closest possible curve found
                        disc_en_bdd = [disc_en_bdd(:,1:bdd_steps(3,K-1)+added_lines),new_en_bdd_point,disc_en_bdd(:,bdd_steps(3,K)+added_lines:end)];
                        char_approx = {char_approx{1:K-1+added_lines},Mchar,char_approx{K+added_lines:end}};
                        new_point_dist = {new_point_dist{1:K-1+added_lines},point_dist,new_point_dist{K+added_lines:end}};
                        added_lines = added_lines+1;
                        Lchar = Mchar;
                        Mchar = char_approx{K+added_lines};
                        break
                    end  
                end
                i = i+1;
            end
        end
    end
    
    new_chars = char_approx;
    new_en_bdd = disc_en_bdd;
end