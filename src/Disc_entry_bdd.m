% Discretize the entry boundary. It will try to find partition such that the distance of
% points on the entry boundary is bounded by constant 'tol'
function [disc_en_bdd, bdd_steps] = Disc_entry_bdd(I,tol,eps)
    format long
    bdd_steps1 = [1;0;1];  % creating array for saving various coefficients used for getting new point on disc_en_bdd
    disc_en_bdd1 = [Entry_bdd(I(1))];
    l=2;
    j = 250;

    for k=2:length(I)  % iterating over interval parametrizing curve which defines entry boundary
        p1 = Entry_bdd(I(k));
        p2 = disc_en_bdd1(:,end);

        if norm(p2 - p1) > tol  % points are too far --> tries to find closer ones
            while norm(p2 - p1) > tol
                i = 1;
                while norm(p2 - p1) > tol
                    p1 = Entry_bdd(I(k)-(i*(tol/eps)));
                    i = i+1;
                    if i == j
                        break
                    end
                end

                j = i-1;
                if i == 2
                    break
                end

                bdd_steps1 = [bdd_steps1,[k;eps-i;l]];
                l=l+1;
                disc_en_bdd1 = [disc_en_bdd1, p1];
                p2 = p1;
                p1 = Entry_bdd(I(k));
            end

            if ismembertol(p1,Entry_bdd(I(end)))
                bdd_steps1 = [bdd_steps1,[k;0;l]];
                l=l+1;
                disc_en_bdd1 = [disc_en_bdd1, p1];
            end
        
        else  % points already had good distance --> add the new one and move on
            bdd_steps1 = [bdd_steps1,[k;0;l]];
            l=l+1;
            disc_en_bdd1 = [disc_en_bdd1, p1];
        end 
    end
    disc_en_bdd = disc_en_bdd1;
    bdd_steps = bdd_steps1;
end