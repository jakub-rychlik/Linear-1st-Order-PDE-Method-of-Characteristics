% Returns true/false if point X lies/doesn't lie on boundary of Omega.
function bool = Set_bdd(X)
    format long
    global set_bdd set;
    x = X(1);
    y = X(2);
    bool = false;

    for k=1:length(set)  % test whether X is out of closure of Omega
        if (set{k}.less(x,y) > set{k}.func(x,y)) || (set{k}.func(x,y) > set{k}.great(x,y))
            bool = false;
            return
        end
    end

    for k=1:length(set_bdd)  % test whether X lies on boundary of Omega
        if ismembertol(set_bdd{k}.func(x,y),set_bdd{k}.equals(x,y))
            bool = true;
            break
        end
    end
end