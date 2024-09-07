% Returns true/false if X lies/doesn't lie inside of the Omega.
function bool = Set_def(X)     
    format long
    global set;
    x = X(1);
    y = X(2);
    bool = true;
    
    for k=1:length(set)  % test whether X lies outside of Omega
        if ~(set{k}.less(x,y) < set{k}.func(x,y)) || ~(set{k}.func(x,y) < set{k}.great(x,y))
            bool = false;
            break
        end
    end
end