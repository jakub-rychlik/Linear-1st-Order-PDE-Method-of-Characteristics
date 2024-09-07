% Function coefficients next to derivatives of u.
function res = B(x,y)
    format long
    global b1 b2;
    res = [b1(x,y); b2(x,y)]; 
end