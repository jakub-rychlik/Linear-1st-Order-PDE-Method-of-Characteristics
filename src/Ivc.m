% Initial value condition which returns value of point (t,s).
function ivc = Ivc(t,s)
    format long
    global g;
    ivc = g(t,s);
end