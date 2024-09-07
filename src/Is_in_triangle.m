% Test whether P lies in triangle with vertices [ver1,ver2,ver3]
function bool = Is_in_triangle(ver1,ver2,ver3,P)
    bool = false;
    a_1 = ver2 - ver1;
    a_2 = ver3 - ver1;
    a_3 = ver3 - ver2;
    S = abs(det([a_1,a_2]))/2;
    s_1 = abs(det([a_1,P-ver1]))/2;
    s_2 = abs(det([a_2,P-ver1]))/2;
    s_3 = abs(det([a_3,P-ver2]))/2;
    if ismembertol(s_1+s_2+s_3,S,10^(-5))
        bool = true;
    end
end