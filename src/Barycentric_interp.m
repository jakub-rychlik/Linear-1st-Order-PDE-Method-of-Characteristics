%Returns interpolated value of P inside triangle [ver1,ver2,ver3].
%Input variable 'Values' represents values of the vertices.
function Result = Barycentric_interp(Values,ver1,ver2,ver3,P)
    a_1 = ver2 - ver1;
    a_2 = ver3 - ver1;
    a_3 = ver3 - ver2;
    S = abs(det([a_1,a_2]))/2;  %calculates volume of the main triangle 
    s_1 = abs(det([a_1,P-ver1]))/2;  %calculates volumes of the little ones
    s_2 = abs(det([a_2,P-ver1]))/2;
    s_3 = abs(det([a_3,P-ver2]))/2;    
    Result = (s_1/S)*Values(3) + (s_2/S)*Values(2) + (s_3/S)*Values(1);
end