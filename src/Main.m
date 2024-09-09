% Main function for my bachelor thesis. It initializes all the data and
% then calls function MoC, which will find approximate solution to the equation.

format long
global b1 b2 c f g set set_bdd Act_en_bdd S;
b1 = inputdlg("b1(x,y) = ");  % initialize b1
b1 = str2func(convertCharsToStrings(append('@(x,y)',b1{1})));

b2 = inputdlg("b2(x,y) = ");  % initialize b2
b2 = str2func(convertCharsToStrings(append('@(x,y)',b2{1})));

c = inputdlg("c(x,y) = ");  % initialize c
if c{1} == '0'
    c = 0;
else
    c = str2func(convertCharsToStrings(append('@(x,y)',c{1})));
end

f = inputdlg("f(x,y) = ");  % initialize f
if f{1} == '0'
    f = 0;
else
    f = str2func(convertCharsToStrings(append('@(x,y)',f{1})));
end

k = inputdlg("How many set defining functions?");  % initialize defining equations of Omega
set = cell(str2double(convertCharsToStrings(k{1})));
set_bdd = cell(str2double(convertCharsToStrings(k{1})));
for j=1:length(set)
    temp = inputdlg(num2str(j));
    ineq = inputdlg("for g < f < h write 'g,h'; if unbounded then Inf or -Inf):");
    ineq = split(convertCharsToStrings(ineq{1}),",");
    less = str2func(append("@(x,y)",ineq(1)));
    great = str2func(append("@(x,y)",ineq(2)));
    set{j}= struct('func',str2func(convertCharsToStrings(append('@(x,y)',temp{1}))),'less',less,'great',great);
    if ineq(1) == "Inf" || ineq(1) == "-Inf"
        set_bdd{j} = struct('func',str2func(convertCharsToStrings(append('@(x,y)',temp{1}))),'equals',great);
    else
        set_bdd{j} = struct('func',str2func(convertCharsToStrings(append('@(x,y)',temp{1}))),'equals',less);
    end
end

k = inputdlg("How many boundary functions?");  % initialize parametrized boundary of Omega
S=cell(1,str2double(convertCharsToStrings(k{1})));
for j=1:length(S)
    temp = inputdlg("(x,y)=(v(s),w(s)), write 'v,w'");
    temp = split(convertCharsToStrings(temp{1}),",");
    interval = inputdlg("Interval [a,b], write as 'a,b'");
    interval = split(convertCharsToStrings(interval{1}),",");
    I = [str2num(interval(1));str2num(interval(2))];
    v = str2func(append("@(s)",temp(1)));
    w = str2func(append("@(s)",temp(2)));
    S{j} = struct('xfunc',v,'yfunc',w,'int',I);
end

k = inputdlg("How many ENTRY boundary functions?");  % initialize parametrized entry boundary of Omega
Q=cell(1,str2double(convertCharsToStrings(k{1})));
for j=1:length(Q)
    temp = inputdlg("(x,y)=(v(s),w(s)), write 'v,w'");
    temp = split(convertCharsToStrings(temp{1}),",");
    interval = inputdlg("Interval [a,b], write as 'a,b'");
    interval = split(convertCharsToStrings(interval{1}),",");
    I = [str2num(interval(1));str2num(interval(2))];
    v = str2func(append("@(s)",temp(1)));
    w = str2func(append("@(s)",temp(2)));
    Q{j} = struct('xfunc',v,'yfunc',w,'int',I);
end

g = inputdlg("Initial value condition: g(x,y) = ");  % initial value condition g
g = str2func(convertCharsToStrings(append('@(x,y)',g{1})));

paraview = inputdlg("Do you want to create .vtk file for Paraview? (y/n)");  % permition for .vtk file generation
if paraview{1} == 'y'
    vtk = true;
else
    vtk = false;
end

adj_tol = inputdlg("Do you want to adjust toleration constant? (y/n)");  % possible adjustment of constant 'tol'
if adj_tol{1} == 'y'
    tol = inputdlg("tol = ");
    tol = str2num(convertCharsToStrings(tol{1}))+10^(-10);
else
    tol = 0.15+10^(-10);  % +10^(-10) because of double-precision in Matlab
end

tic
for k=1:length(Q)
    Act_en_bdd = Q{k};
    MoC(Q{k}.int,vtk,tol,k);
end
toc
