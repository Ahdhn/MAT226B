clc;
clear;
clf;

format long e;

matrix_a = load('HW2_P1_a.mat');
matrix_b = load('HW2_P1_b.mat');


[Ja,Ia] = extract(matrix_a.A);

fprintf("\n Matrix HW2_P1_a.mat J(100000)= %d\n ", Ja(100000));
fprintf("\n Matrix HW2_P1_a.mat J(150000)= %d\n ", Ja(150000));
fprintf("\n Matrix HW2_P1_a.mat J(200000)= %d\n ", Ja(200000));
fprintf("\n Matrix HW2_P1_a.mat J(250000)= %d\n ", Ja(250000));
fprintf("\n Matrix HW2_P1_a.mat J(300000)= %d\n ", Ja(300000));

fprintf("\n Matrix HW2_P1_a.mat I(10000)= %d\n ", Ia(10000));
fprintf("\n Matrix HW2_P1_a.mat I(20000)= %d\n ", Ia(20000));
fprintf("\n Matrix HW2_P1_a.mat I(30000)= %d\n ", Ia(30000));
fprintf("\n Matrix HW2_P1_a.mat I(40000)= %d\n ", Ia(40000));
fprintf("\n Matrix HW2_P1_a.mat I(50000)= %d\n ", Ia(50000));

[Jb,Ib] = extract(matrix_b.A);

fprintf("\n\n Matrix HW2_P1_b.mat Vector I")
Ib

fprintf("\n Matrix HW2_P1_b.mat Vector J")
Jb


function [J,I] = extract(mat)
    [J,K]= find(mat);
    KK = [1; K(2:length(K)) ~= K(1:length(K)-1)];
    I = find(KK);
end