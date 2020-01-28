clc;
clear;
clf;
format long e;

small_A= load('small_ex.mat');
med_A1= load('medium_ex1.mat');
med_A2= load('medium_ex2.mat');

nnzmax = 100000000000;
tlimit = 1000000000000;
[L, flag, k] = spcholesky(small_A.A, nnzmax, tlimit);

fprintf("\n nnz(L) = %d\n",nnz(L) );

fprintf("\n nonzero entries in L at column = %d\n", 4);
L(:,4)

fprintf("\n nonzero entries in L at column = %d\n", 7);
L(:,7)

fprintf("\n nonzero entries in L at column = %d\n", 10);
L(:,10)

fprintf("\n nonzero entries in L at column = %d\n", 13);
L(:,13)

fprintf("\n nonzero entries in L at column = %d\n", 16);
L(:,16)

err = normest(small_A.A - L*transpose(L))/normest(L);
fprintf("\n relative error = %e\n", err);