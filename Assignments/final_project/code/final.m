set(0,'DefaultFigureWindowStyle','docked');
clc;
clear;
close all;
clf;
format long e;

%mat_vec = @(x) A*x;
%mat_trans_vec = @(x) transpose(A)*x;

%mat_vec_prime = @(x) A_prime*x;
%mat_trans_vec_prime = @(x) transpose(A_prime)*x;

load('FP_Ex1.mat');

%s0 =50000;
s0 = 500000;

W = A - s0.*E;

%TODO use 'vector'????
[L, U, P, Q, D] = lu(W);
inv_D = inv(D);

%M = inv(W)*E;
%test = diag(magic(length(b)));
%q = Mv(L, U, P, Q, inv_D, E, test); 
%norm(q - M*test)
%qt = transposeMv(L, U, P, Q, inv_D, E, test); 
%norm(qt - transpose(M)*test)

%mu = computeMoments(L, U, P, Q, inv_D, E, c, b, 5);
%r = -inv(W)*b;
%transpose(c)*r
%transpose(c)*M*r
%transpose(c)*M*M*r

[alph, beta] = textbookAlgo(L, U, P, Q, inv_D, E, c, b, 4);