set(0,'DefaultFigureWindowStyle','docked');
clc;
clear;
close all;
clf;
format long e;

global L  U P Q inv_D E;
 
mat_vec = @(x) Mv(L, U, P, Q, inv_D, E,x);
trans_mat_vec = @(x) transposeMv(L, U, P, Q, inv_D, E, x);

load('FP_Ex1.mat');

%s0 =50000;
s0 = 500000;

W = A - s0.*E;

%% LU factorization 
%TODO use 'vector'????
[L, U, P, Q, D] = lu(W);
inv_D = inv(D);

%% Compute r efficiently 
s = L\(-P*inv_D*b);    
d = U\s;
r = Q*d;
    
%% Test Mv and transposeMv
%{
M = inv(W)*E;
test = diag(magic(length(b)));
q = mat_vec(test);
norm(q - M*test)
qt = trans_mat_vec(test); 
norm(qt - transpose(M)*test)
%}

%% Test computeMoments 
%{
M = inv(W)*E;
mu = computeMoments(L, U, P, Q, inv_D, E, c, b, 5);
r = -inv(W)*b;
transpose(c)*r
transpose(c)*M*r
transpose(c)*M*M*r
%transpose(c)*M*M*M*r
%transpose(c)*M*M*M*M*r
%}

%% Test textbookAlgo 
%[alpha, beta] = textbookAlgo(L, U, P, Q, inv_D, E, c, b, 5);

%% Test zkViaLanczos
zkViaLanczos(mat_vec, trans_mat_vec, r, c, 5);