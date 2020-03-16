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


s0 = 100^1 + 2.0*pi*1i*5.5*10^8;
s = 2*pi*1i*10^8;
k = 5;

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
test = diag(1./magic(length(b)));
q = mat_vec(test);
norm(q - M*test)
qt = trans_mat_vec(test); 
norm(qt - transpose(M)*test)
%}

%% Test computeMoments 
%{
mu = computeMoments(L, U, P, Q, inv_D, E, c, r, k);

M = inv(W)*E;
rr = -inv(W)*b;
transpose(c)*rr
transpose(c)*M*rr
transpose(c)*M*M*rr
transpose(c)*M*M*M*rr
transpose(c)*M*M*M*M*rr
%}

%% Test textbookAlgo 
mu = textbookAlgo(L, U, P, Q, inv_D, E, c, r, k);

M = inv(W)*E;
rr = -inv(W)*b;
transpose(c)*rr
transpose(c)*M*rr
transpose(c)*M*M*rr
transpose(c)*M*M*M*rr
transpose(c)*M*M*M*M*rr



%% Test zkViaLanczos

zkViaLanczos(mat_vec, trans_mat_vec, r, c, k, s, s0);