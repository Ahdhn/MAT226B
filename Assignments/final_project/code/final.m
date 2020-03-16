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


s0 = 1000^1 + 2.0*pi*1i*5.5*10^8;
s = 2*pi*1i*10^8;
k = 5;

W = A - s0.*E;

%% LU factorization 
%TODO use 'vector'????
[L, U, P, Q, D] = lu(W);
inv_D = inv(D);

%% Compute r efficiently 
r = computeR( L, U, P, Q, inv_D, b);

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

Zs = transpose(c)*(inv(s0.*E-A)*b)
Zk_textbook = polyval(flip(mu),s-s0)

%% Test zkViaLanczos
Zk_lanczos = zkViaLanczos(mat_vec, trans_mat_vec, r, c, k, s, s0);
Zk_lanczos 