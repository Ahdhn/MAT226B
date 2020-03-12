set(0,'DefaultFigureWindowStyle','docked');
clc;
clear;
close all;
clf;
format long e;

global A U D1 L L_t U_t D0;

mat_vec = @(x) A*x;
mat_trans_vec = @(x) transpose(A)*x;

%a)
fprintf('\n ########## Part a):\n');
load('HW4_P4a.mat');
T = nonsymmetric_lanczos(mat_vec,mat_trans_vec, r, c, 5);
fprintf('\n T eigenvalues for K = 5');
eig(T)
T = nonsymmetric_lanczos(mat_vec,mat_trans_vec, r, c, 10);
fprintf('\n T eigenvalues for K = 10');
eig(T)
T = nonsymmetric_lanczos(mat_vec,mat_trans_vec, r, c, 20);
fprintf('\n T eigenvalues for K = 20');
eig(T)
fprintf('\n norm(eig(T20)-eig(A)) = %e\n',...
    norm(sort(abs(eig(T))) - sort(abs(eig(A)))));

%b)
fprintf('\n\n ########## Part b):\n');
%load('HW4_P4b.mat');


%preconditioner setup
D0 = diag(diag(A));
F = -tril(A,-1);
G = -triu(A,1);
L = D0 - F;
U = D0 - G;
D1 = D0 - 2*D0;

L_t = transpose(D0) - transpose(F);
U_t = transpose(D0) - transpose(G);

rr = mat_vec_prime(r);

A_prime_trans = transpose(A_prime);

rr = mat_transpose_vec_prime(r);
rrr = A_prime_trans*r;

function out = mat_vec_prime(z)
    global U L D0;
    us = U\z;
	d = diag(-D0).*us;
	v = z + d;
	ls = L\v;
	sum = ls + us;
	out = D0*sum;            
end 

function out = mat_transpose_vec_prime(z)    
    global U_t L_t D0 D1
    Q1 = U_t\z;
    Q2 = L_t\z;
    Q3 = transpose(D1)*Q2;
    Q4 = U_t\z;
    out = transpose(D0)*(Q1 + (Q2 + Q4));
end