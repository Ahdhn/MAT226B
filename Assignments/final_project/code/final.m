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


s0 = 1000 + 2.0*pi*1i*5.5*10^8;

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


f = 2*pi*1i.*linspace(10^8,10^9,1001)';
Zk_textbook = zeros(1001,1);
Zk_lanczos = zeros(1001,1);
Zs = zeros(1001,1);

k = 10;
[alpha, beta] = textbookAlgo(L, U, P, Q, inv_D, E, r, c, k);
Tk = nonsymmetricLanczos(mat_vec, trans_mat_vec, r, c, k);    
for ff =1:length(f)
    s = f(ff);
    
    %Exact
    Zs(ff) = transpose(c)*(inv(s.*E-A)*b);    
    
    %Textbook
    Zk_textbook(ff) = polyval(flip(alpha), s-s0)/polyval(flip(beta),s-s0);    
    
    %Lanczos
    Zk_lanczos(ff) = zkViaLanczos(Tk,c,r, s, s0);    
end




