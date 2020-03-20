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

%s0 = 1000 + 2.0*pi*1i*5.5*10^8;

s0 = 2.0*pi*10^8;


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

num_data = 501;
f = linspace(10^8,8*10^9,num_data)';
s_vector = 2*pi*1i.*f;
Zk_textbook = zeros(num_data,1);
Zk_lanczos = zeros(num_data,1);
Zs = zeros(num_data,1);

k = 20;
[alpha, beta] = textbookAlgo(L, U, P, Q, inv_D, E, r, c, k);
Tk = nonsymmetricLanczos(mat_vec, trans_mat_vec, r, c, k);   
I =speye(k);
e1 = zeros(k,1);
e1(1) = 1;
for ss =1:length(s_vector)
    s = s_vector(ss);       
    
    %Exact
    Zs(ss) = log10(abs(transpose(c)*(inv(s.*E-A)*b)));    
    
    %Textbook
    Zk_textbook(ss) = log10(abs(polyval(flip(alpha), s-s0)/...
        polyval(flip(beta),s-s0)));    
    
    %Lanczos
    Zk_lanczos(ss) = log10(abs(zkViaLanczos(Tk, I, e1, c,r, s, s0)));      
    
    %semilogy(f, Zs);
    %hold on;
    %semilogy(f(1:ss), Zk_textbook(1:ss), '*');
    %hold on;
    %semilogy(f(1:ss), Zk_lanczos(1:ss), 'o');
    %hold off;
end


semilogy(f, Zs);
hold on;
semilogy(f, Zk_textbook, '*');
hold on;
semilogy(f, Zk_lanczos, 'o');
hold off;

