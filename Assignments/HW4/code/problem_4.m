set(0,'DefaultFigureWindowStyle','docked');
clc;
clear;
close all;
clf;
format long e;

global A U L D0 A_prime;

mat_vec = @(x) A*x;
mat_trans_vec = @(x) transpose(A)*x;

mat_vec_prime = @(x) A_prime*x;
mat_trans_vec_prime = @(x) transpose(A_prime)*x;

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
load('HW4_P4b.mat');


%preconditioner setup
D0 = diag(diag(A));
F = -tril(A,-1);
G = -triu(A,1);
L = D0 - F;
U = D0 - G;
I = speye(length(r));
inv_L = L\I;
inv_U = U\I;

A_prime = inv_L*A*inv(inv(D0)*U);

eig_a_prime = eig(full(A_prime));

k_vector = [20 60 100 200 300];
for i = 1:k_vector
    k = k_vector(i);
    plot(eig_a_prime,'*');
    T = nonsymmetric_lanczos(mat_vec_prime,mat_trans_vec_prime, r, c, k);
    [V, H] = arnoldi_process(A_prime, r, k);
    eig_t = eig(T);
    eig_h = eig(H(1:k,1:k));
    hold on
    plot(eig_t,'o');
    hold on
    HH = plot(eig_h,'d');
    legend('eig(A'')', 'Nonsymmetric Lanczos eig(T)','Arnoldi eig(H)','Location','best');
    hold off
    fileName = sprintf('p4b_k_%d',k,'.png');
    saveas(HH, fileName);
end 

