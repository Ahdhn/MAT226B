set(0,'DefaultFigureWindowStyle','docked');
clc;
clear;
close all;
clf;
format long e;

global A A_prime U L D0 I;


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

A_prime_setup();
precond_setup();
eig_a_prime = eig(full(A_prime));
k_vector = int32([20, 60, 100, 200, 300 500]);
for i = 1:length(k_vector)
    k = k_vector(i);
    plot(eig_a_prime,'*');
    T = nonsymmetric_lanczos(@mat_vec_prime_eff,@mat_trans_vec_prime_eff,...
        r, c, k);
    [V, H] = arnoldi_process(A_prime, r, k);
    eig_t = eig(T);
    eig_h = eig(H(1:k,1:k));
    hold on
    plot(eig_t,'o');
    hold on
    HH = plot(eig_h,'d');
    legend('eig(A'')', 'Nonsymmetric Lanczos eig(T)','Arnoldi eig(H)','Location','best');
    hold off
    fileName = sprintf('p4b_k_%d%s',k,'.png');
    saveas(HH, fileName);
end 


%b)
%fprintf('\n\n ########## Part c):\n');
%load('HW4_P4c.mat');
load('HW4_P4b.mat');

precond_setup();

figure;
k_vector = int32([20 60 100 200 300 500 1000]);

for i =1:length(k_vector)
    k = k_vector(i);
	T = nonsymmetric_lanczos(@mat_vec_prime_eff, @mat_trans_vec_prime_eff,...
     r, c, k);
	eig_t = sort(eig(T)); 
    
    hold on;    
    if i == 1
       HH = plot(eig_t,'o');
    end 
    if i == 2
        HH = plot(eig_t,'*');
    end
    if i == 3
        HH = plot(eig_t,'+');
    end
    if i == 4
        HH = plot(eig_t,'d');
    end
    if i == 5
        HH = plot(eig_t,'^');
    end
    if i == 6
        HH = plot(eig_t,'h');
    end
    if i == 7
        HH = plot(eig_t,'p');
    end   
    
   if i > 1     
       fprintf('\n k= %d, difference = %f', k,...
           norm(eig_t(1:length(eig_prv)) - eig_prv));
   end
   eig_prv = eig_t;
end

legend('K=20', 'K=60', 'K=100', 'K=200', 'K=300', 'K=500', ...
    'K=1000','Location','best');

saveas(HH, 'p4c.png');

function out = mat_vec_prime_eff(z)  
    global U L D0;
    us = U\z;
    d = diag(-D0).*us;
    v = z + d;
    ls = L\v;
    sum = ls + us;
    out = D0*sum;            
end 

function out = mat_trans_vec_prime_eff(z)
    global I;
    trans_z = transpose(z);
    a_pr = mat_vec_prime_eff(I);
    sol = trans_z*a_pr;
    out = transpose(sol);
end

function A_prime_setup()
    global A A_prime;
    n = size(A);
    D0 = diag(diag(A));
    F = -tril(A,-1);
    G = -triu(A,1);
    L = D0 - F;
    U = D0 - G;
    I = speye(n(1));
    inv_L = L\I;	
    inv_D0 = D0\I;
    Q = inv_D0*U;
    S = Q\I;
    A_prime = inv_L*A*S;
end 

function precond_setup()
    global A L U I D0
    n = size(A);
    F = -tril(A,-1);
    G = -triu(A,1);
    D0 = diag(diag(A));
    L = D0 - F;
    U = D0 - G;
    I = speye(n(1));
end