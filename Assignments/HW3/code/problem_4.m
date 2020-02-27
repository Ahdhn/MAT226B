set(0,'DefaultFigureWindowStyle','docked');
clc;
clear;
close all;
clf;



format long e;

matrix_1 = load('HW3_P4_1.mat');
run(matrix_1.A, matrix_1.b, 'Mat1_');

matrix_2 = load('HW3_P4_2.mat');
run(matrix_2.A, matrix_2.b, 'Mat2_');

function run(A, b, filePrefix)
    n = length(b);
    x0 = ones(n,1);
    tol = 1e-8;
    close all;
    clf;
        
    %########## (i) GMRES no precond     
    figure(1);
    fprintf('\n *** (i) GMRES Without Preconditioning ***');    
    numMatVec =0;
    [rho, x, iter] = run_gmres(@my_mat_vec, b, x0, tol, [], [], [],...
        'No Precond', '(i) GMRES Without Preconditioning',...
        sprintf('%s%d',filePrefix,1));    
    fprintf('\n Total number matrix-vector products = %d\n', numMatVec);  
    
    
    %########## (ii) Restart GMRES no precond   
    fprintf('\n *** (ii) Restart GMRES Without Preconditioning ***'); 
    figure(2);
    k=[5,10,20];
    for i=1:length(k)                
         numMatVec =0;
         [rho, x, iter] = run_gmres(@my_mat_vec, b, x0, tol, k(i),[], [],...
             'Restart No Precond','(ii) Restart GMRES Without Preconditioning', ...
             sprintf('%s%d',filePrefix,2));        
        fprintf('\n Total number matrix-vector products = %d\n', numMatVec);
    end   
    
    
    %########## (iii) GMRES with diag precond 
    figure(3);
    fprintf('\n *** (iii) GMRES With Diagonal Preconditioning ***'); 
    D0 = diag(diag(A));
    I = speye(n);
    numMatVec =0;
    [rho, x, iter] = run_gmres(@my_mat_vec, b, x0, tol, [], I, D0,...
        'Diag Precond','(iii) GMRES With Diagonal Preconditioning',...
        sprintf('%s%d',filePrefix,3));    
    fprintf('\n Total number matrix-vector products = %d\n', numMatVec);
    
    
    %########## (iv) Restart GMRES with diag precond 
    figure(4);
    fprintf('\n *** (iv) Restart GMRES With Diagonal Preconditioning ***'); 
    for i=1:length(k)                
        numMatVec =0;
        [rho, x, iter] = run_gmres(@my_mat_vec, b, x0, tol, k(i),I, D0,...
            'Restart Diag Precond',...
            '(iv) Restart GMRES With Diagonal Preconditioning',...
            sprintf('%s%d',filePrefix,4));            
        fprintf('\n Total number matrix-vector products = %d\n', numMatVec);
    end

    
    %########## (v) GMRES precond from 3 with D = D0
    figure(5);
    fprintf('\n *** (v) GMRES With Preconditioning from Problem 3, D=D0 ***');     
    numMatVec =0;
    [D, D1, L, U, b_prime, x0_prime] = precond_setup(D0);    
    [rho, x, iter] = run_gmres(@my_mat_vec_prime, b_prime, x0_prime, tol,...
        [],[],[], 'Precond From 3, D=D0',...
        '(v) GMRES With Preconditioning from Problem 3, D=D0',...
        sprintf('%s%d',filePrefix,5));   
    fprintf('\n Total number matrix-vector products = %d\n', numMatVec);
    
    
    %########## (vi) Restart GMRES with prcond from 3 with D = D0
    fprintf('\n *** (vi) Restart GMRES With Preconditioning from Problem 3, D=D0 ***'); 
    figure(6);    
    for i=1:length(k)                
        numMatVec =0;
        [rho, x, iter] = run_gmres(@my_mat_vec_prime, b_prime, x0_prime,...
            tol, k(i),[] , [] , 'Restart Precond From 3, D=D0',...
            '(vi) Restart GMRES With Preconditioning from Problem 3, D=D0',...
            sprintf('%s%d',filePrefix,6));            
        fprintf('\n Total number matrix-vector products = %d\n', numMatVec);
    end

    
    %########## (v) GMRES precond from 3 with D = 10I
    fprintf('\n *** (v) GMRES With Preconditioning from Problem 3, D=10I ***'); 
    figure(7);    
    numMatVec = 0;
    [D, D1, L, U, b_prime, x0_prime] = precond_setup(10*I);
    [rho, x, iter] = run_gmres(@my_mat_vec_prime, b_prime, x0_prime, tol,...
        [],[],[], 'Precond From 3, D=10I', ...
        '(v) GMRES With Preconditioning from Problem 3, D=10I',...
        sprintf('%s%d',filePrefix,7));   
    fprintf('\n Total number matrix-vector products = %d\n', numMatVec);
        
    
    %########## (vi) Restart GMRES with prcond from 3 with D = 10I
    fprintf('\n *** (vi) Restart GMRES With Preconditioning from Problem 3, D=10I ***'); 
    figure(8);    
    for i=1:length(k)                
        numMatVec =0;
        [rho, x, iter] = run_gmres(@my_mat_vec_prime, b_prime, x0_prime,...
            tol, k(i),[] , [] , 'Restart Precond From 3, D=10I', ...
            '(vi) Restart GMRES With Preconditioning from Problem 3, D=10I',...
            sprintf('%s%d',filePrefix,8));            
        fprintf('\n Total number matrix-vector products = %d\n', numMatVec);
    end        
    
    %%%%%%%%%%%%%%%%%%%% Helper Function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function out = my_mat_vec(b)
        out = A*b;
        numMatVec = numMatVec +1;
    end 
    
    function out = my_mat_vec_prime(z)
        us = U\z;
        d = diag(D1).*us;
        v = z + d;
        ls = L\v;
        sum = ls + us;
        out = D*sum;        
        numMatVec = numMatVec +1;
    end 

    function [D, D1, L, U, b_prime, x0_prime] = precond_setup(D)
        F = -tril(A,-1);
        G = -triu(A,1);
        D1 = D0 - 2*D;
        L = D - F;
        U = D - G;
        M1 = L*D^(-1);
        M2=U;
        b_prime = M1\b;
        x0_prime = M2*x0;
    end    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end 



