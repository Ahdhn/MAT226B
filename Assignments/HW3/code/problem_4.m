clc;
clear;
clf;

format long e;

matrix_1 = load('HW3_P4_1.mat');
%matrix_2 = load('HW3_P4_2.mat');


run(matrix_1.A, matrix_1.b);

function run(A, b)
    n = length(b);
    x0 = ones(n,1);
    tol = 1e-8;
    
    %i GMRES no precond     
    fprintf('\n *** (i) GMRES Without Preconditioning ***');
     numMatVec =0;
    [rho, x, iter] = gmres_no_precond(@my_mat_vec, b, x0, tol, 'No Precond');  
    fprintf('\n Total number matrix-vector products = %d\n', numMatVec);  
    
    
    %ii Restart GMRES no precond   
    k=[5,10,20];
    for i=1:length(k)
        fprintf('\n *** (ii) Restart GMRES Without Preconditioning ***'); 
         numMatVec =0;
        [rho, x, iter] = restart_gmres_no_precond(@my_mat_vec, b,x0, tol, k(i),'Restart');        
        fprintf('\n Total number matrix-vector products = %d\n', numMatVec);
    end
    
    %iii GMRES with diag precond 
    fprintf('\n *** (iii) GMRES With Diagonal Preconditioning ***'); 
    D0 = diag(diag(A));
    I = speye(n);
    numMatVec =0;
    [rho, x, iter] = gmres_precond(@my_mat_vec, b,x0, tol, I, D0, 'Diag Precond');
    fprintf('\n Total number matrix-vector products = %d\n', numMatVec);
    
    
    %iv Restart GMRES with diag precond 
    fprintf('\n *** (iv) Restart GMRES With Diagonal Preconditioning ***'); 
    for i=1:length(k)        
        fprintf('\n *** Restart GMRES Without Preconditioning ***');  
        numMatVec =0;
        [rho, x, iter] = restart_gmres_precond(@my_mat_vec, b,x0, tol,  k(i), I, D0, 'Restart Diag Precond');
        fprintf('\n Total number matrix-vector products = %d\n', numMatVec);
    end
    
    %v GMRES precond from 3 with D = D0
    fprintf('\n *** (v) GMRES With Preconditioning from Problem 3 and D=D0 ***');     
    numMatVec =0;
    
    
    %v GMRES precond from 3 with D = 10I
    fprintf('\n *** (v) GMRES With Preconditioning from Problem 3 and D=10I ***'); 
    
    %vi Restart GMRES with prcond from 3 with D = D0
    fprintf('\n *** (vi) Restart GMRES With Preconditioning from Problem 3 and D=D0 ***'); 
    
    
    %vi Restart GMRES with prcond from 3 with D = 10I
    fprintf('\n *** (vi) Restart GMRES With Preconditioning from Problem 3 and D=10I ***'); 
    numMatVec =0;
    
    function out = my_mat_vec(b)
        out = A*b;
        numMatVec = numMatVec +1;
    end 
    
    function out = my_mat_vec_prime(b)
        
    end 

    xlabel('Iteration Number');
    ylabel('Log Relative Residual');
    legend;   
end 



