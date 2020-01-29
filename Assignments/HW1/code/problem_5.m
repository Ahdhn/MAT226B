clc;
clear;
clf;
format long e;

small_A= load('small_ex.mat');
med_A1= load('medium_ex1.mat');
med_A2= load('medium_ex2.mat');

nnzmax = 100000000000;
tlimit = 1000000000000;
[L, flag, k] = spcholesky(small_A.A, nnzmax, tlimit);

fprintf("\n nnz(L) = %d\n",nnz(L) );

fprintf("\n nonzero entries in L at column = %d\n", 4);
L(:,4)

fprintf("\n nonzero entries in L at column = %d\n", 7);
L(:,7)

fprintf("\n nonzero entries in L at column = %d\n", 10);
L(:,10)

fprintf("\n nonzero entries in L at column = %d\n", 13);
L(:,13)

fprintf("\n nonzero entries in L at column = %d\n", 16);
L(:,16)

err = normest(small_A.A - L*transpose(L))/normest(L);
fprintf("\n relative error = %e\n\n\n", err);

clc;

fprintf("\n         $$$$$$ medium_ex1 $$$$$$$\n ");
run_reorder(med_A1.A, 30, 200);

fprintf("\n\n\n         $$$$$$ medium_ex2 $$$$$$$\n ");
run_reorder(med_A2.A, 30, 200);

function run_reorder(A, nnzmax_factor, tlimit)
    run_chol(A, nnzmax_factor*nnz(A), tlimit,...
        '** No reorder **');
    fprintf("\n");
    
    symamd_P = symamd(A);
    run_chol(A(symamd_P,symamd_P), nnzmax_factor*nnz(A), tlimit, ...
        '** SYMAMD **');
    fprintf("\n");

    colamd_P = colamd(A);
    run_chol(A(colamd_P, colamd_P), nnzmax_factor*nnz(A), tlimit, ...
        '** COLAMD **');
    fprintf("\n");
    
    symrcm_P = symrcm(A);
    run_chol(A(symrcm_P, symrcm_P), nnzmax_factor*nnz(A), tlimit, ...
        '** SYMRCM **');
    fprintf("\n");
    
    colperm_P = colperm(A);
    run_chol(A(colperm_P, colperm_P), nnzmax_factor*nnz(A), tlimit, ...
        '** COLPERM **');
    fprintf("\n");
end 

function run_chol (A, nnzmax, tlimit, msg)
    [L, flag, k] = spcholesky(A, nnzmax, tlimit);
    fprintf('%s \n nnz(L) = %d, flag = %d, k= %d\n', msg, nnz(L), flag,k);
    L_large = abs(full(L(:,1000)));
    [L_large, L_large_row] = sort(L_large,'descend');
    fprintf('Ten largest entries in column 1000\n');    
    for id =1:10
        fprintf("\n L(%d, %d) = %e", L_large_row(id), 1000,L_large(id));
    end     
end



