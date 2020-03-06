set(0,'DefaultFigureWindowStyle','docked');
clc;
clear;
close all;
clf;
format long e;

matrix = load('HW4_P2.mat');

fprintf('\n A eigenvalues: \n')
A_eigs = eigs(matrix.A)

ks = [5, 10, 20];
for k = 1:length(ks)
    [V, H] = arnoldi_process(matrix.A, matrix.r, ks(k));
    H_k = size(H);
    H_k(1) = H_k(1)-1;
    H_eigs = eigs(H(1:H_k,1:H_k));  
    fprintf('\n **** K= %d',ks(k));
    fprintf('\n H_k eigenvalues: \n')
    H_eigs
    fprintf('\n=================\n');
    if ks(k) == 20
        fprintf('\n || H_20_eig - A_eig|| =  \n');
        norm(H_eigs-A_eigs)
    end
end



function [V, H] = arnoldi_process(A, r, kmax)
    n = length(r);
    beta = norm(r);
    V = zeros(n, kmax+1);
    V(:,1) = r/beta;
    H = zeros(kmax+1, kmax);
    for k = 1:kmax
        q = A*V(:,k);
        for j=1:k
            H(j,k) =V(:,j)'*q;
            q = q - H(j,k)*V(:,j);
        end
        H(k+1,k) = norm(q);
        if H(k+1,k)==0
            break;
        end
        V(:,k+1)=q/H(k+1,k);
    end 
end 