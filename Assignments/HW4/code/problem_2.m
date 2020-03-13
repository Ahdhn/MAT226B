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
    H_eigs = eig(H(1:H_k,1:H_k));  
    fprintf('\n **** K= %d',ks(k));
    fprintf('\n H_k eigenvalues: \n')
    H_eigs
    fprintf('\n=================\n');
    if ks(k) == 20
        fprintf('\n || H_20_eig - A_eig|| =  \n');
        norm( eigs(H(1:H_k,1:H_k))-A_eigs)
    end
end



