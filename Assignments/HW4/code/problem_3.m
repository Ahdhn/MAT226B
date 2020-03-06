set(0,'DefaultFigureWindowStyle','docked');
clc;
clear;
close all;
clf;
format long e;


%a)
A = make_3d_laplacian(3);
load('HW4_P3a.mat');

T = hermitian_lanczos(A,r,7);
fprintf('\n ***** Part a)\n Approximate Eigenvalues:');
approx_eigs = sort(eig(T))

fprintf('\n Exact Eigenvalues:');
exact_eigs = laplacian_exact_eig(3);
exact_eigs = sort(exact_eigs(:))

%b) 
fprintf('\n\n ***** Part b)');
load('HW4_P3b.mat');
A = make_3d_laplacian(64);
ks = [100, 200, 400, 800, 1200, 1600, 2000];
for k=1:length(ks)
    T = hermitian_lanczos(A,r,ks(k)); 
    fprintf('\n ## K= %d', ks(k));
    print_eigs(eig(T));
end

exact_eigs = laplacian_exact_eig(64);
fprintf('\n ## Exact Eigenvalues');
print_eigs(exact_eigs(:));


%%%%%%%%%%%%%%%%% Helper functions 
function print_eigs(a_eig)
    a_eig= sort(a_eig);
    fprintf('\n Smallest Eigenvalues');
    a_eig(1:10)
    fprintf('\n Largest Eigenvalues');
    a_eig(end-9:end)
     fprintf('\n');
end 

function lap_eig = laplacian_exact_eig(grid_size)
    [I,J,L] = meshgrid(1:grid_size,1:grid_size,1:grid_size);
    f = grid_size+1;
    lap_eig = 2.0.*(3.0 - cos(I.*pi/f) - cos(J.*pi/f) - cos(L.*pi/f));

end