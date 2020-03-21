set(0,'DefaultFigureWindowStyle','docked');
clc;
clear;
close all;
clf;
format long e;

global A E b c;

load('FP_Ex1.mat');

f_min = 10^8;
f_max = 8*10^9;
num_data = 1001;
s0 = 10^5 + 2.0*pi*1i*5.5*10^8;
%s0 =  2.0*pi*5.5*10^8;
k = 100;

[Zk_lanczos, time_lanczos, Zk_textbook, time_textbook, Zs] =...
    run(s0, k, f_min, f_max, num_data, true, false);

time_lanczos
time_textbook

f = linspace(f_min,f_max,num_data)';
%semilogy(f, Zs);
%hold on;
semilogy(f, Zk_textbook, '*');
hold on;
semilogy(f, Zk_lanczos, 'LineWidth', 3.0);
hold off;


function [Zk_lanczos, time_lanczos, Zk_textbook, time_textbook, Zs] = ...
    run(s0, k, f_min, f_max, num_data, run_textbook, run_exact)

    global L  U P Q inv_D  A E b c;
    
    mat_vec = @(x) Mv(L, U, P, Q, inv_D, E,x);
    trans_mat_vec = @(x) transposeMv(L, U, P, Q, inv_D, E, x);
   
    f = linspace(f_min,f_max,num_data)';
    s_vector = 2*pi*1i.*f;
    
    W = A - s0.*E;
    
    %TODO use 'vector'????
    [L, U, P, Q, D] = lu(W);
    inv_D = inv(D);
    r = computeR( L, U, P, Q, inv_D, b);
    
    t = cputime;
    Zk_lanczos = runLanczos(num_data, k, s0, r, c, s_vector, mat_vec,...
        trans_mat_vec );
    time_lanczos = cputime-t;
    
    Zk_textbook = 0;
    if run_textbook
        t = cputime;
        Zk_textbook = runTextbook(num_data,k, s0, r, c, s_vector);
        time_textbook = cputime-t;
    end 
    
    Zs = 0;
    if run_exact
        Zs = runExact(num_data, A, E, b, c, s_vector);
    end
end


