set(0,'DefaultFigureWindowStyle','docked');
clc;
clear;
close all;
clf;
format long e;

global A E b c;

load('FP_Ex1.mat');

Table_2();

Table_1();

%f = linspace(f_min,f_max,num_data)';
%semilogy(f, Zs);
%hold on;
%semilogy(f, Zk_textbook, '*');
%hold on;
%semilogy(f, Zk_lanczos, 'LineWidth', 3.0);
%hold off;

function Table_1()
    fprintf('\n Problem 5:');
    fprintf('\n Generating Table 1 data: \n');
    f_min = 10^8;
    f_max = 8*10^9; 
    num_data = 1001;

    s0 = 10^5 + 2.0*pi*1i*5.5*10^8;
    
    for k = 2:30
        [Zk_lanczos, time_lanczos, Zk_textbook, time_textbook, Zs] =...
            run(s0, k, f_min, f_max, num_data, true, false);
        diff = Zk_lanczos - Zk_textbook;
        fprintf('\n K= %d, Norm Diff = %e, Max Diff= %e\n', k, norm(diff)^2,...
            max(abs(diff)));
    end 
end

function Table_2()
    fprintf('\n Problem 5:');
    fprintf('\n Generating Table 2 data: \n');
    
    f_min = 10^8;
    f_max = 8*10^9; 
    num_data = 1001;   
    
    
    %compute Zs one time 
     s0 = 10^5;
    [Zk_lanczos, time_lanczos, Zk_textbook, time_textbook, Zs] =...
            run(s0, 2, f_min, f_max, num_data, false, true);
    s0_vector = [10^5 + 2.0*pi*1i*5.5*10^8, 10^5 + 2.0*pi*1i*f_min,...
        10^5 + 2.0*pi*1i*f_max, 10^1, 10^5, 10^10] ;   
    for ss = 1:length(s0_vector)
        %s0 = 10^5 + 2.0*pi*1i*5.5*10^8; 
        s0 = s0_vector(ss);
        
        for k = 200:1000
            [Zk_lanczos, time_lanczos] =...
                run(s0, k, f_min, f_max, num_data, false, false);
            diff = Zk_lanczos - Zs; 
            
            if norm(diff)^2 < 10^-5
                figure(ss);
                
                f = linspace(f_min,f_max,num_data)';
                semilogy(f, Zs,'LineWidth', 2.0);
                hold on;       
                semilogy(f, Zk_lanczos, '*');
                hold off;
                fprintf('\n s0= %e + %ei, K= %d, Norm Diff = %e, LanczosTime= %e\n',...
                     real(s0), imag(s0), k, norm(diff)^2, time_lanczos);
                 legend('Exact', 'Lanczos Approach');
                break;
            end       
        end
    end    
     
end 

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
    time_textbook = 0;
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


