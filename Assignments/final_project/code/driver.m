set(0,'DefaultFigureWindowStyle','docked');
clc;
clear;
close all;
clf;
format long e;

global A E b c;

load('FP_Ex1.mat');

Figure_1();

Table_1();

Table_2();

load('FP_Ex2.mat');
Figure_2();

function Figure_1()
    
    f_min = 10^8;
    f_max = 8*10^9; 
    num_data = 1001;
    s0 = 10^5 + 2.0*pi*1i*5.5*10^8;
    
    [Zk_lanczos, time_lanczos, Zk_textbook, time_textbook, Zs] =...
            run(s0, 10 , f_min, f_max, num_data, true, true);
    figure;
    f = linspace(f_min,f_max,num_data)';
    semilogy(f, Zs,'LineWidth', 3.0);
    hold on;
    semilogy(f, Zk_textbook);
    hold on;
    
    [Zk_lanczos, time_lanczos, Zk_textbook, time_textbook, Zs] =...
            run(s0, 100, f_min, f_max, num_data, false, false);
        
    H = semilogy(f, Zk_lanczos, 'LineWidth', 1.5);    
    legend('Exact', 'Textbook Algo.', 'Lanczos-based Algo.', 'Location','best');
    xlabel(gca, 'f');
    ylabel(gca,'log_{10}|Z(s)|');
    
    saveas(H, 'figure1.png');
end

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
        10^5 + 2.0*pi*1i*f_max, 10^10, 10^9];    
    for ss = 1:length(s0_vector)       
        s0 = s0_vector(ss);        
        for k = 200:1000
            [Zk_lanczos, time_lanczos] =...
                run(s0, k, f_min, f_max, num_data, false, false);
            diff = Zk_lanczos - Zs; 
            
            if norm(diff)^2 < 10^-5
                %figure(ss);                
                %f = linspace(f_min,f_max,num_data)';
                %semilogy(f, Zs,'LineWidth', 2.0);
                %hold on;       
                %semilogy(f, Zk_lanczos, '*');
                %hold off;
                fprintf('\n s0= %e + %ei, K= %d, Norm Diff = %e, LanczosTime= %e\n',...
                     real(s0), imag(s0), k, norm(diff)^2, time_lanczos);                 
                break;
            end       
        end
    end    
     
end 

function Figure_2()
    f_min = 10^9;
    f_max = 4*10^10; 
    num_data = 1001;
    
    s0 = 10^5 + 2.0*pi*1i*2.05*10^10;
    k = 100;
    
    [Zk_lanczos, time_lanczos] =...
            run(s0, k , f_min, f_max, num_data, false, false);
        
    figure;
    f = linspace(f_min,f_max,num_data)'; 
    H = loglog(f, Zk_lanczos, 'LineWidth', 1.5);
    legend('Lanczos-based Algo.', 'Location','best');
    xlabel(gca, 'f');
    ylabel(gca,'log_{10}|Z_{k}(s)|');
    
    saveas(H, 'figure2.png');
end

function [Zk_lanczos, time_lanczos, Zk_textbook, time_textbook, Zs] = ...
    run(s0, k, f_min, f_max, num_data, run_textbook, run_exact)

    global A E b c;
   
    f = linspace(f_min,f_max,num_data)';
    s_vector = 2*pi*1i.*f;
    
    W = A - s0.*E;
    
    [L, U, P, Q, D] = lu(W);  
    inv_D = inv(D);
    r = computeR( L, U, P, Q, inv_D, b);
    
    t = cputime;
    Zk_lanczos = runLanczos(num_data, k, s0, r, c, s_vector, L, U, P, Q, inv_D, E);
    time_lanczos = cputime-t;
    
    Zk_textbook = 0;
    time_textbook = 0;
    if run_textbook
        t = cputime;
        Zk_textbook = runTextbook(num_data,k, s0, r, c, s_vector, L, U, P, Q, inv_D, E);
        time_textbook = cputime-t;
    end 
    
    Zs = 0;
    if run_exact
        Zs = runExact(num_data, A, E, b, c, s_vector);
    end
end


