clc;
clear;
clf;
format long e;

fprintf('\n****** Test Case (i)\n');
PoissonSolver(60, 0, 1, 1, 'Case_i');

fprintf('\n****** Test Case (ii)\n');
PoissonSolver(200, 1, 2, 1,'Case_ii');


fprintf('\n****** Test Case (iii)\n');
PoissonSolver(150, 2, 2, 3,'Case_iii');

fprintf('\n****** Test Case (iv)\n');
PoissonSolver(80, 5, 2, 3,'Case_iv');

fprintf('\n****** Test Case (v)\n');
PoissonSolver(200, 5, 3, 5,'Case_v');

function PoissonSolver(m, alpha, beta, gamma, plot_name)
    h = 1/(m+1);
    x = h:h:1-h;
    y = h:h:1-h;
    [X,Y] = meshgrid(x,y);
    
    %RHS
    %d2v_dy2 = -(alpha^2 * pi^2).*(X.^alpha).*cos((beta*pi).*X).*sin((alpha*pi).*Y);
    %t0 = ( (alpha*(alpha-1)).*X.^(alpha-2) + (beta^2*pi^2).*X.^alpha).*cos((beta*pi).*X);
    %t1 = -(2*alpha*beta*pi).*x.^(alpha-1).*sin((beta*pi).*X);
    %d2v_dx2 = (t0+t1).*sin((gamma*pi).*Y);
    %f = -d2v_dx2 -d2v_dy2;     
    f = (beta^2*pi^2).*X.^alpha.*sin((pi*gamma).*Y).*cos((pi*beta).*X)... 
    - (alpha*(alpha-1)).*(X.^(alpha - 2)).*sin((pi*gamma).*Y).*cos((pi*beta).*X)...
    + (gamma^2*pi^2).*(X.^alpha).*sin((pi*gamma).*Y).*cos((pi*beta).*X)...
    + (2*alpha*beta*pi).*(X.^(alpha - 1)).*sin((pi*gamma).*Y).*sin((pi*beta).*X);
    
    b0 = zeros(1,m);
    b1 = sin(gamma*pi).*(x.^alpha).*cos((beta*pi).*x);
    c0 = (0^alpha).*sin((gamma*pi).*y);
    c1 = cos(beta*pi).*sin((gamma*pi).*y);
    
    V_computed = FFTSolver(m,f, b0,b1,c0,c1);
    V_exact = (X.^alpha).*cos( (beta*pi.*X)).*sin((gamma*pi).*Y); 
    
    V_diff = abs(V_computed - V_exact);
    [error, error_id] = max(V_diff(:));
    [error_x, error_y] = ind2sub(size(V_diff),error_id);
    fprintf('\n m = %i, alpha = %f, beta = %f, gamma= %f, error= %f\n\n', ...
        m, alpha, beta,gamma,error);
    
    clf;
    hold on;
    subplot(1,3,1);
    mesh(X,Y,f);    
    xlabel('x');
    ylabel('y');
    title('RHS f');
   
    subplot(1,3,2);
    mesh(X,Y,V_computed);
    xlabel('x');
    ylabel('y');
    title('Computed Solution');
    
    subplot(1,3,3);
    mm_plot = mesh(X,Y,V_exact);
    xlabel('x');
    ylabel('y');
    title('Exact Solution');
    
    filename = sprintf('%s',plot_name,'.png');    
    saveas(mm_plot, filename);
end


