clc;
clear;
clf;
format long e;

A = [1 0 1 0 0 0 1 1 1;
     0 1 0 1 0 0 0 0 0;
     1 0 1 1 0 0 1 0 0;
     0 1 1 1 1 0 0 0 0;
     0 0 0 1 1 0 1 0 0;
     0 0 0 0 0 1 0 1 1;
     1 0 1 0 1 0 1 1 0;
     1 0 0 0 0 1 1 1 0;
     1 0 0 0 0 1 0 0 1];
 
title = sprintf('graph G(A)');
plot_graph(0, A, title);

%Step 1 (Remove Node 2)
A(2,:) = zeros(1,length(A));
A(:,2) = zeros(length(A),1);
title = sprintf('Step 1 - Remove #2');
plot_graph(1, A, title);

%Step 2 (Remove Node 4)
A(4,:) = zeros(1,length(A));
A(:,4) = zeros(length(A),1);
A(3,5) = 1; A(5,3) = 1;
title = sprintf('Step 2 - Remove #4 & Add Edge (3,5)');
plot_graph(2, A, title);

%Step 3 (Remove Node 5)
A(5,:) = zeros(1,length(A));
A(:,5) = zeros(length(A),1);
title = sprintf('Step 3 - Remove #5');
plot_graph(3, A, title);


%Step 4 (Remove Node 3)
A(3,:) = zeros(1,length(A));
A(:,3) = zeros(length(A),1);
title = sprintf('Step 4 - Remove #3');
plot_graph(4, A, title);


%Step 5 (Remove Node 6)
A(6,:) = zeros(1,length(A));
A(:,6) = zeros(length(A),1);
A(8,9) = 1; A(9,8) = 1;
title = sprintf('Step 5 - Remove #6 & Add Edge (8,9)');
plot_graph(5, A, title);

%Step 6 (Remove Node 7)
A(7,:) = zeros(1,length(A));
A(:,7) = zeros(length(A),1);
title = sprintf('Step 6 - Remove #7');
plot_graph(6, A, title);

%Step 7 (Remove Node 1)
A(1,:) = zeros(1,length(A));
A(:,1) = zeros(length(A),1);
title = sprintf('Step 7 - Remove #1');
plot_graph(7, A, title);

%Step 8 (Remove Node 8)
A(8,:) = zeros(1,length(A));
A(:,8) = zeros(length(A),1);
title = sprintf('Step 8 - Remove #8');
plot_graph(8, A, title);

P = [0 0 0 0 0 0 1 0 0;
     1 0 0 0 0 0 0 0 0;
     0 0 0 1 0 0 0 0 0;
     0 1 0 0 0 0 0 0 0;
     0 0 1 0 0 0 0 0 0;
     0 0 0 0 1 0 0 0 0;
     0 0 0 0 0 1 0 0 0;
     0 0 0 0 0 0 0 1 0;
     0 0 0 0 0 0 0 0 1]
syms x;
Asyms = [x 0 x 0 0 0 x x x;
         0 x 0 x 0 0 0 0 0;
         x 0 x x 0 0 x 0 0;
         0 x x x x 0 0 0 0;
         0 0 0 x x 0 x 0 0;
         0 0 0 0 0 x 0 x x;
         x 0 x 0 x 0 x x 0;
         x 0 0 0 0 x x x 0;
         x 0 0 0 0 x 0 0 x];
PT_times_A_times_P = transpose(P)*Asyms*P;
%Then we do Cholesky on this one PT_times_A_times_P

%or add the fill-in's in A 
Asyms_fillin = [x 0 x 0 0 0 x x x;
                0 x 0 x 0 0 0 0 0;
                x 0 x x x 0 x 0 0;
                0 x x x x 0 0 0 0;
                0 0 x x x 0 x 0 0;
                0 0 0 0 0 x 0 x x;
                x 0 x 0 x 0 x x 0;
                x 0 0 0 0 x x x x;
                x 0 0 0 0 x 0 x x];
%Then compute transpose(P)*Asyms_fillin*P
transpose(P)*Asyms_fillin*P
%and take the lower triangle of it 

function plot_graph(no, mat, plot_title)
    mat_draw = mat - eye(size(mat));
    G = graph(mat_draw,'OmitSelfLoops');
    H= plot(G, 'LineWidth',3, 'Layout', 'circle');
    set(gca,'XColor', 'none','YColor','none');
    title(plot_title);
    
    highlight(H,[1:9],'NodeColor','k', 'MarkerSize', 8)
    
    nl = H.NodeLabel;
    H.NodeLabel = '';
    xd = get(H, 'XData');
    yd = get(H, 'YData');
    text(xd, yd, nl, 'FontSize',15, 'FontWeight','bold', 'HorizontalAlignment','left', 'VerticalAlignment','top')

    filename = sprintf('%d%s',no,'.png');
    saveas(H, filename);
end 