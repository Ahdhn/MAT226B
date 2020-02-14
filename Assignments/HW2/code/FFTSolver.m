function V = FFTSolver(m, f, b0, b1, c0, c1)
    h = 1/(m+1);
    factor = 1/(h^2);
    %adjust F
    f(1,:) = f(1,:) + factor*b0;
    f(m,:) = f(m,:) + factor*b1;
    f(:,1) = f(:,1) + factor*c0';
    f(:,m) = f(:,m) + factor*c1';
    
    f = applyFFT(f, h, m);
    f = applyFFT(f', h, m)';
    
    [X,Y] = meshgrid(1:m, 1:m);
    
    L = 2*(1-cos(pi*h.*X) +1 -cos(pi*h.*Y));
    V = ((h^2)*f./L);
    
    V = applyFFT(V, h, m);
    V = applyFFT(V', h, m)';    
end


function out = applyFFT(F,h, m)        
    F_t = [zeros(1,m);F;zeros(m+1,m)];
    w_t = fft(F_t);
    w_h = w_t(2:m+1,:);
    out = -sqrt(2*h)*imag(w_h);
end