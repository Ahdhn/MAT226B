function mu = computeMoments(L, U, P, Q, inv_D, E, c, r, k)
   
    f = r;    
    mu = zeros(2*k,1);
    mu(1) = transpose(c)*f;
    for j = 2:2*k
        f = Mv(L, U, P, Q, inv_D, E, f);
        mu(j) = transpose(c)*f;
    end    
end