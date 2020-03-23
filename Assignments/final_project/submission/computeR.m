function r = computeR(L, U, P, Q, inv_D, b)    
    s = L\(-P*inv_D*b);    
    d = U\s;
    r = Q*d;
end