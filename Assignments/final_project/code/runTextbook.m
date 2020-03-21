function Zk_textbook = runTextbook(num_data,k, s0, r, c, s_vector)
     global L  U P Q inv_D  E;
    Zk_textbook = zeros(num_data,1);
    [alpha, beta] = textbookAlgo(L, U, P, Q, inv_D, E, r, c, k);
    
    for ss =1:length(s_vector)
       s = s_vector(ss);       
       Zk_textbook(ss) = log10(abs(polyval(flip(alpha), s-s0)/...
           polyval(flip(beta),s-s0)));    
    end    
end 