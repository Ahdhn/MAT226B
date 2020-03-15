function [alph, beta] = textbookAlgo(L, U, P, Q, inv_D, E, c, b, k)
  mu = computeMoments(L, U, P, Q, inv_D, E, c, b, k);
  
  M2 = toeplitz( mu(5:8), mu(1:4));
  
end