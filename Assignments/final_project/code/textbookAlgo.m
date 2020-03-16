function [alpha, beta] = textbookAlgo(L, U, P, Q, inv_D, E, c, b, k)
  mu = computeMoments(L, U, P, Q, inv_D, E, c, b, k);
  
  M1 = zeros(k, k);
  for i = 1:k
      M1 = M1 + diag(repmat(mu(i), [1 k-i]),-i);
  end      
  M2 = toeplitz(mu(k:2*k-1), flip(mu(1:k)));
  
  beta = zeros(k+1,1); 
  beta(1) = 1;
  beta(2:end) = M2\mu(k+1:end);
  alpha = mu(1:k) + M1*beta(2:end);
end