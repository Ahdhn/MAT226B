function mu = zkViaLanczos(matVecFunc, matVecTransFunc, r, c, k, s, s0)
  Tk = nonsymmetricLanczos(matVecFunc, matVecTransFunc, r, c, k);
  n = length(r);
  I =speye(n);
  e1 = zeros(n,1);
  e1(1) = 1;
  x = (I - (s-s0).*Tk)\e1;
  %mu = (transpose(c)*r)*transpose(e1)*inv(I - (s-s0).*Tk)*e1;
  mu = (transpose(c)*r)*transpose(e1)*x;
end