function zk = zkViaLanczos(Tk, c, r, s, s0) 
  k = size(Tk,1);
  I =speye(k);
  e1 = zeros(k,1);
  e1(1) = 1;
  x = (I - (s-s0).*Tk)\e1;
  zk = (transpose(c)*r)*transpose(e1)*x;
end 