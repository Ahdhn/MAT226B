function zk = zkViaLanczos(Tk,I, e1, c, r, s, s0)    
  x = (I - (s-s0).*Tk)\e1;
  zk = (transpose(c)*r)*transpose(e1)*x;
end 

