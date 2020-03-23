function q = transposeMv(L, U, P, Q, inv_D, E, v)  
  
    c = transpose(U)\(transpose(Q)*v);
    d = transpose(L)\c;
    g = transpose(P*inv_D)*d;
    q = transpose(E)*g;
end