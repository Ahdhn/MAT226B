function q = Mv(L, U, P, Q, inv_D, E, v)    
    %TODO use 'vector'????
    f = E*v;    
    c = L\(P*inv_D*f); %is this okay?
    d = U\c;
    q = Q*d;
end