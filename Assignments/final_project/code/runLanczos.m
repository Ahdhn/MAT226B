function Zk_lanczos = runLanczos(num_data, k, s0, r, c, s_vector, L, U, P, Q, inv_D, E)

    mat_vec = @(x) Mv(L, U, P, Q, inv_D, E,x);
    trans_mat_vec = @(x) transposeMv(L, U, P, Q, inv_D, E, x);
    
    Zk_lanczos = zeros(num_data,1);
    Tk = nonsymmetricLanczos(mat_vec, trans_mat_vec, r, c, k);   
    I =speye(k);
    e1 = zeros(k,1);
    e1(1) = 1;
    
     for ss =1:length(s_vector)
        s = s_vector(ss);       
        Zk_lanczos(ss) = log10(abs(zkViaLanczos(Tk, I, e1, c,r, s, s0)));           
    end
        
end