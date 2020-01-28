function [L, flag, k] = spcholesky(A, nnzmax, tlimit)
    n = length(A);
    L = tril(A); 
    total_time =0;
   
    for k =1:n
        
        tic; 
        L(k,k) = sqrt(L(k,k));        
        L(k+1:n,k) = (1/L(k,k)).*L(k+1:n,k);
        for j=k+1:n
            if isempty(find(L(j,k), 1)) == 0                
                L(j:n,j) = L(j:n,j) - L(j:n,k).*L(j,k);
            end
        end        
        k_time = toc;           
        
        total_time = total_time + k_time;
        if total_time > tlimit
            flag = 2;
            return;
        end
        
        if nnz(L) > nnzmax
            flag = 1;
            return;
        end
        
    end
    flag =0;
end 