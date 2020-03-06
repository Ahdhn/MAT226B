function [Tk] = hermitian_lanczos(A, r, kmax)
    alpha = zeros(kmax,1);
    beta = zeros(kmax+1,1);
    beta(1) = norm(r);    
    vk = r./beta(1);    
    for k =1:kmax
        q = A*vk;
        if k > 1 
            q = q - beta(k)*vk_1;
        end
        alpha(k) = vk'*q;        
        q = q - alpha(k)*vk;        
        beta(k+1) = norm(q);        
        if beta(k+1) == 0
            break;
        end        
        vk_1 = vk;
        vk = q./beta(k+1);
    end        
    Tk = sparse(diag(alpha(1:k),0)+diag(beta(2:k),-1)+diag(beta(2:k),1));
end