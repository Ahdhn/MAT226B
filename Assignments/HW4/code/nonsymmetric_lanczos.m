function [Tk] = nonsymmetric_lanczos(A, r, c, kmax)
    alpha = zeros(kmax,1);
    gamma = zeros(kmax+1,1);
    beta = zeros(kmax+1,1);
    
    beta(1) = norm(r);    
    gamma(1) = norm(c);
    
    vk = r./beta(1);    
    wk = c./gamma(1);
    
    V = zeros(kmax, kmax);%
    W = zeros(kmax, kmax);%
    
    for k =1:kmax        
        V(:,k) = vk;%
        W(:,k) = wk;%
        
        delta_k = wk'*vk;      
        if abs(delta_k) == 0
            break;
        end
        
        q = A*vk;
        s = transpose(A)*wk; %A'*wk
        
        if k > 1 
            q = q - (gamma(k)*delta_k/delta_k_1).*vk_1;            
            s = s - (beta(k)*delta_k/delta_k_1).*wk_1;            
        end
        
        alpha(k) = wk'*q/delta_k;   
        
        q = q - alpha(k)*vk;        
        s = s - alpha(k)*wk;        
        
        beta(k+1) = norm(q);        
        gamma(k+1) = norm(s);
        
        if beta(k+1) == 0 || gamma(k+1) == 0
            break;
        end        
        vk_1 = vk; 
        wk_1 = wk;
        delta_k_1 = delta_k;
        
        vk = (1/beta(k+1)).*q;
        wk = (1/gamma(k+1)).*s;
    end        
    eta = alpha(1:k).*( gamma(1:k) - gamma(2:k+1) );
    Tk = diag(alpha(1:k),0)+diag(beta(2:k),-1)+diag(eta(2:k),1);
end