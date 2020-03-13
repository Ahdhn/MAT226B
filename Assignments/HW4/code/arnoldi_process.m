function [V, H] = arnoldi_process(A, r, kmax)
    n = length(r);
    beta = norm(r);
    V = zeros(n, kmax+1);
    V(:,1) = r/beta;
    H = zeros(kmax+1, kmax);
    for k = 1:kmax
        q = A*V(:,k);
        for j=1:k
            H(j,k) =V(:,j)'*q;
            q = q - H(j,k)*V(:,j);
        end
        H(k+1,k) = norm(q);
        if H(k+1,k)==0
            break;
        end
        V(:,k+1)=q/H(k+1,k);
    end 
end 