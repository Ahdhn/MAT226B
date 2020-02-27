function [rho, x, iter] = gmres_precond(matVecFun, b,x0, tol, M1, M2, name)
     n = length(b);
    [x,FLAG,relres,iter,rho] = ...
        gmres(matVecFun,b, [], tol, n, M1 , M2, x0);    
    if FLAG ~= 0
        warning('\n WARNING: GMRES with preconditioning finished with flag %d', FLAG);        
    end
    rho = rho./rho(1);
    fprintf('\n Total number of GMRES steps: outer= %d, inner= %d', iter(1),iter(2));
    hold on
    plot( 1:length(rho), log(rho),'DisplayName',name , 'LineWidth',2);
end 




