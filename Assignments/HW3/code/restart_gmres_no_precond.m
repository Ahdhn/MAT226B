function [rho, x, iter] = restart_gmres_no_precond(matVecFun, b,x0, tol, k, name)
     n = length(b);
    [x,FLAG,relres,iter,rho] = ...
        gmres(matVecFun,b, k, tol, n,[], [], x0);    
    if FLAG ~= 0
        warning('\n WARNING: Restart GMRES without preconditioning finished with flag %d', FLAG);        
    end
    rho = rho./rho(1);
    
    fprintf('\n Restart parameter k = %d', k);
    fprintf('\n Total number of GMRES steps: outer= %d, inner= %d', iter(1),iter(2));
    hold on;
    disname = sprintf(' K =%d', k);    
    plot( 1:length(rho), log(rho),'DisplayName',strcat(name,disname), 'LineWidth',2);
end 





