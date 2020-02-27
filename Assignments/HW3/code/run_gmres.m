function [rho, x, iter] = run_gmres(matVecFun, b,x0, tol,  k, M1, M2,...
    displayName, figTitle, fileName)

    n = length(b);
    [x,FLAG,relres,iter,rho] = ...
        gmres(matVecFun,b, k, tol, 1000, M1, M2, x0);    
    if FLAG ~= 0
        warning('\n WARNING: GMRES finished with flag %d', FLAG);        
    end
    rho = rho./rho(1);
    if(k ~= 0)
        fprintf('\n Restart parameter k = %d', k);            
        displayName = strcat(displayName,sprintf(', K =%d', k));
    end
    
    fprintf('\n Total number of GMRES steps: outer= %d, inner= %d, total = %d',...
        iter(1),iter(2), length(rho));
    hold on;    
    
    xylabel = @(fig, t) [xlabel(fig, 'Iteration ');...
        ylabel(fig, 'Relative Residual (Log Scale)');title(t); legend;];
    
    HH = plot( 1:length(rho), log(rho),'DisplayName',displayName, 'LineWidth',2);
    xylabel(gca, figTitle);
    
    hold off;
    fileName = sprintf('%s%s',fileName,'.png');
    saveas(HH, fileName);
end 





