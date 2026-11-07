function theta = fitCopula(uv, id)
switch id
    case 1
        rho = copulafit('Gaussian', uv); 
        theta = rho;
    case 2
        theta = copulafit('Clayton', uv);  
    case 3
        theta = copulafit('Frank', uv);
    case 4
        theta = copulafit('Gumbel', uv);   
    case 5
        [R, nu] = copulafit('t', uv);
        rho_val = R(1,2);        
        theta = [rho_val, nu];     
    case 6
        uv_rot = [1-uv(:,1), uv(:,2)];
        theta = copulafit('Gumbel', uv_rot);
    case 7
        uv_rot = [uv(:,1), 1-uv(:,2)];
        theta = copulafit('Gumbel', uv_rot);
    otherwise
        error('Unknown copula id %d', id);
end
end