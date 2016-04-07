function q_dual = dual_arg(lambda)
    Z = -(Ad' * lambda + Ar'*pr)
    q_dual = sqrt(abs(Z).*(1 ./ r)).*sign(Z)
endfunction

function L = lagrangien(q, lambda)
    L = q' * (r .* q .* abs(q)) + pr' * (Ar * q) + lambda' * (Ad * q - fd)
endfunction

function F = dual(lambda)
    q_dual = dual_arg(lambda)
    F = lagrangien(q_dual, lambda)
endfunction
    
function G = dual_grad(lambda)
    q = dual_arg(lambda)
    G = Ad * q - fd
endfunction

function H = dual_hessian(lambda)
    Z = -(Ad' * lambda + Ar'*pr)
    diag_interm = diag(1/(2 * sqrt(abs(Z) .* r)))
    H = Ad * diag_interm * Ad'
endfunction
    
function [F,G,ind] = OracleDG(lambda,ind)
    select ind
    case 2 then
        F = dual(lambda);
    case 3 then
        G = dual_grad(lambda);
    case 4 then
        F = dual(lambda);
        G = dual_grad(lambda);
    end
endfunction

function [F,G,H,ind] = OracleDH(lambda,ind)
    select ind
    case 2 then
        F = dual(lambda);
    case 3 then
        G = dual_grad(lambda);
    case 4 then
        F = dual(lambda);
        G = dual_grad(lambda);
    case 5 then
        H = dual_hessian(lambda);
    case 6 then
        G = dual_grad(lambda);
        H = dual_hessian(lambda);
    case 7 then
        F = dual(lambda);
        G = dual_grad(lambda);
        H = dual_hessian(lambda);
    end
endfunction

