exec('Wolfe_Skel.sci');

function [fopt, xopt, gopt,log_iter,log_F] = Optim(Oracle, xini, alpha0, iter_max, iter_max_alpha, meth)
    x = xini;
    xp = xini;
    W = ones(3,3);
    log_iter = [];
    log_F = [];
    [F, G] = Oracle(x, 4);
    D = -G;
    for iter = 1:iter_max
        [F, G] = Oracle(x, 4);
        alpha = alpha0;
        
        select meth
        case "GRADF" then
            D = -G;
        case "GRADV" then
            alpha = Wolfe(alpha, x, D, Oracle, iter_max_alpha);
            D = -G;
        case "NEWTN" then
            alpha = Wolfe(alpha, x, D, Oracle, iter_max_alpha);
            D = newton(Oracle, x);
        case "QNEWT" then
            alpha = Wolfe(alpha, x, D, Oracle, iter_max_alpha);
            [D, W] = quasi_newton(Oracle, x);
        case "GRADC" then
            alpha = Wolfe(alpha, x, D, Oracle, iter_max_alpha);
            D = gradient_conjugue(Oracle, xp, x, D);
        end
        xp = x;
        x = xp + alpha * D;
        [Fc, Gc] = Oracle(x, 4);
        log_iter($+1) = iter;
        log_F($+1) = Fc;
    end
    xopt = x;
    [fopt, gopt] = Oracle(xopt, 4);
endfunction
