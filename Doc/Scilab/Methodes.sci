exec('Wolfe_Skel.sci');

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//         RESOLUTION D'UN PROBLEME D'OPTIMISATION SANS CONTRAINTES          //
//                                                                           //
//                Fonction d'optimisation                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

function [fopt, xopt, gopt, log_iter, log_F, log_G, log_P] = Optim(Oracle, xini, xpini, alpha0, iter_max, iter_max_alpha, meth)
    // Initialisation des variables d'optimisation
    x = xini;
    xp = xpini;
    W = eye(dim, dim);
    
    // Evolution de la fonctionnelle, de son gradient et du pas
    
    log_iter = [];
    log_F = [];
    log_G = [];
    log_P = [];
    
    [F, G] = Oracle(x, 4);
    D = -G;
    for iter = 1:iter_max
        [F, G] = Oracle(x, 4);
        alpha = alpha0;
        select meth
        case "GRADF" then
            D = -G;
        case "GRADV" then
            D = -G;
            alpha_init_grad = -2/(G'*D);
            alpha = Wolfe(alpha_init_grad, x, D, Oracle, iter_max_alpha);
        case "NEWTN" then
            D = newton(Oracle, x);
            alpha = Wolfe(1, x, D, Oracle, iter_max_alpha);
        case "QNEWT" then
            [D, W] = quasi_newton(Oracle, xp, x, W);
            alpha = Wolfe(1, x, D, Oracle, iter_max_alpha);
        case "GRADC" then
            D = gradient_conjugue(Oracle, xp, x, D);
            alpha_init_grad = -2/(G'*D);
            alpha = Wolfe(alpha_init_grad, x, D, Oracle, iter_max_alpha);
        end
        // Modification de l'estimasaiton du vecteur minimisant la fonctionnelle
        xp = x;
        x = xp + alpha * D;
        // Remplissage des logs
        [Fc, Gc] = Oracle(x, 4);
        log_iter($+1) = iter;
        log_F($+1) = Fc;
        //if(modulo(iter,50)==0) then
        //    disp(Fc)
        //end
    end
    xopt = x;
    [fopt, gopt] = Oracle(xopt, 4);
endfunction
