///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//         RESOLUTION D'UN PROBLEME D'OPTIMISATION SANS CONTRAINTES          //
//                                                                           //
//                Methode de quasi-Newton                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// Algorithme BFGS d'évaluation d'une approximation de la matrice hessienne

function [Wn] = bfgs(xp, xn, Gp, Gn, Wp)
    delta_x = xn - xp
    delta_g = Gn - Gp
    if norm(delta_g) == 0 then
        Wn = 0
    else
        size_x = size(delta_x, "r");
        dif_1 = eye(size_x,size_x) - (1/(delta_g' * delta_x)).*(delta_x * delta_g');
        dif_2 = eye(size_x,size_x) - (1/(delta_g' * delta_x)).*(delta_g * delta_x');
        Wn = dif_1 * Wp * dif_2 + (1/(delta_g' * delta_x))*(delta_x * delta_x');
    end
endfunction

// Mise à jour de l'estimation de la hessienne, calcul de la direction de descente

function [D, Wn] = quasi_newton(Oracle, xp, xn, Wp)
    ind = 4;
    [Fn,Gn] = Oracle(xn, ind);
    [Fp,Gp] = Oracle(xp, ind);
    // H est une matrice 3x3, on va regarder si elle est inversible
    Wn = bfgs(xp, xn, Gp, Gn, Wp);
    D = -Wn*Gn;
endfunction
