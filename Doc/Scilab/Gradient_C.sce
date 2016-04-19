///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                  RESEAUX DE DISTRIBUTION D'EAUX                           //
//                                                                           //
//                    Methode du gradient conjugue                           //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

function [D] = gradient_conjugue(Oracle, xp, xn, dp)
    ind = 4
    [Fp ,Gp] = Oracle(xp, ind);
    [Fn ,Gn] = Oracle(xn, ind);
    b = Polak(Gp,Gn);
    D = -Gn + b * dp;
endfunction

// Retourne le coefficient d'importance de la direction du pas précédent
function b = Polak(Gp,Gn)
    b = (Gn - Gp)' * Gn;
    b = b / (Gp' * Gp);
endfunction
