function [D] = gradient_conjugue(Oracle, xp, xn, dp)
    ind = 4
    [Fp ,Gp] = Oracle(xp, ind);
    [Fn ,Gn] = Oracle(xn, ind);
    B = Polak(Gp,Gn);
    D = -Gn + B * dp;
    // H est une matrice 3x3, on va regarder si elle est inversible

endfunction

function B = Polak(Gp,Gn)
    // H est une matrice 3x3, on va regarder si elle est inversible
    B = (Gn - Gp)' * Gn;
    B = B / (Gp' * Gp);
endfunction
