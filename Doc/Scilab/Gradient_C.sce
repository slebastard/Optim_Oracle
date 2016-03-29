function [D] = gradient_conjugue(Oracle, xp, xn, dp)
    ind = 4
    [Fp ,Gp] = Oracle(xp, ind);
    [Fn ,Gn] = Oracle(xn, ind);
    b = Polak(Gp,Gn);
    D = -Gn + b * dp;
    // H est une matrice 3x3, on va regarder si elle est inversible

endfunction

function b = Polak(Gp,Gn)
    // H est une matrice 3x3, on va regarder si elle est inversible
    b = (Gn - Gp)' * Gn;
    b = b / (Gp' * Gp);
endfunction
