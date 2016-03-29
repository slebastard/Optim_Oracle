function [Wn] = bfgs(xp, xn, Gp, Gn, Wp)
    delta_x = xn - xp
    delta_g = Gn - Gp
    size_x = size(delta_x, "r");
    dif_1 = ones(size_x,size_x) - (1/(delta_g' * delta_x)).*(delta_x * delta_g');
    dif_2 = ones(size_x,size_x) - (1/(delta_g' * delta_x)).*(delta_g * delta_x');
    Wn = dif_1 * Wp * dif_2 + (1/(delta_g' * delta_x))*(delta_x * delta_x');
endfunction

function [D, Wn] = quasi_newton(Oracle, xp, xn, Wp)
    ind = 4;
    [Fn,Gn] = Oracle(xn, ind);
    [Fp,Gp] = Oracle(xp, ind);
    // H est une matrice 3x3, on va regarder si elle est inversible
    Wn = bfgs(xp, xn, Gp, Gn, Wp);
    D = -Wn*Gn;
endfunction
