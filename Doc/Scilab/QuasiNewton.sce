function [W] = bfgs(xp, xn, F, Gp, Gn, W)
    delta_x = xn - xp
    delta_g = Gn - Gp
    dif_1 = eyes(3) - (1/(delta_g' * delta_x))*(delta_x * delta_g')
    dif_2 = eyes(3) - (1/(delta_g' * delta_x))*(delta_g * delta_x')
    W = dif_1 * W * dif_2 + (1/(delta_g' * delta_x))*(delta_x * delta_x')
endfunction

function [D, W] = quasi_newton(Oracle, xp, xn, W)
    ind = 4
    [Fn,Gn] = Oracle(xn, ind);
    [Fp,Gp] = Oracle(xp, ind);
    // H est une matrice 3x3, on va regarder si elle est inversible
    W = bgfs(xp, xn, Fn, Gp, Gn, W);
    D = -W*Gn;
endfunction
