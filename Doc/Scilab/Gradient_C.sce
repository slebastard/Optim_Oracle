function D = gradient_conjugue(Oracle, x)
    ind = 4
    [F ,G ,H] = Oracle(x, ind);
    // H est une matrice 3x3, on va regarder si elle est inversible
    w = det(H*H);
    if(w ~= 0) then
        W = inv(H*H)
        D = -W * G;
    end
endfunction

function B = Polak(Gp,Gn)
    ind = 4
    [F ,G ,H] = Oracle(x, ind);
    // H est une matrice 3x3, on va regarder si elle est inversible
    w = det(H*H);
    if(w ~= 0) then
        W = inv(H*H)
        D = -W * G;
    end
endfunction
