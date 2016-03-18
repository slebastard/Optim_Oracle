function F = fonctio(qc)
    tmp1 = (q0 + B*qc)' * (r .* (q0 + B*qc) .* abs(q0 + B*qc))
    tmp2 = pr'*(Ar*(q0+B*qc))
    F = ((1/3)*tmp1) + tmp2
endfunction
    
function G = fonctio_grad(qc)
    tmp00 = (q0 + B * qc);
    tmp20 = B' * (r .* tmp00 .* abs(tmp00))
    tmp30 = B' * Ar' * pr
    G = tmp20 + tmp30
endfunction

function H = fonctio_hessian(q0, qc)
    v = r.*abs(q0 + B*qc);
    V = cat(2, v, v, v);
    H = 2*(B' * (V.*B));
endfunction
    
function [F,G,ind] = OraclePG(qc,ind)
    select ind
    case 2 then
        F = fonctio(qc)
    case 3 then
        G = fonctio_grad(qc)
    case 4 then
        F = fonctio(qc)
        G = fonctio_grad(qc)
    end
endfunction

function [F,G,H,ind] = OraclePH(qc,ind)
    select ind
    case 2 then
        F = fonctio(qc)
    case 3 then
        G = fonctio_grad(qc)
    case 4 then
        F = fonctio(qc)
        G = fonctio_grad(qc)
    case 5 then
        H = fonctio_hessian(qc)
    case 6 then
        G = fonctio_grad(qc)
        H = fonctio_hessian(qc)
    case 7 then
        F = fonctio(qc)
        G = fonctio_grad(qc)
        H = fonctio_hessian(qc)
    end
endfunction
