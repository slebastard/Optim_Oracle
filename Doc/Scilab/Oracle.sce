ierr = exec("Probleme_R.sce")
ierr = exec("Structures_R.sce")

function F = fonctio(qc)
    tmp1 = (q0 + B*qc) * (r .* (q0 + B*qc) .* abs(q0 + B*qc))
    tmp2 = pr*(Ar*(q0+B*qc))
    F = 1/3*tmp1 + tmp2
endfunction
    
function G = fonctio_grad(qc)
    tmp00 = (q0 + B * qc);
    tmp10 = (2*r .* abs(tmp00)) * B
    tmp11 = tmp10' * tmp00
    tmp20 = B' * (r .* tmp00 .* abs(tmp00))
    tmp30 = B' * Ar' * pr
    G = tmp11 + tmp20 + tmp30
endfunction

function H = fonctio_hessian(qc)
    // A completer
    // Genere la matrice hessienne de la fonction de cout
endfunction
    
function [F,G,H,ind] = OraclePG(qc,ind = 4)
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
    
function [F#, qc#, G#] = Minimise(Oracle, qc0, eps, seuil_q, seuil_F, iter_max)
    qc = qc0
    iter = 0
    var_q = seuil_q + 1
    var_F = seuil_F + 1
    while (iter < iter_max & var_q > seuil_q & var_F > seuil_F)
        [F,grad] = Oracle(qc)
        qc = qc - eps .* grad
        var_q = eps .* grad
        [F_new] = Oracle(qc)
        var_F = F_new - F
        iter = iter + 1
    end
    qc# = qc
    [F#, G#] = Oracle(qc)
endfunction

function [q#, z#, f#, p#] = ComputeVariables(qc#)
    // A remplir
    // Calcul de toutes les variables utiles au probleme a partir de qc#,
    // valeur optimale de qc
endfunction
