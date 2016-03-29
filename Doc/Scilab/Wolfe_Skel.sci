function [W1_is_OK] = Wolfe1_cond(F,Fd,G,alpha,D,w1)
    W1_is_OK = (Fd <= F + w1*alpha*G'*D);
endfunction

function [W2_is_OK] = Wolfe2_cond(F,G,Gd,alpha,D,w2)
    W2_is_OK = (Gd'*D >= w2*G'*D);
endfunction

function [alphan,ok]=Wolfe(alpha,x,D,Oracle,iter_max)

//////////////////////////////////////////////////////////////
//                                                          //
//   RECHERCHE LINEAIRE SUIVANT LES CONDITIONS DE WOLFE     //
//                                                          //
//                                                          //
//  Arguments en entree                                     //
//  -------------------                                     //
//    alpha  : valeur initiale du pas                       //
//    x      : valeur initiale des variables                //
//    D      : direction de descente                        //
//    Oracle : nom de la fonction Oracle                    //
//                                                          //
//  Arguments en sortie                                     //
//  -------------------                                     //
//    alphan : valeur du pas apres recherche lineaire       //
//    ok     : indicateur de reussite de la recherche       //
//             = 1 : conditions de Wolfe verifiees          //
//             = 2 : indistinguabilite des iteres           //
//                                                          //
//                                                          //
//    omega1 : coefficient pour la 1-ere condition de Wolfe //
//    omega2 : coefficient pour la 2-eme condition de Wolfe //
//                                                          //
//////////////////////////////////////////////////////////////

// -------------------------------------
// Coefficients de la recherche lineaire
// -------------------------------------

   omega1 = 0.1;
   omega2 = 0.9;

   alphamin = 0.0;
   alphamax = %inf;

   ok = 0;
   dltx = 0.00000001;
   iter = 1;

// ---------------------------------
// Algorithme de Fletcher-Lemarechal
// ---------------------------------

   // Appel de l'oracle au point initial
   
   ind = 4;
   [F,G] = Oracle(x, ind);

   // Initialisation de l'algorithme

   alphan = alpha;
   xn = x;

   // Boucle de calcul du pas
   //
   // xn represente le point pour la valeur courante du pas,
   // xp represente le point pour la valeur precedente du pas.

   while (ok == 0 & iter <= iter_max)
      
      xp = xn;
      xn = x + (alphan*D);
      [Fd, Gd] = Oracle(xn, ind);

      // Calcul des conditions de Wolfe

      if (~Wolfe1_cond(F, Fd, G, alphan, D, omega1)) then
          alphamax = alphan;
          alphan = 1/2*(alphamin+alphamax);
      else
          if (~Wolfe2_cond(F, G, Gd, alphan, D, omega2)) then
              if (alphamax == %inf) then
                  alphamin = alphamin * 2;
              else
                  alphamin = (1/2)*(alphamin + alphamax);
              end
          else
              ok = 1;
          end
      end
      

      // Test de la valeur de alphan :
      // - si les deux conditions de Wolfe sont verifiees,
      //   faire ok = 1 : on sort alors de la boucle while
      // - sinon, modifier la valeur de alphan : on reboucle.

      // Test d'indistinguabilite
      if norm(xn-xp) < dltx then
        ok = 2;
      end

   end
endfunction
