///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  MONITEUR D'ENCHAINEMENT POUR LE CALCUL DE L'EQUILIBRE D'UN RESEAU D'EAU  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

// --------------------------------------
// Dimensionnement de l'espace de travail
// --------------------------------------

   stacksize(10000000);

// ------------------------------------------
// Fonctions fournies dans le cadre du projet
// ------------------------------------------

   // Donnees du problemes

   exec('Probleme_R.sce');
   exec('Structures_R.sce');
   
   // Affichage des resultats

   exec('Visualg.sci');
   
   // Verification  des resultats

   exec('HydrauliqueP.sci');
   exec('HydrauliqueD.sci');
   exec('Verification.sci');

// ------------------------------------------
// Fonctions a ecrire dans le cadre du projet
// ------------------------------------------

   // ---> Charger les fonctions  associees a l'oracle du probleme,
   //      aux algorithmes d'optimisation et de recherche lineaire.
   //
   // Exemple : la fonction "optim" de Scilab
   //
   exec('Oracle.sce');
   exec('Oracle_lagrangien.sce');
   exec('Gradient_F.sci');
   exec('Methodes.sci');
   titrgr = "Fonction optim de Scilab sur le probleme primal";

    exec('QuasiNewton.sce');
    exec('Newton.sce');
    exec('Gradient_C.sce');

// ------------------------------
// Initialisation de l'algorithme
// ------------------------------

   // La dimension (dim) est celle du probleme primal
   is_primal = %T;
   
   dim = n - md
   lambdaIni = 0.1 * rand(dim,1);
   lambdaPIni = 0.1 * rand(dim,1);

// ----------------------------
// Minimisation proprement dite
// ----------------------------
    meth = "GRADF";
    iter_max = 1000;
    iter_max_alpha = 1e3;
    alpha0 = 3e-5;

    tic()
    // Calcul du vecteur optimal (primal ou dual selon contexte)
   [fopt, lambdaOpt, gopt, log_iter, log_F, log_G, log_P] = Optim(OraclePH, lambdaIni, lambdaPIni, alpha0, iter_max, iter_max_alpha, meth);
    disp("Temps de calcul de la méthode")
    disp(toc())
    
// ---------------------------------------------------------
// Visualisation de l'évolution des variables d'optimisation
// ---------------------------------------------------------
    disp("Évolution de la fonctionnelle énergie en fct de itération")
    plot(log_iter, log_F);
    
    //disp("Évolution du gradient en fct de itération")
    //plot(log_iter, log_G);
    
    //disp("Évolution du pas en fct de itération")
    //plot(log_iter, log_P);

    if is_primal then
        xopt = lambdaOpt;
        [q,z,f,p] = HydrauliqueP(xopt, %F);
        Verification(q,z,f,p);
    else
        // -----------------
        // Passage au primal
        // -----------------
        xopt = dual_arg(lambdaOpt);
        // --------------------------
        // Verification des resultats
        // --------------------------
        [q,z,f,p] = HydrauliqueP(xopt, %T);
        Verification(q,z,f,p);
    end


//
