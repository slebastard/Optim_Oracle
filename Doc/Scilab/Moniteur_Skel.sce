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
   exec('Gradient_F.sci');
   exec('Wolfe_Skel.sci');
   titrgr = "Fonction optim de Scilab sur le probleme primal";

    exec('QuasiNewton.sce');
    exec('Newton.sce');
    exec('Gradient_C.sce');

// ------------------------------
// Initialisation de l'algorithme
// ------------------------------

   // La dimension (n-md) est celle du probleme primal

   xini = 0.1 * rand(n-md,1);

// ----------------------------
// Minimisation proprement dite
// ----------------------------
    meth = "GRADF";
    iter_max = 10000;
    iter_max_alpha = 10000;
    alpha0 = 0.1;
   // Exemple : la fonction "optim" de Scilab
   //
   [fopt,xopt,gopt,log_iter,log_F] = OptimTest(OraclePH, xini, alpha0, iter_max, iter_max_alpha, meth);
    plot(log_iter, log_F);
   // -----> A completer...

// --------------------------
// Verification des resultats
// --------------------------

   [q,z,f,p] = HydrauliqueP(xopt);

   Verification(q,z,f,p);

//
