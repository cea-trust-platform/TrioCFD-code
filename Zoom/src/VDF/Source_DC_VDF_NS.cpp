/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/
//////////////////////////////////////////////////////////////////////////////
//
// File:        Source_DC_VDF_NS.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/VDF
// Version:     /main/16
//
//////////////////////////////////////////////////////////////////////////////

#include <Dirichlet_entree_fluide_leaves.h>
#include <Navier_Stokes_std.h>
#include <Source_DC_VDF_NS.h>
#include <Zone_VDF.h>

Implemente_instanciable(Source_DC_VDF_NS,"Source_DC_VDF_Face",Source_Correction_Deficitaire);


//// printOn
//

Sortie& Source_DC_VDF_NS::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}


//// readOn
//

Entree& Source_DC_VDF_NS::readOn(Entree& s )
{
  return s ;
}


void Source_DC_VDF_NS::associer_domaines(const Zone_dis& zone_dis,
                                      const Zone_Cl_dis& zone_cl_dis)
{
  const Zone_VDF& zvdf = ref_cast(Zone_VDF,zone_dis.valeur());

  int nb_face = zvdf.nb_faces();
  la_correction.resize(nb_face);
  le_residu.resize(nb_face);
  le_residu=0;
}



//ATTENTION !!!!!!!!!!!!!!!!!!!!!!
//ON NE TRAITE QUE LE CAS OU LA CONVECTION EST NEGLIGEABLE  :-(
//!!!!!!!!!!!!!!!!!!!!
//!!!!!!!!!!!!!!!!!!!!
DoubleTab& Source_DC_VDF_NS::calculer_residu(Connectivites_base& connect, LIST(Prolongement)& P, Equation_base& eqG, Nom& type_eq)
{
  Cerr<<"debut de Source_DC_VDF_NS::calculer_residu"<<finl;

  le_residu = 0.;

  /* On resoud  dU/dt + L U + grad P = F         ou d/dt est une derivee partielle                                      */
  /* le "residu" final doit etre :                                                                                      */
  /*  prol_fin(L_gros U_gros + grad(P_gros) - F_gros) - L_fin(prol_fin(U_gros)) - grad_fin(Prol_fin(P_gros))            */
  /* ie prolongement des operateurs grossiers -  operateur des prolongements - gradients du prolongement de la pression */

  int i;


  if(type_eq == "NS")
    {
      Prolongement& P1 = P(0); //prolongement pour la vitesse sur tout le domaine fin
      Prolongement& P2 = P(1); //prolongement pour la pression
      Prolongement& P3 = P(2); //prolongement pour la vitesse sur les bords du domaine fin

      /* Recuperation des equations de bases et des zones                     */
      Navier_Stokes_std& eqF = ref_cast(Navier_Stokes_std, equation());
      Zone_VDF& le_dom = ref_cast(Zone_VDF, eqG.zone_dis().valeur());
      Zone_VDF& le_dom_fine = ref_cast(Zone_VDF, eqF.zone_dis().valeur());

      /* Recuperation des inconnues vitesses fines et grossieres              */
      DoubleTab& presentG = eqG.inconnue().valeurs();
      DoubleTab& presentF = eqF.inconnue().valeurs();
      int tailleF = presentF.dimension(0); // = nombre de faces fines
      int tailleG = presentG.dimension(0); // = nombre de faces grossieres

      /* on recupere la pression grossiere et fine                            */
      Navier_Stokes_std& eqG_typee = ref_cast(Navier_Stokes_std, eqG);
      DoubleTab& pressionG = eqG_typee.pression().valeurs();
      DoubleTab& pressionF = eqF.pression().valeurs();

      const Bord front_fictive; // Pour les prolongements 0 et 1






      /*----------------------------------------------------------------------*/
      /* Partie 1 : Prolongement des operateurs grossiers                     */
      /*----------------------------------------------------------------------*/

      Cerr << "partie 1" << finl;
      /* Calcul de la somme des operateurs grossiers */
      Cerr << "Calcul de la somme des operateurs grossiers"<<finl;
      DoubleTab op_G;
      op_G.resize(tailleG);
      op_G  = 0.;
      for(i=0; i<eqG.nombre_d_operateurs(); i++)
        {
          eqG.operateur(i).ajouter(presentG,op_G) /*<< finl*/ ;
        }
      //     Cerr << "operateur grossier " << op_G << finl;
      //     Cerr << "Calcul de la somme des operateurs grossiers --->OK"<<finl;

      /* ajout du gradient de pression prolonge aux operateur grossiers */
      Cerr<<"ajout du gradient de pression grossier"<<finl;
      /* calcul du gradient de pression grossier */
      eqG_typee.operateur_gradient().ajouter(pressionG,op_G);
      //     Cerr << "gradient de pression grossier" << grad_pG << finl;

      /* ici on a -somme des operateurs => on multiplie par -1 */
      op_G *= -1.;
      //     Cerr<<"ajout du gradient de pression grossier --->OK"<<finl;


      /* Le prolongement fait intervenir des surfaces grossieres que l'on ne veut pas */
      /* (il nous faut des surfaces fines car on prolonge sur la grille fine)         */
      /* => on divise tout par les surfaces grossiere avant de prolonger              */
      const DoubleVect& surface_G = le_dom.face_surfaces();
      for (i=0; i<tailleG; i++)
        op_G(i)=op_G(i)/surface_G(i);


      /* prolongement de la somme des operateurs grossiers */
      DoubleTab op_G_prol;
      op_G_prol.resize(tailleF);
      op_G_prol = 0.;
      P1.prolonger(le_dom, le_dom_fine,front_fictive, connect.connectivites_elemF_elemG(), op_G, op_G_prol, 1);
      //     Cerr << "prolongement des operateurs grossiers" << op_G_prol << finl;

      Cerr<<"l prolongement des operateurs "<< le_residu << finl;



      DoubleTab deriveeG(presentG);
      DoubleTab derivee(le_residu);
      eqG.derivee_en_temps_inco(deriveeG) ;
      P1.prolonger(le_dom, le_dom_fine,front_fictive, connect.connectivites_elemF_elemG(), deriveeG,derivee , 1);


      Cerr << "Prolongement Derivee en temps inco " << derivee << finl;
      le_residu = derivee;
      le_residu *= -1.;
      Cerr << "le _residu = " << le_residu << finl;

      /*----------------------------------------------------------------------*/
      /* Partie 2 : - Operateurs des prolongements                            */
      /*----------------------------------------------------------------------*/

      Cerr<<"partie 2"<<finl;
      /* Prolongement de l'inconnue grossiere present                         */
      Cerr <<"prolongement de l'inconnue grossiere present" << finl;
      DoubleTab presentG_prol;
      presentG_prol.resize(tailleF);
      presentG_prol = 0.;

      DoubleTab op_presentG_prol;
      op_presentG_prol.resize(tailleF);
      op_presentG_prol = 0.;

      P1.prolonger(le_dom, le_dom_fine,front_fictive, connect.connectivites_elemF_elemG(), presentG, presentG_prol, 1);
      Cerr <<"prolongement de l'inconnue grossiere present ---> ok" << finl;

      //      Cerr<<"vitesse apres prolongement : "<< presentG_prol <<finl;

      /* On ne peut pas prolonger directement car les operateurs fins ont besoin de conditions limites venant de la grille grossiere */
      /* => On construit une condition limite de Dirichlet pour le laplacien du prolongement                                          */
      /*    en fonction des vitesses de la grille grossiere                                                                           */

      /*       Schema_Temps_base& schF_temp = eqF.schema_temps();
               double t_fin = schF_temp.temps_courant();
      */

      /* Recuperation de la zone fine de base */
      Zone_Cl_dis_base& zoneCL_baseF = eqF.zone_Cl_dis().valeur();

      /* recuperation des conditions limites existantes */
      Conds_lim& conds_limF = zoneCL_baseF.les_conditions_limites();

      /* Creation des conditions limites pour le laplacien du prolongement */
      /* en modifiant l'objet deja existant                                */


      /* ---> creation d'une condition limite Frontiere_ouverte_vitesse_imposee */
      Entree_fluide_vitesse_imposee&   CL_vitesse_imposeeF = ref_cast(Entree_fluide_vitesse_imposee,conds_limF[0].valeur());


      const Champ_front& champ_frontF = CL_vitesse_imposeeF.champ_front();
      const Frontiere& frontiereF =champ_frontF.frontiere_dis().frontiere();

      /* ---> On remplit ici les valeurs de la vitesse imposee en prenant celles de la grille grossiere prolongee */
      DoubleTab& vitesse_imposeeF = CL_vitesse_imposeeF.champ_front().valeur().valeurs();
      const int nb_faces_bords    = frontiereF.nb_faces();//eqF.zone_Cl_dis().nb_faces_Cl();
      const int nb_faces_bords_tot    = eqF.zone_Cl_dis().nb_faces_Cl();

      Cerr << "nb_faces_bords :"<<nb_faces_bords<< finl;


      /* vitesse_imposeeF.resize(nb_faces_bords,Objet_U::dimension); */
      /* Cerr << "present G = " << presentG << finl; */

      P3.prolonger(le_dom, le_dom_fine,frontiereF, connect.connectivites_faceF_faceG(), presentG, vitesse_imposeeF, Objet_U::dimension);



      for (i=0; i<10; i++)
        {
          //        vitesse_imposeeF(i,0)=1./40.+i/20.;
          vitesse_imposeeF(i,1)=0;
        }
      for (i=10; i<20; i++)
        {
          vitesse_imposeeF(i,0)=1;
          vitesse_imposeeF(i,1)=0;
        }
      for (i=20; i<30; i++)
        {
          vitesse_imposeeF(i,0)=0;
          vitesse_imposeeF(i,1)=0;
        }
      for (i=30; i<40; i++)
        {
          //        vitesse_imposeeF(i,0)=1./40.+(i-60)/20.;
          vitesse_imposeeF(i,1)=0;
        }
      Cerr << "vitesse_imposee F   = " << vitesse_imposeeF << finl;




      //             zoneCL_baseF.mettre_a_jour(t_fin);
      //             Cerr << "vitesse imposee sur les bords : " << vitesse_imposeeF << finl;




      DoubleTab pressionG_prol;
      pressionG_prol.resize(pressionF.size());


      P2.prolonger(le_dom, le_dom_fine,front_fictive, connect.connectivites_elemF_elemG(), pressionG, pressionG_prol, 1);
      DoubleTab presentF_temp(presentF.size());
      DoubleTab pressionF_temp(pressionF.size());
      presentF_temp = presentF;
      pressionF_temp = pressionF;
      presentF = presentG_prol;
      pressionF = pressionG_prol;


      /* for(i=0; i<presentF.size();i++)
         {
         if (le_dom_fine.orientation(i) == 0)
         presentF(i)=-(le_dom_fine.xv(i,1)*le_dom_fine.xv(i,1)-le_dom_fine.xv(i,1))/10;
         else presentF(i)=0;
         }
      */

      zoneCL_baseF.imposer_cond_lim(eqF.inconnue(),eqF.probleme().presentTime());

      {
        DoubleTab tmp(presentG);
        tmp -= presentF;
        Cerr << "||presentF-presentG|| = " << norme_array(tmp) << finl;
        Cerr << "diff : " << tmp << finl;

        tmp = pressionG_prol;
        tmp -= pressionF;
        Cerr << "||pressionF-pressionG_prol|| = " << norme_array(tmp) << finl;
      }


      /* Cerr << "presentF = " << presentF << finl; */


      DoubleTab deriveeF(presentF.size());
      DoubleTab residu_temp(le_residu.size());
      residu_temp = le_residu;
      le_residu = 0;
      eqF.derivee_en_temps_inco(deriveeF);
      le_residu = residu_temp;

      Cerr << "Derivee en temps inco du prolongement" << deriveeF << finl;
      {
        DoubleTab tmp(deriveeF);
        tmp-=deriveeG;
        Cerr << "||deriveeF-deriveeG|| = " << norme_array(tmp) << finl;
      }

      // ajout au residu
      Cerr<<"ajout au residu"<<finl;
      for(i=0; i<tailleF; i++)
        {
          //le_residu(i) = le_residu(i) + op_presentG_prol(i);
          le_residu(i)+=deriveeF(i);
        }

      const DoubleVect& volentre = le_dom_fine.volumes_entrelaces();
      for(i=0; i<tailleF; i++)
        {
          le_residu(i)*=volentre(i);
        }

      for(i=0; i<nb_faces_bords_tot; i++)
        {
          le_residu(i)=0;
        }

      Cerr<<"ajout au residu --->OK"<<finl;

      Cerr << "le residu = " << le_residu << finl;
      Cerr << "norme residu = " << norme_array(le_residu) << finl;


      /* retour a la normale pour le reste */

      vitesse_imposeeF = 0.0;
      //      zoneCL_baseF.mettre_a_jour(t_fin);
      presentF = presentF_temp;
      pressionF = pressionF_temp;
      zoneCL_baseF.imposer_cond_lim(eqF.inconnue(),eqF.probleme().presentTime());
    }
  else
    {
      Cerr<<"Ce type d'equation n'a pas ete traite !!"<<finl;
      exit();
    }

  //  Cerr<<"voici le residu final : "<< le_residu << finl;
  return le_residu;
}

