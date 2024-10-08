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
// File:        Loi_Paroi_Nu_Impose_VDF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Lois_Paroi/Scal
//
//////////////////////////////////////////////////////////////////////////////

#include <Loi_Paroi_Nu_Impose_VDF.h>
#include <Champ_Uniforme.h>
#include <Champ_Uniforme_Morceaux.h>
#include <Champ_Fonc_Tabule.h>
#include <Domaine_Cl_VDF.h>
#include <Dirichlet_paroi_fixe.h>
#include <Dirichlet_paroi_defilante.h>
#include <Neumann_paroi.h>
#include <Debog.h>
#include <Fluide_base.h>
#include <Convection_Diffusion_std.h>
#include <Probleme_base.h>
#include <EcrFicPartage.h>
#include <Modele_turbulence_scal_base.h>

Implemente_instanciable(Loi_Paroi_Nu_Impose_VDF,"Loi_Paroi_Nu_Impose_VDF",Paroi_scal_hyd_base_VDF);


//     printOn()
/////

Sortie& Loi_Paroi_Nu_Impose_VDF::printOn(Sortie& s) const
{
  return s << que_suis_je() << " " << le_nom();
}

//// readOn
//

Entree& Loi_Paroi_Nu_Impose_VDF::readOn(Entree& s)
{

  // Definition des mots-cles
  Motcles les_mots(2);
  les_mots[0] = "nusselt";
  les_mots[1] = "diam_hydr";
  int nusselt_ok=-1;
  // Lecture et interpretation
  Motcle motlu, accolade_fermee="}", accolade_ouverte="{";
  s >> motlu;
  assert(motlu==accolade_ouverte);
  s >> motlu;
  if(motlu == accolade_ouverte)
    {
      // on passe l'accolade ouvrante
      s >> motlu;
    }
  while (motlu != accolade_fermee)
    {
      int rang=les_mots.search(motlu);
      switch(rang)
        {
        case 0 :   // nusselt expression
          {
            nusselt_ok=0;
            Nom tmp;
            s >> tmp;
            Cerr << "Lecture et interpretation de la fonction " << tmp << " ... ";
            nusselt.setNbVar(5+dimension);
            nusselt.setString(tmp);
            nusselt.addVar("Re");
            nusselt.addVar("Pr");
            nusselt.addVar("t");
            nusselt.addVar("x");
            if (dimension>1)
              nusselt.addVar("y");
            if (dimension>2)
              nusselt.addVar("z");
            nusselt.addVar("Dh");
            nusselt.addVar("Lambda");
            nusselt.parseString();
            break;
          }
        case 1: // diam_hydr
          s >> diam_hydr;
          break;
        default : // non compris
          Cerr << "Mot cle \"" << motlu << "\" non compris lors de la lecture d'un "
               << que_suis_je() << finl;
          exit();
        }
      s >> motlu;
    }

  // Verification de la coherence
  if (nusselt_ok==-1)
    {
      Cerr << "Il faut definir l'expression nusselt(Re,Pr,x,y,z,Dh)" << finl;
      exit();
    }

  if (diam_hydr->nb_comp()!=1)
    {
      Cerr << "Il faut definir le champ diam_hydr a une composante" << finl;
      exit();
    }
  return s;
}


/////////////////////////////////////////////////////////////////////////
//
//    Implementation des fonctions de la classe Loi_Paroi_Nu_Impose_VDF
//
////////////////////////////////////////////////////////////////////////

int  Loi_Paroi_Nu_Impose_VDF::calculer_scal(Champ_Fonc_base& diffusivite_turb)
{
  const Domaine_VDF& domaine_VDF = le_dom_VDF.valeur();
  const IntTab& face_voisins = domaine_VDF.face_voisins();
  const Equation_base& eqn_hydr = mon_modele_turb_scal->equation().probleme().equation(0);
  const DoubleTab& vitesse = eqn_hydr.inconnue().valeurs();
  const Fluide_base& le_fluide = ref_cast(Fluide_base,eqn_hydr.milieu());
  const Champ_Don_base& ch_visco_cin = le_fluide.viscosite_cinematique();
  const IntVect& orientation = domaine_VDF.orientation();
  const DoubleTab& xv=domaine_VDF.xv() ;                   // centres de gravite des faces
  const IntTab& elem_faces = domaine_VDF.elem_faces();

  DoubleVect pos(dimension);

  const DoubleTab& tab_visco = ch_visco_cin.valeurs();
  int l_unif;
  double visco=-1;
  if (sub_type(Champ_Uniforme,ch_visco_cin))
    {
      l_unif = 1;
      visco = std::max(tab_visco(0,0),DMINFLOAT);
    }
  else
    l_unif = 0;

  if ((!l_unif) && (local_min_vect(tab_visco)<DMINFLOAT))
    //   on ne doit pas changer tab_visco ici !
    {
      Cerr<<" visco <=0 ?"<<finl;
      exit();
    }
  //tab_visco+=DMINFLOAT;

  bool dh_constant=sub_type(Champ_Uniforme,diam_hydr.valeur())?true:false;
  double dh_valeur=diam_hydr->valeurs()(0,0);



  int ndeb,nfin;
  int elem;
  double d_visco;

  //  const Convection_Diffusion_std& eqn = mon_modele_turb_scal->equation();
  const Champ_Don_base& alpha = le_fluide.diffusivite();

  // Boucle sur les bords:
  for (int n_bord=0; n_bord<domaine_VDF.nb_front_Cl(); n_bord++)
    {

      // Pour chaque condition limite on regarde son type
      // On applique les lois de paroi thermiques uniquement
      // aux voisinages des parois ou l'on impose la temperature

      const Cond_lim& la_cl = le_dom_Cl_VDF->les_conditions_limites(n_bord);
      if ( (sub_type(Dirichlet_paroi_fixe,la_cl.valeur()))
           || (sub_type(Dirichlet_paroi_defilante,la_cl.valeur())) )
        {

          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          ndeb = le_bord.num_premiere_face();
          nfin = ndeb + le_bord.nb_faces();

          //find the associated boundary
          int boundary_index=-1;
          if (domaine_VDF.front_VF(n_bord).le_nom() == le_bord.le_nom())
            boundary_index=n_bord;
          assert(boundary_index >= 0);

          for (int num_face=ndeb; num_face<nfin; num_face++)
            {
              elem = face_voisins(num_face,0);
              if (elem == -1)
                elem = face_voisins(num_face,1);
              if (l_unif)
                d_visco = visco;
              else
                d_visco = tab_visco[elem];

              double d_alpha=0.;
              if (sub_type(Champ_Uniforme,alpha))
                d_alpha = alpha.valeurs()(0,0);
              else
                {
                  if (alpha.nb_comp()==1)
                    d_alpha = alpha.valeurs()(elem);
                  else
                    d_alpha = alpha.valeurs()(elem,0);
                }
              double Pr = d_visco/d_alpha;

              // On calcule la vitesse debitante comme etant la moyenne des vitesses presentes sur les faces orthogonales a la face courante du bord.
              int ori = orientation(num_face);
              int face1 = elem_faces(elem,(ori+1));
              int face2 = elem_faces(elem,(ori+1+dimension)%(2*dimension));
              double Ud = 0.5*(vitesse(face1)+vitesse(face2));
              Ud *= Ud;
              if (dimension==3)
                {
                  face1 = elem_faces(elem,(ori+2));
                  face2 = elem_faces(elem,(ori+2+dimension)%(2*dimension));
                  double tmp =0.5*(vitesse(face1)+vitesse(face2));
                  Ud += tmp*tmp;
                }

              Ud=sqrt(Ud);

              // Calcul de la position
              for (int i=0; i<dimension; i++)
                pos[i]=xv(num_face,i);

              // Calcul du diametre hydraulique
              if (!dh_constant)
                {
                  dh_valeur=diam_hydr->valeur_a_compo(pos,0);
                }

              double Re=Ud*dh_valeur/d_visco;

              nusselt.setVar("Re",Re);
              nusselt.setVar("Pr",Pr);
              nusselt.setVar("Dh",dh_valeur);
              //nusselt->setVar("Lambda",d_alpha);
              //        nusselt->setVar("t",Pr);
              for (int i=0; i<dimension; i++)
                {
                  nusselt.setVar(3+i,pos[i]);
                }
              double Nu = nusselt.eval();

              // L'expression de d_equiv ne tient pas compte de alpha_t comme en VEF
              // Cela dit, c'est normale car c'est lors du calcul du flux que la
              // turbulence est prise en compte.
              // Ne pas uniformiser l'ecriture avec le VEF, car on tombe sur des problemes
              // au niveau des parois contacts entre plusieurs problemes (alpha_t pas recuperable !).
              int global_face=num_face;
              int local_face=domaine_VDF.front_VF(boundary_index).num_local_face(global_face);
              equivalent_distance_[boundary_index](local_face)=dh_valeur/Nu;
              //Cout << "pos = " << pos[0] << " " << pos[1] << finl;
              //Cout << "Nu = " << Nu << finl;
              //Cout << "diam = " << dh_valeur << finl;
              //Cout << "tab_d_equiv_ = " << tab_d_equiv_[num_face] << finl;
            }
        }
    }
  return 1;
}



void Loi_Paroi_Nu_Impose_VDF::imprimer_nusselt(Sortie&) const
{
  const Domaine_VDF& domaine_VDF = le_dom_VDF.valeur();
  const IntTab& face_voisins = domaine_VDF.face_voisins();
  int ndeb,nfin,elem;
  const Convection_Diffusion_std& eqn = mon_modele_turb_scal->equation();
  const Equation_base& eqn_hydr = eqn.probleme().equation(0);
  const Fluide_base& le_fluide = ref_cast(Fluide_base, eqn_hydr.milieu());
  const Champ_Don_base& conductivite = le_fluide.conductivite();
  const DoubleTab& temperature = eqn.probleme().equation(1).inconnue().valeurs();

  EcrFicPartage Nusselt;
  ouvrir_fichier_partage(Nusselt,"Nusselt");

  bool dh_constant=sub_type(Champ_Uniforme,diam_hydr.valeur())?true:false;
  double dh_valeur=diam_hydr->valeurs()(0,0);
  DoubleVect pos(dimension);

  for (int n_bord=0; n_bord<domaine_VDF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = le_dom_Cl_VDF->les_conditions_limites(n_bord);

      if ( (sub_type(Dirichlet_paroi_fixe,la_cl.valeur())) ||
           (sub_type(Dirichlet_paroi_defilante,la_cl.valeur()) ))
        {
          const Domaine_Cl_VDF& domaine_Cl_VDF_th = ref_cast(Domaine_Cl_VDF, eqn.probleme().equation(1).domaine_Cl_dis());
          const Cond_lim& la_cl_th = domaine_Cl_VDF_th.les_conditions_limites(n_bord);
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());

          //find the associated boundary
          int boundary_index=-1;
          if (domaine_VDF.front_VF(n_bord).le_nom() == le_bord.le_nom())
            boundary_index=n_bord;
          assert(boundary_index >= 0);

          if ( (sub_type(Neumann_paroi,la_cl_th.valeur())))
            {
              const Neumann_paroi& la_cl_neum = ref_cast(Neumann_paroi,la_cl_th.valeur());

              if(je_suis_maitre())
                {
                  Nusselt << finl;
                  Nusselt << "Bord " << le_bord.le_nom() << finl;
                  if (dimension == 2)
                    {
                      Nusselt << "----------------------------------------------------------------------------------------------------------------------------------------------------------------" << finl;
                      Nusselt << "\tFace a\t\t\t\t|" << finl;
                      Nusselt << "----------------------------------------------------------------------------------------------------------------------------------------------------------------" << finl;
                      Nusselt << "X\t\t| Y\t\t\t| dist. carac. (m)\t| Nusselt (local)\t| h (Conv. W/m2/K)\t| Tf cote paroi (K)\t| Tparoi equiv.(K)" << finl;
                      Nusselt << "----------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------" << finl;
                    }
                  if (dimension == 3)
                    {
                      Nusselt << "----------------------------------------------------------------------------------------------------------------------------------------------------------------" << finl;
                      Nusselt << "\tFace a\t\t\t\t\t\t\t|" << finl;
                      Nusselt << "----------------------------------------------------------------------------------------------------------------------------------------------------------------" << finl;
                      Nusselt << "X\t\t| Y\t\t\t| Z\t\t\t| dist. carac. (m)\t| Nusselt (local)\t| h (Conv. W/m2/K)\t| Tf cote paroi (K)\t| Tparoi equiv.(K)" << finl;
                      Nusselt <<        "----------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------" << finl;
                    }
                }
              ndeb = le_bord.num_premiere_face();
              nfin = ndeb + le_bord.nb_faces();
              for (int num_face=ndeb; num_face<nfin; num_face++)
                {
                  // Calcul de la position
                  for (int i=0; i<dimension; i++)
                    pos[i]=domaine_VDF.xv(num_face,i);

                  // Calcul du diametre hydraulique
                  if (!dh_constant)
                    dh_valeur=diam_hydr->valeur_a_compo(pos,0);

                  double x=domaine_VDF.xv(num_face,0);
                  double y=domaine_VDF.xv(num_face,1);
                  double lambda;

                  elem = face_voisins(num_face,0);
                  if (elem == -1)
                    elem = face_voisins(num_face,1);
                  if (sub_type(Champ_Uniforme,conductivite))
                    lambda = conductivite.valeurs()(0,0);
                  else
                    {
                      if (conductivite.nb_comp()==1)
                        lambda = conductivite.valeurs()(elem);
                      else
                        lambda = conductivite.valeurs()(elem,0);
                    }

                  if (dimension == 2)
                    Nusselt << x << "\t| " << y;
                  if (dimension == 3)
                    {
                      double z=domaine_VDF.xv(num_face,2);
                      Nusselt << x << "\t| " << y << "\t| " << z;
                    }

                  int global_face=num_face;
                  int local_face=domaine_VDF.front_VF(boundary_index).num_local_face(global_face);

                  double flux = la_cl_neum.flux_impose(num_face-ndeb);
                  double tparoi = temperature(elem)+flux/lambda*equivalent_distance_[boundary_index](local_face);

                  Nusselt << "\t| " << equivalent_distance_[boundary_index](local_face) << "\t| " << dh_valeur/equivalent_distance_[boundary_index](local_face) << "\t| "
                          << lambda/equivalent_distance_[boundary_index](local_face) << "\t| " << temperature(elem) << "\t|" << tparoi << finl;
                }
            }
          // fin de condition limite flux impose
          else
            {
              if(je_suis_maitre())
                {
                  Nusselt << finl;
                  Nusselt << "Bord " << le_bord.le_nom() << finl;
                  if (dimension == 2)
                    {
                      Nusselt << "----------------------------------------------------------------------------------------------------------------------------------------" << finl;
                      Nusselt << "\tFace a\t\t\t\t|" << finl;
                      Nusselt << "----------------------------------------------------------------------------------------------------------------------------------------" << finl;
                      Nusselt << "X\t\t| Y\t\t\t| dist. carac. (m)\t| Nusselt (local)\t| h (Conv. W/m2/K)\t| Tf cote paroi (K)" << finl;
                      Nusselt << "----------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------" << finl;
                    }
                  if (dimension == 3)
                    {
                      Nusselt << "----------------------------------------------------------------------------------------------------------------------------------------" << finl;
                      Nusselt << "\tFace a\t\t\t\t\t\t\t|" << finl;
                      Nusselt << "----------------------------------------------------------------------------------------------------------------------------------------" << finl;
                      Nusselt << "X\t\t| Y\t\t\t| Z\t\t\t| dist. carac. (m)\t| Nusselt (local)\t| h (Conv. W/m2/K)\t| Tf cote paroi (K)" << finl;
                      Nusselt <<        "----------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------" << finl;
                    }
                }
              ndeb = le_bord.num_premiere_face();
              nfin = ndeb + le_bord.nb_faces();
              for (int num_face=ndeb; num_face<nfin; num_face++)
                {
                  // Calcul de la position
                  for (int i=0; i<dimension; i++)
                    pos[i]=domaine_VDF.xv(num_face,i);

                  // Calcul du diametre hydraulique
                  if (!dh_constant)
                    dh_valeur=diam_hydr->valeur_a_compo(pos,0);

                  double x=domaine_VDF.xv(num_face,0);
                  double y=domaine_VDF.xv(num_face,1);
                  double lambda;
                  elem = face_voisins(num_face,0);
                  if (elem == -1)
                    elem = face_voisins(num_face,1);
                  if (sub_type(Champ_Uniforme,conductivite))
                    lambda = conductivite.valeurs()(0,0);
                  else
                    {
                      if (conductivite.nb_comp()==1)
                        lambda = conductivite.valeurs()(elem);
                      else
                        lambda = conductivite.valeurs()(elem,0);
                    }
                  if (dimension == 2)
                    Nusselt << x << "\t| " << y;
                  if (dimension == 3)
                    {
                      double z=domaine_VDF.xv(num_face,2);
                      Nusselt << x << "\t| " << y << "\t| " << z;
                    }
                  int global_face=num_face;
                  int local_face=domaine_VDF.front_VF(boundary_index).num_local_face(global_face);
                  Nusselt << "\t| " << equivalent_distance_[boundary_index](local_face) << "\t| "
                          << dh_valeur/equivalent_distance_[boundary_index](local_face) << "\t| "
                          << lambda/equivalent_distance_[boundary_index](local_face) << "\t| " << temperature(elem) << finl;
                }
            }
          Nusselt.syncfile();
        }
    }
  if(je_suis_maitre())
    Nusselt << finl << finl;
  Nusselt.syncfile();
}

