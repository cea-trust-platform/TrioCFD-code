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
// File:        Loi_Paroi_Nu_Impose_VEF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Lois_Paroi/Scal
//
//////////////////////////////////////////////////////////////////////////////

#include <Loi_Paroi_Nu_Impose_VEF.h>
#include <Champ_Uniforme.h>
#include <Domaine_Cl_VEF.h>
#include <Dirichlet_paroi_fixe.h>
#include <Dirichlet_paroi_defilante.h>
#include <Symetrie.h>
#include <Neumann_paroi.h>
#include <Fluide_base.h>
#include <Fluide_Quasi_Compressible.h>
#include <Convection_Diffusion_std.h>
#include <Modele_turbulence_scal_base.h>
#include <Probleme_base.h>
#include <EcrFicPartage.h>
#include <SFichier.h>
#include <Modifier_pour_fluide_dilatable.h>

Implemente_instanciable(Loi_Paroi_Nu_Impose_VEF,"Loi_Paroi_Nu_Impose_VEF",Paroi_scal_hyd_base_VEF);


//     printOn()
/////

Sortie& Loi_Paroi_Nu_Impose_VEF::printOn(Sortie& s) const
{
  return s << que_suis_je() << " " << le_nom();
}

//// readOn
//

Entree& Loi_Paroi_Nu_Impose_VEF::readOn(Entree& s)
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
            nusselt.setNbVar(3+dimension);
            nusselt.setString(tmp);
            nusselt.addVar("Re");
            nusselt.addVar("Pr");
            nusselt.addVar("t");
            nusselt.addVar("x");
            if (dimension>1)
              nusselt.addVar("y");
            if (dimension>2)
              nusselt.addVar("z");
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
      Cerr << "Il faut definir l'expression nusselt(Re,Pr)" << finl;
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
//    Implementation des fonctions de la classe Loi_Paroi_Nu_Impose_VEF
//
////////////////////////////////////////////////////////////////////////

int  Loi_Paroi_Nu_Impose_VEF::calculer_scal(Champ_Fonc_base& diffusivite_turb)
{
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  const IntTab& face_voisins = domaine_VEF.face_voisins();
  const Equation_base& eqn_hydr = mon_modele_turb_scal->equation().probleme().equation(0);
  const DoubleTab& vitesse = eqn_hydr.inconnue()->valeurs();
  const Fluide_base& le_fluide = ref_cast(Fluide_base,eqn_hydr.milieu());
  const Champ_Don& ch_visco_cin = le_fluide.viscosite_cinematique();

  const DoubleTab& xv=domaine_VEF.xv() ;                   // centres de gravite des faces
  const IntTab& elem_faces = domaine_VEF.elem_faces();

  ArrOfDouble v(dimension);
  static DoubleVect pos(dimension);

  const DoubleTab& tab_visco = ch_visco_cin->valeurs();
  int l_unif;
  double visco=-1;
  if (sub_type(Champ_Uniforme,ch_visco_cin.valeur()))
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

  int elem;
  double d_visco;
  const Champ_Don& alpha = le_fluide.diffusivite();
  DoubleTab alpha_t = diffusivite_turb.valeurs();

  // Boucle sur les bords:
  for (int n_bord=0; n_bord<domaine_VEF.nb_front_Cl(); n_bord++)
    {

      // Pour chaque condition limite on regarde son type
      // On applique les lois de paroi thermiques uniquement
      // aux voisinages des parois ou l'on impose la temperature

      const Cond_lim& la_cl = le_dom_Cl_VEF->les_conditions_limites(n_bord);
      if ( (sub_type(Dirichlet_paroi_fixe,la_cl.valeur()))
           || (sub_type(Dirichlet_paroi_defilante,la_cl.valeur()))
           || (sub_type(Symetrie,la_cl.valeur())) )
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int size=le_bord.nb_faces();
          for (int ind_face=0; ind_face<size; ind_face++)
            {
              int num_face = le_bord.num_face(ind_face);
              elem = face_voisins(num_face,0);
              if (elem == -1)
                elem = face_voisins(num_face,1);
              if (l_unif)
                d_visco = visco;
              else
                d_visco = tab_visco[elem];

              double d_alpha=0.;
              if (sub_type(Champ_Uniforme,alpha.valeur()))
                d_alpha = alpha->valeurs()(0,0);
              else
                {
                  if (alpha->nb_comp()==1)
                    d_alpha = alpha->valeurs()(elem);
                  else
                    d_alpha = alpha->valeurs()(elem,0);
                }
              double Pr = d_visco/d_alpha;

              // Calcul la moyenne de la vitesse aux faces (autre que bord) dans l'elem, et en calcule le module
              v=0.;
              double Ud = 0.;
              int Compte_face = 0;
              for (int nfa=0; nfa<domaine_VEF.domaine().nb_faces_elem(); nfa++)
                if (domaine_VEF.premiere_face_int()<=elem_faces(elem,nfa))
                  {
                    int face = elem_faces(elem,nfa);
                    for (int dim=0; dim<dimension; dim++)
                      v[dim] += vitesse(face,dim);
                    Compte_face+=1 ;
                  }
              for (int dim=0; dim<dimension; dim++)
                {
                  v[dim] /= Compte_face ;
                  Ud += v[dim]*v[dim];
                }
              Ud=sqrt(Ud) ;

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
              for (int i=0; i<dimension; i++)
                {
                  nusselt.setVar(3+i,pos[i]);
                }
              double Nu = nusselt.eval();
              equivalent_distance_[n_bord](ind_face) = (d_alpha+alpha_t(elem))/d_alpha*dh_valeur/(Nu+DMINFLOAT);
            }
          equivalent_distance_[n_bord].echange_espace_virtuel();
        }
    }
  return 1;
}


void Loi_Paroi_Nu_Impose_VEF::imprimer_nusselt(Sortie&) const
{
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  const IntTab& face_voisins = domaine_VEF.face_voisins();
  int ndeb,nfin,elem;
  const Convection_Diffusion_std& eqn = mon_modele_turb_scal->equation();
  const Equation_base& eqn_hydr = eqn.probleme().equation(0);
  const Fluide_base& le_fluide = ref_cast(Fluide_base, eqn_hydr.milieu());
  const Champ_Don& conductivite = le_fluide.conductivite();
  const DoubleTab& temperature = eqn.probleme().equation(1).inconnue()->valeurs();

  const DoubleTab& conductivite_turbulente =  mon_modele_turb_scal->conductivite_turbulente()->valeurs();

  EcrFicPartage Nusselt;
  ouvrir_fichier_partage(Nusselt,"Nusselt");

  bool dh_constant=sub_type(Champ_Uniforme,diam_hydr.valeur())?true:false;
  double dh_valeur=diam_hydr->valeurs()(0,0);
  DoubleVect pos(dimension);

  for (int n_bord=0; n_bord<domaine_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = le_dom_Cl_VEF->les_conditions_limites(n_bord);
      if ( (sub_type(Dirichlet_paroi_fixe,la_cl.valeur()))
           || (sub_type(Dirichlet_paroi_defilante,la_cl.valeur()))
           || (sub_type(Symetrie,la_cl.valeur())) )
        {
          const Domaine_Cl_VEF& domaine_Cl_VEF_th = ref_cast(Domaine_Cl_VEF, eqn.probleme().equation(1).domaine_Cl_dis());
          const Cond_lim& la_cl_th = domaine_Cl_VEF_th.les_conditions_limites(n_bord);
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());

          if(je_suis_maitre())
            {
              Nusselt << finl;
              Nusselt << "Bord " << le_bord.le_nom() << finl;
              if ( (sub_type(Neumann_paroi,la_cl_th.valeur())))
                {
                  if (dimension == 2)
                    {
                      Nusselt << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << finl;
                      Nusselt << "\tFace a\t\t\t\t|" << finl;
                      Nusselt << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << finl;
                      Nusselt << "X\t\t| Y\t\t\t| dist. carac. (m)\t| Nusselt (local)\t| h (Conv. W/m2/K)\t| Tf cote paroi (K)\t| T face de bord (K)\t| Tparoi equiv.(K) " << finl;
                      Nusselt << "----------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------" << finl;
                    }
                  if (dimension == 3)
                    {
                      Nusselt << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << finl;
                      Nusselt << "\tFace a\t\t\t\t\t\t\t|" << finl;
                      Nusselt << "----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------" << finl;
                      Nusselt << "X\t\t| Y\t\t\t| Z\t\t\t| dist. carac. (m)\t| Nusselt (local)\t| h (Conv. W/m2/K)\t| Tf cote paroi (K)\t| T face de bord (K)\t| Tparoi equiv.(K)" << finl;
                      Nusselt <<        "----------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------" << finl;
                    }
                }
              else
                {
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
                      Nusselt << "----------------------------------------------------------------------------------------------------------------------------------------------------------------" << finl;
                      Nusselt << "\tFace a\t\t\t\t\t\t\t|" << finl;
                      Nusselt << "----------------------------------------------------------------------------------------------------------------------------------------------------------------" << finl;
                      Nusselt << "X\t\t| Y\t\t\t| Z\t\t\t| dist. carac. (m)\t| Nusselt (local)\t| h (Conv. W/m2/K)\t| Tf cote paroi (K)" << finl;
                      Nusselt <<        "----------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------|-----------------------" << finl;
                    }
                }
            }

          ndeb = le_bord.num_premiere_face();
          nfin = ndeb + le_bord.nb_faces();
          for (int num_face=ndeb; num_face<nfin; num_face++)
            {
              // Calcul de la position
              for (int i=0; i<dimension; i++)
                pos[i]=domaine_VEF.xv(num_face,i);

              // Calcul du diametre hydraulique
              if (!dh_constant)
                dh_valeur=diam_hydr->valeur_a_compo(pos,0);

              double x=domaine_VEF.xv(num_face,0);
              double y=domaine_VEF.xv(num_face,1);
              double lambda,lambda_t;
              elem = face_voisins(num_face,0);
              if (elem == -1)
                elem = face_voisins(num_face,1);
              if (sub_type(Champ_Uniforme,conductivite.valeur()))
                lambda = conductivite->valeurs()(0,0);
              else
                {
                  if (conductivite->nb_comp()==1)
                    lambda = conductivite->valeurs()(elem);
                  else
                    lambda = conductivite->valeurs()(elem,0);
                }

              lambda_t=conductivite_turbulente(elem);
              if (dimension == 2)
                Nusselt << x << "\t| " << y;
              if (dimension == 3)
                {
                  double z=domaine_VEF.xv(num_face,2);
                  Nusselt << x << "\t| " << y << "\t| " << z;
                }

              // on doit calculer Tfluide premiere maille sans prendre en compte Tparoi
              double tfluide=0.;
              int nb_faces_elem = domaine_VEF.domaine().nb_faces_elem();
              double surface_face=domaine_VEF.face_surfaces(num_face);
              double surface_pond;
              int j;
              for (int i=0; i<nb_faces_elem; i++)
                {
                  if ( (j=domaine_VEF.elem_faces(elem,i)) != num_face )
                    {
                      surface_pond = 0.;
                      for (int kk=0; kk<dimension; kk++)
                        surface_pond -= (domaine_VEF.face_normales(j,kk)*domaine_VEF.oriente_normale(j,elem)*domaine_VEF.face_normales(num_face,kk)*
                                         domaine_VEF.oriente_normale(num_face,elem))/(surface_face*surface_face);
                      tfluide+=temperature(j)*surface_pond;
                    }
                }

              double d_equiv=equivalent_distance(n_bord,num_face-ndeb);

              if ( (sub_type(Neumann_paroi,la_cl_th.valeur())))
                {
                  // dans ce cas on va imprimer Tfluide (moyenne premiere maille), Tface et on Tparoi recalcule avec d_equiv
                  const Neumann_paroi& la_cl_neum = ref_cast(Neumann_paroi,la_cl_th.valeur());
                  double tparoi = temperature(num_face);
                  double flux = la_cl_neum.flux_impose(num_face-ndeb);
                  double tparoi_equiv = tfluide+flux/lambda*d_equiv;
                  Nusselt << "\t| " << d_equiv << "\t| " << (lambda+lambda_t)/lambda*dh_valeur/d_equiv << "\t| " << (lambda+lambda_t)/d_equiv << "\t| " << tfluide << "\t| " << tparoi << "\t| " << tparoi_equiv << finl;
                }
              else
                {
                  // on imprime Tfluide seulement car normalement Tface=Tparoi est connu
                  Nusselt << "\t| " << d_equiv << "\t| " << (lambda+lambda_t)/lambda*dh_valeur/d_equiv << "\t| " << (lambda+lambda_t)/d_equiv << "\t| " << tfluide << finl;
                }
            }
          Nusselt.syncfile();
        }
    }
  if(je_suis_maitre())
    Nusselt << finl << finl;
  Nusselt.syncfile();
}



