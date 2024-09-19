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
// File:        Paroi_loi_WW_hyd_EF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/EF/Lois_Paroi/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Paroi_loi_WW_hyd_EF.h>
#include <Fluide_base.h>
#include <Champ_Uniforme.h>
#include <Dirichlet_paroi_fixe.h>
#include <Dirichlet_paroi_defilante.h>
#include <Domaine.h>
#include <Modele_turbulence_hyd_base.h>
#include <Equation_base.h>
#include <Param.h>

Implemente_instanciable_sans_constructeur(Paroi_loi_WW_hyd_EF,"loi_WW_hydr_EF",Paroi_std_hyd_EF);

// PrintOn
Sortie& Paroi_loi_WW_hyd_EF::printOn(Sortie& s ) const
{
  return s << que_suis_je() << " " << le_nom();
}

// ReadOn
Entree& Paroi_loi_WW_hyd_EF::readOn(Entree& is )
{
  Cerr<<"Reading of data for a "<<que_suis_je()<<" wall law"<<finl;
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);
  return is ;
}

void Paroi_loi_WW_hyd_EF::set_param(Param& param)
{
  Paroi_std_hyd_EF::set_param(param);
  param.ajouter_flag("impr",&impr);
}

/////////////////////////////////////////////////////////////////////
//
//  Implementation des fonctions de la classe Paroi_loi_WW_hyd_EF
//
/////////////////////////////////////////////////////////////////////

// Remplissage de la table
int Paroi_loi_WW_hyd_EF::init_lois_paroi_hydraulique()
{
  Cmu_ = mon_modele_turb_hyd->get_Cmu();
  A= 8.3 ;
  B= 1./7. ;
  Y0= 11.81 ;
  return 1;
}

// On annule les valeurs des grandeurs turbulentes qui
// correspondent aux mailles de paroi
int Paroi_loi_WW_hyd_EF::preparer_calcul_hyd(DoubleTab& tab)
{
  return 1;
}

// calculer_hyd pour le k-epsilon
int Paroi_loi_WW_hyd_EF::calculer_hyd(DoubleTab& tab_k_eps)
{
  Cerr << " Paroi_loi_WW_hyd_EF::calculer_hyd(DoubleTab& tab_k_eps) " << finl;
  Cerr <<  "on ne doit pas entrer dans cette methode" << finl;
  Cerr << " car elle est definie uniquement pour la LES " << finl ;
  return 1 ;
}



int Paroi_loi_WW_hyd_EF::calculer_hyd(DoubleTab& tab_nu_t,DoubleTab& tab_k)
{
  const Domaine_EF& domaine_EF = le_dom_EF.valeur();
  const IntTab& face_voisins = domaine_EF.face_voisins();
  const Equation_base& eqn_hydr = mon_modele_turb_hyd->equation();
  const Fluide_base& le_fluide = ref_cast(Fluide_base, eqn_hydr.milieu());
  const Champ_Don& ch_visco_cin = le_fluide.viscosite_cinematique();
  const DoubleTab& vitesse = eqn_hydr.inconnue().valeurs();
  const DoubleTab& tab_visco = ch_visco_cin->valeurs();
  double visco=-1;
  int l_unif;
  if (sub_type(Champ_Uniforme,ch_visco_cin.valeur()))
    {
      visco = std::max(tab_visco(0,0),DMINFLOAT);
      l_unif = 1;
    }
  else
    l_unif = 0;

  // preparer_calcul_hyd(tab);
  if ((!l_unif) && (tab_visco.local_min_vect()<DMINFLOAT))
    //   on ne doit pas changer tab_visco ici !
    {
      Cerr<<" visco <=0 ?"<<finl;
      exit();
    }
  //tab_visco+=DMINFLOAT;

  int ndeb,nfin;
  double norm_v;
  double dist_G,dist_som;
  double d_visco;

  ArrOfDouble vit(dimension);
  ArrOfDouble val(dimension);

  int nsom = domaine_EF.nb_som_face();
  int nsom_elem = domaine_EF.domaine().nb_som_elem();
  const IntTab& elems=domaine_EF.domaine().les_elems() ;

  ArrOfDouble vit_face(dimension);
  ArrOfInt nodes_face(nsom);
  int nb_nodes_free = nsom_elem - nsom;

  // Boucle sur les bords

  for (int n_bord=0; n_bord<domaine_EF.nb_front_Cl(); n_bord++)
    {

      // pour chaque condition limite on regarde son type
      // On applique les lois de paroi uniquement
      // aux voisinages des parois

      const Cond_lim& la_cl = le_dom_Cl_EF->les_conditions_limites(n_bord);

      if (sub_type(Dirichlet_paroi_fixe,la_cl.valeur()) )
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());

          // Loop on real faces
          ndeb = 0;
          nfin = le_bord.nb_faces_tot();
          for (int ind_face=ndeb; ind_face<nfin; ind_face++)
            {
              int num_face=le_bord.num_face(ind_face);
              int elem = face_voisins(num_face,0);

              // vitesse face CL
              vit_face=0.;
              nodes_face=0;
              for(int jsom=0; jsom<nsom; jsom++)
                {
                  int num_som = domaine_EF.face_sommets(num_face, jsom);
                  nodes_face[jsom] = num_som;
                  for(int comp=0; comp<dimension; comp++) vit_face[comp]+=vitesse(num_som,comp)/nsom;
                }

              vit=0.;
              // Loop on nodes : vitesse moyenne des noeuds n'appartenant pas a la face CL
              for (int i=0; i<nsom_elem; i++)
                {
                  int node=elems(elem,i);
                  int IOK = 1;
                  for(int jsom=0; jsom<nsom; jsom++)
                    if (nodes_face[jsom] == node) IOK=0;
                  // Le noeud contribue
                  if (IOK)
                    for (int j=0; j<dimension; j++)
                      vit[j]+=(vitesse(node,j)-vit_face[j])/nb_nodes_free; // permet de soustraire la vitesse de glissement eventuelle
                }
              norm_v = norm_vit_lp(vit,num_face,domaine_EF,val);

              if (l_unif)
                d_visco = visco;
              else
                d_visco = tab_visco[elem];

              dist_G = distance_face_elem(num_face,elem,domaine_EF);
              dist_som = 2.0*dist_G;

              // Calcul de u* et des grandeurs turbulentes:
              // a partir de de norm_v
              calculer_local(d_visco,tab_nu_t,tab_k,norm_v,dist_G,dist_som,elem,num_face);

              // Calcul de la contrainte tangentielle:
              // Signe du cisaillement : on le donne dans l OpDift
              double vit_frot = tab_u_star(num_face)*tab_u_star(num_face);
              for (int j=0; j<dimension; j++) Cisaillement_paroi_(num_face,j) = vit_frot*val[j];
            }
        }  //  fin de Dirichlet_paroi_fixe


      // ATTENTION : CODER PAROI DEFILANTE
      else if (sub_type(Dirichlet_paroi_defilante,la_cl.valeur()) )
        {
          Cerr << "pour l instant Werner et Wengle n'est pas code pour la condition de paroi defilante!!" << finl;
        } // fin de Dirichlet_paroi_defilante
    }  // fin du for bord CL

  Cisaillement_paroi_.echange_espace_virtuel();
  tab_nu_t.echange_espace_virtuel();
  tab_k.echange_espace_virtuel();
  return 1;
}  // fin du calcul_hyd (nu-t)


int Paroi_loi_WW_hyd_EF::calculer_local(double d_visco,
                                        DoubleTab& tab_nu_t,DoubleTab& tab_k,double norm_vit,
                                        double dist_G,double dist_som,int elem,int num_face)
{
  // C est la hauteur de la premiere maille en int qui doit etre testee dist_som
  double up_lim = d_visco*Y0*Y0/dist_som;

  if (norm_vit <= up_lim)
    {
      calculer_u_star_sous_couche_visq(norm_vit,d_visco,dist_som,num_face);
      //calculer_sous_couche_visq(tab_nu_t,tab_k,elem);
      tab_nu_t(elem) = 0.;
      tab_k(elem) = 0.;

      if (impr==1)  Cerr << "Domaine lineaire" << finl;
    }
  else
    {
      calculer_u_star_couche_puissance(norm_vit,d_visco,dist_som,num_face);
      calculer_couche_puissance(tab_nu_t,tab_k,dist_som,elem,num_face);
      if (impr==1)  Cerr << "Domaine en puissance" << finl;
    }
  return 1;
}


int Paroi_loi_WW_hyd_EF::calculer_u_star_sous_couche_visq(double norm_vit,
                                                          double d_visco,double dist,
                                                          int face)
{
  // Dans la sous couche visqueuse:  u* = sqrt(u*nu/d)

  tab_u_star_(face) = sqrt(norm_vit*d_visco/dist);
  tab_d_plus_(face) = (dist*tab_u_star_(face))/d_visco;
  uplus_(face) = norm_vit/tab_u_star_(face);

  return 1;
}


int Paroi_loi_WW_hyd_EF::calculer_u_star_couche_puissance(double norm_vit,double d_visco,
                                                          double dist, int face)
{
  // Implementation Werner, Wengle 1993
  // double part1, part2 ;
  // static double Apuiss = pow(A,(1+B)/(1-B));

  // part1= ( (B+1) * pow(d_visco,B) * norm_vit ) / (A *pow(2.*dist,B) );
  // part2= ( (1-B) * pow(d_visco,B+1) * Apuiss ) / (2 * pow(2.*dist,B+1) ) ;

  // tab_u_star_(face)= pow(part1+part2, 1/(B+1) );

  // Implementation Wilhelm, Jacob, Sagaut 2018
  double part1bis;
  part1bis= ( pow(d_visco,B) * norm_vit ) / (A *pow(dist,B) );
  tab_u_star_(face)= pow(part1bis, 1/(B+1) );
  tab_d_plus_(face) = (dist*tab_u_star_(face))/d_visco;
  uplus_(face) = norm_vit/tab_u_star_(face);

  return 1;
}


int Paroi_loi_WW_hyd_EF::calculer_couche_puissance(DoubleTab& nu_t,DoubleTab& tab_k,
                                                   double dist,int elem,int face)
{
  //  nu_t = Cmu*k*k/eps
  //
  //                          2                       3
  //  En utilisant  k =     u*/sqrt(Cmu_)  et eps = u* / Kd
  //
  //  on calcule nu_t en fonction de u*

  double u_star = tab_u_star(face);

  tab_k(elem) = u_star*u_star/sqrt(Cmu_);
  nu_t(elem) =  u_star*Kappa_*dist ;

  return 1;
}
