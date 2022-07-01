/****************************************************************************
* Copyright (c) 2017, CEA
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
// File:        Paroi_std_hyd_EF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/EF/Lois_Paroi/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Paroi_std_hyd_EF.h>
#include <Champ_Q1_EF.h>
#include <Fluide_base.h>
#include <Champ_Uniforme.h>
#include <Dirichlet_paroi_fixe.h>
#include <Dirichlet_paroi_defilante.h>
#include <Periodique.h>
#include <Static_Int_Lists.h>
#include <Debog.h>
#include <TRUSTList.h>
#include <EcrFicPartage.h>
#include <Mod_turb_hyd_RANS_0_eq.h>
#include <Neumann_sortie_libre.h>
#include <Mod_turb_hyd_combin.h>
#include <Param.h>
#include <Paroi_rugueuse.h>
#include <SFichier.h>
#include <Paroi_decalee_Robin.h>
#include <Schema_Temps.h>
#include <communications.h>
#include <Equation_base.h>


Implemente_instanciable_sans_constructeur(Paroi_std_hyd_EF,"loi_standard_hydr_EF",Paroi_hyd_base_EF);
Implemente_instanciable(Loi_expert_hydr_EF,"Loi_expert_hydr_EF",Paroi_std_hyd_EF);

Paroi_std_hyd_EF::Paroi_std_hyd_EF(): is_u_star_impose_ (0) {}

Sortie& Paroi_std_hyd_EF::printOn(Sortie& s) const
{
  return s << que_suis_je() << " " << le_nom();
}

Sortie& Loi_expert_hydr_EF::printOn(Sortie& s) const
{
  return s << que_suis_je() << " " << le_nom();
}
//// readOn
//

Entree& Paroi_std_hyd_EF::readOn(Entree& s)
{
  return Paroi_hyd_base_EF::readOn(s);
}

Entree& Loi_expert_hydr_EF::readOn(Entree& s)
{

  Param param(que_suis_je());
  Paroi_std_hyd_EF::set_param(param);

  param.ajouter("u_star_impose",&u_star_impose_);
  param.lire_avec_accolades_depuis(s);

  Nom mot_test;
  mot_test = "u_star_impose";

  if (param.get_list_mots_lus().rang(Motcle(mot_test))!=-1)
    is_u_star_impose_ = 1;

  return s ;
}

void Paroi_std_hyd_EF::set_param(Param& param)
{
  Paroi_hyd_base_EF::set_param(param);
  Paroi_log_QDM::set_param(param);
}

/////////////////////////////////////////////////////////////////////
//
//  Implementation des fonctions de la classe Paroi_std_hyd_EF
//
/////////////////////////////////////////////////////////////////////

int Paroi_std_hyd_EF::init_lois_paroi()
{
  tab_u_star_.resize(la_zone_EF->nb_faces_tot());
  tab_d_plus_.resize(la_zone_EF->nb_faces_tot());
  uplus_.resize(la_zone_EF->nb_faces_tot());
  if (!Cisaillement_paroi_.get_md_vector().non_nul())
    {
      Cisaillement_paroi_.resize(0, dimension);
      la_zone_EF.valeur().creer_tableau_faces(Cisaillement_paroi_);
    }
  seuil_LP.resize(la_zone_EF->nb_faces_tot());
  iterations_LP.resize(la_zone_EF->nb_faces_tot());

  return init_lois_paroi_hydraulique();

}

// Remplissage de la table

int Paroi_std_hyd_EF::init_lois_paroi_hydraulique()
{
  Cmu = mon_modele_turb_hyd->get_Cmu();
  init_lois_paroi_hydraulique_();
  return 1;
}

int Paroi_std_hyd_EF::calculer_hyd_BiK(DoubleTab& tab_k,DoubleTab& tab_eps)
{
  Cerr << " Paroi_std_hyd_EF::calculer_hyd_BiK(DoubleTab& tab_k_eps,DoubleTab& tab_eps) " << finl;
  Cerr <<  "Methode non definie en discretisation EF " << finl ;
  return 1 ;
}

// calculer_hyd pour le k-epsilon
int Paroi_std_hyd_EF::calculer_hyd(DoubleTab& tab_k_eps)
{
  Cerr << " Paroi_std_hyd_EF::calculer_hyd(DoubleTab& tab_k_eps) " << finl;
  Cerr <<  "on ne doit pas entrer dans cette methode" << finl;
  Cerr << " car elle est definie uniquement pour la LES " << finl ;
  return 1 ;
}

int Paroi_std_hyd_EF::calculer_hyd(DoubleTab& tab_nu_t,DoubleTab& tab_k)
{
  const Zone_EF& zone_EF = la_zone_EF.valeur();
  const IntTab& face_voisins = zone_EF.face_voisins();
  const Equation_base& eqn_hydr = mon_modele_turb_hyd->equation();
  const Fluide_base& le_fluide = ref_cast(Fluide_base, eqn_hydr.milieu());
  const Champ_Don& ch_visco_cin = le_fluide.viscosite_cinematique();
  const DoubleTab& vitesse = eqn_hydr.inconnue().valeurs();
  const DoubleTab& tab_visco = ch_visco_cin->valeurs();
  int nsom = zone_EF.nb_som_face();
  int nsom_elem = zone_EF.zone().nb_som_elem();
  const IntTab& elems=zone_EF.zone().les_elems() ;

  ArrOfInt nodes_face(nsom);
  int nb_nodes_free = nsom_elem - nsom;

  double visco=-1;
  int l_unif;
  if (sub_type(Champ_Uniforme,ch_visco_cin.valeur()))
    {
      visco = std::max(tab_visco(0,0),DMINFLOAT);
      l_unif = 1;
    }
  else
    l_unif = 0;
  if ((!l_unif) && (tab_visco.local_min_vect()<DMINFLOAT))
    //   on ne doit pas changer tab_visco ici !
    {
      Cerr<<" visco <=0 ?"<<finl;
      exit();
    }
  //tab_visco+=DMINFLOAT;

  double dist=-1,d_visco=-1;
  double u_plus_d_plus,u_plus,d_plus,u_star;
  double k,eps;
  ArrOfDouble val(dimension);
  ArrOfDouble vit(dimension);
  Cisaillement_paroi_=0;

  double dist_corr=2.;

  bool LM   =(sub_type(Mod_turb_hyd_RANS_0_eq,mon_modele_turb_hyd.valeur()) ? 1 : 0); // Longueur de Melange
  bool COMB =(sub_type(Mod_turb_hyd_combin,mon_modele_turb_hyd.valeur()) ? 1 : 0);  //Modele Combinaison (fonction analytique et (ou) dependance a des champs sources)

  ArrOfDouble vit_face(dimension);

  // Loop on boundaries
  int nb_bords=zone_EF.nb_front_Cl();
  for (int n_bord=0; n_bord<nb_bords; n_bord++)
    {
      const Cond_lim& la_cl = la_zone_Cl_EF->les_conditions_limites(n_bord);

      // Only Dirichlet conditions:
      if (sub_type(Dirichlet_paroi_fixe,la_cl.valeur()) || (sub_type(Dirichlet_paroi_defilante,la_cl.valeur())))
        {
          // Recuperation de la valeur Erugu
          double erugu=Erugu;
          if (sub_type(Paroi_rugueuse,la_cl.valeur()))
            erugu=ref_cast(Paroi_rugueuse,la_cl.valeur()).get_erugu();

          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());

          // Loop on real faces
          ArrOfDouble vit_face_loc(dimension);
          int ndeb = 0;
          int nfin = le_bord.nb_faces_tot();
          for (int ind_face=ndeb; ind_face<nfin; ind_face++)
            {
              int num_face=le_bord.num_face(ind_face);
              int elem = face_voisins(num_face,0);

              // vitesse face CL
              vit_face=0.;
              nodes_face=0;
              for(int jsom=0; jsom<nsom; jsom++)
                {
                  int num_som = zone_EF.face_sommets(num_face, jsom);
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
              double norm_v = norm_vit_lp(vit,num_face,zone_EF,val);

              dist   = distance_face_elem(num_face,elem,zone_EF);
              dist  *= dist_corr;

              if (l_unif)
                d_visco = visco;
              else
                d_visco = (tab_visco.nb_dim()==1 ? tab_visco(elem) : tab_visco(elem,0));

              u_plus_d_plus = norm_v*dist/d_visco;

              u_plus = calculer_u_plus(ind_face,u_plus_d_plus,erugu);

              if (!is_u_star_impose_)
                {
                  if(u_plus)
                    {
                      u_star = norm_v/u_plus ;
                      d_plus = u_plus_d_plus/u_plus ;
                    }
                  else
                    {
                      u_star = 0.;
                      d_plus = 0.;
                    }
                }
              else
                {
                  u_star = u_star_impose_;
                  d_plus = 0.;
                }

              calculer_k_eps(k,eps,d_plus,u_star,d_visco,dist);

              // Calcul de la contrainte tangentielle
              for (int j=0; j<dimension; j++)
                Cisaillement_paroi_(num_face,j) = u_star*u_star*val[j];

              // Remplissage des tableaux (dans le cas de Longueur de melange on laisse la viscosite telle quelle)
              tab_k(elem) = k;

              if((!LM) && (!COMB)) tab_nu_t(elem) = Cmu*k*k/(eps+DMINFLOAT);

              uplus_(num_face) = u_plus;
              tab_d_plus_(num_face) = d_plus;
              tab_u_star_(num_face) = u_star;

              // Modification de nu_t (et par consequent lambda_t) pour exploiter la valeur de nu_t (lambda_t) en y=deq_lam.
              // La valeur de dist_corr n est valable que dans le cas particuler ou nu_t est fonction lineaire de y
              if (COMB)
                {
                  Mod_turb_hyd_combin& modele_turb = ref_cast(Mod_turb_hyd_combin,mon_modele_turb_hyd.valeur());
                  if (modele_turb.nombre_sources()==0)
                    tab_nu_t(elem) *= dist_corr;
                }

            } // End loop on real faces

        } // End Dirichlet conditions

    } // End loop on boundaries

  Cisaillement_paroi_.echange_espace_virtuel();
  tab_nu_t.echange_espace_virtuel();
  tab_k.echange_espace_virtuel();
  Debog::verifier("Paroi_std_hyd_EF::calculer_hyd k",tab_k);
  Debog::verifier("Paroi_std_hyd_EF::calculer_hyd tab_nu_t",tab_nu_t);
  Debog::verifier("Paroi_std_hyd_EF::calculer_hyd Cisaillement_paroi_",Cisaillement_paroi_);
  return 1;
}  // fin du calcul_hyd (nu-t)

int Paroi_std_hyd_EF::calculer_k_eps(double& k, double& eps , double yp, double u_star,
                                     double d_visco, double dist)
{
  // PQ : 05/04/07 : formulation continue de k et epsilon
  //  assurant le bon comportement asymptotique
  double u_star_carre = u_star * u_star;
  k    = 0.07*yp*yp*(exp(-yp/9.))+(1./(sqrt(Cmu)))*pow((1.-exp(-yp/20.)),2);  // k_plus
  k   *= u_star_carre;
  // PL: 50625=15^4 on evite d'utiliser pow car lent
  eps  = (1./(Kappa*pow(yp*yp*yp*yp+50625,0.25)));  // eps_plus
  eps *= u_star_carre*u_star_carre/d_visco;

  return 1;
}

double Paroi_std_hyd_EF::calculer_u_plus(const int ind_face,const double u_plus_d_plus,const double erugu)
{
  // PQ : 05/04/07 : formulation continue de la loi de paroi
  // construite d'apres la loi de Reichardt

  if(u_plus_d_plus<1.) return sqrt(u_plus_d_plus); // pour eviter de faire tourner la procedure iterative


  double u_plus = u_plus_d_plus/100.;
  double r      = 1.;
  double seuil  = 0.001;
  int    iter   = 0;
  int    itmax  = 25;

  double A = (1/Kappa)*log(erugu/Kappa) ; //  (=7.44, contre 7.8 dans la loi d'origine)
  //  permettant d'avoir en l'infini la loi de
  //  Reichardt se calant sur : u+ = (1/Kappa).ln(Erugu.y+)
  double d_plus,u_plus2;

  while((iter++<itmax) && (r>seuil))
    {
      d_plus  = u_plus_d_plus/ u_plus ;
      u_plus2 = ((1./Kappa)*log(1.+Kappa*d_plus))+A*(1.-exp(-d_plus/11.)-exp(-d_plus/3.)*d_plus/11.); // Equation de Reichardt
      u_plus  = 0.5*(u_plus+u_plus2);
      r       = std::fabs(u_plus-u_plus2)/u_plus;
    }

  seuil_LP(ind_face) = r;
  iterations_LP(ind_face) = iter;

  if (iter >= itmax) erreur_non_convergence();

  return u_plus;
}

void Paroi_std_hyd_EF::imprimer_ustar(Sortie& os) const
{
  const Zone_EF& zone_EF = la_zone_EF.valeur();
  int ndeb,nfin;
  double upmoy,dpmoy,utaumoy;
  double seuil_moy,iter_moy;
  double norme_L2=0.;

  upmoy=0.;
  dpmoy=0.;
  utaumoy=0.;
  seuil_moy=0.;
  iter_moy=0.;

  int compt=0;

  EcrFicPartage Ustar;
  ouvrir_fichier_partage(Ustar,"Ustar");

  for (int n_bord=0; n_bord<zone_EF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = la_zone_Cl_EF->les_conditions_limites(n_bord);
      if ( (sub_type(Dirichlet_paroi_fixe,la_cl.valeur())) ||
           (sub_type(Dirichlet_paroi_defilante,la_cl.valeur())) ||
           (sub_type(Paroi_decalee_Robin,la_cl.valeur()) ) ||
           (la_cl.valeur().que_suis_je() == "Entree_fluide_vitesse_imposee_ALE"))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          if(je_suis_maitre())
            {
              Ustar << finl;
              Ustar << "Bord " << le_bord.le_nom() << finl;
              if (dimension == 2)
                {
                  Ustar << "-------------------------------------------------------------------------------------------";
                  Ustar << "--------------------------------------------------------------------------------------------" << finl;
                  Ustar << "\tFace a\t\t\t\t|\t\t\t\t\t\t\t\t\t| TAU=Nu.Grad(Ut) [m2/s2]" << finl;
                  Ustar << "----------------------------------------|--------------------------------------------------";
                  Ustar << "---------------------|----------------------------------------------------------------------" << finl;
                  Ustar << "X\t\t| Y\t\t\t| u+\t\t\t| d+\t\t\t| u*\t\t\t| ||TAU||_2\t\t| |TAUx|\t\t| |TAUy|" << finl;
                  Ustar << "----------------|-----------------------|-----------------------|-----------------------|--";
                  Ustar << "---------------------|-----------------------|-----------------------|----------------------" << finl;
                }
              if (dimension == 3)
                {
                  Ustar << "-----------------------------------------------------------------------------------------------------------------";
                  Ustar << "-----------------------------------------------------------------------------------------------------------------" << finl;
                  Ustar << "\tFace a\t\t\t\t\t\t\t|\t\t\t\t\t\t\t\t\t| TAU=Nu.Grad(Ut) [m2/s2]" << finl;
                  Ustar << "----------------------------------------------------------------|------------------------------------------------";
                  Ustar << "-----------------------|-----------------------------------------------------------------------------------------" << finl;
                  Ustar << "X\t\t| Y\t\t\t| Z\t\t\t| u+\t\t\t| d+\t\t\t| u*\t\t\t| ||TAU||_2\t\t| |TAUx|\t\t| |TAUy|\t\t| |TAUz|" << finl;
                  Ustar << "----------------|-----------------------|-----------------------|-----------------------|-----------------------|";
                  Ustar << "-----------------------|-----------------------|-----------------------|-----------------------|-----------------" << finl;
                }
            }
          ndeb = le_bord.num_premiere_face();
          nfin = ndeb + le_bord.nb_faces();
          for (int num_face=ndeb; num_face<nfin; num_face++)
            {
              double x=zone_EF.xv(num_face,0);
              double y=zone_EF.xv(num_face,1);
              norme_L2= Cisaillement_paroi_(num_face,0)*Cisaillement_paroi_(num_face,0) + Cisaillement_paroi_(num_face,1)*Cisaillement_paroi_(num_face,1);
              if (dimension == 2)
                Ustar << x << "\t| " << y;
              if (dimension == 3)
                {
                  double z=zone_EF.xv(num_face,2);
                  Ustar << x << "\t| " << y << "\t| " << z;
                  norme_L2+= Cisaillement_paroi_(num_face,2)*Cisaillement_paroi_(num_face,2);
                }
              norme_L2=sqrt(norme_L2);
              Ustar << "\t| " << uplus_(num_face) << "\t| " << tab_d_plus(num_face) << "\t| " << tab_u_star(num_face);
              Ustar << "\t| " << norme_L2 << "\t| " << Cisaillement_paroi_(num_face,0) << "\t| " << Cisaillement_paroi_(num_face,1) ;
              if (dimension == 3)
                Ustar << "\t| " << Cisaillement_paroi_(num_face,2) << finl;
              else
                Ustar << finl;

              // PQ : 03/03 : Calcul des valeurs moyennes (en supposant maillage regulier)

              upmoy +=uplus_(num_face);
              dpmoy +=tab_d_plus(num_face);
              utaumoy +=tab_u_star(num_face);
              seuil_moy += seuil_LP(num_face);
              iter_moy += iterations_LP(num_face);
              compt +=1;
            }
          Ustar.syncfile();
        }
    }
  /* Reduce 6 mp_sum to 1 by using mp_sum_for_each_item:
  upmoy = mp_sum(upmoy);
  dpmoy = mp_sum(dpmoy);
  utaumoy = mp_sum(utaumoy);
  seuil_moy = mp_sum(seuil_moy);
  iter_moy = mp_sum(iter_moy);
  compt = mp_sum(compt);
  */

  ArrOfDouble array(6);
  array[0]=upmoy;
  array[1]=dpmoy;
  array[2]=utaumoy;
  array[3]=seuil_moy;
  array[4]=iter_moy;
  array[5]=compt;
  mp_sum_for_each_item(array);
  upmoy=array[0];
  dpmoy=array[1];
  utaumoy=array[2];
  seuil_moy=array[3];
  iter_moy=array[4];
  compt=(int)array[5];
  if (je_suis_maitre())
    {
      if (compt)
        {
          Ustar << finl;
          Ustar << "-------------------------------------------------------------" << finl;
          Ustar << "Calcul des valeurs moyennes (en supposant maillage regulier):" << finl;
          Ustar << "<u+>= " << upmoy/compt << " <d+>= " << dpmoy/compt << " <u*>= " << utaumoy/compt << " seuil_LP= " << seuil_moy/compt << " iterations_LP= " << iter_moy/compt << finl;
        }
      Ustar << finl << finl;
    }
  Ustar.syncfile();
}


