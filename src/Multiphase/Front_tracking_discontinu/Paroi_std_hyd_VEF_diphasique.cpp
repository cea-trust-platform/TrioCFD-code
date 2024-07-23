/****************************************************************************
* Copyright (c) 2019, CEA
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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Paroi_std_hyd_VEF_diphasique.cpp
// Directory : $TrioCFD_ROOT/Front_tracking_discontinu/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Champ_P1NC.h>
#include <Fluide_Incompressible.h>
#include <Fluide_Diphasique.h>
#include <Dirichlet_paroi_defilante.h>
#include <Debog.h>
#include <EcrFicPartage.h>
#include <Modele_turbulence_hyd_Longueur_Melange_base.h>
#include <Modele_turbulence_hyd_combinaison.h>
#include <Paroi_rugueuse.h>
#include <Paroi_decalee_Robin.h>
#include <Paroi_std_hyd_VEF_diphasique.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <Probleme_base.h>
#include <Champ_Uniforme.h>

Implemente_instanciable( Paroi_std_hyd_VEF_diphasique, "loi_standard_hydr_diphasique_VEF", Paroi_std_hyd_VEF ) ;

Sortie& Paroi_std_hyd_VEF_diphasique::printOn( Sortie& os ) const
{
  return os << que_suis_je() << " " << le_nom();
}

Entree& Paroi_std_hyd_VEF_diphasique::readOn( Entree& is )
{
  Paroi_std_hyd_VEF::readOn(is);
  return is;
}

extern double norm_vit_lp(const ArrOfDouble& vit,int face,const Domaine_VEF& domaine,ArrOfDouble& val);
/* Codee classe mere
double norm_vit_lp(const ArrOfDouble& vit,int face,const Domaine_VEF& domaine,ArrOfDouble& val)
{
  // A reverser dans VEF/Domaine (?)

  const DoubleTab& face_normale = domaine.face_normales();
  int dim=Objet_U::dimension;
  ArrOfDouble r(dim);
  double psc,norm_vit;

  for(int i=0; i<dim; i++) r[i]=face_normale(face,i);

  r/=norme_array(r);
  psc = dotproduct_array(r,vit);

  if(dim==3) norm_vit = vitesse_tangentielle(vit[0],vit[1],vit[2],r[0],r[1],r[2]);
  else             norm_vit = vitesse_tangentielle(vit[0],vit[1],r[0],r[1]);

  for(int i=0; i<dim; i++) val[i]=(vit[i]-psc*r[i])/(norm_vit+DMINFLOAT);

  return norm_vit;
} */


int Paroi_std_hyd_VEF_diphasique::calculer_hyd(DoubleTab& tab_nu_t,DoubleTab& tab_k)
{
  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  const IntTab& face_voisins = domaine_VEF.face_voisins();
  const Equation_base& eqn_hydr = mon_modele_turb_hyd->equation();
  const DoubleTab& vitesse = eqn_hydr.inconnue().valeurs();
  // Physical properties of both phases
  const Fluide_Diphasique& le_fluide = ref_cast(Fluide_Diphasique, eqn_hydr.milieu());
  const Fluide_Incompressible& phase_1 = le_fluide.fluide_phase(1);
  const Fluide_Incompressible& phase_0 = le_fluide.fluide_phase(0);
  const Champ_Don& ch_visco_cin_ph1 = phase_1.viscosite_cinematique();
  const Champ_Don& ch_visco_cin_ph0 = phase_0.viscosite_cinematique();
  const DoubleTab& tab_visco_ph1 = phase_1.viscosite_cinematique()->valeurs();
  const DoubleTab& tab_visco_ph0 = phase_0.viscosite_cinematique()->valeurs();
  const double delta_nu = tab_visco_ph1(0,0) - tab_visco_ph0(0,0);

  // One way to get the Transport equation to pass the indicator DoubleTab
  const Domaine_Cl_dis_base& domaine_Cl_dis_base = eqn_hydr.domaine_Cl_dis().valeur();
  const Equation_base& eqn_trans = domaine_Cl_dis_base.equation().probleme().equation("Transport_Interfaces_FT_Disc");
  const Transport_Interfaces_FT_Disc& eqn_interf = ref_cast(Transport_Interfaces_FT_Disc, eqn_trans);
  const DoubleTab& indic = eqn_interf.inconnue().valeurs();

  const Domaine& domaine = domaine_VEF.domaine();
  int nfac = domaine.nb_faces_elem();

  double visco_ph0=-1;
  int l_unif;

  if (sub_type(Champ_Uniforme,ch_visco_cin_ph1.valeur()) && sub_type(Champ_Uniforme,ch_visco_cin_ph0.valeur()) )
    {
      visco_ph0 = std::max(tab_visco_ph0(0,0),DMINFLOAT);
      l_unif = 1;
    }
  else
    l_unif = 0;
  if ((!l_unif) && ((tab_visco_ph1.local_min_vect()<DMINFLOAT) || (tab_visco_ph0.local_min_vect()<DMINFLOAT) ))
    {
      Cerr << "Negative viscosity !!!" << finl;
      Process::exit();
    }

  double dist=-1,d_visco=-1;
  double u_plus_d_plus,u_plus,d_plus,u_star;
  double k,eps;
  ArrOfDouble val(dimension);
  ArrOfDouble vit(dimension);
  Cisaillement_paroi_=0;

  double dist_corr=1.;
  double coef_vit=nfac;
  if (sub_type(Champ_P1NC,eqn_hydr.inconnue().valeur()))
    {
      dist_corr=double(dimension+1)/double(dimension);
      coef_vit=nfac-1;
    }

  bool LM   =(sub_type(Modele_turbulence_hyd_Longueur_Melange_base,mon_modele_turb_hyd.valeur()) ? 1 : 0); // Longueur de Melange
  bool COMB =(sub_type(Modele_turbulence_hyd_combinaison,mon_modele_turb_hyd.valeur()) ? 1 : 0);  //Modele Combinaison (fonction analytique et (ou) dependance a des champs sources)

  // Loop on boundaries
  int nb_bords=domaine_VEF.nb_front_Cl();
  for (int n_bord=0; n_bord<nb_bords; n_bord++)
    {
      const Cond_lim& la_cl = le_dom_Cl_VEF->les_conditions_limites(n_bord);

      // Only Dirichlet conditions:
      if (sub_type(Dirichlet_paroi_fixe,la_cl.valeur()) || (sub_type(Dirichlet_paroi_defilante,la_cl.valeur())))
        {
          // Recuperation de la valeur Erugu
          double erugu=Erugu;
          if (sub_type(Paroi_rugueuse,la_cl.valeur()))
            erugu=ref_cast(Paroi_rugueuse,la_cl.valeur()).get_erugu();

          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          const IntTab& elem_faces = domaine_VEF.elem_faces();

          // Loop on real faces
          int ndeb = 0;
          int nfin = le_bord.nb_faces_tot();

          for (int ind_face=ndeb; ind_face<nfin; ind_face++)
            {
              int num_face=le_bord.num_face(ind_face);
              int elem = face_voisins(num_face,0);

              vit=0.;
              for (int i=0; i<nfac; i++)
                {
                  int face=elem_faces(elem,i);
                  for (int j=0; j<dimension; j++)
                    vit[j]+=(vitesse(face,j)-vitesse(num_face,j)); // permet de soustraire la vitesse de glissement eventuelle
                }
              vit   /= coef_vit;
              dist   = distance_face_elem(num_face,elem,domaine_VEF);
              dist  *= dist_corr; // pour passer du centre de gravite au milieu des faces en P1NC

              double norm_v = norm_vit_lp(vit,num_face,domaine_VEF,val);

              if (l_unif)
                d_visco = visco_ph0 + indic(elem) * delta_nu;
              else
                d_visco = (tab_visco_ph0.nb_dim()==1 ? (tab_visco_ph0(elem) + indic(elem) * delta_nu) : (tab_visco_ph0(elem,0) + indic(elem) * delta_nu));

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
                  Modele_turbulence_hyd_combinaison& modele_turb = ref_cast(Modele_turbulence_hyd_combinaison,mon_modele_turb_hyd.valeur());
                  if (modele_turb.nombre_sources()==0)
                    tab_nu_t(elem) *= dist_corr;
                }

            } // End loop on real faces

        } // End Dirichlet conditions

      // Robin condition:
      else if (sub_type(Paroi_decalee_Robin,la_cl.valeur()))
        {
          // Recuperation de la valeur Erugu
          double erugu=Erugu;

          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          const Paroi_decalee_Robin& Paroi = ref_cast(Paroi_decalee_Robin,la_cl.valeur());
          const DoubleTab& normales = domaine_VEF.face_normales();
          double delta = Paroi.get_delta();

          // Loop on real faces
          int ndeb = 0;
          int nfin = le_bord.nb_faces_tot();
          for (int ind_face=ndeb; ind_face<nfin; ind_face++)
            {
              int num_face = le_bord.num_face(ind_face);
              int elem = face_voisins(num_face,0);

              double psc=0, norm=0;
              double norm_v=0;

              for(int comp=0; comp<dimension; comp++)
                {
                  psc += vitesse(num_face,comp)*normales(num_face,comp);
                  norm += normales(num_face,comp)*normales(num_face,comp);
                }
              // psc /= norm; // Fixed bug: Arithmetic exception
              if (std::fabs(norm)>=DMINFLOAT) psc/=norm;

              for(int comp=0; comp<dimension; comp++)
                {
                  val[comp]=vitesse(num_face,comp)-psc*normales(num_face,comp);
                  norm_v += val[comp]*val[comp];
                }

              norm_v = sqrt(norm_v);
              val /= norm_v;
              dist = delta;

              // Common to Dirichlet

              if (l_unif)
                d_visco = visco_ph0 + indic(elem) * delta_nu;
              else
                d_visco = (tab_visco_ph0.nb_dim()==1 ? (tab_visco_ph0(elem) + indic(elem) * delta_nu) : (tab_visco_ph0(elem,0) + indic(elem) * delta_nu));

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
                  Modele_turbulence_hyd_combinaison& modele_turb = ref_cast(Modele_turbulence_hyd_combinaison,mon_modele_turb_hyd.valeur());
                  if (modele_turb.nombre_sources()==0)
                    tab_nu_t(elem) *= dist_corr;
                }

              // End common to Dirichlet

            } // End loop on real faces

        } // End Robin condition

    } // End loop on boundaries

  Cisaillement_paroi_.echange_espace_virtuel();
  tab_nu_t.echange_espace_virtuel();
  tab_k.echange_espace_virtuel();
  Debog::verifier("Paroi_std_hyd_VEF::calculer_hyd k",tab_k);
  Debog::verifier("Paroi_std_hyd_VEF::calculer_hyd tab_nu_t",tab_nu_t);
  Debog::verifier("Paroi_std_hyd_VEF::calculer_hyd Cisaillement_paroi_",Cisaillement_paroi_);

  return 1;
}  // fin du calcul_hyd (nu-t)

int Paroi_std_hyd_VEF_diphasique::calculer_hyd(DoubleTab& tab_k_eps)
{
  Cerr << " Paroi_std_hyd_VEF_diphasique::calculer_hyd(DoubleTab& tab_k_eps) " << finl;
  Cerr <<  "on ne doit pas entrer dans cette methode" << finl;
  Cerr << " car elle est definie uniquement pour la LES " << finl ;
  Process::exit();
  return 1 ;
} // fin de calcul_hyd (K-eps)


