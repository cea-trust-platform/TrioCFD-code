/****************************************************************************
* Copyright (c) 2021, CEA
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
// File:        Dirichlet_loi_paroi_QDM.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Incompressible/Cond_Lim
// Version:     /main/28
//
//////////////////////////////////////////////////////////////////////////////

#include <Dirichlet_loi_paroi_QDM.h>

#include <Viscosite_turbulente_k_omega.h>
#include <Viscosite_turbulente_k_tau.h>
#include <Loi_paroi_adaptative.h>
#include <Operateur_Diff_base.h>
#include <Champ_Face_base.h>
#include <Pb_Multiphase.h>
#include <Domaine_VF.h>
#include <math.h>

Implemente_instanciable(Dirichlet_loi_paroi_QDM,"Dirichlet_loi_paroi_QDM",Dirichlet_loi_paroi);

Sortie& Dirichlet_loi_paroi_QDM::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

Entree& Dirichlet_loi_paroi_QDM::readOn(Entree& s )
{
  Param param(que_suis_je());
  param.ajouter("beta_omega", &beta_omega);
  param.ajouter("beta_k", &beta_k);
  param.ajouter("von_karman", &von_karman_);
  param.ajouter("y_p_prod_k", &y_p_prod_k_);
  param.ajouter("fac_prod_k", &fac_prod_k_);
  param.ajouter("y_p_prod_k_grand", &y_p_prod_k_grand_);
  param.ajouter("fac_prod_k_grand", &fac_prod_k_grand_);
  param.lire_avec_accolades_depuis(s);

  le_champ_front.typer("Champ_front_vide");

  Dirichlet_loi_paroi::readOn(s);
  return s;
}

void Dirichlet_loi_paroi_QDM::completer()
{
  if (!ref_cast(Operateur_Diff_base, domaine_Cl_dis().equation().operateur(0).l_op_base()).is_turb()) Process::exit(que_suis_je() + " : diffusion operator must be turbulent !");
  if sub_type(Viscosite_turbulente_k_tau, (*ref_cast(Operateur_Diff_base, domaine_Cl_dis().equation().operateur(0).l_op_base()).correlation_viscosite_turbulente()).valeur())
    {
      if (fac_prod_k_<-1.e7) fac_prod_k_ = 1.2;
      if (y_p_prod_k_<-1.e7) y_p_prod_k_ =  4.;
      if (fac_prod_k_grand_<-1.e7) fac_prod_k_grand_ = .2;
      if (y_p_prod_k_grand_<-1.e7) y_p_prod_k_grand_ = 150.;
    }
  else if sub_type(Viscosite_turbulente_k_omega, (*ref_cast(Operateur_Diff_base, domaine_Cl_dis().equation().operateur(0).l_op_base()).correlation_viscosite_turbulente()).valeur())
    {
      if (fac_prod_k_<-1.e7) fac_prod_k_ = 1.0;
      if (y_p_prod_k_<-1.e7) y_p_prod_k_ =  4.;
      if (fac_prod_k_grand_<-1.e7) fac_prod_k_grand_ = .6;
      if (y_p_prod_k_grand_<-1.e7) y_p_prod_k_grand_ = 120.;
    }
  else
    {
      if (fac_prod_k_<-1.e7) fac_prod_k_ =  0.;
      if (y_p_prod_k_<-1.e7) y_p_prod_k_ =  4.;
      if (fac_prod_k_grand_<-1.e7) fac_prod_k_grand_ =  0.;
      if (y_p_prod_k_grand_<-1.e7) y_p_prod_k_grand_ = 120.;
    }
}

int Dirichlet_loi_paroi_QDM::compatible_avec_eqn(const Equation_base& eqn) const
{
  Motcle dom_app=eqn.domaine_application();
  Motcle Turbulence="Hydraulique";

  if (dom_app==Turbulence) return 1;
  else err_pas_compatible(eqn); return 0;
}

void Dirichlet_loi_paroi_QDM::mettre_a_jour(double tps)
{
  if (mon_temps!=tps) {me_calculer() ; mon_temps=tps;}
}

int Dirichlet_loi_paroi_QDM::initialiser(double temps)
{
  d_.resize(0,domaine_Cl_dis().equation().inconnue().valeurs().line_size()*dimension);
  la_frontiere_dis.valeur().frontiere().creer_tableau_faces(d_);

  correlation_loi_paroi_ = ref_cast(Pb_Multiphase, domaine_Cl_dis().equation().probleme()).get_correlation("Loi_paroi");

  return 1;
}

void Dirichlet_loi_paroi_QDM::me_calculer()
{
  Loi_paroi_adaptative& corr_loi_paroi = ref_cast(Loi_paroi_adaptative, correlation_loi_paroi_.valeur().valeur());
  const Domaine_VF& domaine = ref_cast(Domaine_VF, domaine_Cl_dis().equation().probleme().domaine_dis().valeur());

  const DoubleTab& u_tau = corr_loi_paroi.get_tab("u_tau"), y = corr_loi_paroi.get_tab("y"); // ces tables sont numerotees selon les faces du domaine
  const DoubleTab& visc  = ref_cast(Navier_Stokes_std, domaine_Cl_dis().equation().probleme().equation(0)).diffusivite_pour_pas_de_temps().valeurs();
  const DoubleTab& vit   = domaine_Cl_dis().equation().probleme().get_champ("vitesse").valeurs(),
                   &rho = domaine_Cl_dis().equation().probleme().get_champ("masse_volumique").valeurs(),
                    *alp = sub_type(Pb_Multiphase, domaine_Cl_dis().equation().probleme()) ? &domaine_Cl_dis().equation().probleme().get_champ("alpha").valeurs() : nullptr;

  int nf = la_frontiere_dis.valeur().frontiere().nb_faces(), nf_tot = domaine.nb_faces_tot(), f1 = la_frontiere_dis.valeur().frontiere().num_premiere_face();
  int N = domaine_Cl_dis().equation().inconnue().valeurs().line_size(), D = dimension;
  int Nv = domaine_Cl_dis().equation().probleme().equation(0).inconnue().valeurs().line_size();

  const DoubleTab& n_f = domaine.face_normales();
  const DoubleVect& fs = domaine.face_surfaces();
  const IntTab& f_e = domaine.face_voisins();

  // Recuperation des vitesses aux elements
  DoubleTab pvit_elem(0, Nv * dimension);
  if (nf_tot == vit.dimension_tot(0))
  {
  const Champ_Face_base& ch = ref_cast(Champ_Face_base, domaine_Cl_dis().equation().probleme().equation(0).inconnue().valeur());
  domaine.domaine().creer_tableau_elements(pvit_elem);
  ch.get_elem_vector_field(pvit_elem, true);
  }

  // On recupere la diffusivite totale
  DoubleTab mu = ref_cast(Navier_Stokes_std, domaine_Cl_dis().equation().probleme().equation(0)).diffusivite_pour_transport().passe() ;  // Copie expres !!!
  if (ref_cast(Operateur_Diff_base, domaine_Cl_dis().equation().operateur(0).l_op_base()).is_turb())
  {
  const Viscosite_turbulente_base& corr_visc = ref_cast(Viscosite_turbulente_base, (*ref_cast(Operateur_Diff_base, domaine_Cl_dis().equation().operateur(0).l_op_base()).correlation_viscosite_turbulente()).valeur());
  corr_visc.modifier_mu(mu) ;
  }

  int n = 0 ; // la phase turbulente frotte
  for (int f =0 ; f < nf ; f++)
    {
      int f_domaine = f + f1; // number of the face in the domaine
      int e_domaine = (f_e(f_domaine,0)>=0) ? f_e(f_domaine,0) : f_e(f_domaine,1) ; // Make orientation vdf-proof

      double u_orth = 0 ;
      DoubleTrav u_parallel(D);
      if (nf_tot == vit.dimension_tot(0)) // PolyMAC
      {
        for (int d = 0; d <D ; d++) u_orth -= pvit_elem(e_domaine, Nv*d+n)*n_f(f_domaine,d)/fs(f_domaine); // ! n_f pointe vers la face 1 donc vers l'exterieur de l'element, d'ou le -
        for (int d = 0 ; d < D ; d++) u_parallel(d) = pvit_elem(e_domaine, Nv*d+n) - u_orth*(-n_f(f_domaine,d))/fs(f_domaine) ; // ! n_f pointe vers la face 1 donc vers l'exterieur de l'element, d'ou le -
      }
      else // VDF
      {
        for (int d = 0; d <D ; d++) u_orth -= vit(nf_tot + e_domaine * D+d, n)*n_f(f_domaine,d)/fs(f_domaine); // ! n_f pointe vers la face 1 donc vers l'exterieur de l'element, d'ou le -
        for (int d = 0 ; d < D ; d++) u_parallel(d) = vit(nf_tot + e_domaine * D + d, n) - u_orth*(-n_f(f_domaine,d))/fs(f_domaine) ; // ! n_f pointe vers la face 1 donc vers l'exterieur de l'element, d'ou le -
      }
      double norm_u_parallel = std::sqrt(domaine.dot(&u_parallel(0), &u_parallel(0)));

      double y_loc = y(f_domaine,  n);
      double y_plus_loc = y_loc * u_tau(f_domaine, n)/ visc(e_domaine, n) ;
//      double fac_coeff_grad_ = 1 + fac_prod_k_ * std::tanh( (y_plus_loc/y_p_prod_k_)*(y_plus_loc/y_p_prod_k_) ) - fac_prod_k_grand_ * std::tanh( (y_plus_loc/y_p_prod_k_grand_)*(y_plus_loc/y_p_prod_k_grand_)*(y_plus_loc/y_p_prod_k_grand_) );
      if (y_plus_loc>1)
        {
          double coeff_frottement = (alp ? (*alp)(e_domaine, n) : 1) * rho(e_domaine, n) * u_tau(f_domaine, n)*u_tau(f_domaine, n)/norm_u_parallel; // f_tau = - alpha_k rho_k u_tau**2 n_par, coeff = u_tau**2 /u_par
          for (int d=0 ; d<D ; d++) d_(f, N*d + n) = u_parallel(d) * (1.- coeff_frottement*y_loc/mu(e_domaine, n)) ;
        }
      else for (int d=0 ; d<D ; d++) d_(f, N*d + n) = 0; // les phases non turbulentes sont non porteuses : dirichlet vitesse paroi = vitesse fluide
    }

  for (n=1 ; n<N ; n++)
   for (int d=0 ; d<D ; d++)
    for (int f =0 ; f < nf ; f++) d_(f, N*d + n) = 0; // les phases non turbulentes sont non porteuses : dirichlet vitesse paroi = vitesse fluide

  d_.echange_espace_virtuel();
}