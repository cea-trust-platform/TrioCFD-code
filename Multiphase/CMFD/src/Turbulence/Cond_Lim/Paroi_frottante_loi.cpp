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
// File:        Paroi_frottante_loi.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Incompressible/Cond_Lim
// Version:     /main/28
//
//////////////////////////////////////////////////////////////////////////////

#include <Paroi_frottante_loi.h>
#include <Motcle.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <Convection_Diffusion_Concentration.h>
#include <Loi_paroi_adaptative.h>
#include <Frontiere_dis_base.h>
#include <Frontiere.h>
#include <Pb_Multiphase.h>
#include <Navier_Stokes_std.h>
#include <Domaine_Poly_base.h>
#include <Op_Diff_PolyMAC_base.h>
#include <Op_Diff_PolyMAC_P0_base.h>
#include <Op_Diff_Turbulent_PolyMAC_P0_Face.h>
#include <Viscosite_turbulente_k_omega.h>
#include <Viscosite_turbulente_k_tau.h>

#include <math.h>

Implemente_instanciable(Paroi_frottante_loi,"Paroi_frottante_loi", Frottement_global_impose);
// XD Paroi_frottante_loi condlim_base Paroi_frottante_loi 1 Adaptive wall-law boundary condition for velocity

Sortie& Paroi_frottante_loi::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

Entree& Paroi_frottante_loi::readOn(Entree& s )
{
  if (app_domains.size() == 0) app_domains = { Motcle("turbulence") };

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

  return s;
}

void Paroi_frottante_loi::completer()
{
  if (!sub_type(Op_Diff_Turbulent_PolyMAC_P0_Face, domaine_Cl_dis().equation().operateur(0).l_op_base())) Process::exit(que_suis_je() + " : diffusion operator must be turbulent !");
  if sub_type(Viscosite_turbulente_k_tau, ref_cast(Op_Diff_Turbulent_PolyMAC_P0_Face, domaine_Cl_dis().equation().operateur(0).l_op_base()).correlation())
    {
      if (fac_prod_k_<-1.e7) fac_prod_k_ = 1.2;
      if (y_p_prod_k_<-1.e7) y_p_prod_k_ =  4.;
      if (fac_prod_k_grand_<-1.e7) fac_prod_k_grand_ = .2;
      if (y_p_prod_k_grand_<-1.e7) y_p_prod_k_grand_ = 150.;
    }
  else if sub_type(Viscosite_turbulente_k_omega, ref_cast(Op_Diff_Turbulent_PolyMAC_P0_Face, domaine_Cl_dis().equation().operateur(0).l_op_base()).correlation())
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


void Paroi_frottante_loi::liste_faces_loi_paroi(IntTab& tab)
{
  int nf = la_frontiere_dis.valeur().frontiere().nb_faces(), f1 = la_frontiere_dis.valeur().frontiere().num_premiere_face();
  int N = tab.line_size();

  for (int f =0 ; f < nf ; f++)
    for (int n = 0 ; n<N ; n++)
      tab(f + f1, n) |= 1;
}

double Paroi_frottante_loi::coefficient_frottement(int i) const
{
  return valeurs_coeff_(i,0);
}

double Paroi_frottante_loi::coefficient_frottement(int i,int j) const
{
  return valeurs_coeff_(i,j);
}

double Paroi_frottante_loi::coefficient_frottement_grad(int i) const
{
  return valeurs_coeff_grad_(i,0);
}

double Paroi_frottante_loi::coefficient_frottement_grad(int i,int j) const
{
  return valeurs_coeff_grad_(i,j);
}

int Paroi_frottante_loi::initialiser(double temps)
{
  valeurs_coeff_.resize(0, ref_cast(Pb_Multiphase, domaine_Cl_dis().equation().probleme()).nb_phases());
  la_frontiere_dis.valeur().frontiere().creer_tableau_faces(valeurs_coeff_);
  valeurs_coeff_ = 0 ;

  valeurs_coeff_grad_.resize(0, ref_cast(Pb_Multiphase, domaine_Cl_dis().equation().probleme()).nb_phases());
  la_frontiere_dis.valeur().frontiere().creer_tableau_faces(valeurs_coeff_grad_);
  valeurs_coeff_grad_ = 0 ;

  correlation_loi_paroi_ = ref_cast(Pb_Multiphase, domaine_Cl_dis().equation().probleme()).get_correlation("Loi_paroi");
  return 1;
}

void Paroi_frottante_loi::mettre_a_jour(double tps)
{
  if (mon_temps < tps) {me_calculer() ; mon_temps=tps;}
}

void Paroi_frottante_loi::me_calculer()
{
  Loi_paroi_adaptative& corr_loi_paroi = ref_cast(Loi_paroi_adaptative, correlation_loi_paroi_.valeur().valeur());
  const Domaine_Poly_base& domaine = ref_cast(Domaine_Poly_base, domaine_Cl_dis().equation().probleme().domaine_dis().valeur());

  const DoubleTab& u_tau = corr_loi_paroi.get_tab("u_tau"); // y_p est numerote selon les faces du domaine
  const DoubleTab& visc  = ref_cast(Navier_Stokes_std, domaine_Cl_dis().equation().probleme().equation(0)).diffusivite_pour_pas_de_temps().valeurs();
  const DoubleTab& vit   = domaine_Cl_dis().equation().probleme().get_champ("vitesse").valeurs(),
                   &rho = domaine_Cl_dis().equation().probleme().get_champ("masse_volumique").valeurs(),
                    *alp = sub_type(Pb_Multiphase, domaine_Cl_dis().equation().probleme()) ? &domaine_Cl_dis().equation().probleme().get_champ("alpha").valeurs() : NULL,
                     &mu = sub_type(Op_Diff_PolyMAC_base   , domaine_Cl_dis().equation().operateur(0).l_op_base()) ? ref_cast(Op_Diff_PolyMAC_base, domaine_Cl_dis().equation().operateur(0).l_op_base()).nu() :
                           ref_cast(Op_Diff_PolyMAC_P0_base, domaine_Cl_dis().equation().operateur(0).l_op_base()).nu();

  int nf = la_frontiere_dis.valeur().frontiere().nb_faces(), nf_tot = domaine.nb_faces_tot(), f1 = la_frontiere_dis.valeur().frontiere().num_premiere_face();
  int N = domaine_Cl_dis().equation().inconnue().valeurs().line_size(), D = dimension;

  const DoubleTab& n_f = domaine.face_normales();
  const DoubleVect& fs = domaine.face_surfaces();
  const IntTab& f_e = domaine.face_voisins();

  int n = 0 ; // la phase turbulente frotte
  for (int f =0 ; f < nf ; f++)
    {
      int f_domaine = f + f1; // number of the face in the domaine
      if (f_e(f_domaine, 1) >= 0) Process::exit("Paroi_frottante_loi : error in the definition of the boundary faces for wall laws");
      int e = f_e(f_domaine,0);

      double u_orth = 0 ;
      for (int d = 0; d <D ; d++) u_orth -= vit(nf_tot + e * D+d, n)*n_f(f_domaine,d)/fs(f_domaine); // ! n_f pointe vers la face 1 donc vers l'exterieur de l'element, d'ou le -

      DoubleTrav u_parallel(D);
      for (int d = 0 ; d < D ; d++) u_parallel(d) = vit(nf_tot + e * D + d, n) + n_f(f_domaine,d) * u_orth/fs(f_domaine) ; // ! + car on a mis - au-dessus
      double residu = 0 ;
      for (int d = 0; d <D ; d++) residu += u_parallel(d)*n_f(f_domaine,d)/fs(f_domaine);
      if (residu > 1e-8) Process::exit("Paroi_frottante_loi : Error in the calculation of the parallel velocity for wall laws");
      double norm_u_parallel = std::sqrt(domaine.dot(&u_parallel(0), &u_parallel(0)));

      double y_loc = domaine.dist_face_elem0(f_domaine,  e);
      double y_plus_loc = y_loc * u_tau(f_domaine, n)/ visc(e, n) ;
      double fac_coeff_grad_ = 1 + fac_prod_k_ * std::tanh( (y_plus_loc/y_p_prod_k_)*(y_plus_loc/y_p_prod_k_) ) - fac_prod_k_grand_ * std::tanh( (y_plus_loc/y_p_prod_k_grand_)*(y_plus_loc/y_p_prod_k_grand_)*(y_plus_loc/y_p_prod_k_grand_) );
      if (y_plus_loc>1)
        {
          valeurs_coeff_(f, n) = (alp ? (*alp)(e, n) : 1) * rho(e, n) * u_tau(f_domaine, n)*u_tau(f_domaine, n)/norm_u_parallel; // f_tau = - alpha_k rho_k u_tau**2 n_par, coeff = u_tau**2 /u_par
          valeurs_coeff_grad_(f, n) = fac_coeff_grad_ * 1/mu(e,n) * (alp ? (*alp)(e, n) : 1) * rho(e, n) * u_tau(f_domaine, n)*u_tau(f_domaine, n)/norm_u_parallel; // f_tau = - alpha_k rho_k u_tau**2 n_par, coeff = u_tau**2 /u_par
//          valeurs_coeff_grad_(f, n) = (alp ? (*alp)(e, n) : 1) * 1./y_loc ; // dirichlet for calculation of gradient
        }
      else
        {
          valeurs_coeff_(f, n) = (alp ? (*alp)(e, n) : 1) * rho(e, n) * visc(e, n)/y_loc; // viscous case : if u_tau is small
          if (norm_u_parallel > 1.e-6) valeurs_coeff_grad_(f, n) = fac_coeff_grad_ * (alp ? (*alp)(e, n) : 1) * 1./y_loc ; // dirichlet for calculation of gradient
          else                         valeurs_coeff_grad_(f, n) = fac_coeff_grad_ * (alp ? (*alp)(e, n) : 1) * 1/mu(e,n) * rho(e, n) * visc(e, n)/y_loc; // viscous case : if u_tau is small
        }
    }

  for (n=1 ; n<N ; n++)
    for (int f =0 ; f < nf ; f++)
      {
        valeurs_coeff_(f, n) = 0; // les phases non turbulentes sont non porteuses : pas de contact paroi => des symmetries
        valeurs_coeff_grad_(f, n) = 0; // les phases non turbulentes sont non porteuses : pas de contact paroi => des symmetries

        /*
                // Test : faire frotter un peu quand meme
                int f_domaine = f + f1; // number of the face in the domaine
                int e = f_e(f_domaine,0);
                valeurs_coeff_(f, n)      = valeurs_coeff_(f, n)      * mu(e,n) / mu(e,0) ;
                valeurs_coeff_grad_(f, n) = valeurs_coeff_grad_(f, n) ;
        */
      }

  valeurs_coeff_.echange_espace_virtuel();
  valeurs_coeff_grad_.echange_espace_virtuel();
}

