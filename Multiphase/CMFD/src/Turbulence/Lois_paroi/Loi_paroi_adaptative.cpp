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
// File:        Loi_paroi_adaptative.cpp
// Directory:   $TRUST_ROOT/src/Turbulence/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Loi_paroi_adaptative.h>
#include <Correlation_base.h>
#include <Pb_Multiphase.h>
#include <Navier_Stokes_std.h>
#include <Zone_CoviMAC.h>
#include <DoubleTrav.h>
#include <Neumann_loi_paroi.h>
#include <Paroi_frottante_loi.h>
#include <Echange_impose_base.h>
#include <Cond_lim_base.h>
#include <math.h>
#include <Nom.h>
#include <Motcle.h>
#include <Champ_base.h>

Implemente_instanciable(Loi_paroi_adaptative, "Loi_paroi_adaptative", Loi_paroi_base);

Sortie& Loi_paroi_adaptative::printOn(Sortie& os) const
{
  return os;
}

Entree& Loi_paroi_adaptative::readOn(Entree& is)
{
  return is;
}

void Loi_paroi_adaptative::completer()
{
  const DoubleTab& vit = pb_.valeur().get_champ("vitesse").valeurs() ;
  Zone_CoviMAC& zone = ref_cast(Zone_CoviMAC, pb_.valeur().domaine_dis().zone_dis(0).valeur());
  int nf_tot = zone.nb_faces_tot();

  valeurs_loi_paroi_["y_plus"] = DoubleTab(0,1); // pour l'instant, turbulence dans seulement une phase
  valeurs_loi_paroi_["y"] = DoubleTab(0,1);
  valeurs_loi_paroi_["u_plus"] = DoubleTab(0,1);
  valeurs_loi_paroi_["dyp_u_plus"] = DoubleTab(0,1);
  valeurs_loi_paroi_["u_tau"] = DoubleTab(0,1);
  Faces_a_calculer_ = IntTab(nf_tot, 1);

  MD_Vector_tools::creer_tableau_distribue(vit.get_md_vector(), valeurs_loi_paroi_["y_plus"]);
  MD_Vector_tools::creer_tableau_distribue(vit.get_md_vector(), valeurs_loi_paroi_["y"]);
  MD_Vector_tools::creer_tableau_distribue(vit.get_md_vector(), valeurs_loi_paroi_["u_plus"]);
  MD_Vector_tools::creer_tableau_distribue(vit.get_md_vector(), valeurs_loi_paroi_["dyp_u_plus"]);
  MD_Vector_tools::creer_tableau_distribue(vit.get_md_vector(), valeurs_loi_paroi_["u_tau"]);

  for (int i = 0 ; i <pb_.valeur().nombre_d_equations() ; i++)
    for (int j = 0 ; j<pb_.valeur().equation(i).zone_Cl_dis()->nb_cond_lim(); j++)
      {
        Cond_lim cond_lim_loc = pb_.valeur().equation(i).zone_Cl_dis()->les_conditions_limites(j);
        if sub_type(Neumann_loi_paroi, cond_lim_loc.valeur())
          ref_cast(Neumann_loi_paroi, cond_lim_loc.valeur()).liste_faces_loi_paroi(Faces_a_calculer_);  // met des 1 si doit remplir la table
        else if sub_type(Paroi_frottante_loi, cond_lim_loc.valeur())
          ref_cast(Paroi_frottante_loi, cond_lim_loc.valeur()).liste_faces_loi_paroi(Faces_a_calculer_);  // met des 1 si doit remplir la table
        else if sub_type(Echange_impose_base, cond_lim_loc.valeur())
          ref_cast(Echange_impose_base, cond_lim_loc.valeur()).liste_faces_loi_paroi(Faces_a_calculer_);  // met des 1 si doit remplir la table
      }
}

void Loi_paroi_adaptative::mettre_a_jour(double temps)
{
  const DoubleTab& vit = pb_->get_champ("vitesse").valeurs();
  const DoubleTab& nu  = ref_cast(Navier_Stokes_std, pb_->equation(0)).diffusivite_pour_pas_de_temps().valeurs();
  calc_u_tau_y_plus(vit, nu);
  valeurs_loi_paroi_["y_plus"].echange_espace_virtuel();
  valeurs_loi_paroi_["y"].echange_espace_virtuel();
  valeurs_loi_paroi_["u_plus"].echange_espace_virtuel();
  valeurs_loi_paroi_["dyp_u_plus"].echange_espace_virtuel();
  valeurs_loi_paroi_["u_tau"].echange_espace_virtuel();
}


void Loi_paroi_adaptative::calc_u_tau_y_plus(const DoubleTab& vit, const DoubleTab& nu_visc)
{
  Zone_CoviMAC& zone = ref_cast(Zone_CoviMAC, pb_.valeur().domaine_dis().zone_dis(0).valeur());
  DoubleTab& u_t = valeurs_loi_paroi_["u_tau"], &y_p = valeurs_loi_paroi_["y_plus"], &y = valeurs_loi_paroi_["y"], &u_p = valeurs_loi_paroi_["u_plus"], &d_u_p = valeurs_loi_paroi_["dyp_u_plus"];
  const DoubleTab& n_f = zone.face_normales();
  const DoubleVect& fs = zone.face_surfaces();
  const IntTab& f_e = zone.face_voisins();

  int nf = zone.nb_faces(), nf_tot = zone.nb_faces_tot(), N = vit.line_size(), D = dimension;

  for (int f = 0 ; f < nf ; f ++) for (int n = 0 ; n < N ; n++) if (Faces_a_calculer_(f,0)==1)
        {
          if (f_e(f, 1) >= 0) Process::exit("Error in the definition of the boundary conditions for wall laws");
          int e = f_e(f,0);
          double u_orth = - zone.dot(&vit(nf_tot + e * D, n), &n_f(f,0))/fs(f); // ! n_f pointe vers la face 1 donc vers l'exterieur de l'element, d'ou le -
          DoubleTrav u_parallel(D);
          for (int d = 0 ; d < D ; d++) u_parallel(d) = vit(nf_tot + e * D + d, n) - u_orth*(-n_f(f,d))/fs(f) ; // ! n_f pointe vers la face 1 donc vers l'exterieur de l'element, d'ou le -
          if (zone.dot(&u_parallel(0), &n_f(f,0))/fs(f) > 1e-8) Process::exit("Error in the calculation of the parallel velocity for wall laws");
          double norm_u_parallel = std::sqrt(zone.dot(&u_parallel(0), &u_parallel(0)));
          double y_loc = zone.dist_face_elem0(f,  e);
          u_t(f, n) = calc_u_tau_loc(norm_u_parallel, nu_visc(e, n), y_loc);
          y_p(f, n) = y_loc*u_t(f, n)/nu_visc(e, n);
          y(f,n) = y_loc;
          if ( abs(norm_u_parallel/u_t(f, n) - u_plus_de_y_plus(y_p(f, n))) > 1e-4) Process::exit(Nom("No convergence on the Dichotomic algorithm ; u_t=") + Nom(u_t(f, n)) + Nom("u_parr=") + Nom(norm_u_parallel) +Nom("u_plus=") + Nom(u_plus_de_y_plus(y_p(f, n))));
          u_p(f, n) = norm_u_parallel/u_t(f, n);
          d_u_p(f,n)= deriv_u_plus_de_y_plus(y_p(f, n));
        }
}

double Loi_paroi_adaptative::calc_u_tau_loc(double u_par, double nu, double y)
{
  double eps = 1.e-6;
  int iter_max = 50;
  int n_iter = 0;

  double u_tau_1 = 1. ;
  while (to_zero(u_tau_1, u_par, nu, y) < 0) u_tau_1 *= 10 ;
  double u_tau_0 = 0.1 ;
  while (to_zero(u_tau_0, u_par, nu, y) > 0) u_tau_0 /= 10 ;
  if (u_tau_1 == 1.) u_tau_1 = u_tau_0*10 ;
  else if (u_tau_0 == 0.1) u_tau_0 = u_tau_1/10;

  eps *= u_tau_0;

  /* Implementing Dichotomic method */
  while ((abs(u_tau_1-u_tau_0)>eps) and (n_iter <= iter_max))
    {
      double f_0 = to_zero( (u_tau_0+u_tau_1)/2, u_par, nu, y);
      if (f_0 > 0) u_tau_1 -= (u_tau_1-u_tau_0)/2 ;
      else u_tau_0 += (u_tau_1-u_tau_0)/2;
      n_iter+=1;
    }
  if (n_iter > iter_max) Process::exit("Wall law has not converged !");

  return u_tau_0;
}

double Loi_paroi_adaptative::to_zero(double u_tau, double u_par, double nu, double y) // fonction for which we are looking for the root
{
  double y_p = y*u_tau/nu;
  return u_plus_de_y_plus(y_p) - u_par/u_tau;
}

double Loi_paroi_adaptative::u_plus_de_y_plus(double y_p) // Blended Reichardt model
{
  double reichardt = log(1+0.4*y_p)/von_karman_;
  reichardt += 7.8;
  reichardt += -7.8*exp(-y_p/11);
  reichardt += -7.8*y_p/11*exp(-y_p/3);

  double log_law = log(y_p)/von_karman_ + 5.1;

  double blending = tanh( y_p/27*y_p/27*y_p/27*y_p/27);

  return (1-blending)*reichardt + blending*log_law;
}

double Loi_paroi_adaptative::deriv_u_plus_de_y_plus(double y_p)
{
  double reichardt = log(1+0.4*y_p)/von_karman_ + 7.8 -7.8*exp(-y_p/11) -7.8*y_p/11*exp(-y_p/3);
  double log_law = log(y_p)/von_karman_ + 5.1;
  double blending = tanh( y_p/27*y_p/27*y_p/27*y_p/27);

  double d_reichardt = 0.4/(1+0.4*y_p)*1/von_karman_;
  d_reichardt += 7.8/11*exp(-y_p/11);
  d_reichardt += -7.8/11*exp(-y_p/3) + 7.8*y_p/33*exp(-y_p/3) ;
  double d_log_law = 1/(y_p+limiteur_y_p);
  double d_blending = 4/27*(y_p/27*y_p/27*y_p/27)*(1-blending*blending);

  return (1-blending)*d_reichardt - reichardt*d_blending + blending*d_log_law + d_blending*log_law ;
}

