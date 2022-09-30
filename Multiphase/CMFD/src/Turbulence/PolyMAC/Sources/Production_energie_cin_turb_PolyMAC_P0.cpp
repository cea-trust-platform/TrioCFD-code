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
// File:        Production_energie_cin_turb_PolyMAC_P0.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/PolyMAC_P0
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Production_energie_cin_turb_PolyMAC_P0.h>

#include <Zone_PolyMAC_P0.h>
#include <Champ_Elem_PolyMAC_P0.h>
#include <Matrix_tools.h>
#include <Probleme_base.h>
#include <grad_Champ_Face_PolyMAC_P0.h>
#include <Champ_Uniforme.h>
#include <Flux_interfacial_base.h>
#include <Milieu_composite.h>
#include <Operateur_Diff.h>
#include <Op_Diff_Turbulent_PolyMAC_P0_Face.h>
#include <Navier_Stokes_std.h>
#include <Viscosite_turbulente_base.h>
#include <TRUSTTab_parts.h>
#include <Loi_paroi_adaptative.h>
#include <Pb_Multiphase.h>
#include <Echelle_temporelle_turbulente.h>
#include <Taux_dissipation_turbulent.h>

Implemente_instanciable(Production_energie_cin_turb_PolyMAC_P0,"Production_energie_cin_turb_Elem_PolyMAC_P0", Source_base);
// XD Production_energie_cin_turb source_base Production_energie_cin_turb 0 Production source term for the TKE equation


Sortie& Production_energie_cin_turb_PolyMAC_P0::printOn(Sortie& os) const
{
  return os;
}

Entree& Production_energie_cin_turb_PolyMAC_P0::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("limiter_production", &limiter_prod_);
  param.lire_avec_accolades_depuis(is);

  return is;
}

void Production_energie_cin_turb_PolyMAC_P0::completer()
{
  const Zone_PolyMAC_P0&                      zone = ref_cast(Zone_PolyMAC_P0, equation().zone_dis().valeur());
  fac_.resize(zone.nb_elem_tot(), equation().inconnue().valeur().valeurs().line_size(), 2); // 3rd column : derivative
  fac_ = 1. ;

  for (int e = 0 ; e < fac_.dimension_tot(0) ; e++)
    for (int n = 0 ; n < fac_.dimension_tot(1) ; n++) fac_(e, n, 1) = 0. ; // 0 derivative by default

  if (ref_cast(Pb_Multiphase, equation().probleme()).has_correlation("Loi_paroi"))
    correlation_loi_paroi_ = ref_cast(Pb_Multiphase, equation().probleme()).get_correlation("Loi_paroi");
}

void Production_energie_cin_turb_PolyMAC_P0::mettre_a_jour(double tps)
{
  if (tps != mon_temps_) {calculer_fac() ; mon_temps_ = tps;}
}
void Production_energie_cin_turb_PolyMAC_P0::calculer_fac()
{
  if (correlation_loi_paroi_.non_nul())
    {
      const Zone_PolyMAC_P0&                      zone = ref_cast(Zone_PolyMAC_P0, equation().zone_dis().valeur());
      const IntTab& f_e = zone.face_voisins();

      Loi_paroi_adaptative& corr_loi_paroi = ref_cast(Loi_paroi_adaptative, correlation_loi_paroi_.valeur().valeur());
      DoubleTab& u_tau = corr_loi_paroi.get_tab("u_tau");
      const DoubleTab& tab_k = equation().probleme().get_champ("k").passe();

      int nf = zone.nb_faces(), N = tab_k.line_size() ;

      fac_ = 1.;
      for (int e = 0 ; e < fac_.dimension_tot(0) ; e++)
        for (int n = 0 ; n < fac_.dimension_tot(1) ; n++) fac_(e, n, 1) = 0. ; // 0 derivative by default

      for(int f = 0 ; f < nf ; f++)
        for (int n = 0 ; n < N ; n++)
          if (corr_loi_paroi.a_calculer(f))
            if (u_tau(f, n) > 1.e-6)
              {
                int e = f_e(f, 0) ;

                double kp = tab_k(e, n) / (u_tau(f, n)*u_tau(f, n));
                if ( (kp-3.33) > 1. ) fac_(e, n, 0) = 1. - limiter_prod_ ;
                else if ( (kp-3.33) < 0. ) fac_(e, n, 0) = 1. ;
                else fac_(e, n, 0) = 1. - limiter_prod_*(kp-3.33), fac_(e, n, 1) = - limiter_prod_  / (u_tau(f, n)*u_tau(f, n)) ;
              }
    }
}

void Production_energie_cin_turb_PolyMAC_P0::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
// empty : no derivative for the turbulent kinetic energy production to add in the blocks

  const Zone_PolyMAC_P0& zone = ref_cast(Zone_PolyMAC_P0, equation().zone_dis().valeur());
  const DoubleTab& k 	 = ref_cast(Champ_Elem_PolyMAC_P0, equation().inconnue().valeur()).valeurs();
  const int ne = zone.nb_elem(), ne_tot = zone.nb_elem_tot(), Nk = k.line_size();

  std::string Type_diss = ""; // omega or tau dissipation
  for (int i = 0 ; i < equation().probleme().nombre_d_equations() ; i++)
    {
      if sub_type(Echelle_temporelle_turbulente, equation().probleme().equation(i)) Type_diss = "tau";
      else if sub_type(Taux_dissipation_turbulent, equation().probleme().equation(i)) Type_diss = "omega";
    }
  if (Type_diss == "") return;

  assert(Nk == 1); // si plus d'une phase turbulente, il vaut mieux iterer sur les id_composites des phases turbulentes modelisees par un modele k-tau
  if (Type_diss == "tau") assert(ref_cast(Champ_Elem_PolyMAC_P0,equation().probleme().get_champ("tau")).valeurs().line_size() == 1);
  if (Type_diss == "omega") assert(ref_cast(Champ_Elem_PolyMAC_P0,equation().probleme().get_champ("omega")).valeurs().line_size() == 1);

  for (auto &&n_m : matrices)
	  if (n_m.first == "alpha" || n_m.first == "tau" || n_m.first == "omega" || n_m.first == "temperature" || n_m.first == "pression")
	  {
		  Matrice_Morse& mat = *n_m.second, mat2;
		  const DoubleTab& dep = equation().probleme().get_champ(n_m.first.c_str()).valeurs();
		  int nc = dep.dimension_tot(0),
		      M  = dep.line_size();
		  IntTrav sten(0, 2);
		  sten.set_smart_resize(1);
		  if (n_m.first == "alpha" || n_m.first == "temperature" || n_m.first == "tau"|| n_m.first == "omega")
			  for (int e = 0; e < ne; e++)
				  for (int n = 0; n < Nk; n++)
					  if (n < M) sten.append_line(Nk * e + n, M * e + n);
		  if (n_m.first == "pression" )
			  for (int e = 0; e < ne; e++)
				  for (int n = 0, m = 0; n < Nk; n++, m+=(M>1)) sten.append_line(Nk * e + n, M * e + m);
		  Matrix_tools::allocate_morse_matrix(Nk * ne_tot, M * nc, sten, mat2);
		  mat.nb_colonnes() ? mat += mat2 : mat = mat2;
	  }
}

void Production_energie_cin_turb_PolyMAC_P0::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
	const Zone_PolyMAC_P0&                   zone = ref_cast(Zone_PolyMAC_P0, equation().zone_dis().valeur());
	const Probleme_base&                       pb = ref_cast(Probleme_base, equation().probleme());
	const Navier_Stokes_std&               eq_qdm = ref_cast(Navier_Stokes_std, pb.equation(0));
	const DoubleTab&                     tab_grad = pb.get_champ("gradient_vitesse").passe();
	const Op_Diff_Turbulent_PolyMAC_P0_Face& Op_diff = ref_cast(Op_Diff_Turbulent_PolyMAC_P0_Face, eq_qdm.operateur(0).l_op_base());
	const Viscosite_turbulente_base&    visc_turb = ref_cast(Viscosite_turbulente_base, Op_diff.correlation().valeur());
	const DoubleVect& pe = equation().milieu().porosite_elem(), &ve = zone.volumes();

  std::string Type_diss = ""; // omega or tau dissipation
  if sub_type(Echelle_temporelle_turbulente, equation()) Type_diss = "tau";
  else if sub_type(Taux_dissipation_turbulent, equation()) Type_diss = "omega";
  double limiter_ = visc_turb.limiteur();
  double nut_l, fac;

  const DoubleTab& tab_rho = equation().probleme().get_champ("masse_volumique").passe(),
                   &tab_alp = equation().probleme().get_champ("alpha").passe(),
                    &nu =  equation().probleme().get_champ("viscosite_cinematique").passe(),
                     &k = equation().probleme().get_champ("k").valeurs(),
                      *diss = equation().probleme().has_champ(Type_diss) ? &equation().probleme().get_champ(Type_diss).valeurs() : nullptr ;

  const Champ_base&   ch_alpha_rho = sub_type(Pb_Multiphase,equation().probleme()) ? ref_cast(Pb_Multiphase,equation().probleme()).eq_masse.champ_conserve() : equation().milieu().masse_volumique().valeur();
  const DoubleTab&       alpha_rho = ch_alpha_rho.valeurs();
  const tabs_t&      der_alpha_rho = ref_cast(Champ_Inc_base, ch_alpha_rho).derivees(); // dictionnaire des derivees

  int Nph = pb.get_champ("vitesse").valeurs().dimension(1), nb_elem = zone.nb_elem(), D = dimension, nf_tot = zone.nb_faces_tot() ;
  int N = equation().inconnue()->valeurs().line_size(),
      Np = equation().probleme().get_champ("pression").valeurs().line_size(),
      Nt = equation().probleme().get_champ("temperature").valeurs().line_size(),
      Na = sub_type(Pb_Multiphase,equation().probleme()) ? equation().probleme().get_champ("alpha").valeurs().line_size() : 1;
  int e, n, mp;

//  DoubleTrav Rij(0, Nph, D, D);
// MD_Vector_tools::creer_tableau_distribue(eq_qdm.pression()->valeurs().get_md_vector(), Rij); //Necessary to compare size in reynolds_stress()
//  visc_turb.reynolds_stress(Rij);

  if (Type_diss == "")
    {
      DoubleTrav nut(0, Nph);
      MD_Vector_tools::creer_tableau_distribue(eq_qdm.pression()->valeurs().get_md_vector(), nut); //Necessary to compare size in reynolds_stress()
      visc_turb.eddy_viscosity(nut);

      Matrice_Morse *mat = matrices.count("k") ? matrices.at("k") : nullptr;

      for( e = 0 ; e < nb_elem ; e++)
        for( n = 0; n<N ; n++)
          {
            double secmem_en = 0.;
            for (int d_U = 0; d_U < D; d_U++)
              for (int d_X = 0; d_X < D; d_X++)
                secmem_en += ( tab_grad(nf_tot + d_U + e * D , D * n + d_X) + tab_grad(nf_tot + d_X + e * D , D * n + d_U) ) * tab_grad(nf_tot + d_X + e * D , D * n + d_U) ;
//            secmem_en += Rij(e, n, d_X, d_U) * tab_grad(nf_tot + d_X + e * D , D * n + d_U) ;
            secmem_en *= pe(e) * ve(e) * tab_alp(e, n) * tab_rho(e, n) * nut(e, n) ;
//        secmem_en *= (-1) * pe(e) * ve(e) * tab_alp(e, n) * tab_rho(e, n) ;

//        secmem(e, n) += fac_(e, n, 0) * std::max(secmem_en, 0.);
            secmem(e, n) += std::max(secmem_en, 0.);

            if (mat) (*mat)(N * e + n, N * e + n) -= 0.;//fac_(e, n, 1) * std::max(secmem_en, 0.);
          }
    }

  else
    {
      for( e = 0 ; e < nb_elem ; e++)
        for(n = 0, mp = 0; n<N ; n++, mp += (Np > 1))
          {
            double grad_grad = 0.;
            for (int d_U = 0; d_U < D; d_U++)
              for (int d_X = 0; d_X < D; d_X++)
                grad_grad += ( tab_grad(nf_tot + d_U + e * D , D * n + d_X) + tab_grad(nf_tot + d_X + e * D , D * n + d_U) ) * tab_grad(nf_tot + d_X + e * D , D * n + d_U) ;
            fac = std::max(grad_grad, 0.) * pe(e) * ve(e) ;
            if (Type_diss == "tau")        nut_l =                         std::max(k(e, n) * (*diss)(e, n), limiter_ * nu(e, n)) ;
            else if (Type_diss == "omega") nut_l = ( ((*diss)(e,n) > 0.) ? std::max(k(e, n) / (*diss)(e, n), limiter_ * nu(e, n)) : limiter_ * nu(e, n) );
            else Process::exit(que_suis_je() + " : ajouter_blocs : probleme !!!") ;

            secmem(e, n) += fac * alpha_rho(e, n) * nut_l ;
            for (auto &&i_m : matrices)
              {
                Matrice_Morse& mat = *i_m.second;
                if (i_m.first == "alpha") 	    mat(N * e + n, Na * e + n) -= fac * nut_l * (der_alpha_rho.count("alpha") ? der_alpha_rho.at("alpha")(e, n) : 0 );			      // derivee par rapport au taux de vide
                if (i_m.first == "temperature") mat(N * e + n, Nt * e + n) -= fac * nut_l * (der_alpha_rho.count("temperature") ? der_alpha_rho.at("temperature")(e, n) : 0 );// derivee par rapport a la temperature
                if (i_m.first == "pression")    mat(N * e + n, Np * e + mp)-= fac * nut_l * (der_alpha_rho.count("pression") ? der_alpha_rho.at("pression")(e, n) : 0 );		  // derivee par rapport a la pression
              }
            if (Type_diss == "tau")
              for (auto &&i_m : matrices)
                {
                  Matrice_Morse& mat = *i_m.second;
                  if (i_m.first == "k")           mat(N * e + n,  N * e + n) -= fac * alpha_rho(e, n) * (                       k(e, n)*(*diss)(e, n)>limiter_*nu(e, n) ? (*diss)(e,n) : 0. ) ;
                  if (i_m.first == "tau")         mat(N * e + n,  N * e + n) -= fac * alpha_rho(e, n) * (                       k(e, n)*(*diss)(e, n)>limiter_*nu(e, n) ?       k(e,n) : 0. ) ;
                }
            if (Type_diss == "omega")
              for (auto &&i_m : matrices)
                {
                  Matrice_Morse& mat = *i_m.second;
                  if (i_m.first == "k")           mat(N * e + n,  N * e + n) -= fac * alpha_rho(e, n) * ( ((*diss)(e,n) > 0.) ?(k(e, n)/(*diss)(e, n)>limiter_*nu(e, n) ? 1/(*diss)(e, n)                    : 0.) : 0. ) ;
                  if (i_m.first == "omega")       mat(N * e + n,  N * e + n) -= fac * alpha_rho(e, n) * ( ((*diss)(e,n) > 0.) ?(k(e, n)/(*diss)(e, n)>limiter_*nu(e, n) ?-k(e,n)/((*diss)(e, n)*(*diss)(e, n)) : 0.) : 0. ) ;
                }
          }
    }
}
