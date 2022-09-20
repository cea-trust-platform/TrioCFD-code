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
// File:        Rupture_Yao_Morel.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Rupture_bulles_1groupe_PolyMAC_P0.h>
#include <Pb_Multiphase.h>
#include <Champ_Elem_PolyMAC_P0.h>
#include <Milieu_composite.h>
#include <Op_Diff_Turbulent_PolyMAC_P0_Face.h>
#include <Viscosite_turbulente_base.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <Rupture_bulles_1groupe_base.h>
#include <math.h>

Implemente_instanciable(Rupture_bulles_1groupe_PolyMAC_P0, "Rupture_bulles_1groupe_PolyMAC_P0", Source_base);

Sortie& Rupture_bulles_1groupe_PolyMAC_P0::printOn(Sortie& os) const
{
  return os;
}

Entree& Rupture_bulles_1groupe_PolyMAC_P0::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("beta_k", &beta_k_);
  param.lire_avec_accolades_depuis(is);


  const Pb_Multiphase *pbm = sub_type(Pb_Multiphase, equation().probleme()) ? &ref_cast(Pb_Multiphase, equation().probleme()) : NULL;

  if (!pbm || pbm->nb_phases() == 1) Process::exit(que_suis_je() + " : not needed for single-phase flow!");
  for (int n = 0; n < pbm->nb_phases(); n++) //recherche de n_l, n_g : phase {liquide,gaz}_continu en priorite
    if (pbm->nom_phase(n).debute_par("liquide") && (n_l < 0 || pbm->nom_phase(n).finit_par("continu")))  n_l = n;
  if (n_l < 0) Process::exit(que_suis_je() + " : liquid phase not found!");

  if (pbm->has_correlation("Rupture_bulles_1groupe")) correlation_ = pbm->get_correlation("Rupture_bulles_1groupe"); //correlation fournie par le bloc correlation
  else correlation_.typer_lire(*pbm, "Rupture_bulles_1groupe", is); //sinon -> on la lit

  return is;
}

void Rupture_bulles_1groupe_PolyMAC_P0::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const Zone_PolyMAC_P0& zone = ref_cast(Zone_PolyMAC_P0, equation().zone_dis().valeur());
  const int ne = zone.nb_elem(), ne_tot = zone.nb_elem_tot(), N = equation().inconnue().valeurs().line_size();

  for (auto &&n_m : matrices)
    if (n_m.first == "alpha" || n_m.first == "k" || n_m.first == "tau" || n_m.first == "omega" || n_m.first == "interfacial_area")
      {
        Matrice_Morse& mat = *n_m.second, mat2;
        const DoubleTab& dep = equation().probleme().get_champ(n_m.first.c_str()).valeurs();
        int nc = dep.dimension_tot(0),
            M  = dep.line_size();
        IntTrav sten(0, 2);
        sten.set_smart_resize(1);
        if (n_m.first == "alpha")
          for (int e = 0; e < ne; e++)
            for (int n = 0; n < N; n++)
              {
                sten.append_line(N * e + n, N * e + n_l);
                if (n != n_l) sten.append_line(N * e + n, N * e + n);
              }
        if (n_m.first == "k" || n_m.first == "tau" || n_m.first == "omega") // N <= M
          for (int e = 0; e < ne; e++)
            for (int n = 0; n < N; n++) sten.append_line(N * e + n, M * e +n_l);
        Matrix_tools::allocate_morse_matrix(N * ne_tot, M * nc, sten, mat2);
        mat.nb_colonnes() ? mat += mat2 : mat = mat2;
      }
}

void Rupture_bulles_1groupe_PolyMAC_P0::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Zone_PolyMAC_P0& zone = ref_cast(Zone_PolyMAC_P0, equation().zone_dis().valeur());
  const DoubleVect& pe = zone.porosite_elem(), &ve = zone.volumes();

  const DoubleTab& inco = equation().inconnue().valeurs(),
                   &d_bulles_p = equation().probleme().get_champ("diametre_bulles").passe(),
                    &alpha = ref_cast(Pb_Multiphase, equation().probleme()).eq_masse.inconnue().valeurs(),
                     &alpha_p = ref_cast(Pb_Multiphase, equation().probleme()).eq_masse.inconnue().passe(),
                      &press_p = ref_cast(Pb_Multiphase, equation().probleme()).eq_qdm.pression().passe(),
                       &temp_p  = ref_cast(Pb_Multiphase, equation().probleme()).eq_energie.inconnue().passe(),
                        &rho_p   = equation().milieu().masse_volumique().passe(),
                         &nu_p = equation().probleme().get_champ("viscosite_cinematique").passe(),
                          *tab_k_p = equation().probleme().has_champ("k") ? &equation().probleme().get_champ("k").passe() : NULL,
                           *tab_k = equation().probleme().has_champ("k") ? &equation().probleme().get_champ("k").valeurs() : NULL,
                            *tau = equation().probleme().has_champ("tau") ? &equation().probleme().get_champ("tau").valeurs() : NULL,
                             *omega = equation().probleme().has_champ("omega") ? &equation().probleme().get_champ("omega").valeurs() : NULL ;

  const Milieu_composite& milc = ref_cast(Milieu_composite, equation().milieu());
  int N = ref_cast(Pb_Multiphase, equation().probleme()).nb_phases(), Nk = (tab_k) ? (*tab_k).line_size() : -1, Np = equation().probleme().get_champ("pression").valeurs().line_size();

  std::string Type_diss = "other"; // omega, tau or other dissipation
  if (tau) Type_diss = "tau";
  else if (omega) Type_diss = "omega";

  DoubleTrav epsilon(alpha);
  const Op_Diff_Turbulent_PolyMAC_P0_Face& op_diff 		= ref_cast(Op_Diff_Turbulent_PolyMAC_P0_Face, equation().probleme().equation(0).operateur(0).l_op_base());
  const Viscosite_turbulente_base&   	visc_turb 		= ref_cast(Viscosite_turbulente_base, op_diff.correlation().valeur());
  visc_turb.eps(epsilon); // Epsilon is in the past
  double limiter = visc_turb.limiteur();

  Matrice_Morse  *Ma = matrices.count("alpha") ? matrices.at("alpha") : NULL,
                  *Mk = matrices.count("k") ? matrices.at("k") : NULL,
                   *Mtau = matrices.count("tau") ? matrices.at("tau") : NULL,
                    *Momega = matrices.count("omega") ? matrices.at("omega") : NULL,
                     *Mai = matrices.count("interfacial_area") ? matrices.at("interfacial_area") : NULL;

  int cR = (rho_p.dimension_tot(0) == 1), cM = (nu_p.dimension_tot(0) == 1), n, k, e;
  DoubleTrav a_l(N), p_l(N), T_l(N), rho_l(N), nu_l(N), sigma_l(N,N), dv(N, N), d_bulles_l(N), eps_l(Nk), k_l(Nk), coeff(N, N); //arguments pour coeff
  double dh;
  const Rupture_bulles_1groupe_base& correlation_rupt = ref_cast(Rupture_bulles_1groupe_base, correlation_.valeur());

  /* elements */
  for (e = 0; e < zone.nb_elem(); e++)
    {

      for (n = 0; n < N; n++) a_l(n)   = alpha_p(e, n);
      for (n = 0; n < N; n++) p_l(n)   = press_p(e, n * (Np > 1));
      for (n = 0; n < N; n++) T_l(n)   =  temp_p(e, n);
      for (n = 0; n < N; n++) rho_l(n) =   rho_p(!cR * e, n);
      for (n = 0; n < N; n++) nu_l(n)  =    nu_p(!cM * e, n);
      for (n = 0; n < N; n++)
        {
          for (k = 0; k < N; k++)
            if(milc.has_interface(n, k))
              {
                Interface_base& sat = milc.get_interface(n, k);
                sigma_l(n,k) = sat.sigma(temp_p(e,n), press_p(e,n * (Np > 1)));
              }
        }
      for (n = 0; n < N; n++) d_bulles_l(n) = d_bulles_p(e,n);
      for (n = 0; n <Nk; n++) eps_l(n) =epsilon(e, n) ;
      for (n = 0; n <Nk; n++) k_l(n)   = (tab_k_p)   ? (*tab_k_p)(e,0) : 0;

      correlation_rupt.coefficient(a_l, p_l, T_l, rho_l, nu_l, sigma_l, dh, dv, d_bulles_l, eps_l, k_l, coeff); // Explicit coeff


      for (k = 0 ; k<N ; k++)
        {

          double eps_valeurs ;

          if (Type_diss == "tau")        eps_valeurs = beta_k_ * ((*tab_k)(e, n_l)>1.e-8 ? (*tab_k)(e, n_l)*(*tab_k)(e, n_l)/ std::max((*tab_k)(e, n_l) * (*tau)(e, n_l), limiter * nu_p(e, n_l)) : 0 );
          else if (Type_diss == "omega") eps_valeurs = beta_k_ * ((*tab_k)(e, n_l)*(*omega)(e, n_l)) ;
          else eps_valeurs = epsilon(e, n_l);

          double fac = pe(e) * ve(e) * M_PI /(3*std::pow(6, 5./3.)) * coeff(k, n_l);

          secmem(e , k) += fac * std::pow(alpha(e, k), -2./3.) * alpha(e, n_l) * std::pow(inco(e, k), 5./3.) * std::pow(eps_valeurs, 1./3.) ;

          if (Ma)  (*Ma)(N * e + k , N * e + k) -= fac * -2./3.* std::pow(alpha(e, k), -5./3.) * alpha(e, n_l) * std::pow(inco(e, k), 5./3.) * std::pow(eps_valeurs, 1./3.) ;
          if (Ma)  (*Ma)(N * e + k , N * e +n_l)-= fac * std::pow(alpha(e, k), -2./3.) * std::pow(inco(e, k), 5./3.) * std::pow(eps_valeurs, 1./3.) ;
          if (Mai)(*Mai)(N * e + k , N * e + k) -= fac * std::pow(alpha(e, k), -2./3.) * alpha(e, n_l) * 5./3. * std::pow(inco(e, k), 2./3.) * std::pow(eps_valeurs, 1./3.) ;
          if (Type_diss == "tau")
            {
              if ((*tab_k)(e, n_l) * (*tau)(e, n_l) > limiter * nu_p(e, n_l)) // derivee en k ; depend de l'activation ou non du limiteur
                {
                  if (Mk) (*Mk)(N * e + k, Nk * e + n_l)   -=
                      fac * std::pow(alpha(e, k), -2./3.) * alpha(e, n_l) * std::pow(inco(e, k), 5./3.) * 1./3. * std::pow(beta_k_, 1./3.) * std::pow((*tab_k)(e, n_l),-2./3.) /std::pow((*tau)(e, n_l), 1./3.);
                  if (Mtau)(*Mtau)(N * e + k, Nk * e + n_l)-=
                      fac * std::pow(alpha(e, k), -2./3.) * alpha(e, n_l) * std::pow(inco(e, k), 5./3.) *-1./3. * std::pow(beta_k_, 1./3.) * std::pow((*tab_k)(e, n_l), 1./3.) *std::pow((*tau)(e, n_l),-4./3.);
                }
              else if (Mk) (*Mk)(N * e + k, Nk * e + n_l)   -=
                  fac * std::pow(alpha(e, k), -2./3.) * alpha(e, n_l) * std::pow(inco(e, k), 5./3.) * 1./3. * std::pow((*tab_k)(e, n_l),-1./3.) /std::pow(limiter * nu_p(e, n_l), 1./3.);
            }
          if (Type_diss == "omega")
            {
              if (Momega)(*Momega)(N * e + k , Nk * e + n_l) -=
                  fac * std::pow(alpha(e, k), -2./3.) * alpha(e, n_l) * std::pow(inco(e, k), 5./3.) * 1./3. * std::pow((*tab_k)(e, n_l), 1./3.) * std::pow((*omega)(e, n_l),-2./3.);
              if (Mk)        (*Mk)(N * e + k , Nk * e + n_l) -=
                  fac * std::pow(alpha(e, k), -2./3.) * alpha(e, n_l) * std::pow(inco(e, k), 5./3.) * 1./3. * std::pow((*tab_k)(e, n_l),-2./3.) * std::pow((*omega)(e, n_l), 1./3.);
            }

        }
    }
}
