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
// File:        Coalescence_Yao_Morel.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Coalescence_Yao_Morel_PolyMAC_V2.h>
#include <Pb_Multiphase.h>
#include <Champ_P0_PolyMAC_V2.h>
#include <Milieu_composite.h>
#include <Op_Diff_Turbulent_PolyMAC_V2_Face.h>
#include <Viscosite_turbulente_base.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <math.h>

Implemente_instanciable(Coalescence_Yao_Morel, "Coalescence_Yao_Morel_P0_PolyMAC_V2", Source_base);

Sortie& Coalescence_Yao_Morel::printOn(Sortie& os) const
{
  return os;
}

Entree& Coalescence_Yao_Morel::readOn(Entree& is)
{
  return is;
}

void Coalescence_Yao_Morel::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const Zone_PolyMAC_V2& zone = ref_cast(Zone_PolyMAC_V2, equation().zone_dis().valeur());
  const int ne = zone.nb_elem(), ne_tot = zone.nb_elem_tot(), N = equation().inconnue().valeurs().line_size();

  for (auto &&n_m : matrices) if (n_m.first == "alpha" || n_m.first == "k" || n_m.first == "tau" || n_m.first == "omega")
      {
        Matrice_Morse& mat = *n_m.second, mat2;
        const DoubleTab& dep = equation().probleme().get_champ(n_m.first.c_str()).valeurs();
        int nc = dep.dimension_tot(0),
            M  = dep.line_size();
        IntTrav sten(0, 2);
        sten.set_smart_resize(1);
        if (n_m.first == "alpha")
          for (int e = 0; e < ne; e++) for (int n = 0; n < N; n++) sten.append_line(N * e + n, N * e + n);
        if (n_m.first == "k" || n_m.first == "tau" || n_m.first == "omega") // N <= M
          for (int e = 0; e < ne; e++) for (int n = 0; n < N; n++) for (int n_l = 0; n_l < M; n_l++) sten.append_line(N * e + n, M * e +n_l);
        Matrix_tools::allocate_morse_matrix(N * ne_tot, M * nc, sten, mat2);
        mat.nb_colonnes() ? mat += mat2 : mat = mat2;
      }
}

void Coalescence_Yao_Morel::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Zone_PolyMAC_V2& zone = ref_cast(Zone_PolyMAC_V2, equation().zone_dis().valeur());
  const DoubleVect& pe = zone.porosite_elem(), &ve = zone.volumes();

  const DoubleTab& inco = equation().inconnue().valeurs(),
                   &d_bulles_p = equation().probleme().get_champ("diametre_bulles").passe(),
                    &alpha = ref_cast(Pb_Multiphase, equation().probleme()).eq_masse.inconnue().valeurs(),
                     &alpha_p = ref_cast(Pb_Multiphase, equation().probleme()).eq_masse.inconnue().passe(),
                      &press = ref_cast(Pb_Multiphase, equation().probleme()).eq_qdm.pression().passe(),
                       &temp  = ref_cast(Pb_Multiphase, equation().probleme()).eq_energie.inconnue().passe(),
                        &rho   = equation().milieu().masse_volumique().passe(),
                         &nu = equation().probleme().get_champ("viscosite_cinematique").passe(),
                          *tab_k = equation().probleme().has_champ("k") ? &equation().probleme().get_champ("k").valeurs() : NULL,
                           *tau = equation().probleme().has_champ("tau") ? &equation().probleme().get_champ("tau").valeurs() : NULL,
                            *omega = equation().probleme().has_champ("omega") ? &equation().probleme().get_champ("omega").valeurs() : NULL ;

  const Milieu_composite& milc = ref_cast(Milieu_composite, equation().milieu());
  int N = ref_cast(Pb_Multiphase, equation().probleme()).nb_phases(), Nk = (tab_k) ? (*tab_k).line_size() : -1, Np = equation().probleme().get_champ("pression").valeurs().line_size();

  std::string Type_diss = "other"; // omega, tau or other dissipation
  if (tau) Type_diss = "tau";
  else if (omega) Type_diss = "omega";

  DoubleTrav epsilon(alpha);
  const Op_Diff_Turbulent_PolyMAC_V2_Face& op_diff 		= ref_cast(Op_Diff_Turbulent_PolyMAC_V2_Face, equation().probleme().equation(0).operateur(0).l_op_base());
  const Viscosite_turbulente_base&   	visc_turb 		= ref_cast(Viscosite_turbulente_base, op_diff.correlation().valeur());
  visc_turb.eps(epsilon);
  double limiter = visc_turb.limiteur();

  Matrice_Morse  *Ma = matrices.count("alpha") ? matrices.at("alpha") : NULL,
                  *Mk = matrices.count("k") ? matrices.at("k") : NULL,
                   *Mtau = matrices.count("tau") ? matrices.at("tau") : NULL,
                    *Momega = matrices.count("omega") ? matrices.at("omega") : NULL,
                     *Mai = matrices.count("interfacial_area") ? matrices.at("interfacial_area") : NULL;

  /* elements */
  int n_l = 0 ; // phase porteuse
  for (int e = 0; e < zone.nb_elem(); e++)
    for (int k = 0 ; k<N ; k++) if (k != n_l) //phase gazeuse
        if (alpha(e, k) > 1.e-6)
          {
            Interface_base& interface = milc.get_interface(n_l, k) ;
            double sigma = interface.sigma_(temp(e,n_l),press(e,n_l * (Np > 1)));
            double We = 2 * rho(e, n_l) * std::pow(epsilon(e, n_l)*d_bulles_p(e,k), 2./3.) * d_bulles_p(e,k) / sigma ;
            double g_alpha = (alpha_max_1_3 - std::pow(alpha_p(e,k), 1./3.) ) / alpha_max_1_3 ;
            double eps_loc ;

            if (Type_diss == "tau") eps_loc = (*tab_k)(e, n_l)>1.e-8 ? (*tab_k)(e, n_l)*(*tab_k)(e, n_l)/ std::max((*tab_k)(e, n_l) * (*tau)(e, n_l), limiter * nu(e, n_l)) : 0 ;
            else if (Type_diss == "omega") eps_loc = (*tab_k)(e, n_l)*(*omega)(e, n_l) ;
            else eps_loc = epsilon(e, n_l);

            double fac = pe(e) * ve(e) * M_PI /(3*std::pow(6, 5./3.)) * Kc1 *1/(g_alpha+Kc2 *alpha_p(e, k)*std::sqrt(We/We_cr))*std::exp(-Kc3*std::sqrt(We/We_cr));

            secmem(e , k) += - fac * std::pow(alpha(e, k), 1./3.) * std::pow(inco(e, k), 5./3.) * std::pow(eps_loc, 1./3.) ;

            if (Ma)  (*Ma)(N * e + k , N * e + k) -= - fac * 1./3. * std::pow(alpha(e, k),-2./3.) * std::pow(inco(e, k), 5./3.) * std::pow(eps_loc, 1./3.) ;
            if (Mai)(*Mai)(N * e + k , N * e + k) -= - fac * std::pow(alpha(e, k), 1./3.) * 5./3. * std::pow(inco(e, k), 2./3.) * std::pow(eps_loc, 1./3.) ;
            if (Type_diss == "tau")
              {
                if ((*tab_k)(e, k) * (*tau)(e, k) > limiter * nu(e, k)) // derivee en k ; depend de l'activation ou non du limiteur
                  {
                    if (Mk) (*Mk)(N * e + k, Nk * e + n_l)   -=
                        - fac * std::pow(alpha(e, k), 1./3.) * std::pow(inco(e, k), 5./3.) * 1./3. * std::pow((*tab_k)(e, n_l),-2./3.) /std::pow((*tau)(e, n_l), 1./3.);
                    if (Mtau)(*Mtau)(N * e + k, Nk * e + n_l)-=
                        - fac * std::pow(alpha(e, k), 1./3.) * std::pow(inco(e, k), 5./3.) *-1./3. * std::pow((*tab_k)(e, n_l), 1./3.) *std::pow((*tau)(e, n_l),-4./3.);
                  }
                else if (Mk) (*Mk)(N * e + k, Nk * e + n_l)   -=
                    - fac * std::pow(alpha(e, k), 1./3.) * std::pow(inco(e, k), 5./3.) * 1./3. * std::pow((*tab_k)(e, n_l),-1./3.) /std::pow(limiter * nu(e, n_l), 1./3.);
              }
            if (Type_diss == "omega")
              {
                if (Momega)(*Momega)(N * e + k , Nk * e + n_l) -= - fac * std::pow(alpha(e, k), 1./3.) * std::pow(inco(e, k), 5./3.) * 1./3. * std::pow((*tab_k)(e, n_l), 1./3.) * std::pow((*omega)(e, n_l),-2./3.);
                if (Mk)        (*Mk)(N * e + k , Nk * e + n_l) -= - fac * std::pow(alpha(e, k), 1./3.) * std::pow(inco(e, k), 5./3.) * 1./3. * std::pow((*tab_k)(e, n_l),-2./3.) * std::pow((*omega)(e, n_l), 1./3.);
              }
          }
}