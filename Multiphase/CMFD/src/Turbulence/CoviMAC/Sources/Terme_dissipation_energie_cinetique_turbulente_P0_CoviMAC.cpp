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
// File:        Travail_pression_CoviMAC.cpp
// Directory:   $TRUST_ROOT/src/Turbulence/CoviMAC/Sources
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Terme_dissipation_energie_cinetique_turbulente_P0_CoviMAC.h>
#include <Zone_CoviMAC.h>
#include <Champ_P0_CoviMAC.h>
#include <Equation_base.h>
#include <Pb_Multiphase.h>
#include <Milieu_composite.h>
#include <Array_tools.h>
#include <Matrix_tools.h>
#include <Op_Diff_Turbulent_CoviMAC_Face.h>
#include <Viscosite_turbulente_k_tau.h>

Implemente_instanciable(Terme_dissipation_energie_cinetique_turbulente_P0_CoviMAC,"Terme_dissipation_energie_cinetique_turbulente_P0_CoviMAC", Source_base);
// XD Terme_dissipation_energie_cinetique_turbulente source_base Terme_dissipation_energie_cinetique_turbulente 0 Source term which corresponds to the dissipation source term that appears in the turbulent kinetic energy equation

Sortie& Terme_dissipation_energie_cinetique_turbulente_P0_CoviMAC::printOn(Sortie& os) const
{
  return os;
}

Entree& Terme_dissipation_energie_cinetique_turbulente_P0_CoviMAC::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("beta_k", &beta_k);
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Terme_dissipation_energie_cinetique_turbulente_P0_CoviMAC::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const Zone_CoviMAC& zone = ref_cast(Zone_CoviMAC, equation().zone_dis().valeur());
  const int ne = zone.nb_elem(), ne_tot = zone.nb_elem_tot(), Nk = equation().inconnue().valeurs().line_size(); //, Nt = equation().probleme().get_champ("temperature").valeurs().line_size(), Ntau = equation().probleme().get_champ("tau").valeurs().line_size(), Np = equation().probleme().get_champ("pression").valeurs().line_size();
//  const int Nphases = equation().probleme().get_champ("alpha").valeurs().line_size();

//  for (auto &&i_m : matrices) if (i_m.first == "alpha" || i_m.first == "temperature" || i_m.first == "pression" || i_m.first == "tau")
//      {
//        Matrice_Morse mat;
//        IntTrav stencil(0, 2);
//        stencil.set_smart_resize(1);
//        if (i_m.first == "alpha") for (int e = 0; e < ne; e++) for(int n = 0; n < Nk ; n++) stencil.append_line(Ntau * e + n, Nphases * e + n);
//        if (i_m.first == "tau") for (int e = 0; e < ne; e++) for(int n = 0, mtau = 0; n < Nk ; n++, mtau += (Ntau > 1)) stencil.append_line(Ntau * e + n, Ntau * e + mtau);
//        if (i_m.first == "temperature") for (int e = 0; e < ne; e++) for(int n = 0; n < Nk ; n++) stencil.append_line(Ntau * e + n, Nt * e + n);
//        if (i_m.first == "pression") for (int e = 0; e < ne; e++) for(int n = 0, mp = 0; n < Nk ; n++, mp += (Np > 1)) stencil.append_line(Ntau * e + n, Np * e + mp); // attention reduction en pression
//        tableau_trier_retirer_doublons(stencil);
//        Matrix_tools::allocate_morse_matrix(ne_tot, ne_tot, stencil, mat);
//        i_m.second->nb_colonnes() ? *i_m.second += mat : *i_m.second = mat;
//      }

  for (auto &&n_m : matrices) if (n_m.first == "alpha" || n_m.first == "tau" || n_m.first == "temperature" || n_m.first == "pression")
      {
        Matrice_Morse& mat = *n_m.second, mat2;
        const DoubleTab& dep = equation().probleme().get_champ(n_m.first.c_str()).valeurs();
        int nc = dep.dimension_tot(0),
            M  = dep.line_size();
        IntTrav sten(0, 2);
        sten.set_smart_resize(1);
        if (n_m.first == "alpha" || n_m.first == "temperature")
          for (int e = 0; e < ne; e++) for (int n = 0; n < Nk; n++) sten.append_line(Nk * e + n, M * e + n);
        if (n_m.first == "pression" || n_m.first == "tau")
          for (int e = 0; e < ne; e++) for (int n = 0, m = 0; n < Nk; n++, m+=(M>1)) sten.append_line(Nk * e + n, M * e + m);
        Matrix_tools::allocate_morse_matrix(Nk * ne_tot, M * nc, sten, mat2);
        mat.nb_colonnes() ? mat += mat2 : mat = mat2;
      }
}

void Terme_dissipation_energie_cinetique_turbulente_P0_CoviMAC::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl)  const
{
  const Zone_CoviMAC& 					zone 			= ref_cast(Zone_CoviMAC, equation().zone_dis().valeur());
  const Champ_P0_CoviMAC& 				ch_k 			= ref_cast(Champ_P0_CoviMAC, equation().inconnue().valeur());		// Champ k
  const DoubleTab& 						k 				= ch_k.valeurs();
  const Champ_P0_CoviMAC& 				ch_tau 			= ref_cast(Champ_P0_CoviMAC,equation().probleme().get_champ("tau")); // Champ tau
  const DoubleTab& 						tau 			= ch_tau.valeurs() ; //== NULL ? ch_tau.valeurs() : NULL;
  const Champ_Inc_base& 				ch_alpha_rho_k 	= equation().champ_conserve();
  const DoubleTab& 						alpha_rho_k		= ch_alpha_rho_k.valeurs();
  const tabs_t& 						der_alpha_rho_k = ref_cast(Champ_Inc_base, ch_alpha_rho_k).derivees(); // dictionnaire des derivees
  const Navier_Stokes_std&     			eq_qdm 			= ref_cast(Navier_Stokes_std, equation().probleme().equation(0));
  const Op_Diff_Turbulent_CoviMAC_Face& op_diff 		= ref_cast(Op_Diff_Turbulent_CoviMAC_Face, eq_qdm.operateur(0).l_op_base());
  const Viscosite_turbulente_k_tau&   	visc_turb 		= ref_cast(Viscosite_turbulente_k_tau, op_diff.corr.valeur());
  const DoubleTab&                      nu 				= equation().probleme().get_champ("viscosite_cinematique").passe();

  const int Nk = k.line_size(), Ntau = tau.line_size(), Np = equation().probleme().get_champ("pression").valeurs().line_size(), Na = equation().probleme().get_champ("alpha").valeurs().line_size(), Nt = equation().probleme().get_champ("temperature").valeurs().line_size(), nb_elem = zone.nb_elem();

  Matrice_Morse *Ma = matrices.count("alpha") ? matrices.at("alpha") : NULL,
                 *Mk = matrices.count(ch_k.le_nom().getString()) ? matrices.at(ch_k.le_nom().getString()) : NULL,
                  *Mtau = matrices.count("tau") ? matrices.at("tau") : NULL,
                   *Mp = matrices.count("pression") ? matrices.at("pression") : NULL,
                    *Mt	= matrices.count("temperature") ? matrices.at("temperature") : NULL;

  double inv_tau = 0;

  for (int e = 0; e < nb_elem; e++) for (int mk = 0, mtau = 0, mp = 0; mk < Nk; mk++, mtau += (Ntau > 1), mp += (Np > 1))
      {
        inv_tau = k(e, mk) / max(k(e, mk) * tau(e, mtau), visc_turb.limiteur() * nu(e, mk));
        secmem(e, mk) -= beta_k * alpha_rho_k(e,mk) * inv_tau;
        if (Ma) 	(*Ma)(Nk * e + mk, Na * e + mk)   	  += beta_k * (der_alpha_rho_k.count("alpha") ? der_alpha_rho_k.at("alpha")(e,mk) : NULL ) * inv_tau;	// derivee en alpha
        if (Mk) 	(*Mk)(Nk * e + mk, Nk * e + mk)       += beta_k * (der_alpha_rho_k.count("k") ? der_alpha_rho_k.at("k")(e,mk) : NULL ) * inv_tau; // derivee en k
        if (Mtau) 	(*Mtau)(Nk * e + mk, Ntau * e + mtau) += beta_k * alpha_rho_k(e, mk) * (-pow(inv_tau,2)); // derivee en tau
        if (Mt) 	(*Mt)(Nk * e + mk, Nt * e + mk)       += beta_k * (der_alpha_rho_k.count("temperature") ? der_alpha_rho_k.at("temperature")(e, mk) : NULL ) * inv_tau;	// derivee par rapport a la temperature
        if (Mp) 	(*Mp)(Nk * e + mk, Np * e + mp)       += beta_k * (der_alpha_rho_k.count("pression") ? der_alpha_rho_k.at("pression")(e, mk) : NULL ) * inv_tau;		// derivee par rapport a la pression
      }
}

