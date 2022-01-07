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

#include <Terme_dissipation_echelle_temporelle_turbulente_P0_CoviMAC.h>
#include <Zone_CoviMAC.h>
#include <Champ_P0_CoviMAC.h>
#include <Equation_base.h>
#include <Pb_Multiphase.h>
#include <Milieu_composite.h>
#include <Array_tools.h>
#include <Matrix_tools.h>

Implemente_instanciable(Terme_dissipation_echelle_temporelle_turbulente_P0_CoviMAC,"Terme_dissipation_echelle_temporelle_turbulente_P0_CoviMAC", Source_base);
// XD Terme_dissipation_echelle_temporelle_turbulente_P0_CoviMAC source_base Terme_dissipation_echelle_temporelle_turbulente_P0_CoviMAC 0 Source term which corresponds to the dissipation source term that appears in the transport equation for tau (in the k-tau turbulence model)

Sortie& Terme_dissipation_echelle_temporelle_turbulente_P0_CoviMAC::printOn(Sortie& os) const
{
  return os;
}

Entree& Terme_dissipation_echelle_temporelle_turbulente_P0_CoviMAC::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("beta_omega", &beta_omega);
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Terme_dissipation_echelle_temporelle_turbulente_P0_CoviMAC::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const Zone_CoviMAC& zone = ref_cast(Zone_CoviMAC, equation().zone_dis().valeur());
  const int ne = zone.nb_elem(), ne_tot = zone.nb_elem_tot(), Ntau = equation().inconnue().valeurs().line_size(); //, Nt = equation().probleme().get_champ("temperature").valeurs().line_size(), Np = equation().probleme().get_champ("pression").valeurs().line_size();
  //const int Nphases = equation().probleme().get_champ("alpha").valeurs().line_size();

//  for (auto &&i_m : matrices) if (i_m.first == "alpha" || i_m.first == "temperature" || i_m.first == "pression")
//      {
//        Matrice_Morse mat;
//        IntTrav stencil(0, 2);
//        stencil.set_smart_resize(1);
//        if (i_m.first == "alpha") for (int e = 0; e < ne; e++) for(int n = 0; n < Ntau ; n++) stencil.append_line(Ntau * e + n, Nphases * e + n);
//        if (i_m.first == "temperature") for (int e = 0; e < ne; e++) for(int n = 0; n < Ntau ; n++) stencil.append_line(Ntau * e + n, Nt * e + n);
//        if (i_m.first == "pression") for (int e = 0; e < ne; e++) for(int n = 0, mp = 0; n < Ntau ; n++, mp += (Np > 1)) stencil.append_line(Ntau * e + n, Np * e + mp); // attention reduction en pression
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
          for (int e = 0; e < ne; e++) for (int n = 0; n < Ntau; n++) sten.append_line(Ntau * e + n, M * e + n);
        if (n_m.first == "pression" )
          for (int e = 0; e < ne; e++) for (int n = 0, m = 0; n < Ntau; n++, m+=(M>1)) sten.append_line(Ntau * e + n, M * e + m);
        Matrix_tools::allocate_morse_matrix(Ntau * ne_tot, M * nc, sten, mat2);
        mat.nb_colonnes() ? mat += mat2 : mat = mat2;
      }
}

void Terme_dissipation_echelle_temporelle_turbulente_P0_CoviMAC::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl)  const
{
  const Zone_CoviMAC& 		zone 			= ref_cast(Zone_CoviMAC, equation().zone_dis().valeur());
  const Champ_P0_CoviMAC& 	ch_tau 			= ref_cast(Champ_P0_CoviMAC,equation().inconnue().valeur()); // Champ tau
  const DoubleTab& 			tau 			= ch_tau.valeurs() ;
  const Champ_base& 		ch_alpha_rho 	= sub_type(Pb_Multiphase,equation().probleme()) ? ref_cast(Pb_Multiphase,equation().probleme()).eq_masse.champ_conserve() : equation().milieu().masse_volumique();
  const DoubleTab& 			alpha_rho		= ch_alpha_rho.valeurs();
  const tabs_t& 			der_alpha_rho 	= ref_cast(Champ_Inc_base, ch_alpha_rho).derivees(); // dictionnaire des derivees

  const int nb_elem = zone.nb_elem(), Ntau = tau.line_size(), Np = equation().probleme().get_champ("pression").valeurs().line_size(), Nt = equation().probleme().get_champ("temperature").valeurs().line_size();
  const int Na = sub_type(Pb_Multiphase,equation().probleme()) ? equation().probleme().get_champ("alpha").valeurs().line_size() : 1;

  int mtau, mp;

  for (int e = 0; e < nb_elem; e++) for (mtau = 0, mp = 0; mtau < Ntau; mtau++, mp += (Np > 1))
      {
        secmem(e, mtau) -= beta_omega * alpha_rho(e, mtau) ;
        for (auto &&i_m : matrices)
          {
            Matrice_Morse& mat = *i_m.second;
            if (i_m.first == "alpha") 		mat(Ntau * e + mtau, Na * e + mtau) += beta_omega * (der_alpha_rho.count("alpha") ? der_alpha_rho.at("alpha")(e, mtau) : NULL );
            if (i_m.first == "temperature") mat(Ntau * e + mtau, Nt * e + mtau) += beta_omega * (der_alpha_rho.count("temperature") ? der_alpha_rho.at("temperature")(e, mtau) : NULL );	// derivee par rapport a la temperature
            if (i_m.first == "pression") 	mat(Ntau * e + mtau, Np * e + mp) 	+= beta_omega * (der_alpha_rho.count("pression") ? der_alpha_rho.at("pression")(e, mtau) : NULL );		// derivee par rapport a la pression
          }
      }
}


