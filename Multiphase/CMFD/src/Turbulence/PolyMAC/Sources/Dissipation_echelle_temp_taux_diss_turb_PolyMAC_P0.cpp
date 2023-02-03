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
// File:        Dissipation_echelle_temp_taux_diss_turb_PolyMAC_P0.cpp
// Directory:   $TRUST_ROOT/src/Turbulence/PolyMAC_P0/Sources
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Dissipation_echelle_temp_taux_diss_turb_PolyMAC_P0.h>
#include <Domaine_PolyMAC_P0.h>
#include <Champ_Elem_PolyMAC_P0.h>
#include <Equation_base.h>
#include <Pb_Multiphase.h>
#include <Milieu_composite.h>
#include <Array_tools.h>
#include <Matrix_tools.h>
#include <Echelle_temporelle_turbulente.h>
#include <Taux_dissipation_turbulent.h>

Implemente_instanciable(Dissipation_echelle_temp_taux_diss_turb_PolyMAC_P0,"Dissipation_echelle_temp_taux_diss_turb_Elem_PolyMAC_P0", Source_base);
// XD Terme_dissipation_echelle_temporelle_turbulente_Elem_PolyMAC_P0 source_base Terme_dissipation_echelle_temporelle_turbulente_Elem_PolyMAC_P0 0 Source term which corresponds to the dissipation source term that appears in the transport equation for tau (in the k-tau turbulence model)

Sortie& Dissipation_echelle_temp_taux_diss_turb_PolyMAC_P0::printOn(Sortie& os) const
{
  return os;
}

Entree& Dissipation_echelle_temp_taux_diss_turb_PolyMAC_P0::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("beta_omega", &beta_omega);
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Dissipation_echelle_temp_taux_diss_turb_PolyMAC_P0::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const Domaine_PolyMAC_P0& domaine = ref_cast(Domaine_PolyMAC_P0, equation().domaine_dis().valeur());
  const int ne = domaine.nb_elem(), ne_tot = domaine.nb_elem_tot(), N = equation().inconnue().valeurs().line_size();

  assert( N == 1 ); // si Ntau > 1 il vaut mieux iterer sur les id_composites des phases turbulentes
  for (auto &&n_m : matrices)
    if (n_m.first == "alpha" || n_m.first == "temperature" || n_m.first == "pression")
      {
        Matrice_Morse& mat = *n_m.second, mat2;
        const DoubleTab& dep = equation().probleme().get_champ(n_m.first.c_str()).valeurs();
        int nc = dep.dimension_tot(0),
            M  = dep.line_size();
        IntTrav sten(0, 2);
        sten.set_smart_resize(1);
        if (n_m.first == "alpha" || n_m.first == "temperature") // N <= M
          for (int e = 0; e < ne; e++)
            for (int n = 0; n < N; n++) sten.append_line(N * e + n, M * e + n);
        if (n_m.first == "pression" )
          for (int e = 0; e < ne; e++)
            for (int n = 0, m = 0; n < N; n++, m+=(M>1)) sten.append_line(N * e + n, M * e + m);
        Matrix_tools::allocate_morse_matrix(N * ne_tot, M * nc, sten, mat2);
        mat.nb_colonnes() ? mat += mat2 : mat = mat2;
      }
}

void Dissipation_echelle_temp_taux_diss_turb_PolyMAC_P0::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl)  const
{
  const Domaine_PolyMAC_P0& 		domaine 			= ref_cast(Domaine_PolyMAC_P0, equation().domaine_dis().valeur());
  const Champ_Elem_PolyMAC_P0& 	ch_diss			= ref_cast(Champ_Elem_PolyMAC_P0,equation().inconnue().valeur()); // Champ tau ou omega
  const DoubleTab& 			diss 			= ch_diss.valeurs() ;
  const DoubleTab& 			pdiss 			= ch_diss.passe() ;
  const Champ_base& 		ch_alpha_rho 	= sub_type(Pb_Multiphase,equation().probleme()) ? ref_cast(Pb_Multiphase,equation().probleme()).eq_masse.champ_conserve() : equation().milieu().masse_volumique().valeur();
  const DoubleTab& 			alpha_rho		= ch_alpha_rho.valeurs();
  const tabs_t& 			der_alpha_rho 	= ref_cast(Champ_Inc_base, ch_alpha_rho).derivees(); // dictionnaire des derivees
  const DoubleVect& pe = equation().milieu().porosite_elem(), &ve = domaine.volumes();

  const int nb_elem = domaine.nb_elem(), N = diss.line_size(), Np = equation().probleme().get_champ("pression").valeurs().line_size(), Nt = equation().probleme().get_champ("temperature").valeurs().line_size();
  const int Na = sub_type(Pb_Multiphase,equation().probleme()) ? equation().probleme().get_champ("alpha").valeurs().line_size() : 1;

  std::string Type_diss = ""; // omega or tau dissipation
  if sub_type(Echelle_temporelle_turbulente, equation()) Type_diss = "tau";
  else if sub_type(Taux_dissipation_turbulent, equation()) Type_diss = "omega";
  if (Type_diss == "") abort();

  int m, mp;

  assert( N == 1 ); // si Ntau > 1 il vaut mieux iterer sur les id_composites des phases turbulentes

  for (int e = 0; e < nb_elem; e++)
    for (m = 0, mp = 0; m < N; m++, mp += (Np > 1))
      {
        if (Type_diss == "tau")
          {
            double secmem_en  = pe(e) * ve(e) * beta_omega * alpha_rho(e, m) ;
            secmem(e, m) += secmem_en ;
            for (auto &&i_m : matrices)
              {
                Matrice_Morse& mat = *i_m.second;
                if (i_m.first == "alpha") 		mat(N * e + m, Na * e + m)   -= pe(e) * ve(e) * beta_omega * (der_alpha_rho.count("alpha") ? der_alpha_rho.at("alpha")(e, m) : 0 );			// derivee par rapport au taux de vide
                if (i_m.first == "temperature") mat(N * e + m, Nt * e + m) -= pe(e) * ve(e) * beta_omega * (der_alpha_rho.count("temperature") ? der_alpha_rho.at("temperature")(e, m) : 0 );// derivee par rapport a la temperature
                if (i_m.first == "pression") 	mat(N * e + m, Np * e + mp)  -= pe(e) * ve(e) * beta_omega * (der_alpha_rho.count("pression") ? der_alpha_rho.at("pression")(e, m) : 0 );		// derivee par rapport a la pression
              }
          }
        else if (Type_diss == "omega")
          {
            double secmem_en  = - pe(e) * ve(e) * beta_omega * alpha_rho(e, m) * pdiss(e,m) * (  pdiss(e,m) + 2 * (diss(e,m) - pdiss(e,m) )  );
            secmem(e, m) += secmem_en ;
            for (auto &&i_m : matrices)
              {
                Matrice_Morse& mat = *i_m.second;
                if (i_m.first == "alpha") 		mat(N * e + m, Na * e + m)  += pe(e) * ve(e) * beta_omega * (der_alpha_rho.count("alpha") ? der_alpha_rho.at("alpha")(e, m) : 0 ) * pdiss(e,m) * (  pdiss(e,m) + 2 * (diss(e,m) - pdiss(e,m) )  );			// derivee par rapport au taux de vide
                if (i_m.first == "temperature") mat(N * e + m, Nt * e + m)+= pe(e) * ve(e) * beta_omega * (der_alpha_rho.count("temperature") ? der_alpha_rho.at("temperature")(e, m) : 0 ) * pdiss(e,m) * (  pdiss(e,m) + 2 * (diss(e,m) - pdiss(e,m) )  );// derivee par rapport a la temperature
                if (i_m.first == "pression") 	mat(N * e + m, Np * e + mp) += pe(e) * ve(e) * beta_omega * (der_alpha_rho.count("pression") ? der_alpha_rho.at("pression")(e, m) : 0 ) * pdiss(e,m) * (  pdiss(e,m) + 2 * (diss(e,m) - pdiss(e,m) )  );		// derivee par rapport a la pression
                if (i_m.first == "omega")   mat(N * e + m, N * e + m)   	+= pe(e) * ve(e) * beta_omega * alpha_rho(e, m) *2* pdiss(e,m) ;
              }
          }
      }
}


