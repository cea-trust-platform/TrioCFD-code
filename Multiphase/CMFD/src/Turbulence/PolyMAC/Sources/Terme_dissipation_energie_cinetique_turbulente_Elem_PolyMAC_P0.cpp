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
// File:        Travail_pression_PolyMAC_P0.cpp
// Directory:   $TRUST_ROOT/src/Turbulence/PolyMAC_P0/Sources
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Terme_dissipation_energie_cinetique_turbulente_Elem_PolyMAC_P0.h>
#include <Zone_PolyMAC_P0.h>
#include <Champ_Elem_PolyMAC_P0.h>
#include <Equation_base.h>
#include <Pb_Multiphase.h>
#include <Milieu_composite.h>
#include <Array_tools.h>
#include <Matrix_tools.h>
#include <Op_Diff_Turbulent_PolyMAC_P0_Face.h>
#include <Viscosite_turbulente_k_tau.h>
#include <Viscosite_turbulente_k_omega.h>
#include <Echelle_temporelle_turbulente.h>
#include <Taux_dissipation_turbulent.h>

Implemente_instanciable(Terme_dissipation_energie_cinetique_turbulente_Elem_PolyMAC_P0,"Terme_dissipation_energie_cinetique_turbulente_Elem_PolyMAC_P0", Source_base);
// XD Terme_dissipation_energie_cinetique_turbulente source_base Terme_dissipation_energie_cinetique_turbulente 0 Source term which corresponds to the dissipation source term that appears in the turbulent kinetic energy equation

Sortie& Terme_dissipation_energie_cinetique_turbulente_Elem_PolyMAC_P0::printOn(Sortie& os) const
{
  return os;
}

Entree& Terme_dissipation_energie_cinetique_turbulente_Elem_PolyMAC_P0::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("beta_k", &beta_k);
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Terme_dissipation_energie_cinetique_turbulente_Elem_PolyMAC_P0::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const Zone_PolyMAC_P0& zone = ref_cast(Zone_PolyMAC_P0, equation().zone_dis().valeur());
  const DoubleTab& k 	 = ref_cast(Champ_Elem_PolyMAC_P0, equation().inconnue().valeur()).valeurs();
  const int ne = zone.nb_elem(), ne_tot = zone.nb_elem_tot(), Nk = k.line_size();

  std::string Type_diss = ""; // omega or tau dissipation
  for (int i = 0 ; i < equation().probleme().nombre_d_equations() ; i++)
    {
      if sub_type(Echelle_temporelle_turbulente, equation().probleme().equation(i)) Type_diss = "tau";
      else if sub_type(Taux_dissipation_turbulent, equation().probleme().equation(i)) Type_diss = "omega";
    }
  if (Type_diss == "") abort();

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

void Terme_dissipation_energie_cinetique_turbulente_Elem_PolyMAC_P0::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl)  const
{
  const Zone_PolyMAC_P0& 					zone 		= ref_cast(Zone_PolyMAC_P0, equation().zone_dis().valeur());
  const Champ_Elem_PolyMAC_P0& 				ch_k	= ref_cast(Champ_Elem_PolyMAC_P0, equation().inconnue().valeur());		// Champ k
  const DoubleTab& 						      k 	= ch_k.valeurs();
  const Champ_Inc_base& ch_alpha_rho_k 	= equation().champ_conserve();
  const DoubleTab& 				alpha_rho_k		= ch_alpha_rho_k.valeurs();
  const tabs_t& 				der_alpha_rho_k = ref_cast(Champ_Inc_base, ch_alpha_rho_k).derivees(); // dictionnaire des derivees
  const Navier_Stokes_std&     eq_qdm 	= ref_cast(Navier_Stokes_std, equation().probleme().equation(0));
  const Op_Diff_Turbulent_PolyMAC_P0_Face& op_diff 		= ref_cast(Op_Diff_Turbulent_PolyMAC_P0_Face, eq_qdm.operateur(0).l_op_base());
  const Viscosite_turbulente_base&   	visc_turb 		= ref_cast(Viscosite_turbulente_base, op_diff.correlation().valeur());
  const DoubleTab&                      nu 		  		= equation().probleme().get_champ("viscosite_cinematique").passe();
  const DoubleVect& pe = zone.porosite_elem(), &ve = zone.volumes();

  std::string Type_diss = ""; // omega or tau dissipation
  for (int i = 0 ; i < equation().probleme().nombre_d_equations() ; i++)
    {
      if sub_type(Echelle_temporelle_turbulente, equation().probleme().equation(i)) Type_diss = "tau";
      else if sub_type(Taux_dissipation_turbulent, equation().probleme().equation(i)) Type_diss = "omega";
    }
  if (Type_diss == "") abort();
  const Champ_Elem_PolyMAC_P0& 			ch_diss 	= ref_cast(Champ_Elem_PolyMAC_P0,equation().probleme().get_champ(Nom(Type_diss.c_str()))); // Champ tau ou omega
  const DoubleTab& 					      	diss	= ch_diss.valeurs() ;


  const int Nk = k.line_size(), Np = equation().probleme().get_champ("pression").valeurs().line_size(), Na = equation().probleme().get_champ("alpha").valeurs().line_size(), Nt = equation().probleme().get_champ("temperature").valeurs().line_size(), nb_elem = zone.nb_elem();

  Matrice_Morse *Ma = matrices.count("alpha") ? matrices.at("alpha") : nullptr,
                 *Mk = matrices.count(ch_k.le_nom().getString()) ? matrices.at(ch_k.le_nom().getString()) : nullptr,
                  *Mdiss = matrices.count(Type_diss) ? matrices.at(Type_diss) : nullptr,
                   *Mp = matrices.count("pression") ? matrices.at("pression") : nullptr,
                    *Mt	= matrices.count("temperature") ? matrices.at("temperature") : nullptr;

  for (int e = 0; e < nb_elem; e++)
    for (int mk = 0, mp = 0; mk < Nk; mk++, mp += (Np > 1))
      {
        if (Type_diss == "tau")
          {
            double inv_tau =  k(e, mk) / std::max(k(e, mk) * diss(e, mk), visc_turb.limiteur() * nu(e, mk));
            secmem(e, mk) -= pe(e) * ve(e) * beta_k * alpha_rho_k(e,mk) * inv_tau;
            if (!(Ma==nullptr)) 	(*Ma)(Nk * e + mk, Na * e + mk)   	  += pe(e) * ve(e) * beta_k * (der_alpha_rho_k.count("alpha") ? der_alpha_rho_k.at("alpha")(e,mk) : 0 ) * inv_tau;	// derivee en alpha
            if (!(Mt==nullptr)) 	(*Mt)(Nk * e + mk, Nt * e + mk)       += pe(e) * ve(e) * beta_k * (der_alpha_rho_k.count("temperature") ? der_alpha_rho_k.at("temperature")(e, mk) : 0 ) * inv_tau;	// derivee par rapport a la temperature
            if (!(Mp==nullptr)) 	(*Mp)(Nk * e + mk, Np * e + mp)       += pe(e) * ve(e) * beta_k * (der_alpha_rho_k.count("pression") ? der_alpha_rho_k.at("pression")(e, mk) : 0 ) * inv_tau;		// derivee par rapport a la pression
            if (!(Mk==nullptr))
              {
                if (k(e, mk) * diss(e, mk) > visc_turb.limiteur() * nu(e, mk))
                  (*Mk)(Nk * e + mk, Nk * e + mk)       += pe(e) * ve(e) * beta_k * (der_alpha_rho_k.count("k") ? der_alpha_rho_k.at("k")(e,mk) : 0 ) / diss(e, mk); // derivee en k ; depend de l'activation ou non du limiteur
                else
                  (*Mk)(Nk * e + mk, Nk * e + mk)       += pe(e) * ve(e) * 2 * beta_k * alpha_rho_k(e, mk) / (visc_turb.limiteur() * nu(e, mk)); // derivee en k
              }
            if (!(Mdiss==nullptr))
              {
                if ( k(e, mk) * diss(e, mk) > visc_turb.limiteur() * nu(e, mk))
                  (*Mdiss)(Nk * e + mk, Nk * e + mk)       += pe(e) * ve(e) * beta_k * alpha_rho_k(e, mk) * (-1)/(diss(e,mk)*diss(e,mk)); // derivee en tau  ; depend de l'activation ou non du limiteur
                else
                  (*Mdiss)(Nk * e + mk, Nk * e + mk)       += 0;
              }
          }
        else if (Type_diss == "omega")
          {
            secmem(e, mk) -= pe(e) * ve(e) * beta_k * alpha_rho_k(e,mk) * diss(e, mk);
            if (!(Ma==nullptr)) 	(*Ma)(Nk * e + mk, Na * e + mk)   	  += pe(e) * ve(e) * beta_k * (der_alpha_rho_k.count("alpha") ? der_alpha_rho_k.at("alpha")(e,mk) : 0 ) * diss(e, mk);	// derivee en alpha
            if (!(Mt==nullptr)) 	(*Mt)(Nk * e + mk, Nt * e + mk)       += pe(e) * ve(e) * beta_k * (der_alpha_rho_k.count("temperature") ? der_alpha_rho_k.at("temperature")(e, mk) : 0 ) * diss(e, mk);	// derivee par rapport a la temperature
            if (!(Mp==nullptr)) 	(*Mp)(Nk * e + mk, Np * e + mp)       += pe(e) * ve(e) * beta_k * (der_alpha_rho_k.count("pression") ? der_alpha_rho_k.at("pression")(e, mk) : 0 ) * diss(e, mk);		// derivee par rapport a la pression
            if (!(Mk==nullptr))   (*Mk)(Nk * e + mk, Nk * e + mk)       += pe(e) * ve(e) * beta_k * (der_alpha_rho_k.count("k") ? der_alpha_rho_k.at("k")(e,mk) : 0 ) * diss(e, mk); // derivee en k
            if (!(Mdiss==nullptr)) (*Mdiss)(Nk * e + mk, Nk * e + mk)   += pe(e) * ve(e) * beta_k * alpha_rho_k(e, mk); // derivee en omega
          }
      }
}

