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
// File:        Source_Dissipation_energie_cin_turb.cpp
// Directory:   $TRUST_ROOT/src/Turbulence/Sources
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Dissipation_energie_cin_turb.h>

#include <Echelle_temporelle_turbulente.h>
#include <Taux_dissipation_turbulent.h>
#include <Viscosite_turbulente_base.h>
#include <Pb_Multiphase.h>
#include <Array_tools.h>
#include <Matrix_tools.h>

Implemente_base(Source_Dissipation_energie_cin_turb,"Source_Dissipation_energie_cin_turb", Sources_Multiphase_base);
// XD Terme_dissipation_energie_cinetique_turbulente source_base Terme_dissipation_energie_cinetique_turbulente 1 Dissipation source term used in the TKE equation
// XD attr beta_k floattant beta_k 1 Constant for the used model

// XD Production_echelle_temp_taux_diss_turb source_base Production_echelle_temp_taux_diss_turb -1 Production source term used in the tau and omega equations
// XD attr alpha_omega floattant alpha_omega 1 Constant for the used model

// XD Dissipation_echelle_temp_taux_diss_turb source_base Dissipation_echelle_temp_taux_diss_turb -1 Dissipation source term used in the tau and omega equations
// XD attr beta_omega floattant beta_omega 1 Constant for the used model

// XD Diffusion_croisee_echelle_temp_taux_diss_turb source_base Diffusion_croisee_echelle_temp_taux_diss_turb -1 Cross-diffusion source term used in the tau and omega equations
// XD attr sigma_d floattant sigma_d 1 Constant for the used model


Sortie& Source_Dissipation_energie_cin_turb::printOn(Sortie& os) const
{
  return os;
}

Entree& Source_Dissipation_energie_cin_turb::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("beta_k", &beta_k);
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Source_Dissipation_energie_cin_turb::completer()
{
  const Navier_Stokes_std&     eq_qdm 	= ref_cast(Navier_Stokes_std, equation().probleme().equation(0));
  if (ref_cast(Operateur_Diff_base, eq_qdm.operateur(0).l_op_base()).correlation_viscosite_turbulente()==nullptr) Process::exit(que_suis_je() + " : the momentum diffusion must be turbulent !");
}

void Source_Dissipation_energie_cin_turb::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const Domaine_VF& domaine = ref_cast(Domaine_VF, equation().domaine_dis().valeur());
  const DoubleTab& k 	 = equation().inconnue().valeur().valeurs();
  const int ne = domaine.nb_elem(), ne_tot = domaine.nb_elem_tot(), Nk = k.line_size();

  std::string Type_diss = ""; // omega or tau dissipation
  for (int i = 0 ; i < equation().probleme().nombre_d_equations() ; i++)
    {
      if sub_type(Echelle_temporelle_turbulente, equation().probleme().equation(i)) Type_diss = "tau";
      else if sub_type(Taux_dissipation_turbulent, equation().probleme().equation(i)) Type_diss = "omega";
    }
  if (Type_diss == "") abort();

  assert(Nk == 1); // si plus d'une phase turbulente, il vaut mieux iterer sur les id_composites des phases turbulentes modelisees par un modele k-tau
  if (Type_diss == "tau") assert(equation().probleme().get_champ("tau").valeurs().line_size() == 1);
  if (Type_diss == "omega") assert(equation().probleme().get_champ("omega").valeurs().line_size() == 1);

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

void Source_Dissipation_energie_cin_turb::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl)  const
{
  const Domaine_VF&             domaine = ref_cast(Domaine_VF, equation().domaine_dis().valeur());
  const DoubleTab&                    k = equation().inconnue().valeur().valeurs();
  const Champ_Inc_base&  ch_alpha_rho_k = equation().champ_conserve();
  const DoubleTab&          alpha_rho_k = ch_alpha_rho_k.passe();
  const tabs_t&         der_alpha_rho_k = ref_cast(Champ_Inc_base, ch_alpha_rho_k).derivees(); // dictionnaire des derivees
  const Navier_Stokes_std&       eq_qdm = ref_cast(Navier_Stokes_std, equation().probleme().equation(0));
  const Viscosite_turbulente_base& visc_turb = ref_cast(Viscosite_turbulente_base, (*ref_cast(Operateur_Diff_base, eq_qdm.operateur(0).l_op_base()).correlation_viscosite_turbulente()).valeur());
  const DoubleTab& nu = equation().probleme().get_champ("viscosite_cinematique").passe();
  const DoubleVect& pe = equation().milieu().porosite_elem(), &ve = domaine.volumes();
  double dt = equation().schema_temps().pas_de_temps();

  std::string Type_diss = ""; // omega or tau dissipation
  for (int i = 0 ; i < equation().probleme().nombre_d_equations() ; i++)
    {
      if sub_type(Echelle_temporelle_turbulente, equation().probleme().equation(i)) Type_diss = "tau";
      else if sub_type(Taux_dissipation_turbulent, equation().probleme().equation(i)) Type_diss = "omega";
    }
  if (Type_diss == "") abort();
  const DoubleTab&                diss = equation().probleme().get_champ(Nom(Type_diss.c_str())).passe() ;
  const DoubleTab&               pdiss = equation().probleme().get_champ(Nom(Type_diss.c_str())).passe() ;

  const int Nk = k.line_size(), Np = equation().probleme().get_champ("pression").valeurs().line_size(), Na = equation().probleme().get_champ("alpha").valeurs().line_size(), Nt = equation().probleme().get_champ("temperature").valeurs().line_size(), nb_elem = domaine.nb_elem();

  Matrice_Morse *Ma = matrices.count("alpha") ? matrices.at("alpha") : nullptr,
                 *Mk = matrices.count("k") ? matrices.at("k") : nullptr,
                  *Mdiss = matrices.count(Type_diss) ? matrices.at(Type_diss) : nullptr,
                   *Mp = matrices.count("pression") ? matrices.at("pression") : nullptr,
                    *Mt	= matrices.count("temperature") ? matrices.at("temperature") : nullptr;

  for (int e = 0; e < nb_elem; e++)
    for (int mk = 0, mp = 0; mk < Nk; mk++, mp += (Np > 1))
      {
        if (Type_diss == "tau")
          {
            double inv_tau = (k(e, mk) * diss(e, mk) > visc_turb.limiteur() * nu(e, mk))
                             ? 1./diss(e,mk)
                             : k(e, mk) / (visc_turb.limiteur() * nu(e, mk)) ;
            secmem(e, mk) -= pe(e) * ve(e) * beta_k * alpha_rho_k(e,mk) * inv_tau;
            if (!(Ma==nullptr)) 	(*Ma)(Nk * e + mk, Na * e + mk)   	  += pe(e) * ve(e) * beta_k * (der_alpha_rho_k.count("alpha") ?       der_alpha_rho_k.at("alpha")(e,mk) : 0 )        * inv_tau;	// derivee en alpha
            if (!(Mt==nullptr)) 	(*Mt)(Nk * e + mk, Nt * e + mk)       += pe(e) * ve(e) * beta_k * (der_alpha_rho_k.count("temperature") ? der_alpha_rho_k.at("temperature")(e, mk) : 0 ) * inv_tau;	// derivee par rapport a la temperature
            if (!(Mp==nullptr)) 	(*Mp)(Nk * e + mk, Np * e + mp)       += pe(e) * ve(e) * beta_k * (der_alpha_rho_k.count("pression") ?    der_alpha_rho_k.at("pression")(e, mp) : 0 )    * inv_tau;		// derivee par rapport a la pression
            if (!(Mk==nullptr))
              {
                if (k(e, mk) * diss(e,mk) > visc_turb.limiteur() * nu(e, mk))
                  (*Mk)(Nk * e + mk, Nk * e + mk)       += pe(e) * ve(e) * beta_k * (der_alpha_rho_k.count("k") ? der_alpha_rho_k.at("k")(e,mk) : 0 ) * inv_tau; // derivee en k ; depend de l'activation ou non du limiteur
                else
                  (*Mk)(Nk * e + mk, Nk * e + mk)       += pe(e) * ve(e) * 2 * beta_k * alpha_rho_k(e, mk) / (visc_turb.limiteur() * nu(e, mk)); // derivee en k
              }
            if (!(Mdiss==nullptr))
              {
                if ( k(e, mk) * diss(e,mk) > visc_turb.limiteur() * nu(e, mk))
                  (*Mdiss)(Nk * e + mk, Nk * e + mk)       += pe(e) * ve(e) * beta_k * alpha_rho_k(e, mk) * (-1)/(diss(e,mk)*diss(e,mk)); // derivee en tau  ; depend de l'activation ou non du limiteur
                else
                  (*Mdiss)(Nk * e + mk, Nk * e + mk)       += 0*pdiss(e, mk);
              }
          }
        else if (Type_diss == "omega")
          {
            secmem(e, mk) -= pe(e) * ve(e) * beta_k * (alpha_rho_k(e,mk)*0. + k(e,mk)) * diss(e, mk) * dt/dt ;
            if (Ma) 	(*Ma)(Nk * e + mk, Na * e + mk)   	  += pe(e) * ve(e) * beta_k * (der_alpha_rho_k.count("alpha") ?       der_alpha_rho_k.at("alpha")(e,mk) : 0 )        * diss(e, mk)*0;	// derivee en alpha
            if (Mt) 	(*Mt)(Nk * e + mk, Nt * e + mk)       += pe(e) * ve(e) * beta_k * (der_alpha_rho_k.count("temperature") ? der_alpha_rho_k.at("temperature")(e, mk) : 0 ) * diss(e, mk)*0;	// derivee par rapport a la temperature
            if (Mp) 	(*Mp)(Nk * e + mk, Np * e + mp)       += pe(e) * ve(e) * beta_k * (der_alpha_rho_k.count("pression") ?    der_alpha_rho_k.at("pression")(e, mp) : 0 )    * diss(e, mk)*0;		// derivee par rapport a la pression
            if (Mk)   (*Mk)(Nk * e + mk, Nk * e + mk)       += pe(e) * ve(e) * beta_k * (der_alpha_rho_k.count("k") ?           der_alpha_rho_k.at("k")(e,mk) : 0 )            * diss(e, mk)*0; // derivee en k
            if (Mk)   (*Mk)(Nk * e + mk, Nk * e + mk)       += pe(e) * ve(e) * beta_k * diss(e, mk) ; // derivee en k
            if (Mdiss) (*Mdiss)(Nk * e + mk, Nk * e + mk)   += pe(e) * ve(e) * beta_k * k(e, mk) *0.; // derivee en omega
          }
      }
}

