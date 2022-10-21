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
// File:        Production_echelle_temp_taux_diss_turb_PolyMAC_P0.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/PolyMAC_P0
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Production_echelle_temp_taux_diss_turb_PolyMAC_P0.h>

#include <Domaine_PolyMAC_P0.h>
#include <Champ_Elem_PolyMAC_P0.h>
#include <Probleme_base.h>
#include <grad_Champ_Face_PolyMAC_P0.h>
#include <Champ_Uniforme.h>
#include <Flux_interfacial_base.h>
#include <Milieu_composite.h>
#include <Operateur_Diff.h>
#include <Op_Diff_Turbulent_PolyMAC_P0_Face.h>
#include <Navier_Stokes_std.h>
#include <Viscosite_turbulente_k_tau.h>
#include <Energie_cinetique_turbulente.h>
#include <Array_tools.h>
#include <Matrix_tools.h>
#include <Echelle_temporelle_turbulente.h>
#include <Taux_dissipation_turbulent.h>
#include <Pb_Multiphase.h>


Implemente_instanciable(Production_echelle_temp_taux_diss_turb_PolyMAC_P0,"Production_echelle_temp_taux_diss_turb_Elem_PolyMAC_P0", Source_base);

Sortie& Production_echelle_temp_taux_diss_turb_PolyMAC_P0::printOn(Sortie& os) const
{
  return os;
}

Entree& Production_echelle_temp_taux_diss_turb_PolyMAC_P0::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("alpha_omega", &alpha_omega_, Param::REQUIRED);
  param.lire_avec_accolades_depuis(is);
  Cout << "alpha_omega = " << alpha_omega_ << finl ;

  return is;
}

void Production_echelle_temp_taux_diss_turb_PolyMAC_P0::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const Domaine_PolyMAC_P0&       domaine = ref_cast(Domaine_PolyMAC_P0, equation().domaine_dis().valeur());
  int ne = domaine.nb_elem(), ne_tot = domaine.nb_elem_tot(), Nk = equation().inconnue().valeurs().line_size();
  assert(Nk == 1); // si plus d'une phase turbulente, il vaut mieux iterer sur les id_composites des phases turbulentes modelisees par un modele k-tau

  for (auto &&n_m : matrices)
    if (n_m.first == "alpha" || n_m.first == "temperature" || n_m.first == "pression")
      {
        Matrice_Morse& mat = *n_m.second, mat2;
        const DoubleTab& dep = equation().probleme().get_champ(n_m.first.c_str()).valeurs();
        int nc = dep.dimension_tot(0),
            M  = dep.line_size();
        IntTrav sten(0, 2);
        sten.set_smart_resize(1);
        if (n_m.first == "alpha" || n_m.first == "temperature")
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

void Production_echelle_temp_taux_diss_turb_PolyMAC_P0::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Domaine_PolyMAC_P0&                   domaine = ref_cast(Domaine_PolyMAC_P0, equation().domaine_dis().valeur());
  const Probleme_base&                       pb = ref_cast(Probleme_base, equation().probleme());
  const Navier_Stokes_std&               eq_qdm = ref_cast(Navier_Stokes_std, pb.equation(0));
  const grad_Champ_Face_PolyMAC_P0&        grad = ref_cast(grad_Champ_Face_PolyMAC_P0, eq_qdm.get_champ("gradient_vitesse"));
  const DoubleTab&                     tab_grad = grad.valeurs();
  const DoubleTab&                     tab_diss = ref_cast(Champ_Elem_PolyMAC_P0, equation().inconnue().valeur()).valeurs(); // tau ou omega selon l'equation
  const DoubleTab&                    tab_pdiss = ref_cast(Champ_Elem_PolyMAC_P0, equation().inconnue().valeur()).passe(); // tau ou omega selon l'equation
  const DoubleVect& pe = equation().milieu().porosite_elem(), &ve = domaine.volumes();

  const Champ_base&   ch_alpha_rho = sub_type(Pb_Multiphase,equation().probleme()) ? ref_cast(Pb_Multiphase,equation().probleme()).eq_masse.champ_conserve() : equation().milieu().masse_volumique().valeur();
  const DoubleTab&       alpha_rho = ch_alpha_rho.valeurs();
  const tabs_t&      der_alpha_rho = ref_cast(Champ_Inc_base, ch_alpha_rho).derivees(); // dictionnaire des derivees

  int nf_tot = domaine.nb_faces_tot(), ne = domaine.nb_elem(), D = dimension ;
  int N = equation().inconnue()->valeurs().line_size(),
      Np = equation().probleme().get_champ("pression").valeurs().line_size(),
      Nt = equation().probleme().get_champ("temperature").valeurs().line_size(),
      Na = sub_type(Pb_Multiphase,equation().probleme()) ? equation().probleme().get_champ("alpha").valeurs().line_size() : 1;
  int e, n, mp;

  std::string Type_diss = ""; // omega or tau dissipation
  if sub_type(Echelle_temporelle_turbulente, equation()) Type_diss = "tau";
  else if sub_type(Taux_dissipation_turbulent, equation()) Type_diss = "omega";
  if (Type_diss == "") Process::exit(que_suis_je() + " : you must have tau or omega dissipation ! ");

  // Second membre
  for(e = 0 ; e < ne ; e++)
    for(n = 0, mp = 0; n<N ; n++, mp += (Np > 1))
      {
        double grad_grad = 0.;
        for (int d_U = 0; d_U < D; d_U++)
          for (int d_X = 0; d_X < D; d_X++)
            grad_grad += ( tab_grad(nf_tot + d_U + e * D , D * n + d_X) + tab_grad(nf_tot + d_X + e * D , D * n + d_U) ) * tab_grad(nf_tot + d_X + e * D , D * n + d_U) ;

        double fac = std::max(grad_grad, 0.) * pe(e) * ve(e) * alpha_omega_ ;

        if (Type_diss == "tau")
          {
            secmem(e, n) -= fac * alpha_rho(e, n) * (2*tab_diss(e, n)-tab_pdiss(e, n)) * tab_pdiss(e, n) ; // tau has negative production
            for (auto &&i_m : matrices)
              {
                Matrice_Morse& mat = *i_m.second;
                if (i_m.first == "tau")         mat(N * e + n, N  * e + n) += fac *  2 *                                 tab_pdiss(e, n) * alpha_rho(e, n) ;
                if (i_m.first == "alpha") 	    mat(N * e + n, Na * e + n) += fac * (2*tab_diss(e, n)-tab_pdiss(e, n)) * tab_pdiss(e, n) * (der_alpha_rho.count("alpha") ?       der_alpha_rho.at("alpha")(e, n) : 0 );		   // derivee par rapport au taux de vide
                if (i_m.first == "temperature") mat(N * e + n, Nt * e + n) += fac * (2*tab_diss(e, n)-tab_pdiss(e, n)) * tab_pdiss(e, n) * (der_alpha_rho.count("temperature") ? der_alpha_rho.at("temperature")(e, n) : 0 );// derivee par rapport a la temperature
                if (i_m.first == "pression")    mat(N * e + n, Np * e + mp)+= fac * (2*tab_diss(e, n)-tab_pdiss(e, n)) * tab_pdiss(e, n) * (der_alpha_rho.count("pression") ?    der_alpha_rho.at("pression")(e, mp) : 0 );	 // derivee par rapport a la pression
              }
          }
        else if (Type_diss == "omega")
          {
            secmem(e, n) +=   fac * alpha_rho(e, n);
            for (auto &&i_m : matrices)
              {
                Matrice_Morse& mat = *i_m.second;
                if (i_m.first == "alpha") 	    mat(N * e + n, Na * e + n) += fac * (der_alpha_rho.count("alpha") ?       der_alpha_rho.at("alpha")(e, n) : 0 );		   // derivee par rapport au taux de vide
                if (i_m.first == "temperature") mat(N * e + n, Nt * e + n) += fac * (der_alpha_rho.count("temperature") ? der_alpha_rho.at("temperature")(e, n) : 0 );// derivee par rapport a la temperature
                if (i_m.first == "pression")    mat(N * e + n, Np * e + mp)+= fac * (der_alpha_rho.count("pression") ?    der_alpha_rho.at("pression")(e, mp) : 0 );	 // derivee par rapport a la pression
              }
          }
      }
}
