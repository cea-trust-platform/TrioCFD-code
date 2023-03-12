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
// File:        Production_echelle_temp_taux_diss_turb_VDF.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/PolyMAC_P0
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Production_echelle_temp_taux_diss_turb_VDF.h>

#include <Echelle_temporelle_turbulente.h>
#include <Taux_dissipation_turbulent.h>
#include <Navier_Stokes_std.h>
#include <Domaine_Cl_VDF.h>
#include <Champ_Face_VDF.h>
#include <Pb_Multiphase.h>
#include <Domaine_VF.h>

Implemente_instanciable(Production_echelle_temp_taux_diss_turb_VDF,"Production_echelle_temp_taux_diss_turb_VDF_P0_VDF", Source_Production_echelle_temp_taux_diss_turb);

Sortie& Production_echelle_temp_taux_diss_turb_VDF::printOn(Sortie& os) const { return Source_Production_echelle_temp_taux_diss_turb::printOn(os); }
Entree& Production_echelle_temp_taux_diss_turb_VDF::readOn(Entree& is) { return Source_Production_echelle_temp_taux_diss_turb::readOn(is);}

void Production_echelle_temp_taux_diss_turb_VDF::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Domaine_VF&                     domaine = ref_cast(Domaine_VF, equation().domaine_dis().valeur());
  const Probleme_base&                       pb = ref_cast(Probleme_base, equation().probleme());
  const Navier_Stokes_std&               eq_qdm = ref_cast(Navier_Stokes_std, pb.equation(0));
  const DoubleTab&                     tab_diss = equation().inconnue()->valeurs(); // tau ou omega selon l'equation
  const DoubleTab&                    tab_pdiss = equation().inconnue()->passe(); // tau ou omega selon l'equation
  const DoubleVect& pe = equation().milieu().porosite_elem(), &ve = domaine.volumes();

  const Champ_base&   ch_alpha_rho = sub_type(Pb_Multiphase,equation().probleme()) ? ref_cast(Pb_Multiphase,equation().probleme()).eq_masse.champ_conserve() : equation().milieu().masse_volumique().valeur();
  const DoubleTab&       alpha_rho = ch_alpha_rho.valeurs();
  const tabs_t&      der_alpha_rho = ref_cast(Champ_Inc_base, ch_alpha_rho).derivees(); // dictionnaire des derivees

  int ne = domaine.nb_elem(), D = dimension, ne_tot = domaine.nb_elem_tot() ;
  int N = equation().inconnue()->valeurs().line_size(),
      Np = equation().probleme().get_champ("pression").valeurs().line_size(),
      Nt = equation().probleme().get_champ("temperature").valeurs().line_size(),
      Na = sub_type(Pb_Multiphase,equation().probleme()) ? equation().probleme().get_champ("alpha").valeurs().line_size() : 1;
  int e, n, mp;

  std::string Type_diss = ""; // omega or tau dissipation
  if sub_type(Echelle_temporelle_turbulente, equation()) Type_diss = "tau";
  else if sub_type(Taux_dissipation_turbulent, equation()) Type_diss = "omega";
  if (Type_diss == "") Process::exit(que_suis_je() + " : you must have tau or omega dissipation ! ");

  const Champ_Face_VDF& ch_vit = ref_cast(Champ_Face_VDF, eq_qdm.inconnue().valeur());
  DoubleTrav    tab_grad(ne_tot, D, D, ch_vit.valeurs().line_size());//, tab_vit_liq(nf_tot);
//  boucle remplit  tab_vit_liq;
  ch_vit.calcul_duidxj(ch_vit.passe(), tab_grad, ref_cast(Domaine_Cl_VDF, ch_vit.domaine_Cl_dis().valeur()));

  // Second membre
  for(e = 0 ; e < ne ; e++)
    for(n = 0, mp = 0; n<N ; n++, mp += (Np > 1))
      {
        double grad_grad = 0.;
        for (int d_U = 0; d_U < D; d_U++)
          for (int d_X = 0; d_X < D; d_X++)
            grad_grad += ( tab_grad( e, d_U, d_X) + tab_grad( e, d_X, d_U) ) * tab_grad( e, d_U, d_X) ;

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
            secmem(e, n) +=   fac ;//* alpha_rho(e, n);
            /*            for (auto &&i_m : matrices)
                          {
                            Matrice_Morse& mat = *i_m.second;
                            if (i_m.first == "alpha") 	    mat(N * e + n, Na * e + n) += fac * (der_alpha_rho.count("alpha") ?       der_alpha_rho.at("alpha")(e, n) : 0 );		   // derivee par rapport au taux de vide
                            if (i_m.first == "temperature") mat(N * e + n, Nt * e + n) += fac * (der_alpha_rho.count("temperature") ? der_alpha_rho.at("temperature")(e, n) : 0 );// derivee par rapport a la temperature
                            if (i_m.first == "pression")    mat(N * e + n, Np * e + mp)+= fac * (der_alpha_rho.count("pression") ?    der_alpha_rho.at("pression")(e, mp) : 0 );	 // derivee par rapport a la pression
                          }
            */
          }
      }
}
