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

#include <Echelle_temporelle_turbulente.h>
#include <Taux_dissipation_turbulent.h>
#include <grad_Champ_Face_PolyMAC_P0.h>
#include <Champ_Elem_PolyMAC_P0.h>
#include <Navier_Stokes_std.h>
#include <Pb_Multiphase.h>
#include <Domaine_VF.h>

Implemente_instanciable(Production_echelle_temp_taux_diss_turb_PolyMAC_P0,"Production_echelle_temp_taux_diss_turb_Elem_PolyMAC_P0", Source_Production_echelle_temp_taux_diss_turb);

Sortie& Production_echelle_temp_taux_diss_turb_PolyMAC_P0::printOn(Sortie& os) const { return Source_Production_echelle_temp_taux_diss_turb::printOn(os); }
Entree& Production_echelle_temp_taux_diss_turb_PolyMAC_P0::readOn(Entree& is) { return Source_Production_echelle_temp_taux_diss_turb::readOn(is);}

void Production_echelle_temp_taux_diss_turb_PolyMAC_P0::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Domaine_VF&                     domaine = ref_cast(Domaine_PolyMAC_P0, equation().domaine_dis());
  const DoubleTab&                     tab_diss = ref_cast(Champ_Elem_PolyMAC_P0, equation().inconnue()).valeurs(); // tau ou omega selon l'equation
  const DoubleTab&                    tab_pdiss = ref_cast(Champ_Elem_PolyMAC_P0, equation().inconnue()).passe(); // tau ou omega selon l'equation
  const DoubleTab&                     tab_grad = equation().probleme().get_champ("gradient_vitesse").passe();
  const DoubleTab&                          alp = equation().probleme().get_champ("alpha").valeurs();
  const DoubleVect& pe = equation().milieu().porosite_elem(), &ve = domaine.volumes();

  int nf_tot = domaine.nb_faces_tot(), ne = domaine.nb_elem(), D = dimension ;
  int N = equation().inconnue().valeurs().line_size();
  int Na = sub_type(Pb_Multiphase,equation().probleme()) ? equation().probleme().get_champ("alpha").valeurs().line_size() : 1;
  int e, n;

  std::string Type_diss = ""; // omega or tau dissipation
  if sub_type(Echelle_temporelle_turbulente, equation()) Type_diss = "tau";
  else if sub_type(Taux_dissipation_turbulent, equation()) Type_diss = "omega";
  if (Type_diss == "") Process::exit(que_suis_je() + " : you must have tau or omega dissipation ! ");

  // Second membre
  for(e = 0 ; e < ne ; e++)
    for(n = 0; n<N ; n++)
      {
        double grad_grad = 0.;
        for (int d_U = 0; d_U < D; d_U++)
          for (int d_X = 0; d_X < D; d_X++)
            grad_grad += ( tab_grad(nf_tot + d_U + e * D , D * n + d_X) + tab_grad(nf_tot + d_X + e * D , D * n + d_U) ) * tab_grad(nf_tot + d_X + e * D , D * n + d_U) ;

        double fac = std::max(grad_grad, 0.) * pe(e) * ve(e) * alpha_omega_ ;

        if (Type_diss == "tau")
          {
            secmem(e, n) -= fac * (2*tab_diss(e, n)-tab_pdiss(e, n)) * tab_pdiss(e, n) * alp(e,n) ; // tau has negative production
            for (auto &&i_m : matrices)
              {
                if (i_m.first == "tau")    (*i_m.second)(N*e+n, N*e+n) += fac *  2 * tab_pdiss(e, n) * alp(e,n) ;
                if (i_m.first == "alpha")  (*i_m.second)(N*e+n, Na*e+n) += fac * (2*tab_diss(e, n)-tab_pdiss(e, n)) * tab_pdiss(e, n) ;
              }
          }
        else if (Type_diss == "omega")
          {
            secmem(e, n) +=   fac * alp(e,n);
            for (auto &&i_m : matrices)
              if (i_m.first == "alpha")  (*i_m.second)(N*e+n, Na*e+n) += fac ;
          }
      }
}
