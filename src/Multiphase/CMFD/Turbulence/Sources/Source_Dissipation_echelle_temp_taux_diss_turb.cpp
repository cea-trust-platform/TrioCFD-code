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
// File:        Source_Dissipation_echelle_temp_taux_diss_turb.cpp
// Directory:   $TRUST_ROOT/src/Turbulence/Sources
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Dissipation_echelle_temp_taux_diss_turb.h>

#include <Echelle_temporelle_turbulente.h>
#include <Taux_dissipation_turbulent.h>
#include <Milieu_composite.h>
#include <Equation_base.h>
#include <Pb_Multiphase.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <Domaine_VF.h>

Implemente_base(Source_Dissipation_echelle_temp_taux_diss_turb,"Source_Dissipation_echelle_temp_taux_diss_turb", Sources_Multiphase_base);
// XD Source_Dissipation_echelle_temp_taux_diss_turb source_base Source_Dissipation_echelle_temp_taux_diss_turb 0 Source term which corresponds to the dissipation source term that appears in the transport equation for tau (in the k-tau turbulence model)

Sortie& Source_Dissipation_echelle_temp_taux_diss_turb::printOn(Sortie& os) const {   return os;}

Entree& Source_Dissipation_echelle_temp_taux_diss_turb::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("beta_omega", &beta_omega);
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Source_Dissipation_echelle_temp_taux_diss_turb::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl)  const
{
  const Domaine_VF& domaine = ref_cast(Domaine_VF, equation().domaine_dis());
  const DoubleTab& diss     = equation().inconnue().valeurs() ;
  const DoubleTab& pdiss    = equation().inconnue().passe() ;
  const DoubleVect& pe      = equation().milieu().porosite_elem(), &ve = domaine.volumes();

  const int nb_elem = domaine.nb_elem(), N = diss.line_size();

  std::string Type_diss = ""; // omega or tau dissipation
  if sub_type(Echelle_temporelle_turbulente, equation()) Type_diss = "tau";
  else if sub_type(Taux_dissipation_turbulent, equation()) Type_diss = "omega";
  if (Type_diss == "") abort();

  int n;

  assert( N == 1 );

  for (int e = 0; e < nb_elem; e++)
    for (n = 0; n<N; n++)
      {
        if (Type_diss == "tau")
          {
            double secmem_en  = pe(e) * ve(e) * beta_omega ;
            secmem(e, n) += secmem_en ;
          }
        else if (Type_diss == "omega")
          {
            double secmem_en  = - pe(e) * ve(e) * beta_omega * pdiss(e,n) * (  pdiss(e,n) + 2 * (diss(e,n) - pdiss(e,n) ) );
            secmem(e, n) += secmem_en ;
            for (auto &&i_m : matrices)
              if (i_m.first == "omega") (*i_m.second)(N*e+n, N*e+n) += pe(e) * ve(e) * beta_omega * 2* pdiss(e,n) ;
          }
      }
}


