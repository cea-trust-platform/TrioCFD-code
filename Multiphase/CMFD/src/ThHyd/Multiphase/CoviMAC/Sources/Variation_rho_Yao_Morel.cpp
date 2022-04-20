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
// File:        Variation_rho_Yao_Morel.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Correlations
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Variation_rho_Yao_Morel.h>
#include <Pb_Multiphase.h>
#include <Champ_P0_CoviMAC.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <math.h>

Implemente_instanciable(Variation_rho_Yao_Morel, "Variation_rho_Yao_Morel_CoviMAC", Source_base);

Sortie& Variation_rho_Yao_Morel::printOn(Sortie& os) const
{
  return os;
}

Entree& Variation_rho_Yao_Morel::readOn(Entree& is)
{
  return is;
}

void Variation_rho_Yao_Morel::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const Champ_base& ch_rho = equation().milieu().masse_volumique();
  const Champ_Inc_base *pch_rho = sub_type(Champ_Inc_base, ch_rho) ? &ref_cast(Champ_Inc_base, ch_rho) : NULL;

  if (!(pch_rho)) return; //rien a faire : le terme est nul

  const Zone_CoviMAC& zone = ref_cast(Zone_CoviMAC, equation().zone_dis().valeur());
  const int ne = zone.nb_elem(), ne_tot = zone.nb_elem_tot(), N = equation().inconnue().valeurs().line_size();

  assert( N == 1 ); // si Ntau > 1 il vaut mieux iterer sur les id_composites des phases turbulentes
  for (auto &&n_m : matrices) if (n_m.first == "interfacialarea" || n_m.first == "temperature" || n_m.first == "pression")
      {
        Matrice_Morse& mat = *n_m.second, mat2;
        const DoubleTab& dep = equation().probleme().get_champ(n_m.first.c_str()).valeurs();
        int nc = dep.dimension_tot(0),
            M  = dep.line_size();
        IntTrav sten(0, 2);
        sten.set_smart_resize(1);
        if (n_m.first == "interfacialarea" || n_m.first == "temperature") // N <= M
          for (int e = 0; e < ne; e++) for (int n = 0; n < N; n++) sten.append_line(N * e + n, M * e + n);
        if (n_m.first == "pression" )
          for (int e = 0; e < ne; e++) for (int n = 0, m = 0; n < N; n++, m+=(M>1)) sten.append_line(N * e + n, M * e + m);
        Matrix_tools::allocate_morse_matrix(N * ne_tot, M * nc, sten, mat2);
        mat.nb_colonnes() ? mat += mat2 : mat = mat2;
      }
}

void Variation_rho_Yao_Morel::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Champ_base& ch_rho = equation().milieu().masse_volumique();
  const Champ_Inc_base *pch_rho = sub_type(Champ_Inc_base, ch_rho) ? &ref_cast(Champ_Inc_base, ch_rho) : NULL;

  if (!(pch_rho)) return; //rien a faire : le terme est nul

  const Zone_CoviMAC& zone = ref_cast(Zone_CoviMAC, equation().zone_dis().valeur());
  const DoubleVect& pe = zone.porosite_elem(), &ve = zone.volumes();

  const DoubleTab& inco = equation().inconnue().valeurs(),
                   &rho = (*pch_rho).valeurs(),
                    &rho_p = (*pch_rho).passe(),
                     &dP_rho = pch_rho->derivees().at("pression"),
                      &dT_rho = pch_rho->derivees().at("temperature");

  double pas_tps = equation().probleme().schema_temps().pas_de_temps();
  int N = ref_cast(Pb_Multiphase, equation().probleme()).nb_phases(), Np = equation().probleme().get_champ("pression").valeurs().line_size();

  Matrice_Morse *Mp = matrices.count("pression")    ? matrices.at("pression")    : NULL,
                 *Mt = matrices.count("temperature") ? matrices.at("temperature") : NULL,
                  *Mai = matrices.count("interfacialarea") ? matrices.at("interfacialarea") : NULL;

  /* elements */
  int n_l = 0 ; // phase porteuse
  for (int e = 0; e < zone.nb_elem(); e++)
    for (int k = 0, mp = 0 ; k<N ; k++ , mp += (Np > 1)) if (k != n_l) //phase gazeuse
        {
          double fac = 2./3.*1/pas_tps * pe(e) * ve(e) ;
          secmem(e , k) += fac * inco(e, k) * ( 1 - rho_p(e, k)/rho(e, k));
          if (Mp) (*Mp)(N * e + k , e) -= fac * inco(e, k) * rho_p(e, k) * -dP_rho(e, k)/(rho(e, k)*rho(e, k));
          if (Mt) (*Mp)(N * e + k , N * e + k) -= fac * inco(e, k) * rho_p(e, k) * -dT_rho(e, k)/(rho(e, k)*rho(e, k));
          if (Mai)(*Mai)(N * e + k, N * e + k) -= fac * ( 1 - rho_p(e, k)/rho(e, k));
        }
}