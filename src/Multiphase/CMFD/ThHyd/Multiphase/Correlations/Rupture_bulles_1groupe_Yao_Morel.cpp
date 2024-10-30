/****************************************************************************
* Copyright (c) 2022, CEA
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

#include <Rupture_bulles_1groupe_Yao_Morel.h>
#include <Pb_Multiphase.h>

Implemente_instanciable(Rupture_bulles_1groupe_Yao_Morel, "Rupture_bulles_1groupe_Yao_Morel", Rupture_bulles_1groupe_base);

Sortie& Rupture_bulles_1groupe_Yao_Morel::printOn(Sortie& os) const
{
  return os;
}

Entree& Rupture_bulles_1groupe_Yao_Morel::readOn(Entree& is)
{
  const Pb_Multiphase *pbm = sub_type(Pb_Multiphase, pb_.valeur()) ? &ref_cast(Pb_Multiphase, pb_.valeur()) : nullptr;

  if (!pbm || pbm->nb_phases() == 1) Process::exit(que_suis_je() + " : not needed for single-phase flow!");
  for (int n = 0; n < pbm->nb_phases(); n++) //recherche de n_l, n_g : phase {liquide,gaz}_continu en priorite
    if (pbm->nom_phase(n).debute_par("liquide") && (n_l < 0 || pbm->nom_phase(n).finit_par("continu")))  n_l = n;
  if (n_l < 0) Process::exit(que_suis_je() + " : liquid phase not found!");

  return is;
}

void Rupture_bulles_1groupe_Yao_Morel::coefficient(const DoubleTab& alpha, const DoubleTab& p, const DoubleTab& T,
                                                   const DoubleTab& rho, const DoubleTab& nu, const DoubleTab& sigma, double Dh,
                                                   const DoubleTab& ndv, const DoubleTab& d_bulles,
                                                   const DoubleTab& eps, const DoubleTab& k_turb,
                                                   DoubleTab& coeff) const
{
  int N = alpha.dimension(0);

  for (int k = 0 ; k<N ; k++)
    if (k != n_l) //phase gazeuse
      if (alpha(k) > 1.e-6)
        {

          double We = 2 * rho(n_l) * std::pow(eps(n_l)*d_bulles(k), 2./3.) * d_bulles(k) / sigma(k, n_l) ;

          coeff(k, n_l) = Kb1 *1/(1+Kb2 *alpha(n_l)*std::sqrt(We/We_cr))*std::exp(-std::sqrt(We_cr/We));
          coeff(n_l, k) = coeff(k, n_l);
        }
}

