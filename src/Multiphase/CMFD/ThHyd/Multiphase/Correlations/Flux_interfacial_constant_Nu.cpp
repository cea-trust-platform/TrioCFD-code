/****************************************************************************
* Copyright (c) 2023, CEA
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

#include <Flux_interfacial_constant_Nu.h>
#include <Pb_Multiphase.h>

Implemente_instanciable(Flux_interfacial_constant_Nu, "Flux_interfacial_constant_Nu", Flux_interfacial_base);

Sortie& Flux_interfacial_constant_Nu::printOn(Sortie& os) const
{
  return os;
}

Entree& Flux_interfacial_constant_Nu::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("Nu", &Nu_, Param::REQUIRED);
  param.lire_avec_accolades_depuis(is);

  return is;
}

void Flux_interfacial_constant_Nu::completer()
{
  Pb_Multiphase *pbm = sub_type(Pb_Multiphase, pb_.valeur()) ? &ref_cast(Pb_Multiphase, pb_.valeur()) : nullptr;

  if (!pbm || pbm->nb_phases() == 1) Process::exit(que_suis_je() + " : not needed for single-phase flow!");
  for (int n = 0; n < pbm->nb_phases(); n++) //recherche de n_l, n_g : phase {liquide,gaz}_continu en priorite
    if (pbm->nom_phase(n).debute_par("liquide") && (n_l < 0 || pbm->nom_phase(n).finit_par("continu")))  n_l = n;

  if (n_l < 0) Process::exit(que_suis_je() + " : liquid phase not found!");
}

void Flux_interfacial_constant_Nu::coeffs(const input_t& in, output_t& out) const
{
  int k, N = out.hi.dimension(0);
  for (k = 0; k < N; k++)
    if (k != n_l)
      {
        // Jusqu'au 19/02/2024
        // int    ind_trav = (k>n_l) ? (n_l*(N-1)-(n_l-1)*(n_l)/2) + (k-n_l-1) : (k*(N-1)-(k-1)*(k)/2) + (n_l-k-1);

        // double ur_ishii_zuber_loc = std::pow(4.*9.81*in.sigma[ind_trav]*(in.rho[n_l]-in.rho[k])/(in.rho[n_l]*in.rho[n_l]),.25);
        // double db_WeberGrav = 2*in.sigma[ind_trav]/(in.rho[n_l]*ur_ishii_zuber_loc*ur_ishii_zuber_loc);
        // double Nu_ = 15.*std::pow(std::max(0.2, .7-in.alpha[k]),-2);
        // double da_Nu_mano = 0.2 > .7-in.alpha[k] ? 0. : 15.*2.*std::pow(.7-in.alpha[k],-3);

        // out.hi(n_l, k) = Nu_ * in.lambda[n_l] / db_WeberGrav * 6 * std::max(in.alpha[k], 1e-3) / db_WeberGrav  ; // std::max() pour que le flux interfacial sont non nul
        // out.da_hi(n_l, k, k) = in.alpha[k] > 1e-3 ?
        //                        Nu_ * in.lambda[n_l] * 6. / (db_WeberGrav *db_WeberGrav ) + da_Nu_mano * in.lambda[n_l] / db_WeberGrav * 6 * in.alpha[k] / db_WeberGrav
        //                        : da_Nu_mano * in.lambda[n_l] / db_WeberGrav * 6 * 1e-3 / db_WeberGrav;

        // Apres 19/02/2024
        int    ind_trav = (k>n_l) ? (n_l*(N-1)-(n_l-1)*(n_l)/2) + (k-n_l-1) : (k*(N-1)-(k-1)*(k)/2) + (n_l-k-1);

        double L_cap = std::sqrt( in.sigma[ind_trav] / (9.81*(in.rho[n_l]-in.rho[k])) );

        out.hi(n_l, k) = Nu_ * in.lambda[n_l] / L_cap * 6 * std::max(in.alpha[k], 1.e-3) / L_cap  ; // std::max() pour que le flux interfacial sont non nul quand le taux de vide tombe
        out.da_hi(n_l, k, k)  = in.alpha[k] > 1.e-3 ? Nu_ * in.lambda[n_l] / L_cap * 6 / L_cap : 0. ; // std::max() pour que le flux interfacial sont non nul quand le taux de vide tombe
        out.hi(k, n_l) = 1e10;
      }
}
