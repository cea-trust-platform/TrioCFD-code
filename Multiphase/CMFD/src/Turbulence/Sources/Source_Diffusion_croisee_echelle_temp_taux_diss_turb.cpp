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
// File:        Source_Diffusion_croisee_echelle_temp_taux_diss_turb.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/Sources
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Diffusion_croisee_echelle_temp_taux_diss_turb.h>

#include <Neumann_loi_paroi_faible_tau_omega.h>
#include <Neumann_loi_paroi_faible_k.h>
#include <Domaine_Cl_PolyMAC.h>
#include <Pb_Multiphase.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <Domaine_VF.h>

Implemente_base(Source_Diffusion_croisee_echelle_temp_taux_diss_turb,"Source_Diffusion_croisee_echelle_temp_taux_diss_turb", Sources_Multiphase_base);

Sortie& Source_Diffusion_croisee_echelle_temp_taux_diss_turb::printOn(Sortie& os) const
{
  return os;
}

Entree& Source_Diffusion_croisee_echelle_temp_taux_diss_turb::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("sigma_d", &sigma_d);
  param.lire_avec_accolades_depuis(is);
  return is;
}

void Source_Diffusion_croisee_echelle_temp_taux_diss_turb::completer()
{
  const Pb_Multiphase& pb = ref_cast(Pb_Multiphase,  equation().probleme());

  for (int i = 0 ; i <pb.nombre_d_equations() ; i++)
    for (int j = 0 ; j<pb.equation(i).domaine_Cl_dis()->nb_cond_lim(); j++)
      {
        const Cond_lim& cond_lim_loc = pb.equation(i).domaine_Cl_dis()->les_conditions_limites(j);
        if      sub_type(Neumann_loi_paroi_faible_k, cond_lim_loc.valeur())         f_grad_k_fixe = 0;
        else if sub_type(Neumann_loi_paroi_faible_tau_omega, cond_lim_loc.valeur()) f_grad_tau_omega_fixe = 0;
      }

}

void Source_Diffusion_croisee_echelle_temp_taux_diss_turb::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const Domaine_VF& 		domaine 		= ref_cast(Domaine_VF, equation().domaine_dis().valeur());
  const int N = equation().inconnue().valeur().valeurs().line_size(), nb_elem = domaine.nb_elem();
  int e, n;

  assert(N == 1); // si N > 1 il vaut mieux iterer sur les id_composites des phases turbulentes

  for (auto &&n_m : matrices)
    if (n_m.first == "alpha" || n_m.first == "temperature" || n_m.first == "pression")
      {
        Matrice_Morse& mat = *n_m.second, mat2;
        const DoubleTab& dep = equation().probleme().get_champ(n_m.first.c_str()).valeurs();
        int m,
            nc = dep.dimension_tot(0),	// nombre d'elements total
            M  = dep.line_size();		// nombre de composantes
        IntTrav sten(0, 2);
        sten.set_smart_resize(1);
        if (n_m.first == "alpha" || n_m.first == "temperature")	// N <= M
          for (e = 0; e < nb_elem; e++)
            for (n = 0; n < N; n++) sten.append_line(N * e + n, M * e + n);
        if (n_m.first == "pression")
          for (e = 0; e < nb_elem; e++)
            for (n = 0; n < N; n++)
              for (m = 0; m<M; m++) sten.append_line(N * e + n, M * e + m);
        Matrix_tools::allocate_morse_matrix(N * domaine.nb_elem_tot(), M * nc, sten, mat2);
        mat.nb_colonnes() ? mat += mat2 : mat = mat2;
      }
}
