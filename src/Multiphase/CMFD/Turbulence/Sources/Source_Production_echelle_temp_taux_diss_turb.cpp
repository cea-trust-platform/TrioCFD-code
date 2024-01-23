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
// File:        Source_Production_echelle_temp_taux_diss_turb.cpp
// Directory:   $TRUST_ROOT/src/Turbulence/Sources
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Production_echelle_temp_taux_diss_turb.h>

#include <Echelle_temporelle_turbulente.h>
#include <Taux_dissipation_turbulent.h>
#include <Pb_Multiphase.h>
#include <Matrix_tools.h>
#include <Array_tools.h>
#include <Domaine_VF.h>

Implemente_base(Source_Production_echelle_temp_taux_diss_turb,"Source_Production_echelle_temp_taux_diss_turb", Sources_Multiphase_base);

Sortie& Source_Production_echelle_temp_taux_diss_turb::printOn(Sortie& os) const
{
  return os;
}

Entree& Source_Production_echelle_temp_taux_diss_turb::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("alpha_omega", &alpha_omega_, Param::REQUIRED);
  param.lire_avec_accolades_depuis(is);
  Cout << "alpha_omega = " << alpha_omega_ << finl ;

  equation().probleme().creer_champ("gradient_vitesse"); // Besoin du gradient de vitesse

  return is;
}

void Source_Production_echelle_temp_taux_diss_turb::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
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
    if (n_m.first == "alpha")
      {
        Matrice_Morse& mat = *n_m.second, mat2;
        const DoubleTab& dep = equation().probleme().get_champ(n_m.first.c_str()).valeurs();
        int nc = dep.dimension_tot(0),
            M  = dep.line_size();
        IntTrav sten(0, 2);
        sten.set_smart_resize(1);
        for (int e = 0; e < ne; e++)
          for (int n = 0; n < Nk; n++)
            if (n < M) sten.append_line(Nk * e + n, M * e + n);
        Matrix_tools::allocate_morse_matrix(Nk * ne_tot, M * nc, sten, mat2);
        mat.nb_colonnes() ? mat += mat2 : mat = mat2;
      }
}
