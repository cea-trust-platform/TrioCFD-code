/****************************************************************************
* Copyright (c) 2024, CEA
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
// File:        Modele_turbulence_hyd_LES_1elt_VEF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Modeles_Turbulence/LES/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_hyd_LES_1elt_VEF.h>
#include <Schema_Temps_base.h>
#include <Equation_base.h>
#include <Domaine_VEF.h>
#include <Debog.h>

Implemente_instanciable(Modele_turbulence_hyd_LES_1elt_VEF, "Modele_turbulence_hyd_sous_maille_1elt_VEF", Modele_turbulence_hyd_LES_VEF_base);

Sortie& Modele_turbulence_hyd_LES_1elt_VEF::printOn(Sortie& s) const { return s << que_suis_je() << " " << le_nom(); }

Entree& Modele_turbulence_hyd_LES_1elt_VEF::readOn(Entree& s) { return Modele_turbulence_hyd_LES_VEF_base::readOn(s); }

Champ_Fonc& Modele_turbulence_hyd_LES_1elt_VEF::calculer_viscosite_turbulente()
{
  const Domaine_VEF& domaine_VEF = ref_cast(Domaine_VEF, le_dom_VF_.valeur());
  double temps = mon_equation_->inconnue().temps();
  DoubleTab& visco_turb = la_viscosite_turbulente_->valeurs();
  const int nb_elem = domaine_VEF.nb_elem();
  int num_elem;

  F2_.resize(nb_elem);

  calculer_fonction_structure();

  if (visco_turb.size() != nb_elem)
    {
      Cerr << "Size error for the array containing the values of the turbulent viscosity." << finl;
      exit();
    }
  visco_turb = 0.;

  Debog::verifier("Modele_turbulence_hyd_LES_1elt_VEF::calculer_viscosite_turbulente visco_turb 0", visco_turb);
  for (num_elem = 0; num_elem < nb_elem; num_elem++)
    visco_turb(num_elem) = Csm1_ * l_[num_elem] * sqrt(F2_(num_elem));

  Debog::verifier("Modele_turbulence_hyd_LES_1elt_VEF::calculer_viscosite_turbulente visco_turb 1", visco_turb);

  la_viscosite_turbulente_->changer_temps(temps);
  return la_viscosite_turbulente_;
}

void Modele_turbulence_hyd_LES_1elt_VEF::calculer_fonction_structure()
{
  const DoubleTab& la_vitesse = mon_equation_->inconnue().valeurs();
  const Domaine_VEF& domaine_VEF = ref_cast(Domaine_VEF, le_dom_VF_.valeur());
  const int nb_elem = domaine_VEF.nb_elem();
  const IntTab& elem_faces = domaine_VEF.elem_faces();
  const DoubleTab& xv = domaine_VEF.xv();
  const DoubleTab& xp = domaine_VEF.xp();
  const Domaine& domaine = domaine_VEF.domaine();
  int nfac = domaine.nb_faces_elem();

  int num_elem, num_face, i;
  double vit_x_elem, vit_y_elem, vit_z_elem;
  double dist, F2_int1, F2_int2, delta_c_val, dist1;
  static double un_tiers = 1. / 3.;
  double delta_c_calc;

  for (num_elem = 0; num_elem < nb_elem; num_elem++)
    {
      delta_c_val = l_[num_elem];
      // Vitesse au centre de gravite
      vit_x_elem = 0.;
      vit_y_elem = 0.;
      vit_z_elem = 0.;
      dist1 = 0;
      delta_c_calc = 0.;
      for (i = 0; i < nfac; i++)
        {
          num_face = elem_faces(num_elem, i);
          vit_x_elem += la_vitesse(num_face, 0);
          vit_y_elem += la_vitesse(num_face, 1);
          vit_z_elem += la_vitesse(num_face, 2);
          dist1 = (xp(num_elem, 0) - xv(num_face, 0)) * (xp(num_elem, 0) - xv(num_face, 0));
          dist1 += (xp(num_elem, 1) - xv(num_face, 1)) * (xp(num_elem, 1) - xv(num_face, 1));
          dist1 += (xp(num_elem, 2) - xv(num_face, 2)) * (xp(num_elem, 2) - xv(num_face, 2));
          delta_c_calc = std::min(delta_c_calc, dist1);
        }
      vit_x_elem /= (nfac * 1.);
      vit_y_elem /= (nfac * 1.);
      vit_z_elem /= (nfac * 1.);

      // fonction de structure
      F2_int2 = 0.;
      for (i = 0; i < nfac; i++)
        {
          num_face = elem_faces(num_elem, i);

          dist = (xp(num_elem, 0) - xv(num_face, 0)) * (xp(num_elem, 0) - xv(num_face, 0));
          dist += (xp(num_elem, 1) - xv(num_face, 1)) * (xp(num_elem, 1) - xv(num_face, 1));
          dist += (xp(num_elem, 2) - xv(num_face, 2)) * (xp(num_elem, 2) - xv(num_face, 2));

          dist = 1. / dist;
          dist *= (delta_c_val * delta_c_val);
          dist = exp(un_tiers * log(dist));

          F2_int1 = (vit_x_elem - la_vitesse(num_face, 0)) * (vit_x_elem - la_vitesse(num_face, 0));
          F2_int1 += (vit_y_elem - la_vitesse(num_face, 1)) * (vit_y_elem - la_vitesse(num_face, 1));
          F2_int1 += (vit_z_elem - la_vitesse(num_face, 2)) * (vit_z_elem - la_vitesse(num_face, 2));
          F2_int2 += F2_int1 * dist;
        }
      F2_(num_elem) = F2_int2 / nfac;
    }
}
