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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Source_Transport_Realisable_VEF_Face_base.cpp
// Directory : $TURBULENCE_ROOT/src/Specializations/VEF/Sources/RANS/Base
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_Realisable_VEF_Face_base.h>
#include <Modele_Shih_Zhu_Lumley_VEF.h>
#include <Champ_Uniforme.h>
#include <Equation_base.h>
#include <Fluide_base.h>
#include <Champ_P1NC.h>
#include <Domaine_VEF.h>
#include <Debog.h>

Implemente_base_sans_constructeur( Source_Transport_Realisable_VEF_Face_base, "Source_Transport_Realisable_VEF_Face_base", Source_Transport_VEF_Face_base ) ;

Sortie& Source_Transport_Realisable_VEF_Face_base::printOn( Sortie& os ) const { return os << que_suis_je(); }
Entree& Source_Transport_Realisable_VEF_Face_base::readOn(Entree& is) { return Source_Transport_proto::readOn_real(is,que_suis_je()); }

DoubleTab& Source_Transport_Realisable_VEF_Face_base::ajouter_keps_real(DoubleTab& resu) const
{
  Debog::verifier("Source_Transport_Realisable_VEF_Face_base::ajouter resu 0", resu);

  const DoubleTab& visco_turb = get_visc_turb(); // voir les classes filles
  const Modele_Fonc_Realisable_base& mon_modele_fonc = get_modele_fonc(); // voir les classes filles
  const Fluide_base& fluide = ref_cast(Fluide_base, eq_hydraulique->milieu());
  const Champ_Don_base& ch_visco_cin = fluide.viscosite_cinematique();
  const DoubleTab& tab_visco = ch_visco_cin.valeurs();
  int is_visco_const = sub_type(Champ_Uniforme, ch_visco_cin);

  double visco = -1;
  if (is_visco_const) visco = std::max(tab_visco(0, 0), DMINFLOAT);

  const DoubleTab& vit = eq_hydraulique->inconnue().valeurs();
  const DoubleVect& vol_ent = le_dom_VEF->volumes_entrelaces();

  DoubleTab vitesse_filtree(vit);
  ref_cast(Champ_P1NC,eq_hydraulique->inconnue()).filtrer_L2(vitesse_filtree);

  const int nb_faces = le_dom_VEF->nb_faces();
  DoubleTrav P(nb_faces);

  const DoubleTab& CC1 = mon_modele_fonc.get_C1(), &S = mon_modele_fonc.get_S();

  calculer_terme_production_real(vitesse_filtree, visco_turb, P); // voir les classes filles

  Debog::verifier("Source_Transport_Realisable_VEF_Face_base::ajouter P 0", P);
  Debog::verifier("Source_Transport_Realisable_VEF_Face_base::ajouter C1 0", CC1);

  for (int num_face = 0; num_face < nb_faces; num_face++)
    {
      if (!is_visco_const)
        {
          const int elem0 = le_dom_VEF->face_voisins(num_face, 0), elem1 = le_dom_VEF->face_voisins(num_face, 1);
          if (elem1 != -1)
            {
              visco = tab_visco(elem0) * le_dom_VEF->volumes(elem0) + tab_visco(elem1) * le_dom_VEF->volumes(elem1);
              visco /= le_dom_VEF->volumes(elem0) + le_dom_VEF->volumes(elem1);
            }
          else visco = tab_visco(elem0);
        }

      assert(visco > 0.);
      fill_resu_real(num_face, vol_ent, P, CC1, S, visco, resu); // voir les classes filles
    }

  return resu;
}

void Source_Transport_Realisable_VEF_Face_base::contribuer_a_avec(const DoubleTab&, Matrice_Morse& matrice) const
{
  const int size = get_size_k_eps(); // voir les classes filles
  const Fluide_base& fluide = ref_cast(Fluide_base, eq_hydraulique->milieu());
  const Champ_Don_base& ch_visco_cin = fluide.viscosite_cinematique();
  const DoubleTab& tab_visco = ch_visco_cin.valeurs();
  const int is_visco_const = sub_type(Champ_Uniforme, ch_visco_cin);
  double visco = -1;
  if (is_visco_const) visco = std::max(tab_visco(0, 0), DMINFLOAT);

  const Domaine_VEF& domaine_VEF = le_dom_VEF.valeur();
  const DoubleVect& porosite_face = fluide.porosite_face(), &volumes_entrelaces = domaine_VEF.volumes_entrelaces();
  for (int face = 0; face < size; face++)
    {
      if (!is_visco_const)
        {
          const int elem0 = domaine_VEF.face_voisins(face, 0), elem1 = domaine_VEF.face_voisins(face, 1);
          if (elem1 != -1)
            {
              visco = tab_visco(elem0) * domaine_VEF.volumes(elem0) + tab_visco(elem1) * domaine_VEF.volumes(elem1);
              visco /= domaine_VEF.volumes(elem0) + domaine_VEF.volumes(elem1);
            }
          else visco = tab_visco(elem0);
        }

      assert(visco > 0.);
      fill_coeff_matrice(face,porosite_face,volumes_entrelaces,visco,matrice); // voir les classes filles
    }
}
