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
// File      : Source_Transport_Eps_aniso_concen_VEF_Face.cpp
// Directory : $TURBULENCE_ROOT/src/Specializations/VEF/Sources/RANS
//
/////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_Eps_aniso_concen_VEF_Face.h>
#include <Modele_turbulence_hyd_K_Eps_Bicephale.h>
#include <Zone_VEF.h>

Implemente_instanciable_sans_constructeur(Source_Transport_Eps_aniso_concen_VEF_Face,"Source_Transport_Eps_aniso_concen_VEF_P1NC",Source_Transport_Eps_VEF_Face);

Sortie& Source_Transport_Eps_aniso_concen_VEF_Face::printOn(Sortie& s) const { return s << que_suis_je() ; }

Entree& Source_Transport_Eps_aniso_concen_VEF_Face::readOn(Entree& is)
{
  Source_Transport_Eps_VEF_Face::verifier_pb_keps_concen(mon_equation->probleme(),que_suis_je());
  return Source_Transport_Eps_VEF_Face::readOn_concen(is,que_suis_je());
}

void Source_Transport_Eps_aniso_concen_VEF_Face::associer_pb(const Probleme_base& pb)
{
  Source_Transport_Eps_VEF_Face::verifier_milieu_concen(pb,que_suis_je());
  Source_Transport_Eps_VEF_Face::associer_pb(pb);
  Source_Transport_Eps_VEF_Face::associer_pb_concen(pb);
}

void Source_Transport_Eps_aniso_concen_VEF_Face::fill_resu_concen(const DoubleTrav& G, const DoubleVect& volumes_entrelaces, DoubleTab& resu) const
{
  const DoubleTab& K = mon_eq_transport_K->inconnue().valeurs(), &Eps = mon_eq_transport_Eps->inconnue().valeurs();
  double C3_loc, LeK_MIN = mon_eq_transport_Eps->modele_turbulence().get_LeK_MIN();

  for (int face = 0; face < la_zone_VEF->nb_faces(); face++)
    if (K(face) >= LeK_MIN)
      {
        C3_loc = G(face) > 0. ? 0. : C3;
        resu(face) += C1 * (1 - C3_loc) * G(face) * volumes_entrelaces(face) * Eps(face) / K(face);
      }
}

DoubleTab& Source_Transport_Eps_aniso_concen_VEF_Face::ajouter(DoubleTab& resu) const
{
  Source_Transport_Eps_VEF_Face::ajouter(resu); // VB : plutot que de calculer P on appelle ajouter de la classe mere
  return Source_Transport_Eps_VEF_Face::ajouter_concen(resu);
}
