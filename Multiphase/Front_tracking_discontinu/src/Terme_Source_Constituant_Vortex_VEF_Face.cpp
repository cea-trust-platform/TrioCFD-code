/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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
// File:        Terme_Source_Constituant_Vortex_VEF_Face.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/7
//
//////////////////////////////////////////////////////////////////////////////
#include <Terme_Source_Constituant_Vortex_VEF_Face.h>
#include <Param.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <communications.h>
#include <Probleme_base.h>

Implemente_instanciable(Terme_Source_Constituant_Vortex_VEF_Face, "Source_Constituant_Vortex_VEF_P1NC", Terme_Source_Constituant_VEF_Face);

Entree& Terme_Source_Constituant_Vortex_VEF_Face::readOn(Entree& is)
{
  Cerr << "Reading Terme_Source_Constituant_Vortex_VEF_Face\n"
       << " (parameters for expression are (x,y,z,u,v,w) where u,v,w is the interface sensor position)" << finl;
  senseur_.associer_pb(equation().probleme());
  Param param(que_suis_je());
  integrale_ = -1.;
  param.ajouter("rayon_spot", &rayon_spot_, Param::REQUIRED);
  param.ajouter("delta_spot", &delta_spot_, Param::REQUIRED);
  param.ajouter("integrale", &integrale_, Param::REQUIRED);
  param.ajouter("debit", &debit_, Param::REQUIRED);
  param.ajouter("senseur_interface", &senseur_, Param::REQUIRED);

  param.lire_avec_accolades_depuis(is);

  const int dim = dimension;

  if (delta_spot_.size_array() != dim || rayon_spot_ <= 0.)
    {
      Cerr << "Error: wrong segments dimension or nb_points_tests" << finl;
      exit();
    }
  return is;
}

Sortie& Terme_Source_Constituant_Vortex_VEF_Face::printOn(Sortie& os) const
{
  return os;
}

void Terme_Source_Constituant_Vortex_VEF_Face::associer_pb (const Probleme_base& pb)
{
}

DoubleTab& Terme_Source_Constituant_Vortex_VEF_Face::calculer(DoubleTab& tab) const
{
  const Equation_base& eq = equation();
  const Zone_VEF& zone_vef = ref_cast(Zone_VEF, eq.zone_dis().valeur());
  const Zone_Cl_VEF& zone_cl_vef = ref_cast(Zone_Cl_VEF, eq.zone_Cl_dis().valeur());
  const DoubleVect& volumes = zone_vef.volumes_entrelaces();
  const DoubleVect& volumes_cl = zone_cl_vef.volumes_entrelaces_Cl();
  const int premiere_face_std = zone_vef.premiere_face_std();
  const DoubleVect& porosites = zone_vef.porosite_face();
  const DoubleTab& xv = zone_vef.xv();  // Coord des centres des faces
  const int nb_faces = zone_vef.nb_faces();
  const int dim = Objet_U::dimension;
  const ArrOfInt& faces_doubles = zone_vef.faces_doubles();

  // Pour chaque colonne de xv, a quelle variable faut-il attribuer cette valeur ?
  // On met dans expr le terme source et on calcule l'integrale
  double integrale = 0.;
  const double omega = 4. / (rayon_spot_ * rayon_spot_);
  for (int i_face = 0; i_face < nb_faces; i_face++)
    {
      double f;
      {
        double x = xv(i_face, 0) - position_spot_[0];
        double y = xv(i_face, 1) - position_spot_[1];
        double z = (dim == 3) ? (xv(i_face, 2) - position_spot_[2]) : 0.;
        f = exp((-x*x-y*y-z*z) * omega);
      }
      double volume = (i_face < premiere_face_std) ? volumes_cl[i_face] : volumes[i_face];
      double por = porosites(i_face);
      double x = f * volume * por;
      tab(i_face) = x;
      double facteur = 1.;
      if (faces_doubles[i_face])
        facteur = 0.5;
      integrale += x * facteur;
    }

  integrale = mp_sum(integrale);
  if (integrale_ > 0. && integrale > 0.)
    {
      const double facteur = integrale_ / integrale;
      tab *= facteur;
    }
  return tab;
}

DoubleTab& Terme_Source_Constituant_Vortex_VEF_Face::ajouter(DoubleTab& tab) const
{
  const int n = tab.dimension(0);
  // On doit de toutes facons creer un tableau temporaire pour pouvoir renormaliser la source.
  // tableau tmp de la bonne taille mais non initialise:
  DoubleTab tmp;
  tmp.set_smart_resize(1);
  tmp.resize(n);
  calculer(tmp);
  for (int i = 0; i < n; i++)
    tab(i) += tmp(i);
  return tab;
}

void Terme_Source_Constituant_Vortex_VEF_Face::mettre_a_jour(double temps)
{
  senseur_.calculer_position(position_spot_);
  position_spot_ += delta_spot_;
  Cerr << "Source_Constituant_Vortex "
       << position_spot_[0] << " " << position_spot_[1] << " "
       << ((position_spot_.size_array()==3)?position_spot_[2]:0.)
       << finl;
}

void Terme_Source_Constituant_Vortex_VEF_Face::ajouter_terme_div_u(DoubleVect& secmem_pression, double dt) const
{
  const Equation_base& eq = equation();
  const Zone_VEF& zone_vef = ref_cast(Zone_VEF, eq.zone_dis().valeur());
  const DoubleVect& volumes = zone_vef.volumes();
  const DoubleTab& xp = zone_vef.xp();  // Centres des elements
  const int nb_elem = zone_vef.zone().nb_elem();
  const int dim = xp.dimension(1);

  // GF mieux comme cela
  const int n = secmem_pression.size_totale();
  ArrOfDouble tmp;
  tmp.set_smart_resize(1);
  tmp.resize_array(n);
  double integrale = 0.;
  const double omega = 4. / (rayon_spot_ * rayon_spot_);
  int elem;
  for (elem = 0; elem < nb_elem; elem++)
    {
      double x = xp(elem, 0) - position_spot_[0];
      double y = xp(elem, 1) - position_spot_[1];
      double z = (dim == 3) ? (xp(elem, 2) - position_spot_[2]) : 0.;
      double f = exp((-x*x-y*y-z*z) * omega) * volumes[elem];
      integrale += f;
      tmp[elem] = f;
    }
  integrale = mp_sum(integrale);
  double facteur = 0.;
  if (integrale > 0.)
    // Le facteur 1./dimension est du a l'operateur div en vef
    facteur = debit_ / Objet_U::dimension / integrale / dt;

  Cerr << "Terme_Source_Constituant_Vortex_VEF_Face::ajouter_terme_div_u " << integrale << finl;
  for (elem = 0; elem < nb_elem; elem++)
    secmem_pression[elem] += tmp[elem] * facteur;
}
