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
// File:        Senseur_Interface.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/1
//
//////////////////////////////////////////////////////////////////////////////
#include <Senseur_Interface.h>
#include <Champ_base.h>
#include <TRUSTTab.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <Domaine.h>
#include <Domaine_VF.h>
#include <Probleme_base.h>
#include <Param.h>

Implemente_instanciable(Senseur_Interface, "Senseur_Interface", Objet_U);

Entree& Senseur_Interface::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("equation_interface", &equation_interface_, Param::REQUIRED);
  param.ajouter("segment_senseur_1", &segment_senseur_1_, Param::REQUIRED);
  param.ajouter("segment_senseur_2", &segment_senseur_2_, Param::REQUIRED);
  param.ajouter("nb_points_tests", &nb_points_tests_, Param::REQUIRED);
  param.lire_avec_accolades_depuis(is);

  const int dim = dimension;

  if (segment_senseur_1_.size_array() != dim
      || segment_senseur_2_.size_array() != dim
      || nb_points_tests_ < 2)
    {
      Cerr << "Error: wrong segments dimension or nb_points_tests" << finl;
      exit();
    }
  equation_ = ref_cast(Transport_Interfaces_FT_Disc, probleme_->get_equation_by_name(equation_interface_));

  return is;
}

Sortie& Senseur_Interface::printOn(Sortie& os) const
{
  return os;
}

int Senseur_Interface::calculer_position(ArrOfDouble& pos) const
{
  int i, j;
  int ok = 0;

  const Transport_Interfaces_FT_Disc& eq = ref_cast(Transport_Interfaces_FT_Disc, equation_.valeur());
  const Champ_base& indicatrice = eq.inconnue();
  // Recherche de la position de l'interface:
  const int dim = segment_senseur_1_.size_array();
  DoubleTab positions(nb_points_tests_, dim);
  const double facteur = 1. / (double)(nb_points_tests_-1.);
  for (j = 0; j < dim; j++)
    {
      const double x1 = segment_senseur_1_[j];
      const double x2 = segment_senseur_2_[j];
      for (i = 0; i < nb_points_tests_; i++)
        {
          double t = ((double) i) * facteur;
          positions(i,j) = x2 * t + x1 * (1. - t);
        }
    }
  IntVect num_element(nb_points_tests_);
  const Domaine& domaine = eq.domaine_dis().domaine();
  domaine.chercher_elements(positions, num_element);
  //DoubleVect dist(nb_points_tests_);
  // Invalide le numero d'element pour les elements virtuels
  const int nb_elem = domaine.nb_elem();
  for (i = 0; i < nb_points_tests_; i++)
    {
      const int elem = num_element[i];
      if (elem >= 0 && elem < nb_elem)
        {
          const double ind = indicatrice.valeurs()(elem);
          if (ind > 0. && ind < 1.) // On a trouve l'interface !
            break;
        }
    }
  double t_injection; // Entre 0 et 1, coordonnee de l'interface le long du segment
  const DoubleTab& xp = ref_cast(Domaine_VF, eq.domaine_dis()).xp();
  const DoubleTab& normale_interf = eq.get_update_normale_interface().valeurs();
  const DoubleTab& distance_interf = eq.get_update_distance_interface().valeurs();
  if (i < nb_points_tests_)
    {
      const int elem = num_element[i];
      double s1s2_scal_n = 0.; // produit scalaire entre la normale unitaire a l'interface et le vecteur "segment"
      double norme2_s1s2 = 0.;
      double s1s2_scal_s1elem = 0.;
      for (j = 0; j < dim; j++)
        {
          double s1s2 = segment_senseur_2_[j] - segment_senseur_1_[j];
          s1s2_scal_n += s1s2 * normale_interf(elem, j);
          norme2_s1s2 += s1s2 * s1s2;
          s1s2_scal_s1elem += s1s2 * (xp(elem, j) - segment_senseur_1_[j]);
        }
      t_injection = s1s2_scal_s1elem / norme2_s1s2 - distance_interf(elem) / s1s2_scal_n;
    }
  else
    {
      t_injection = 1.e30;
    }
  // Quel processeur a trouve une maille avec l'interface ?
  t_injection = mp_min(t_injection);
  if (t_injection < 0.99e30)
    {
      pos.resize_array(dim);
      for (i = 0; i < dim; i++)
        {
          const double x1 = segment_senseur_1_[i];
          const double x2 = segment_senseur_2_[i];
          const double x = x2 * t_injection + x1 * (1. - t_injection);
          pos[i] = x;
          ok = 1;
        }
    }
  else
    {
      pos = segment_senseur_1_;
      ok = 0;
    }
  return ok;
}
