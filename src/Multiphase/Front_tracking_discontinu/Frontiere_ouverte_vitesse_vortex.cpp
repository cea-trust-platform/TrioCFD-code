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
// File:        Frontiere_ouverte_vitesse_vortex.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/7
//
//////////////////////////////////////////////////////////////////////////////
#include <Frontiere_ouverte_vitesse_vortex.h>
#include <EChaine.h>
#include <Param.h>
#include <Sous_Domaine.h>
#include <Domaine_VF.h>
#include <Domaine_Cl_dis_base.h>
#include <Probleme_base.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <Domaine.h>

Implemente_instanciable(Frontiere_ouverte_vitesse_vortex,"Frontiere_ouverte_vitesse_vortex",Entree_fluide_vitesse_imposee);

Entree& Frontiere_ouverte_vitesse_vortex::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("sous_domaine", &nom_sous_domaine_, Param::REQUIRED);
  param.ajouter("equation", &nom_equation_, Param::REQUIRED);
  param.ajouter("integrale_reference", &integrale_reference_, Param::REQUIRED);
  param.ajouter("signe", &signe_, Param::REQUIRED);
  param.ajouter("coeff_vitesse", &coeff_vitesse_, Param::REQUIRED);
  param.lire_avec_accolades_depuis(is);

  if (coeff_vitesse_.size_array() != dimension)
    {
      Cerr << "Frontiere_ouverte_vitesse_vortex: wrong dimension for coeff_vitesse" << finl;
      exit();
    }
  Nom n("Champ_Front_Uniforme ");
  n += Nom(dimension);
  for (int i = 0; i < dimension; i++)
    n += " 0.";
  EChaine s(n);
  s >> le_champ_front;
  return is;
}

Sortie& Frontiere_ouverte_vitesse_vortex::printOn(Sortie& os) const
{
  return os;
}

void Frontiere_ouverte_vitesse_vortex::mettre_a_jour(double temps)
{
  // Calcul de l'integrale du champ sur la sous-domaine
  const Domaine_Cl_dis_base& domaine_Cl_dis_base = domaine_Cl_dis();
  const Equation_base& equation = domaine_Cl_dis_base.equation().probleme().get_equation_by_name(nom_equation_);
  const Transport_Interfaces_FT_Disc& eq = ref_cast(Transport_Interfaces_FT_Disc, equation);
  const DoubleTab& indic = eq.inconnue().valeurs();
  const Sous_Domaine& sous_domaine = domaine_Cl_dis_base.domaine().ss_domaine(nom_sous_domaine_);
  const DoubleVect& volume = ref_cast(Domaine_VF, domaine_Cl_dis_base.domaine_dis().valeur()).volumes();

  double integrale = 0.;
  const int n = sous_domaine.nb_elem_tot();
  int i;
  for (i = 0; i < n; i++)
    {
      const int elem = sous_domaine[i];
      const double x = indic(elem);
      const double v = volume(elem);
      integrale += x * v;
    }
  DoubleTab& val = le_champ_front->valeurs();

  double facteur = integrale - integrale_reference_;
  if (facteur * signe_ < 0.)
    facteur = 0.;

  for (i = 0; i < dimension; i++)
    val(0, i) = coeff_vitesse_[i] * facteur;

  le_champ_front.mettre_a_jour(temps);
  Cerr << "Frontiere_ouverte_vitesse_vortex integrale=" << integrale << " vitesse="
       << val(0,0) << " " << val(0,1) << " " << ((dimension==3)?val(0,2):0.) << finl;
}
