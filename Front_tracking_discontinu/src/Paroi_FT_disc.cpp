/****************************************************************************
* Copyright (c) 2015, CEA
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
// File:        Paroi_FT_disc.cpp
// Directory:   $TRUST_ROOT/src/Front_tracking_discontinu
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#include <Paroi_FT_disc.h>
#include <Motcle.h>
#include <Equation_base.h>

Implemente_instanciable(Paroi_FT_disc,"Paroi_FT_disc",Symetrie);

Sortie& Paroi_FT_disc::printOn(Sortie& os ) const
{
  return os << que_suis_je() << "\n";
}

// Description:
//  Lecture de l'angle impose dans le jeu de donnees.
//  Format :
//   soit   "constant   champ_frontiere_a_une_composante"
//   soit   "hysteresis champ_frontiere_a_deux_composantes"
//  Les valeurs du champ sont un angle en degres.
//  Convention:
//   L'angle specifie est celui entre l'interface et la paroi
//   du cote de la phase 0.
//               /
//    phase 1   /\        phase 0
//             /  | theta
//    ________/___|________
//
Entree& Paroi_FT_disc::readOn(Entree& is)
{
  Motcle motlu;
  is >> motlu;
  Motcles motcles(3);
  motcles[0] = "constant";
  motcles[1] = "hysteresis";
  motcles[2] = "symetrie";
  int rang = motcles.search(motlu);

  switch(rang)
    {
    case 0:
      type_ = CONSTANT;
      break;
    case 1:
      type_ = HYSTERESIS;
      break;
    case 2:
      type_ = SYMETRIE;
      break;
    default :
      Cerr << "Paroi_FT::readOn: expected a contact angle model.\n You say: "
           << motlu << "\n Valid keywords are: " << motcles << finl;
      exit();
    }

  if (type_ != SYMETRIE)
    is >> le_champ_front;
  else
    {
      le_champ_front.typer("Champ_Front_Uniforme");
      le_champ_front.valeur().fixer_nb_comp(0);
    }

  return is;
}

Paroi_FT_disc::Type_modele Paroi_FT_disc::get_type_modele() const
{
  return type_;
}

// Description: Cette condition aux limites ne support que l'equation
//  Transport_Interfaces_FT_disc
int Paroi_FT_disc::compatible_avec_eqn(const Equation_base& eqn) const
{
  if (eqn.que_suis_je() == "Transport_Interfaces_FT_Disc")
    {
      return 1;
    }
  else
    {
      err_pas_compatible(eqn);
      return 0;
    }
}
