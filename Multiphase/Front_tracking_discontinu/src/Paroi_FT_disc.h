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
// File:        Paroi_FT_disc.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Paroi_FT_disc_included
#define Paroi_FT_disc_included

#include <Symetrie.h>

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//   classe Paroi_FT_disc
//   Condition aux limites d'angle de contact pour l'equation
//    "Transport_Interfaces_FT_disc"
//   Le champ frontiere contient une valeur d'angle de contact exprimee
//   en degres (voir methode readOn)
//   Attention: l'heritage de "Symetrie" est necessaire pour
//    - le calcul du gradient de l'indicatrice (nul au bord)
//    - le calcul du champ indicatrice_som du fichier lml pour
//      les cas-tests de non regression
// .SECTION voir aussi
//    Symetrie
//////////////////////////////////////////////////////////////////////////////

class Paroi_FT_disc : public Symetrie
{
  Declare_instanciable(Paroi_FT_disc);
public:
  // Si l'angle est "constant" , le champ frontiere doit avoir une composante
  // Si l'angle est "hysteresis", le champ frontiere doit avoir deux composantes
  enum Type_modele { CONSTANT = 0, HYSTERESIS = 1, SYMETRIE = 2 };

  int         compatible_avec_eqn(const Equation_base&) const override;
  Type_modele    get_type_modele() const;

protected:
  Type_modele type_;
};

#endif
