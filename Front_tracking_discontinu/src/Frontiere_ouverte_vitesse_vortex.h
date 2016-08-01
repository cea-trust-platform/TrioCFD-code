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
// File:        Frontiere_ouverte_vitesse_vortex.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/6
//
//////////////////////////////////////////////////////////////////////////////
#ifndef Frontiere_ouverte_vitesse_vortex_included
#define Frontiere_ouverte_vitesse_vortex_included

#include <Entree_fluide_vitesse_imposee.h>

class Frontiere_ouverte_vitesse_vortex : public Entree_fluide_vitesse_imposee
{
  Declare_instanciable(Frontiere_ouverte_vitesse_vortex);
public:
  void mettre_a_jour(double temps);
protected:
  // Sous-zone sur laquelle on veut calculer "integrale"
  Nom nom_sous_zone_;
  // Equation interface
  Nom nom_equation_;
  // La vitesse de sortie est (integrale-integrale_reference)*coeff_vitesse_
  //  si (integrale-integrale_reference) est du meme signe que "signe_",
  //  sinon la vitesse est nulle
  double integrale_reference_;
  int signe_;
  ArrOfDouble coeff_vitesse_;
};
#endif
