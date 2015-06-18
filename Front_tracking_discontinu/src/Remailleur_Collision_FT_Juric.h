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
// File:        Remailleur_Collision_FT_Juric.h
// Directory:   $TRUST_ROOT/src/Front_tracking_discontinu
// Version:     /main/11
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Remailleur_Collision_FT_Juric_included
#define Remailleur_Collision_FT_Juric_included

#include <Remailleur_Collision_FT_base.h>

class Maillage_FT_Disc;
class Champ_base;

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    classe Remailleur_Collision_FT_Juric
//     Classe implementant un remailleur d'interfaces entrees en collision et
//     utilisant la reconstruction de Juric sur l'interpolation de indicatrice aux sommets
// .SECTION voir aussi
//     Transport_Interfaces_FT_Disc Maillage_FT_Disc Remailleur_Collision_FT_base
//////////////////////////////////////////////////////////////////////////////

class Remailleur_Collision_FT_Juric : public Remailleur_Collision_FT_base
{

  Declare_instanciable(Remailleur_Collision_FT_Juric);
public:
  int traite_RuptureCoalescenceInterfaces(Maillage_FT_Disc& maillage,
                                          int nb_fa7Intersectees,
                                          const IntTab& couples_fa7Intersectees,
                                          const DoubleTab& segmentsInter_fa7Intersectees,
                                          Champ_base& indicatrice);
  int traite_RuptureCoalescenceInterfaces(Maillage_FT_Disc& maillage,
                                          Champ_base& indicatrice) const;

protected:
  enum Source_Isovaleur { INDICATRICE = 0, FONCTION_DISTANCE = 1};
  // Si Source_Isovaleur==INDICATRICE, on recontruit l'interface
  //   a partir de l'indicatrice (plus robuste)
  // Sinon on utilise la fonction distance (plus precis quand la
  //   topologie est simple mais produit parfois de grosses inclusions
  //   parasites.
  Source_Isovaleur source_isovaleur_;
};

#endif
