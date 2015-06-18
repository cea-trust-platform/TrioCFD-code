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
// File:        Remailleur_Collision_FT_base.h
// Directory:   $TRUST_ROOT/src/Front_tracking_discontinu
// Version:     /main/10
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Remailleur_Collision_FT_base_included
#define Remailleur_Collision_FT_base_included

#include <Objet_U.h>
class Maillage_FT_Disc;
class IntTab;
class DoubleTab;
class Champ_base;

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    classe Remailleur_Collision_FT_base
//     Classe de base pour la hierarchie des remailleurs d'interfaces entrees en collision
// .SECTION voir aussi
//     Transport_Interfaces_FT_Disc Maillage_FT_Disc Remailleur_Collision_FT_base
//////////////////////////////////////////////////////////////////////////////

class Remailleur_Collision_FT_base : public Objet_U
{

  Declare_base_sans_constructeur(Remailleur_Collision_FT_base);

public:
  Remailleur_Collision_FT_base();
  virtual int traite_RuptureCoalescenceInterfaces_Conservatif(Maillage_FT_Disc&, Champ_base&);

  virtual int traite_RuptureCoalescenceInterfaces(Maillage_FT_Disc& maillage,
                                                  Champ_base& indicatrice) const=0;

  virtual int traite_RuptureCoalescenceInterfaces(Maillage_FT_Disc& maillage,
                                                  int nb_fa7Intersectees,
                                                  const IntTab& couples_fa7Intersectees,
                                                  const DoubleTab& segmentsInter_fa7Intersectees,
                                                  Champ_base& indicatrice)=0;

  inline int set_impr(int impr);

protected:
  int impr_;
};

inline int Remailleur_Collision_FT_base::set_impr(int impr)
{
  impr_ = impr;
  return impr_;
}

#endif
