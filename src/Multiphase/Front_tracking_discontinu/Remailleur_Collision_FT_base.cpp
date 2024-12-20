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
// File:        Remailleur_Collision_FT_base.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/7
//
//////////////////////////////////////////////////////////////////////////////

#include <Remailleur_Collision_FT_base.h>
#include <SFichier.h>

Implemente_base_sans_constructeur(Remailleur_Collision_FT_base,"Remailleur_Collision_FT_base",Objet_U);



Remailleur_Collision_FT_base::Remailleur_Collision_FT_base() :
  impr_(-1)
{

}

Sortie& Remailleur_Collision_FT_base::printOn(Sortie& os) const
{
  return os;
}

Entree& Remailleur_Collision_FT_base::readOn(Entree& is)
{
  return is;
}

/*! @brief algorithme de remaillage qui tente de conserver le volume.
 *
 * voir class Remailler_Collision_FT_Thomas
 *   Par defaut, appelle le remailleur non conservatif traite_RuptureCoalescenceInterfaces()
 *
 */
int Remailleur_Collision_FT_base::traite_RuptureCoalescenceInterfaces_Conservatif(Maillage_FT_Disc& maillage,
                                                                                  Champ_base& indicatrice)
{
  int ok = traite_RuptureCoalescenceInterfaces(maillage, indicatrice);
  return ok;
}
