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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Entree_fluide_vitesse_imposee.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/ALE/src
// Version:     /main/17
//
//////////////////////////////////////////////////////////////////////////////

#include <Entree_fluide_vitesse_imposee_ALE.h>

Implemente_instanciable(Entree_fluide_vitesse_imposee_ALE,"Frontiere_ouverte_vitesse_imposee_ALE",Entree_fluide_vitesse_imposee);
// XD frontiere_ouverte_vitesse_imposee_ALE dirichlet Frontiere_ouverte_vitesse_imposee_ALE -1 Class for velocity boundary condition on a mobile boundary (ALE framework). NL2 The imposed velocity field is vectorial of type Ch_front_input_ALE, Champ_front_ALE or Champ_front_ALE_Beam.  NL2 Example: frontiere_ouverte_vitesse_imposee_ALE Champ_front_ALE 2 0.5*cos(0.5*t) 0.0
// XD  attr ch front_field_base ch 0 Boundary field type.
Sortie& Entree_fluide_vitesse_imposee_ALE::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

/*! @brief Simple appel a: Cond_lim_base::readOn(Entree& )
 *
 * @param (Entree& s) un flot d'entree
 * @return (Entree& s) le flot d'entree modifie
 */
Entree& Entree_fluide_vitesse_imposee_ALE::readOn(Entree& s)
{
  return Entree_fluide_vitesse_imposee::readOn(s);
}
