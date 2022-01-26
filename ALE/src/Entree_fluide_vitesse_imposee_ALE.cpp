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
// File:        Entree_fluide_vitesse_imposee.cpp
// Directory:   $TRUST_ROOT/src/ThHyd
// Version:     /main/17
//
//////////////////////////////////////////////////////////////////////////////

#include <Entree_fluide_vitesse_imposee_ALE.h>

Implemente_instanciable(Entree_fluide_vitesse_imposee_ALE,"Frontiere_ouverte_vitesse_imposee_ALE",Entree_fluide_vitesse_imposee);
//XD Entree_fluide_vitesse_imposee_ALE dirichlet Frontiere_ouverte_vitesse_imposee_ALE -1 Class for velocity boundary condition on a mobile boundary (ALE framework). NL2 To be used when Reichardt's wall law is applied on a moving boundary. NL2 The imposed velocity field is vectorial of type Ch_front_input_ALE or Champ_front_ALE.  NL2 Example: frontiere_ouverte_vitesse_imposee_ALE Champ_front_ALE 2 0.5*cos(0.5*t) 0.0

// Description:
// Class for velocity boundary condition on a mobile boundary (ALE framework)
// To be used when Reichardt's wall law is applied on a moving boundary.
// Allowing the adaptation of the nonlinear Reichardt wall law for the turbulence models to a mobile grid.
// See the mathematical considerations in the reference:
// Computation of unsteady viscous flows around moving bodies using the (k, eps) turbulence model on unstructured dynamic grids,
// Bruno Koobus, Charbel Farhat, Hai Tran, Comput. Methods Appl. Mech. Engrg. 2000, 1441-1466

Sortie& Entree_fluide_vitesse_imposee_ALE::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

// Description:
//    Simple appel a: Cond_lim_base::readOn(Entree& )
// Precondition:
// Parametre: Entree& s
//    Signification: un flot d'entree
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Entree& s
//    Signification: le flot d'entree modifie
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
Entree& Entree_fluide_vitesse_imposee_ALE::readOn(Entree& s)
{
  return Entree_fluide_vitesse_imposee::readOn(s);
}
