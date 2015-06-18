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
// File:        Source_Con_Phase_field_base.cpp
// Directory:   $TRUST_ROOT/src/Phase_field
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Con_Phase_field_base.h>


Implemente_base(Source_Con_Phase_field_base,"Source_Con_Phase_field_base",Source_base);


// Description:
//    Surcharge Source_base::printOn(Sortie&)
//    Ecriture d'un probleme sur un flot de sortie.
//    !! Attention n'est pas symetrique de la lecture !!
//    On ecrit les equations, le postraitement et le domaine
//    discretise.
// Precondition:
// Parametre: Sortie& os
//    Signification: flot de sortie
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Sortie&
//    Signification: le flot de sortie modifie
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
Sortie& Source_Con_Phase_field_base::printOn(Sortie& os) const
{
  return Source_base::printOn(os);
}

// Description:
//    Lecture d'un probleme dans un flot d'entree, et ouverture du
//    flot de sauvegarde.
//    Format:
//     {
//     nom_equation bloc de lecture d'une equation
//     Postraitement bloc de lecture postraitement
//     reprise | sauvegarde
//     formatte | binaire
//     nom_de_fichier
//     }
// Precondition:
// Parametre: Entree& is
//    Signification: flot d'entree
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Entree&
//    Signification: le flot d'entre modifie
//    Contraintes:
// Exception: pas d'accolade ouvrante en debut de format
// Exception: mot clef "Postraitement" n'est pas la
// Exception: format de sauvegarde doit etre "binaire" ou "formatte"
// Exception: pas d'accolade fermante en fin de jeu de donnee
// Effets de bord:
// Postcondition: Le pb est associe au postraitement.
Entree& Source_Con_Phase_field_base::readOn(Entree& is)
{
  return Source_base::readOn(is);
}
