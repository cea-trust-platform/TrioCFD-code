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
// File:        Connectivites_IndGros.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Geometrie
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Connectivites_IndGros_included
#define Connectivites_IndGros_included


#include <Domaine_VF.h>
#include <Connectivites_base.h>

#include <Domaine.h>
/*! @brief class Connectivites_IndGros
 *
 */

//////////////////////////////////////////////////////////////////////////////////
//
// CLASS: Connectivites_IndGros
// Cette classe etend la classe de Connectivites_base. Elle fournit en plus un tableau
// indiquant si un element grossier se trouve par dessus un element fin
/////////////////////////////////////////////////////////////////////////////////

class Connectivites_IndGros : public Connectivites_base
{
  Declare_instanciable(Connectivites_IndGros);

public:
  void calculer_connectivites(Domaine_VF& domainef, Domaine_VF& domaineg,
                              Domaine& domg) override;
  inline IntVect& indice_gros(); // tableau d'entier de taille NELEM_GROS : 1 si l'elem grossier recouvre un elem fin, 0 sinon.


private:
  IntVect ind_gros; // indique si un element grossier recouvre un element fin

  void calculer_indice_gros(Domaine_VF& domainef, Domaine_VF& domaineg);

};

// Fonctions inline



inline IntVect& Connectivites_IndGros::indice_gros()
{
  return ind_gros;
}


#endif
