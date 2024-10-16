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
// File:        Connectivites_IndGros.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Geometrie
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#include <Connectivites_IndGros.h>

Implemente_instanciable(Connectivites_IndGros,"Connectivites_IndGros",Connectivites_base);

//// printOn
//

Sortie& Connectivites_IndGros::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

//// readOn
//

Entree& Connectivites_IndGros::readOn(Entree& s )
{
  return s ;
}


void Connectivites_IndGros::calculer_indice_gros(Domaine_VF& domaine_vfF,
                                                 Domaine_VF& domaine_vfG)
{

  Domaine& domaineF = domaine_vfF.domaine();
  Domaine& domaineG = domaine_vfG.domaine();

  int nb_elemF = domaineF.nb_elem();
  int nb_elemG = domaineG.nb_elem();

  int nb_elementF;

  ind_gros.resize(nb_elemG);
  ind_gros=-1;

  for(nb_elementF = 0; nb_elementF<nb_elemF; nb_elementF++)
    {
      int ii = connect_elemF_elemG(nb_elementF);
      if (ii >= 0) ind_gros(ii) = 1;
    }

}



/*! @brief Calcul des connectivites entre elements et elements et entre face et face:
 *
 *     (simple appel aux 2 methodes precedentes)
 *
 * @param (Domaine& domaineG, Domaine& domaineF) domaine discretisee grossiere et fine
 * @throws da
 */
void Connectivites_IndGros::calculer_connectivites(Domaine_VF& domainef,
                                                   Domaine_VF& domaineg,
                                                   Domaine& domg)
{
  Connectivites_base::calculer_connectivites(domainef, domaineg, domg);
  calculer_indice_gros(domainef, domaineg);
}

