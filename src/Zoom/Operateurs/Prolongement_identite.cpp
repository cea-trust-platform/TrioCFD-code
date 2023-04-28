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
// File:        Prolongement_identite.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Operateurs
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#include <Prolongement_identite.h>
#include <Domaine.h>
#include <Domaine_VF.h>

Implemente_instanciable(Prolongement_identite,"Prolongement_identite",Prolongement_base);

//// printOn
//

Sortie& Prolongement_identite::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

//// readOn
//

Entree& Prolongement_identite::readOn(Entree& s )
{
  return s ;
}




/*! @brief Prolongement du centre de gravite des faces grossieres au centre de gravite de TOUTES les faces fines
 *
 */
void Prolongement_identite::prolonger(Domaine_VF& domaine_VFG,
                                      Domaine_VF& domaine_VFF,
                                      const Frontiere& frontF,
                                      IntVect& connect,
                                      const DoubleTab& valG, DoubleTab& tab,
                                      int nb_comp)
{

  //Cerr<<"debut de Prolongement_identite::prolonger"<<finl;

  //Pour chaque element fin,
  //on met directement la valeur de l'element
  //grossier dans lequel il se trouve
  int nb_elemF = connect.size_array();
  int nbelemsF;
  int num_elemG;
  //Pour chaque elemF
  for (nbelemsF=0; nbelemsF<nb_elemF; nbelemsF++)
    {
      num_elemG = connect(nbelemsF);
      tab(nbelemsF) = valG(num_elemG);
    }
  //Cerr<<"fin de Prolongement_identite::prolonger"<<finl;
}




//NE FAIT RIEN
void Prolongement_identite::calculer(Domaine_VF& domainef,
                                     Domaine_VF& domaineg,
                                     IntVect& connect_ff)
{
  //ne fait rien mais c'est normal!!!
}
