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
// File:        Associer_pbMG_pbFin.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Noyau
// Version:     /main/6
//
//////////////////////////////////////////////////////////////////////////////

#include <Associer_pbMG_pbFin.h>

Implemente_instanciable(Associer_pbMG_pbFin,"Associer_pbMG_pbFin",Interprete);


/*! @brief Simple appel a: Interprete::printOn(Sortie&)
 *
 *     Imprime l'interprete sur un flot de sortie
 *
 * @param (Sortie& os) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Associer_pbMG_pbFin::printOn(Sortie& os) const
{
  return Interprete::printOn(os);
}


/*! @brief Simple appel a: Interprete::readOn(Entree&)
 *
 *     Ecrit l'interprete sur un flot d'entree.
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Associer_pbMG_pbFin::readOn(Entree& is)
{
  return Interprete::readOn(is);
}

/*! @brief Fonction principale de l'interprete Associer_pbMG_pbFin: associer deux objets.
 *
 *     On associe 1 a 2
 *
 * @param (Entree& is) un flot d'entree, a partir duquel on lit les noms des objets a associer
 * @return (Entree&) le flot d'entree modifie
 * @throws on ne sait pas associer obj1 et obj2
 */
Entree& Associer_pbMG_pbFin::interpreter(Entree& is)
{
  // Acquisition des parametres :
  Nom nom1, nom2;
  is >> nom1 >> nom2;
  /*   if(!is){
       Cerr << "Probleme de lecture dans Associer_pbMG_pbGglobal << finl";
       exit();
       }
  */
  Objet_U& ob1=objet(nom1);
  Objet_U& ob2=objet(nom2);
  Pb_MG& pb1=ref_cast(Pb_MG, ob1);
  Probleme_base& pb2=ref_cast(Probleme_base, ob2);
  // if (!pb1.associer_pbMG_pbGglobal_(pb2))
  //       {
  //         Cerr << "On ne sait pas associer " << pb1.que_suis_je()
  //              << " et " << pb2.que_suis_je()  << finl;
  //         exit();
  //       }
  pb1.associer_pbMG_pbFin_(pb2);
  return is;
}
