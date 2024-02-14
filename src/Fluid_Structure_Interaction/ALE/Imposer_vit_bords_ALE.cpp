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
// File:        Imposer_vit_bords_ALE.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/ALE/src
// Version:     /main/10
//
//////////////////////////////////////////////////////////////////////////////

#include <Imposer_vit_bords_ALE.h>
#include <Domaine_ALE.h>


Implemente_instanciable(Imposer_vit_bords_ALE,"Imposer_vit_bords_ALE",Interprete_geometrique_base);
// XD  imposer_vit_bords_ale interprete imposer_vit_bords_ale 0 For the Arbitrary Lagrangian-Eulerian framework: block to indicate the number of mobile boundaries of the domain and specify the speed that must be imposed on them.
// XD  attr dom ref_domaine dom 0 Name of domain.
// XD attr bloc bloc_lecture bloc 0 between the braces, you must specify the numbers of the mobile borders of the domain then list these mobile borders and indicate the speed which must be imposed on them  NL2 Example:  Imposer_vit_bords_ALE dom_name  { 1  boundary_name  Champ_front_ALE 2 -(y-0.1)*0.01 (x-0.1)*0.01 }





/*! @brief Simple appel a: Interprete::printOn(Sortie&)
 *
 *     Imprime l'interprete sur un flot de sortie
 *
 * @param (Sortie& os) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Imposer_vit_bords_ALE::printOn(Sortie& os) const
{
  return Interprete::printOn(os);
}


/*! @brief Simple appel a: Interprete::readOn(Entree&)
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Imposer_vit_bords_ALE::readOn(Entree& is)
{
  return Interprete::readOn(is);
}


/*! @brief Fonction principale de l'interprete Imposer_deplacements_bords_ALE: discretiser un probleme avec une discretisation
 *
 *     On controle dynamiquement le type du probleme et
 *     de la discretisation que l'on a lu, on verifie
 *     que le schema en temps du probleme a ete lu et on
 *     discretise le probleme avec la discretisation.
 *
 * @throws discretisation non reconnue
 * @throws on ne sait pas discretiser car le probleme n'est pas d'un
 * sous-type reconnu (Probleme,Probleme_base,Probleme_Couple)
 */
Entree& Imposer_vit_bords_ALE::interpreter_(Entree& is)
{
  associer_domaine(is);
  Nom dom_name = domaine().que_suis_je();
  if (!dom_name.contient("ALE"))
    {
      Cerr <<"Mobile domain:  replace  " <<domaine().que_suis_je()<<" with "<< domaine().que_suis_je() <<"_ALE and restart!"<< finl;
      exit();
    }
  Domaine_ALE& dom=ref_cast(Domaine_ALE, domaine());
  dom.reading_vit_bords_ALE(is);
  return is;
}
