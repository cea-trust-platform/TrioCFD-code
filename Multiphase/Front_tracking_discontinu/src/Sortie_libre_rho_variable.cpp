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
// File:        Sortie_libre_rho_variable.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/7
//
//////////////////////////////////////////////////////////////////////////////
#include <Sortie_libre_rho_variable.h>

Implemente_instanciable(Sortie_libre_rho_variable,"Sortie_libre_rho_variable",Sortie_libre_pression_imposee);

Entree& Sortie_libre_rho_variable::readOn(Entree& is)
{
  Sortie_libre_pression_imposee::readOn(is);
  return is;
}

Sortie& Sortie_libre_rho_variable::printOn(Sortie& os) const
{
  Sortie_libre_pression_imposee::printOn(os);
  return os;
}

/*! @brief Renvoie la valeur de la pression imposee.
 *
 * @param (i) numero de la face de la cond_lim ou on cherche la pression Contrainte: 0 <= i <
 */
double Sortie_libre_rho_variable::flux_impose(int i) const
{
  const DoubleTab& tab = le_champ_front.valeurs();
  double pimpose = 0.;

  if (tab.size()==1)    // Champ uniforme ?
    pimpose = tab(0,0);
  else
    pimpose = tab(i,0);
  return pimpose;
}

/*! @brief Idem avec plusieurs composantes.
 *
 * .. Qu'est-ce que ca peut bien vouloir dire ? L'appel a cette methode produit une erreur.
 *
 */
double Sortie_libre_rho_variable::flux_impose(int i, int j) const
{
  Cerr << "Sortie_libre_rho_variable::flux_impose(int i, int j)\n ";
  Cerr << " Qu'est-ce que ca veut dire ?" << finl;
  assert(0);
  exit();
  return 0.;
}

void Sortie_libre_rho_variable::completer()
{
}
