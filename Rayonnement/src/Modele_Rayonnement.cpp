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
// File:        Modele_Rayonnement.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement/src
// Version:     /main/10
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_Rayonnement.h>

Implemente_deriv(Modele_Rayonnement_base);
Implemente_instanciable(Modele_Rayonnement,"Modele_Rayonnement",DERIV(Modele_Rayonnement_base));

Entree& Modele_Rayonnement::readOn(Entree& is)
{

  return is;
}

Sortie& Modele_Rayonnement::printOn(Sortie& os) const
{

  return os;
}

void Modele_Rayonnement::preparer_calcul()
{
  valeur().preparer_calcul();
}

void Modele_Rayonnement::mettre_a_jour(double temps)
{

  valeur().mettre_a_jour(temps);
}

void Modele_Rayonnement::calculer_temperatures()
{
  valeur().calculer_temperatures();
}

void Modele_Rayonnement::calculer_radiosites()
{
  valeur().calculer_radiosites();
}

void Modele_Rayonnement::calculer_flux_radiatifs()
{
  valeur().calculer_flux_radiatifs();
}

void Modele_Rayonnement::imprimer_flux_radiatifs(Sortie& os) const
{
  valeur().imprimer_flux_radiatifs(os);
}

