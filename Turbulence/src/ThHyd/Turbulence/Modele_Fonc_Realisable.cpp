/****************************************************************************
* Copyright (c) 2019, CEA
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
// File:        Modele_Fonc_Realisable.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Turbulence
// Version:     /main/14
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_Fonc_Realisable.h>
#include <Motcle.h>
#include <Equation_base.h>
#include <Discretisation_base.h>

Implemente_deriv(Modele_Fonc_Realisable_base);
Implemente_instanciable(Modele_Fonc_Realisable,"Modele_Fonc_Realisable",DERIV(Modele_Fonc_Realisable_base));

// XD Modele_Fonc_Realisable Modele_Fonc_Realisable_base Modele_Fonc_Realisable -1 Deriv for instanciation of functions necessary to Realizable K-Epsilon Turbulence Model

//
// printOn et readOn

Sortie& Modele_Fonc_Realisable::printOn(Sortie& s ) const
{
  return  s << valeur().que_suis_je() << finl;
}

Entree& Modele_Fonc_Realisable::readOn(Entree& is )
{
  Motcle typ;
  is >> typ;
  Motcle nom1("Modele_");

  nom1 += typ;

  {
    nom1 += "_";
    Cerr << nom1 << finl;
    Nom discr = equation().discretisation().que_suis_je();
    if (discr=="VEFPreP1B") discr = "VEF";
    nom1 += discr;
  }
  DERIV(Modele_Fonc_Realisable_base)::typer(nom1);
  valeur().associer_eqn(equation());
  valeur().associer(equation().zone_dis(), equation().zone_Cl_dis());
  valeur().associer_pb(equation().probleme());
  is >> valeur();
  return is;
}
