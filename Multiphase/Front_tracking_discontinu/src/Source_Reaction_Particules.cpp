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
// File:        Source_Reaction_Particules.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/7
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Reaction_Particules.h>
#include <Transport_Marqueur_FT.h>
#include <Domaine.h>
#include <Domaine_VF.h>
#include <Dirichlet.h>
#include <Dirichlet_homogene.h>
#include <Milieu_base.h>


Implemente_base(Source_Reaction_Particules,"Two_Way_Coupling",Source_base);


Sortie& Source_Reaction_Particules::printOn(Sortie& os) const
{
  return os;
}

Entree& Source_Reaction_Particules::readOn(Entree& is)
{
  return is;
}

void Source_Reaction_Particules::associer_domaines(const Domaine_dis& zdis,const Domaine_Cl_dis& zcldis)
{
  //Ne fait rien
}

void Source_Reaction_Particules::associer_pb(const Probleme_base& pb)
{
  //Ne fait rien
}

void Source_Reaction_Particules::mettre_a_jour(double temps)
{
  //Ne fait rien
}

DoubleTab& Source_Reaction_Particules::calculer(DoubleTab& resu) const
{
  resu = 0;
  return ajouter(resu);
}

void Source_Reaction_Particules::associer_eq_marqueur(const Equation_base& une_equation)
{
  const Transport_Marqueur_FT& eq = ref_cast(Transport_Marqueur_FT,une_equation);
  eq_marqueur_ = eq;
}
