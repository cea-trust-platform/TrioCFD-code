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
// File:        Source_Flottabilite.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/7
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Flottabilite.h>
#include <Probleme_base.h>
#include <Milieu_base.h>

Implemente_instanciable(Source_Flottabilite,"Flottabilite",Source_Action_Particules);


Sortie& Source_Flottabilite::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}


Entree& Source_Flottabilite::readOn(Entree& s )
{
  return s;
}

//L expression codee est precisee dans l entete de la classe
DoubleTab& Source_Flottabilite::ajouter(DoubleTab& resu) const
{
  const DoubleTab& gradient_p = grad_pression();
  ///const DoubleTab& rho_f = rho_fluide();
  const Champ_Don_base& champ_gravite = equation().probleme().milieu().gravite();
  const DoubleTab& gravite = champ_gravite.valeurs();
  const DoubleTab& rho_p = rho_particules();
  const DoubleTab& volume_p = volumes_particules();

  const int dim0 = rho_p .dimension(0);
  const int dim1 =  Objet_U::dimension;
  double masse_particule, terme_pression;

  for (int i=0; i<dim0; i++)
    {
      for (int j=0; j<dim1; j++)
        {
          masse_particule = rho_p(i,0)*volume_p(i,0)*gravite(0,j);
          terme_pression = gradient_p(i,j)*volume_p(i,0);
          resu(i,j) += masse_particule-terme_pression;
        }
    }

  return resu;
}


DoubleTab& Source_Flottabilite::calculer(DoubleTab& resu) const
{
  resu = 0;
  return ajouter(resu);
}
