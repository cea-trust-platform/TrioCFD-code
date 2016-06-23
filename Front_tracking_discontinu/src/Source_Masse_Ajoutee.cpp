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
// File:        Source_Masse_Ajoutee.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/7
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Masse_Ajoutee.h>
#include <Probleme_base.h>
#include <Schema_Temps_base.h>

Implemente_instanciable_sans_constructeur(Source_Masse_Ajoutee,"Masse_Ajoutee",Source_Action_Particules);

Source_Masse_Ajoutee::Source_Masse_Ajoutee()
{
  C_MA = 0.5;
}

Sortie& Source_Masse_Ajoutee::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}


Entree& Source_Masse_Ajoutee::readOn(Entree& s )
{
  return s;
}

//L expression codee est precisee dans l entete de la classe
DoubleTab& Source_Masse_Ajoutee::ajouter(DoubleTab& resu) const
{
  const DoubleTab& vitesse_f = vitesse_fluide();
  const DoubleTab& rho_f = rho_fluide();
  const DoubleTab& vitesse_p = vitesse_particules();
  const DoubleTab& delta_vit = delta_v();
  const DoubleTab& volume_p = volumes_particules();
  const int& dim0 = vitesse_p.dimension(0);
  const int& dim1 = Objet_U::dimension;
  const Transport_Marqueur_FT& eq = ref_cast(Transport_Marqueur_FT,equation());
  //const double& dt =eq.probleme().schema_temps().pas_de_temps();
  const double& dt = eq.dela_t();
  double derivee;
  const int& impl = eq.resol_implicite();

  if (!impl)
    {
      for (int i=0; i<dim0; i++)
        {
          for (int j=0; j<dim1; j++)
            {

              derivee = (vitesse_f(i,j)-vitesse_p(i,j)-delta_vit(i,j))/dt;
              resu(i,j) += C_MA*rho_f(i)*volume_p(i,0)*derivee;
            }
        }
    }
  else
    {
      for (int i=0; i<dim0; i++)
        {
          for (int j=0; j<dim1; j++)
            {

              derivee = (vitesse_f(i,j)-delta_vit(i,j))/dt;
              resu(i,j) += C_MA*rho_f(i)*volume_p(i,0)*derivee;
            }
        }
    }

  return resu;
}


DoubleTab& Source_Masse_Ajoutee::calculer(DoubleTab& resu) const
{
  resu = 0;
  return ajouter(resu);
}

