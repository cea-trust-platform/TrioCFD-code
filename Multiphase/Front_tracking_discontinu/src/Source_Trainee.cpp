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
// File:        Source_Trainee.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Trainee.h>
#include <Probleme_base.h>
#include <Fluide_Diphasique.h>

Implemente_instanciable(Source_Trainee,"Trainee",Source_Action_Particules);


Sortie& Source_Trainee::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}


Entree& Source_Trainee::readOn(Entree& s )
{
  return s;
}

//L expression codee est precisee dans l entete de la classe
DoubleTab& Source_Trainee::ajouter(DoubleTab& resu) const
{

  const DoubleTab& vitesse_f = vitesse_fluide();
  const DoubleTab& rho_f = rho_fluide();
  const DoubleTab& visco_dyn_f = visco_dyn_fluide();
  const DoubleTab& vitesse_p = vitesse_particules();
  const DoubleTab& diametre_p = diametre_particules();
  const DoubleTab& rho_p = rho_particules();

  double grav = 0.;
  double sigma = 0.;
  double Cd_p, Cd_diph;
  Cd_p = 0.;
  Cd_diph = 0.;

  const Equation_base& eq_ns = equation().probleme().equation(0);
  const Milieu_base& mil = eq_ns.milieu();
  if (sub_type(Fluide_Diphasique,mil))
    {
      const Fluide_Diphasique& fluide = ref_cast(Fluide_Diphasique,mil);
      if (fluide.a_gravite())
        {
          const DoubleTab& gravite = fluide.gravite().valeurs();
          assert((gravite.nb_dim() == 2) && (gravite.dimension(1) == dimension));
          for (int j=0; j<dimension; j++)
            grav += gravite(0,j)*gravite(0,j);
          grav = sqrt(grav);
        }

      ////grav=9.81;

      sigma = fluide.sigma();
    }

  const int& dim0 =  vitesse_p.dimension(0);
  const int& dim1 =  vitesse_p.dimension(1);
  assert(vitesse_p.dimension(0)==vitesse_f.dimension(0));

  double norme_delta_v, Reynolds_p, Cd, Surface, pi;
  pi = M_PI;

  DoubleVect deltav(Objet_U::dimension);

  for (int i=0; i<dim0; i++)
    {
      norme_delta_v = 0.;
      for (int j=0; j<dim1; j++)
        {
          deltav[j] = vitesse_p(i,j)-vitesse_f(i,j);
          norme_delta_v += deltav[j]*deltav[j];
        }
      norme_delta_v = sqrt(norme_delta_v);

      Reynolds_p = rho_f(i)*norme_delta_v*diametre_p(i,0)/visco_dyn_f(i);
      Cd_p =24/std::max(1e-12,Reynolds_p);

      if (sub_type(Fluide_Diphasique,mil))
        {
          Cd_diph = std::min((2./3.)*sqrt((diametre_p(i,0)*diametre_p(i,0))*grav*(std::fabs(rho_f(i)-rho_p(i,0)))/sigma),8./3.);
        }

      Cd = std::max(Cd_p,Cd_diph);

      Surface = (dim1==3?pi*diametre_p(i,0)*diametre_p(i,0)/4.:diametre_p(i,0));

      for (int j=0; j<dim1; j++)
        resu(i,j) += -0.5*Cd*rho_f(i)*Surface*norme_delta_v*deltav[j];
    }

  return resu;
}


DoubleTab& Source_Trainee::calculer(DoubleTab& resu) const
{
  resu = 0;
  return ajouter(resu);
}


