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
// File:        Schema_Euler_explicite_ALE.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/ALE/src/New
//
//////////////////////////////////////////////////////////////////////////////

#include <Schema_Euler_explicite_ALE.h>
#include <Equation_base.h>
#include <Debog.h>
#include <Probleme_base.h>
#include <Domaine.h>
#include <Domaine_ALE.h>

Implemente_instanciable(Schema_Euler_explicite_ALE,"Schema_euler_explicite_ALE|Scheme_euler_explicit_ALE",Schema_Euler_explicite);
// XD scheme_euler_explicite_ALE schema_temps_base schema_euler_explicite_ALE -1 This is the Euler explicit scheme used for ALE problems.

Sortie& Schema_Euler_explicite_ALE::printOn(Sortie& s) const
{
  return  Schema_Euler_explicite::printOn(s);
}


Entree& Schema_Euler_explicite_ALE::readOn(Entree& s)
{
  return Schema_Euler_explicite::readOn(s) ;
}



int Schema_Euler_explicite_ALE::faire_un_pas_de_temps_eqn_base(Equation_base& eqn)
{

  DoubleTab& present = eqn.inconnue().valeurs(); // Un
  DoubleTab& futur   = eqn.inconnue().futur();   // Un+1
  DoubleTab dudt(futur);
  Debog::verifier("Schema_Euler_explicite_ALE::faire_un_pas_de_temps_eqn_base -futur avant", futur);

  // Boundary conditions applied on Un+1:
  eqn.domaine_Cl_dis().imposer_cond_lim(eqn.inconnue(),temps_courant()+pas_de_temps());

  // On tourne la roue pour que les operateurs utilisent les champs au temps futur
  eqn.inconnue().avancer();
  eqn.derivee_en_temps_inco(dudt);
  eqn.inconnue().reculer();

  // Un+1=Un+dt_*dU/dt
  futur=dudt;
  futur*=dt_;

  //Adding ALE Jacobians. Jacobians are renewed in Navier_Stokes_std::corriger_derivee_impl().
  // In ALE Un+1=(Jn/Jn+1)*Un+dt_*dU/dt
  Probleme_base& problem=pb_base();
  Domaine_ALE& domaineALE=ref_cast(Domaine_ALE, problem.domaine());
  DoubleTab ALEjacobian_Old=domaineALE.getOldJacobian();
  DoubleTab ALEjacobian_New=domaineALE.getNewJacobian();

  for (int num_face=0; num_face<(futur.size()/dimension); num_face++)
    {
      for (int dim=0; dim<dimension; dim++)
        {
          futur(num_face,dim)+=present(num_face,dim)*(ALEjacobian_Old(num_face,dim)/ALEjacobian_New(num_face,dim));
        }
    }

  futur.echange_espace_virtuel();
  Debog::verifier("Schema_Euler_explicite_ALE::faire_un_pas_de_temps_eqn_base -futur apres", futur);

  eqn.domaine_Cl_dis().imposer_cond_lim(eqn.inconnue(),temps_courant()+pas_de_temps());
  update_critere_statio(dudt, eqn);

  return 1;

}



