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
// File:        Implicite_ALE.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/ALE/src/New
//
//////////////////////////////////////////////////////////////////////////////

#include <Implicite_ALE.h>
#include <Op_Conv_ALE_VEF.h>
#include <Zone_ALE.h>
#include <Probleme_base.h>
#include <Schema_Temps_base.h>
#include <Debog.h>

Implemente_instanciable(Implicite_ALE,"Implicite_ALE",Implicite);
// XD implicite_ALE implicite implicite_ALE -1 Implicite solver used for ALE problem
Sortie& Implicite_ALE::printOn(Sortie& os ) const
{
  return Implicite::printOn(os);
}

Entree& Implicite_ALE::readOn(Entree& is )
{
  Implicite::readOn(is);
  return is;
}


void Implicite_ALE::first_special_treatment(Equation_base& eqn, Navier_Stokes_std& eqnNS, DoubleTab& current, double dt, DoubleTrav& resu)
{
  //Renewing ALE Jacobians
  int TimeStepNr=eqn.probleme().schema_temps().nb_pas_dt();
  Zone_ALE& dom_ale=ref_cast(Zone_ALE, eqn.probleme().domaine());

  DoubleTab New_ALEjacobian_Old=dom_ale.getNewJacobian(); //New  value for ALEjacobian_old
  DoubleTab New_ALEjacobian_New(New_ALEjacobian_Old);
  Op_Conv_ALE_VEF& opALEforJacob=ref_cast(Op_Conv_ALE_VEF, eqnNS.get_terme_convectif().valeur());
  opALEforJacob.calculateALEjacobian(New_ALEjacobian_New); //New  value for ALEjacobian_new
  dom_ale.update_ALEjacobians(New_ALEjacobian_Old, New_ALEjacobian_New, TimeStepNr); // Update new values of ALEjacobian_old and ALEjacobian_new saved in Domaine_ALE
  //End of renewing ALE Jacobians

  Cerr << "Adding ALE contribution..." << finl;
  Op_Conv_ALE& opale=ref_cast(Op_Conv_ALE, eqnNS.get_terme_convectif().valeur());
  DoubleTrav ALE(resu); // copie de la structure, initialise a zero
  opale.ajouterALE(current, ALE);
  ALE.echange_espace_virtuel();
  //solveur_masse.appliquer(ALE); do not need to divide by mass
  //ALE.echange_espace_virtuel();
  resu+=ALE;                        //resu + ALE convection tem
  resu.echange_espace_virtuel();
  Debog::verifier("Implicite_ALE::first_special_treatment resu after adding ALE",resu);
}

void Implicite_ALE::second_special_treatment(Equation_base& eqn,DoubleTab& current, DoubleTrav& resu, Matrice_Morse& matrice)
{
  Zone_ALE& dom_ale=ref_cast(Zone_ALE, eqn.probleme().domaine());
  DoubleVect ALEjacobian_New=dom_ale.getNewJacobian();

  Debog::verifier(" Test ALE Jacobian New", ALEjacobian_New );

  //First step
  // Adding Jacobians to "matrice"  Jn+1[A[Un]]
  int MatriceNbLines=matrice.nb_lignes();
  for (int num_face=0; num_face<MatriceNbLines; num_face++)
    matrice(num_face,num_face)*=ALEjacobian_New[num_face];
  // Adding Jacobians to "matrice" end

  //Second step
  //Adding ALE Jacobians to "resu" start to obtain Jn+1[-gradP + ALE convective term + Sv + Ss]+Jn[(M/dt)Un].
  // Obtaining (M/dt)Un start.
  DoubleTrav deltaJMdtUn(resu); // copie de la structure, initialise a zero. To hold (M/dt)Un.
  double timestep=eqn.probleme().schema_temps().pas_de_temps();
  int pen=0;
  eqn.solv_masse().ajouter_masse(timestep,deltaJMdtUn,eqn.inconnue().passe(),pen);
  // Obtaining (M/dt)Un end.

  DoubleTab ALEjacobian_Old=dom_ale.getOldJacobian();
  DoubleTab ALEjacobian=dom_ale.getNewJacobian();

  int nb_faces = resu.size()/dimension;
  for (int num_face=0; num_face<nb_faces; num_face++) // Obtaining Jn+1[-gradP + ALE convective term + (M/dt)Un + Sv + Ss]+(Jn-Jn+1)[(M/dt)Un]=Jn+1[-gradP + ALE convective term + Sv + Ss]+Jn[(M/dt)Un].
    {
      for (int dim=0; dim<dimension; dim++)
        {
          resu(num_face,dim)=ALEjacobian(num_face,dim)*resu(num_face,dim)+(ALEjacobian_Old(num_face,dim)-ALEjacobian(num_face,dim))*deltaJMdtUn(num_face,dim);
        }
    }
  resu.echange_espace_virtuel();
}
