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
// File:        Navier_Stokes_std_ALE.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/ALE/src/New
//
//////////////////////////////////////////////////////////////////////////////

#include <Navier_Stokes_std_ALE.h>
#include <Probleme_base.h>
#include <Domaine_ALE.h>
#include <Schema_Temps_base.h>
#include <Op_Conv_ALE_VEF.h>
#include <Discretisation_base.h>
#include <DoubleTrav.h>
#include <Debog.h>

Implemente_instanciable(Navier_Stokes_std_ALE,"Navier_Stokes_standard_ALE",Navier_Stokes_std);

Sortie& Navier_Stokes_std_ALE::printOn(Sortie& os ) const
{
  return Navier_Stokes_std::printOn(os);
}

Entree& Navier_Stokes_std_ALE::readOn(Entree& is )
{
  Navier_Stokes_std::readOn(is);
  return is;
}


void Navier_Stokes_std_ALE::renewing_jacobians( DoubleTab& derivee )
{
  //Renewing ALE Jacobians
  int TimeStepNr=probleme().schema_temps().nb_pas_dt();
  Domaine_ALE& dom_ale=ref_cast(Domaine_ALE, probleme().domaine());

  DoubleTab New_ALEjacobian_Old=dom_ale.getNewJacobian(); //New  value for ALEjacobian_old
  DoubleTab New_ALEjacobian_New(New_ALEjacobian_Old);

  Op_Conv_ALE_VEF& opALEforJacob=ref_cast(Op_Conv_ALE_VEF, terme_convectif.valeur());
  opALEforJacob.calculateALEjacobian(New_ALEjacobian_New); //New  value for ALEjacobian_new
  dom_ale.update_ALEjacobians(New_ALEjacobian_Old, New_ALEjacobian_New, TimeStepNr); // Update new values of ALEjacobian_old and ALEjacobian_new saved in Domaine_ALE

  //End of renewing ALE Jacobians

  Nom discr=discretisation().que_suis_je();
  if (discr != "VEFPreP1B")
    {
      Cerr<<"volume_entrelace_Cl used in the mass matrix is wrong for the ALE treatment on the boundaries"<<finl;
      Cerr<<"(vol_entrelace_Cl=0 for some of the boundaries indeed a correction is necessary) :"<<finl;
      Cerr<<"the VEFPreP1B discretization must be used to avoid this problem. "<<finl;
      exit();
    }
  Cerr << "Adding ALE contribution..." << finl;
  Op_Conv_ALE& opale=ref_cast(Op_Conv_ALE, terme_convectif.valeur());
  DoubleTrav ALE(derivee); // copie de la structure, initialise a zero
  opale.ajouterALE(la_vitesse.valeurs(), ALE);
  ALE.echange_espace_virtuel();
  solveur_masse.appliquer(ALE);
  ALE.echange_espace_virtuel();
  derivee+=ALE; // M-1(F + ALEconvectiveTerm - BtP(n))=derivee_withALEconvectiveTerm
  derivee.echange_espace_virtuel();
  //Cerr << "ALE => norme(derivee) = " << mp_norme_vect(derivee) << finl;
  Debog::verifier("derivee_pression Navier_Stokes_std::corriger_derivee_impl",derivee);
}

void Navier_Stokes_std_ALE::div_ale_derivative( DoubleTrav& deriveeALE, double timestep, DoubleTab& derivee, DoubleTrav& secmemP )
{
  Domaine_ALE& dom_ale=ref_cast(Domaine_ALE, probleme().domaine());
  DoubleTab ALEjacobian_Old=dom_ale.getOldJacobian();
  DoubleTab ALEjacobian_New=dom_ale.getNewJacobian();
  DoubleTab& vitesse_faces_ALE= dom_ale.vitesse_faces();
  DoubleTab term_Jacobian_ratio_U_n(la_vitesse);               // For (Jacobian n/Jacobian n+1) * Un / timestep , initialized every iteration with with new la_vitesse values.

  for (int num_face=0; num_face<(vitesse_faces_ALE.size()/dimension); num_face++)
    {
      for (int dim=0; dim<dimension; dim++)
        {
          term_Jacobian_ratio_U_n(num_face,dim)*=(ALEjacobian_Old(num_face,dim)/(ALEjacobian_New(num_face,dim)*timestep));
          deriveeALE(num_face,dim)=term_Jacobian_ratio_U_n(num_face,dim)+derivee(num_face,dim); //(J_{n}/J_{n+1})*(U_{n}/timestep)+derivee_out
        }
    }
  deriveeALE.echange_espace_virtuel();

  divergence.calculer(deriveeALE, secmemP); //Div((J_{n}/J_{n+1})*(U_{n}/timestep)+derivee_out)
  secmemP *= -1; // car div =-B
  secmemP.echange_espace_virtuel();
  //Debog::verifier("secmemP  modifier Navier_Stokes_std::corriger_derivee_impl",secmemP);
  // Correction du second membre d'apres les conditions aux limites :
  assembleur_pression_.modifier_secmem(secmemP);
  secmemP.echange_espace_virtuel();
  Debog::verifier("secmemP Navier_Stokes_std::corriger_derivee_impl",secmemP);
}

void Navier_Stokes_std_ALE::update_pressure_matrix( void )
{
  //In case of zero ALE mesh velocity (..._coeffs()=1), BM-1Bt matrix stays unchanged.
  Domaine_ALE& dom_ale=ref_cast(Domaine_ALE, probleme().domaine());
  if(dom_ale.update_or_not_matrix_coeffs() == 0)
    {
      assembleur_pression_.assembler(matrice_pression_); // Here B M-1 Bt is assembled.
      solveur_pression_->reinit();
    }
}