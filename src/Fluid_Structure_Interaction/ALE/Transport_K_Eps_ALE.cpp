/****************************************************************************
* Copyright (c) 2022, CEA
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
// File:        Transport_K_Eps_ALE.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/ALE/src/
//
//////////////////////////////////////////////////////////////////////////////

#include <Transport_K_Eps_ALE.h>
#include <Pb_Hydraulique_Turbulent_ALE.h>
#include <Op_Conv_ALE.h>
#include <TRUSTTrav.h>
#include <Debog.h>
#include <Discretisation_base.h>

Implemente_instanciable(Transport_K_Eps_ALE,"Transport_K_Eps_ALE",Transport_K_Eps);

Sortie& Transport_K_Eps_ALE::printOn(Sortie& s ) const
{
  return Transport_K_Eps::printOn(s);
}

Entree& Transport_K_Eps_ALE::readOn(Entree& s )
{
  return Transport_K_Eps::readOn(s);
}

void Transport_K_Eps_ALE::corriger_derivee_impl_ALE(DoubleTab& d)
{
  // ALE part start
  if( sub_type(Op_Conv_ALE, terme_convectif.valeur()) )
    {
      Nom discr=discretisation().que_suis_je();
      if (discr != "VEFPreP1B")
        {
          Cerr<<"volume_entrelace_Cl used in the mass matrix is wrong for the ALE treatment on the boundaries"<<finl;
          Cerr<<"(vol_entrelace_Cl=0 for some of the boundaries indeed a correction is necessary) :"<<finl;
          Cerr<<"the VEFPreP1B discretization must be used to avoid this problem. "<<finl;
          exit();
        }
      Cerr << "Adding ALE contribution to K Eps..." << finl;
      Op_Conv_ALE& opale=ref_cast(Op_Conv_ALE, terme_convectif.valeur());
      DoubleTrav ALE(d); // copie de la structure, initialise a zero
      opale.ajouterALE(inconnue().valeurs(), ALE);
      ALE.echange_espace_virtuel();
      solveur_masse.appliquer(ALE);
      ALE.echange_espace_virtuel();
      d+=ALE;
      d.echange_espace_virtuel();
      Debog::verifier("Transport_K_Eps::corriger_derivee_impl",d);
    }
}

