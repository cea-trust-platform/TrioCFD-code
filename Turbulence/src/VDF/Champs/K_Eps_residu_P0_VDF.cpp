/****************************************************************************
* Copyright (c) 2017, CEA
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
// File:        K_Eps_residu_Champ_Face.cpp
//////////////////////////////////////////////////////////////////////////////

#include <K_Eps_residu_P0_VDF.h>
#include <Champ_Face.h>
#include <Zone_VDF.h>
#include <Champ_Inc.h>
#include <Equation_base.h>


Implemente_instanciable(K_Eps_residu_P0_VDF,"K_Eps_residu_P0_VDF",Champ_Fonc_P0_VDF);


//     printOn()
/////

Sortie& K_Eps_residu_P0_VDF::printOn(Sortie& s) const
{
  return s << que_suis_je() << " " << le_nom();
}

//// readOn
//

Entree& K_Eps_residu_P0_VDF::readOn(Entree& s)
{
  return s ;
}

void K_Eps_residu_P0_VDF::me_calculer(double tps)
{ 
  const Zone_VDF& zone_VDF = ref_cast(Zone_VDF,K_Eps_->zone_dis_base());
  const int nb_elem = zone_VDF.nb_elem();
  const DoubleTab K_Eps = K_Eps_.valeur( ).passe( );
  DoubleTab& champ = le_champ( ).valeurs( );

  if( nb_elem != K_Eps.dimension( 0 ) )
    {
      Cerr << "Error in K_Eps_residu_P0_VDF::me_calculer "<<finl;
      Cerr << "There are "<<nb_elem<<" elements and "<<K_Eps.dimension( 0 )<<" values for K_Eps field "<<finl;
      Process::abort( );
    }

  if( tps > 0.0 )
	  champ = K_Eps_.valeur( ).equation( ).get_residuals();

  else
    {
      champ = -10000.0 ;
      Cerr << "[Information] K_Eps_residu_P0_VDF::me_calculer : le residu est mis a -10000.0 au temps initial"<<finl;
    }
}


