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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Operateur_Conv_sensibility_VEF.cpp
// Directory : $BALTIK_COUPLAGE_ROOT/src/New
//
/////////////////////////////////////////////////////////////////////////////

#include <Operateur_Conv_sensibility_VEF.h>
#include <Op_Conv_VEF_base.h>
#include <DoubleTrav.h>
#include <Op_Conv_VEF_Face.h>
#include <Probleme_base.h>
#include <Schema_Temps_base.h>
#include <Porosites_champ.h>
#include <Debog.h>


Implemente_instanciable( Operateur_Conv_sensibility_VEF, "Op_Conv_sensibility_VEF_P1NC",Operateur_Conv_sensibility ) ;


Sortie&  Operateur_Conv_sensibility_VEF::printOn(Sortie& os) const
{
  return os;
}


Entree&  Operateur_Conv_sensibility_VEF::readOn(Entree& is)
{
  Cerr << " Operateur_Conv_sensibility_VEF ::readOn " << finl;
  op_conv.associer_eqn(equation());
  op_conv.associer_vitesse(la_vitesse.valeur());
  op_conv.lire(is);
  Cerr << "Operateur_Conv_sensibility_VEF : " << op_conv.que_suis_je() << finl;
  return is;
}

void  Operateur_Conv_sensibility_VEF::associer (const Zone_dis& zone_dis ,
                                                const Zone_Cl_dis& zone_cl_dis,
                                                const Champ_Inc& inco )
{
  Cerr << " Operateur_Conv_sensibility_VEF::associer" << finl;
  const Zone_VEF& zvef = ref_cast(Zone_VEF,zone_dis.valeur());
  const Zone_Cl_VEF& zclvef = ref_cast(Zone_Cl_VEF,zone_cl_dis.valeur());

  la_zone_vef = zvef;
  la_zcl_vef = zclvef;
  Operateur_Conv_sensibility::associer(zone_dis,zone_cl_dis,inco);
}

DoubleTab& Operateur_Conv_sensibility_VEF::ajouter(const DoubleTab& inco, DoubleTab& resu) const
{
  Cerr << "Operateur_Conv_sensibility_VEF::ajouter" << finl;

  //Check convection scheme discretization type:
  const Op_Conv_VEF_base& opConvVEFbase = ref_cast(Op_Conv_VEF_base, op_conv.valeur());
  const Op_Conv_VEF_Face& opConvVEFFace = ref_cast(Op_Conv_VEF_Face, op_conv.valeur());
  int convectionSchemeDiscrType; // amont=0, muscl=1, centre=2
  opConvVEFFace.get_type_op(convectionSchemeDiscrType);

  if(convectionSchemeDiscrType==0) // Convection scheme discr "amont".
    {
      Cerr << "Operateur_Conv_sensibility_VEF:: Scheme Conv AMONT" << finl;
      ajouter_Lstate_sensibility_Amont(inco, inco, resu);
      ajouter_Lsensibility_state_Amont(inco, inco, resu);
      opConvVEFbase.modifier_flux(*this); // Multiplication by density in case of incompressible Navier Stokes
    }
  else
    {
      Cout << "Sensibility cannot use currently convection scheme " << op_conv.type() <<". Try available Sensibility amont convection scheme." << finl;
      exit();
    }

  resu.echange_espace_virtuel();
  Debog::verifier("resu dansOperateur_Conv_sensibility_VEF::ajouter", resu);
  return resu;
}

void Operateur_Conv_sensibility_VEF::ajouter_Lstate_sensibility_Amont(const DoubleTab& state, const DoubleTab& inco, DoubleTab& resu ) const
{

}
void Operateur_Conv_sensibility_VEF::ajouter_Lsensibility_state_Amont(const DoubleTab& inco, const DoubleTab& state, DoubleTab& resu ) const
{

}
