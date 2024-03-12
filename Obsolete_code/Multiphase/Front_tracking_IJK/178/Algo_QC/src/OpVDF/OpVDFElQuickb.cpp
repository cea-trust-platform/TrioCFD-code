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
/////////////////////////////////////////////////////////////////////////////
//
// File      : OpVDFElQuickb.cpp
// Directory : $IJK_ROOT/src/IJK/OpVDF
//
/////////////////////////////////////////////////////////////////////////////

#include <OpVDFElQuickb.h>

Implemente_instanciable_sans_constructeur(Op_Conv_Quickb_VDF_Elem,"Op_Conv_Quickb_VDF_P0_VDF",Op_Conv_VDF_base);

implemente_It_VDF_Elem(Eval_Quickb_VDF_Elem)

//// printOn
//

Sortie& Op_Conv_Quickb_VDF_Elem::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

//// readOn
//

Entree& Op_Conv_Quickb_VDF_Elem::readOn(Entree& s )
{
  return s ;
}


// Description:
// complete l'iterateur et l'evaluateur
void Op_Conv_Quickb_VDF_Elem::associer(const Zone_dis& zone_dis,
                                       const Zone_Cl_dis& zone_cl_dis,
                                       const Champ_Inc& ch_transporte)
{
  const Zone_VDF& zvdf = ref_cast(Zone_VDF,zone_dis.valeur());
  const Zone_Cl_VDF& zclvdf = ref_cast(Zone_Cl_VDF,zone_cl_dis.valeur());
  const Champ_P0_VDF& inco = ref_cast(Champ_P0_VDF,ch_transporte.valeur());

  iter.associer(zvdf, zclvdf, *this);

  Eval_Quickb_VDF_Elem& eval_conv = (Eval_Quickb_VDF_Elem&) iter.evaluateur();
  eval_conv.associer_zones(zvdf, zclvdf );          // Evaluateur_VDF::associer_zones
  eval_conv.associer_inconnue(inco );        // Eval_VDF_Elem::associer_inconnue
  if (Process::nproc()>1 && zvdf.zone().nb_joints() && zvdf.zone().joint(0).epaisseur()<2)
    {
      Cerr << "Overlapping width (given by larg_joint option) of  " << zvdf.zone().joint(0).epaisseur() << finl;
      Cerr << "is not enough for Quickb scheme in VDF parallel calculation." << finl;
      Cerr << "Please, partition your mesh with an overlapping width of 2 with larg_joint option." << finl;
      Process::exit();
    }
}

const Champ_Inc_base& Op_Conv_Quickb_VDF_Elem::vitesse() const
{
  Eval_Quickb_VDF_Elem& eval_conv = (Eval_Quickb_VDF_Elem&) iter.evaluateur();
  return eval_conv.vitesse();
}

Champ_Inc_base& Op_Conv_Quickb_VDF_Elem::vitesse()
{
  Eval_Quickb_VDF_Elem& eval_conv = (Eval_Quickb_VDF_Elem&) iter.evaluateur();
  return eval_conv.vitesse();
}



// Description:
// constructeur
Op_Conv_Quickb_VDF_Elem::Op_Conv_Quickb_VDF_Elem()
  :Op_Conv_VDF_base(It_VDF_Elem(Eval_Quickb_VDF_Elem)())
{
}
