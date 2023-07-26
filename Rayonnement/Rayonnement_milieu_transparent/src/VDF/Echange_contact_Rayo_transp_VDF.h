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
// File:        Echange_contact_Rayo_transp_VDF.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement/src/VDF
// Version:     /main/15
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Echange_contact_Rayo_transp_VDF_included
#define Echange_contact_Rayo_transp_VDF_included

#include <Echange_contact_VDF.h>
#include <Convection_Diffusion_Temperature.h>
#include <Conduction.h>
#include <Champ_front_calc.h>
#include <Fluide_Incompressible.h>
#include <Front_VF.h>
#include <Cond_Lim_Rayo.h>
#include <TRUST_Ref.h>

class Domaine_Cl_VDF;
class Domaine_VDF;

class Echange_contact_Rayo_transp_VDF : public Cond_Lim_Rayo,public Echange_contact_VDF
{

  Declare_instanciable(Echange_contact_Rayo_transp_VDF);

public :
  int compatible_avec_eqn(const Equation_base&) const override { return 1; }
  void completer() override;
  void mettre_a_jour(double ) override;
  void calculer_Teta_paroi(DoubleTab& tab_p,const DoubleTab& mon_h,const DoubleTab& autre_h,int is_pb_fluide,double temps) override;
  void calculer_Teta_equiv(DoubleTab& Teta_equiv,const DoubleTab& mon_h,const DoubleTab& autre_h,int is_pb_fluide,double temps) override;
  //int verifier_correspondance() const;
protected :
  int num_premiere_face_dans_pb_fluide;
  double alpha_;
};

class Echange_contact_Rayo_transp_sans_relax_VDF : public Echange_contact_Rayo_transp_VDF
{
  Declare_instanciable(Echange_contact_Rayo_transp_sans_relax_VDF);
};


#endif
