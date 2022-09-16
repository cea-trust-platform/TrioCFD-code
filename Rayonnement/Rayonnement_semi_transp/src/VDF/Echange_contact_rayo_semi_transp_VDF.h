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
// File:        Echange_contact_rayo_semi_transp_VDF.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src/VDF
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Echange_contact_rayo_semi_transp_VDF_included
#define Echange_contact_rayo_semi_transp_VDF_included

#include <Cond_Lim_rayo_semi_transp.h>
#include <Echange_contact_VDF.h>

/*! @brief classe Echange_contact_rayo_semi_transp_VDF Cette classe est utilisee pour realiser un couplage entre une
 *
 *     paroi solide et un fluide semi transparent en prenant en compte
 *     le rayonnement de la paroi.
 *     Elle est utilisee uniquement pour une discretisation VDF.
 *
 *
 */
class Echange_contact_rayo_semi_transp_VDF: public Cond_Lim_rayo_semi_transp, public Echange_contact_VDF
{
  Declare_instanciable(Echange_contact_rayo_semi_transp_VDF);

public :
  const Cond_lim_base& la_cl() const override;
  Champ_front& temperature_bord();
  void calculer_temperature_bord(double temps);
  int compatible_avec_eqn(const Equation_base&) const override;
  void mettre_a_jour(double temps) override;

  void completer() override;
  void calculer_Teta_paroi(DoubleTab& Teta_p,const DoubleTab& mon_h,const DoubleTab& autre_h,int i,double temps) override;
  void calculer_Teta_equiv(DoubleTab& Teta_eq,const DoubleTab& mon_h,const DoubleTab& autre_h,int i,double temps) override;
  Echange_contact_rayo_semi_transp_VDF& la_Cl_opposee();

protected :

  Champ_front T_paroi;
  int num_premiere_face_dans_pb_fluide;
};


#endif

