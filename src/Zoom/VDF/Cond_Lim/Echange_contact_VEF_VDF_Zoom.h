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
// File:        Echange_contact_VEF_VDF_Zoom.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/VDF/Cond_Lim
// Version:     /main/11
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Echange_contact_VEF_VDF_Zoom_included
#define Echange_contact_VEF_VDF_Zoom_included

#include <Echange_contact_VDF_Zoom_base.h>
#include <Champ_front_zoom.h>
#include <Pb_MG.h>
#include <Domaine_VEF.h>
#include <Milieu_base.h>
#include <Champ_front_fonc.h>
#include <Champ_Uniforme.h>
#include <Temperature_imposee_paroi.h>

/*! @brief classe : Echange_contact_VEF_VDF_Zoom
 *
 *
 *
 */

class Echange_contact_VEF_VDF_Zoom  : public Temperature_imposee_paroi
{

  Declare_instanciable(Echange_contact_VEF_VDF_Zoom);

public :
  void mettre_a_jour(double ) override ;
  int compatible_avec_eqn(const Equation_base&) const override { return 1; }
};


#endif
