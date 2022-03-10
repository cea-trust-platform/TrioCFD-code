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
// File      : Echange_contact_VDF_FT_Disc_solid.h
// Directory : $FRONT_TRACKING_DISCONTINU_ROOT/src/VDF
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Echange_contact_VDF_FT_Disc_solid_included
#define Echange_contact_VDF_FT_Disc_solid_included

#include <Echange_contact_VDF_FT_Disc.h>

/////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION : class Echange_contact_VDF_FT_Disc_solid
//
// <Description of class Echange_contact_VDF_FT_Disc_solid>
//
/////////////////////////////////////////////////////////////////////////////

class Echange_contact_VDF_FT_Disc_solid : public Echange_contact_VDF_FT_Disc
{

  Declare_instanciable( Echange_contact_VDF_FT_Disc_solid ) ;

public :
  void mettre_a_jour(double ) override;
  int initialiser(double temps) override;
  void completer() override;


  void changer_temps_futur(double temps,int i) override;
  int avancer(double temps) override;
  int reculer(double temps) override;
  inline Champ_front& T_autre_pb() override
  {
    if (numero_T_==0)
      return T_autre_pb_;
    assert(numero_T_==1);
    return T2_autre_pb_;
  };
  inline const Champ_front& T_autre_pb() const override
  {
    if (numero_T_==0)
      return T_autre_pb_;
    assert(numero_T_==1);
    return T2_autre_pb_;
  };
protected :
  Champ_front T2_autre_pb_;
  int numero_T_;
  Nom nom_champ_T2_autre_pb_;
};

#endif /* Echange_contact_VDF_FT_Disc_solid_included */
