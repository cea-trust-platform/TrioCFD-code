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
// File      : Echange_contact_VDF_FT_Disc.h
// Directory : $FRONT_TRACKING_DISCONTINU_ROOT/src/VDF
//
/////////////////////////////////////////////////////////////////////////////

#ifndef Echange_contact_VDF_FT_Disc_included
#define Echange_contact_VDF_FT_Disc_included

#include <Echange_contact_VDF.h>

/*! @brief : class Echange_contact_VDF_FT_Disc
 *
 *  <Description of class Echange_contact_VDF_FT_Disc>
 *
 *
 *
 */
class Echange_contact_VDF_FT_Disc : public Echange_contact_VDF
{

  Declare_instanciable( Echange_contact_VDF_FT_Disc ) ;

public :
  void mettre_a_jour(double ) override;
  int initialiser(double temps) override;
  void completer() override;
  virtual double Ti_wall(int num) const;
  void changer_temps_futur(double temps,int i) override;
  int avancer(double temps) override;
  int reculer(double temps) override;
  void set_temps_defaut(double temps) override;
protected :
  Champ_front indicatrice_;

  // T_wall_ is defined as DoubleTab in Echange_contact_VDF
  // Ti_wall_ temperature at wall-fluid interface
  Champ_front Ti_wall_;
  double indicatrice_ref_;
  Nom nom_champ_indicatrice_, nom_bord_oppose_;
};

#endif /* Echange_contact_VDF_FT_Disc_included */
