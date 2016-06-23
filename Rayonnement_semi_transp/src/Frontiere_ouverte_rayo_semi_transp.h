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
// File:        Frontiere_ouverte_rayo_semi_transp.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Frontiere_ouverte_rayo_semi_transp_included
#define Frontiere_ouverte_rayo_semi_transp_included

#include <Cond_Lim_rayo_semi_transp.h>
class Champ_front;

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    classe Frontiere_ouverte_rayo_semi_transp
//
// .SECTION
//////////////////////////////////////////////////////////////////////////////
class Frontiere_ouverte_rayo_semi_transp: public Cond_Lim_rayo_semi_transp, public Neumann_sortie_libre
{
  Declare_instanciable(Frontiere_ouverte_rayo_semi_transp);

public :

  const Cond_lim_base& la_cl() const;
  void completer();

  inline Champ_front& temperature_bord();
  inline const Champ_front& temperature_bord() const;
  void calculer_temperature_bord(double temps);

protected :

};

inline Champ_front& Frontiere_ouverte_rayo_semi_transp::temperature_bord()
{
  return le_champ_front;
}

inline const Champ_front& Frontiere_ouverte_rayo_semi_transp::temperature_bord() const
{
  return le_champ_front;
}


#endif
