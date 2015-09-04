/****************************************************************************
* Copyright (c) 2015, CEA
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
// File:        Neumann_paroi_rayo_semi_transp_VDF.h
// Directory:   $TRUST_ROOT/src/Rayonnement_semi_transp/VDF
// Version:     /main/11
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Neumann_paroi_rayo_semi_transp_VDF_included
#define Neumann_paroi_rayo_semi_transp_VDF_included

#include <Cond_Lim_rayo_semi_transp.h>
#include <Neumann_paroi.h>

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    classe Neumann_paroi_rayo_semi_transp_VDF
//
// .SECTION
//////////////////////////////////////////////////////////////////////////////
class Neumann_paroi_rayo_semi_transp_VDF: public Cond_Lim_rayo_semi_transp, public Neumann_paroi
{
  Declare_instanciable(Neumann_paroi_rayo_semi_transp_VDF);

public :
  int compatible_avec_eqn(const Equation_base&) const;
  const Cond_lim_base& la_cl() const;
  virtual double flux_impose(int i) const;
  virtual double flux_impose(int i,int j) const;
  inline Champ_front& temperature_bord();
  inline const Champ_front& temperature_bord() const;
  void calculer_temperature_bord(double temps);

  void completer();

protected :
  Champ_front temperature_bord_;
};

inline Champ_front& Neumann_paroi_rayo_semi_transp_VDF::temperature_bord()
{
  return temperature_bord_;
}


inline const Champ_front& Neumann_paroi_rayo_semi_transp_VDF::temperature_bord() const
{
  return temperature_bord_;
}

#endif
