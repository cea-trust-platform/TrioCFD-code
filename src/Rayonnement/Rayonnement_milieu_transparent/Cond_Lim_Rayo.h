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
// File:        Cond_Lim_Rayo.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement/src
// Version:     /main/10
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Cond_Lim_Rayo_included
#define Cond_Lim_Rayo_included

#include <Cond_lim_base.h>
#include <TRUST_Ref.h>

class Modele_Rayonnement_Milieu_Transparent;

class Modele_Rayonnement_base;
class Cond_Lim_Rayo
{

public:
  virtual void completer();
  void preparer_surface(const Frontiere_dis_base& ,const Domaine_Cl_dis_base& );
  virtual void associer_modele_rayo(Modele_Rayonnement_base& );
  //virtual const Cond_lim_base& la_cl() const =0;
  inline virtual double surface(int ) const ;
  inline virtual double teta_i(int )  ;
  inline virtual double teta_i(int ) const;
  virtual inline ~Cond_Lim_Rayo();
  inline Modele_Rayonnement_Milieu_Transparent& modele_rayo();
  inline const Modele_Rayonnement_Milieu_Transparent& modele_rayo() const;

protected :
  DoubleVect Surf_i;
  DoubleVect Teta_i;
  REF(Modele_Rayonnement_Milieu_Transparent) le_modele_rayo;
};

inline Cond_Lim_Rayo::~Cond_Lim_Rayo()
{}

inline Modele_Rayonnement_Milieu_Transparent& Cond_Lim_Rayo::modele_rayo()
{
  return le_modele_rayo.valeur();
}

inline const Modele_Rayonnement_Milieu_Transparent& Cond_Lim_Rayo::modele_rayo() const
{
  return le_modele_rayo.valeur();
}
inline double Cond_Lim_Rayo::surface(int numfa) const
{
  return Surf_i[numfa];
}


inline double Cond_Lim_Rayo::teta_i(int numfa)
{
  return Teta_i[numfa];
}

inline double Cond_Lim_Rayo::teta_i(int numfa) const
{
  return Teta_i[numfa];
}




#endif
