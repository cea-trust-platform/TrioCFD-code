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
// File:        Flux_radiatif.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Flux_radiatif_included
#define Flux_radiatif_included

#include <Flux_radiatif_base.h>
class Zone_Cl_dis_base;

Declare_deriv(Flux_radiatif_base);

//////////////////////////////////////////////////////////////////////////////
//
// .DESCRIPTION
//    classe Flux_radiatif
//    Classe generique de la hierarchie des flux radiatifs
//    Flux_radiatif peut referencer n'importe quel
//    objet derivant de Flux_radiatif_base.
// .SECTION voir aussi
//      Flux_radiatif_base
//////////////////////////////////////////////////////////////////////////////
class Flux_radiatif : public DERIV(Flux_radiatif_base)
{
  Declare_instanciable(Flux_radiatif);

public :

  inline Champ_front& emissivite();
  inline const Champ_front& emissivite() const;
  inline double& A();
  inline const double& A() const;

  inline void calculer_flux_radiatif(const Equation_base& eq_temp);
  inline Champ_front& flux_radiatif();
  inline const Champ_front& flux_radiatif() const;
  inline void completer();

  inline double flux_impose(int i) const;
  inline double flux_impose(int i,int j) const;
  inline Zone_Cl_dis_base& zone_Cl_dis();

  void typer();
  inline void typer(const Nom& type);
};


inline Zone_Cl_dis_base& Flux_radiatif::zone_Cl_dis()
{
  return valeur().zone_Cl_dis();
}


inline Champ_front& Flux_radiatif::emissivite()
{
  return valeur().emissivite();
}

inline const Champ_front& Flux_radiatif::emissivite() const
{
  return valeur().emissivite();
}

inline double& Flux_radiatif::A()
{
  return valeur().A();
}

inline const double& Flux_radiatif::A() const
{
  return valeur().A();
}

inline void Flux_radiatif::calculer_flux_radiatif(const Equation_base& eq_temp)
{
  valeur().calculer_flux_radiatif(eq_temp);
}

inline Champ_front& Flux_radiatif::flux_radiatif()
{
  return valeur().flux_radiatif();
}

inline const Champ_front& Flux_radiatif::flux_radiatif() const
{
  return valeur().flux_radiatif();
}

inline void Flux_radiatif::completer()
{
  valeur().completer();
}

inline double Flux_radiatif::flux_impose(int i) const
{
  return valeur().flux_impose(i);
}

inline double Flux_radiatif::flux_impose(int i,int j) const
{
  return valeur().flux_impose(i,j);
}

inline void Flux_radiatif::typer(const Nom& typ)
{
  DERIV(Flux_radiatif_base)::typer(typ);
}

#endif

