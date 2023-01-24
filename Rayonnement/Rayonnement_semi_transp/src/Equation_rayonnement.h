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
// File:        Equation_rayonnement.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src
// Version:     /main/13
//
//////////////////////////////////////////////////////////////////////////////

//

#ifndef Equation_rayonnement_included
#define Equation_rayonnement_included

#include <Equation_rayonnement_base.h>
#include <TRUST_Deriv.h>


class Equation_rayonnement : public DERIV(Equation_rayonnement_base)
{
  Declare_instanciable(Equation_rayonnement);

public :

  inline void associer_zone_dis(const Zone_dis& Zd);
  inline void associer_milieu_base(const Milieu_base&);
  inline void associer_fluide(const Fluide_base&);
  inline void associer_modele_rayonnement(const Modele_rayo_semi_transp&);
  inline void associer_pb_base(const Probleme_base&);
  inline void associer_sch_tps_base( const Schema_Temps_base&);

  inline Milieu_base& milieu();
  inline const Milieu_base& milieu() const;
  inline const Modele_rayo_semi_transp&  Modele() const;
  inline Modele_rayo_semi_transp& Modele();
  inline int nombre_d_operateurs() const;
  inline const Operateur& operateur(int) const;
  inline Operateur& operateur(int);
  inline const Champ_Inc& inconnue() const;
  inline Champ_Inc& inconnue();
  inline void discretiser();
  inline Fluide_base& fluide();
  inline const Fluide_base& fluide() const;
  inline Operateur_Grad& operateur_gradient();
  inline const Operateur_Grad& operateur_gradient() const;

  inline Zone_Cl_dis& zone_Cl_dis();
  inline const Zone_Cl_dis& zone_Cl_dis() const;
  inline Zone_dis& zone_dis();
  inline const Zone_dis& zone_dis() const;

  inline void mettre_a_jour(double temps);
  inline void modifier_matrice();
  inline void resoudre(double temps);
  inline void assembler_matrice();

  inline void completer();
  inline void typer_op_grad();

  void typer(Nom& type, const Equation_base&);


};


inline void Equation_rayonnement::assembler_matrice()
{
  valeur().assembler_matrice();
}

inline Operateur_Grad& Equation_rayonnement::operateur_gradient()
{
  return valeur().operateur_gradient();
}

inline void Equation_rayonnement::typer_op_grad()
{
  valeur().typer_op_grad();
}

inline const Operateur_Grad& Equation_rayonnement::operateur_gradient() const
{
  return valeur().operateur_gradient();
}

inline void Equation_rayonnement::completer()
{
  valeur().completer();
}

inline Zone_dis& Equation_rayonnement::zone_dis()
{
  return valeur().zone_dis();
}

inline const Zone_dis& Equation_rayonnement::zone_dis() const
{
  return valeur().zone_dis();
}


inline void Equation_rayonnement::associer_sch_tps_base( const Schema_Temps_base& sch)
{
  valeur().associer_sch_tps_base(sch);
}

inline void Equation_rayonnement::associer_pb_base(const Probleme_base& pb)
{
  valeur().associer_pb_base(pb);
}


inline Fluide_base& Equation_rayonnement::fluide()
{
  return valeur().fluide();
}

inline const Fluide_base& Equation_rayonnement::fluide() const
{
  return valeur().fluide();
}


inline void Equation_rayonnement::associer_milieu_base(const Milieu_base& un_milieu)
{
  valeur().associer_milieu_base(un_milieu);
}

inline void Equation_rayonnement::associer_fluide(const Fluide_base& un_fluide)
{
  valeur().associer_fluide(un_fluide);
}

inline void Equation_rayonnement::associer_modele_rayonnement(const Modele_rayo_semi_transp& modele)
{
  valeur().associer_modele_rayonnement(modele);
}

inline Milieu_base& Equation_rayonnement::milieu()
{
  return valeur().milieu();
}

inline const Milieu_base& Equation_rayonnement::milieu() const
{
  return valeur().milieu();
}

inline const Modele_rayo_semi_transp&  Equation_rayonnement::Modele() const
{
  return valeur().Modele();
}

inline Modele_rayo_semi_transp& Equation_rayonnement::Modele()
{
  return valeur().Modele();
}

inline int Equation_rayonnement::nombre_d_operateurs() const
{
  return valeur().nombre_d_operateurs();
}

inline const Operateur& Equation_rayonnement::operateur(int i) const
{
  return valeur().operateur(i);
}

inline Operateur& Equation_rayonnement::operateur(int i)
{
  return valeur().operateur(i);
}

inline const Champ_Inc& Equation_rayonnement::inconnue() const
{
  return valeur().inconnue();
}

inline Champ_Inc& Equation_rayonnement::inconnue()
{
  return valeur().inconnue();
}

inline void Equation_rayonnement::discretiser()
{
  valeur().discretiser();
}

inline void Equation_rayonnement::modifier_matrice()
{
  valeur().modifier_matrice();
}

inline void Equation_rayonnement::resoudre(double temps)
{

  valeur().resoudre(temps);
}

inline void Equation_rayonnement::mettre_a_jour(double temps)
{
  valeur().mettre_a_jour(temps);
}

inline Zone_Cl_dis& Equation_rayonnement::zone_Cl_dis()
{
  return valeur().zone_Cl_dis();
}

inline const Zone_Cl_dis& Equation_rayonnement::zone_Cl_dis() const
{
  return valeur().zone_Cl_dis();
}

inline void Equation_rayonnement::associer_zone_dis(const Zone_dis& Zd)
{
  valeur().associer_zone_dis(Zd);
}

#endif

