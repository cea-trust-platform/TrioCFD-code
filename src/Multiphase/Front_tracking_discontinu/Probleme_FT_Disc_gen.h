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
// File:        Probleme_FT_Disc_gen.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/11
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Probleme_FT_Disc_gen_included
#define Probleme_FT_Disc_gen_included

#include <Pb_Fluide_base.h>
#include <TRUST_Vector.h>
#include <Triple_Line_Model_FT_Disc.h>
#include <TRUST_Ref.h>

class Chimie;
class Equation_base;
class Milieu_base;
class Navier_Stokes_FT_Disc;
class Transport_Interfaces_FT_Disc;

class Probleme_FT_Disc_gen : public Pb_Fluide_base
{
  Declare_instanciable(Probleme_FT_Disc_gen);
public:
  //
  // Methodes virtuelles pures de Probleme_base :
  //
  int nombre_d_equations(void) const override;
  const Equation_base& equation(int i) const override;
  Equation_base& equation(int i) override;

  //
  // Methodes reimplementees de Probleme_base :
  //
  const Equation_base& get_equation_by_name(const Nom& le_nom) const override;
  Equation_base& getset_equation_by_name(const Nom& le_nom) override;
  void   associer_milieu_base(const Milieu_base& milieu) override;
  int associer_(Objet_U& ob) override;
  void   completer(void) override;
  int    verifier(void) override;
  void   preparer_calcul(void) override;
  double calculer_pas_de_temps(void) const override;
  void   mettre_a_jour(double temps) override;
  //
  // Methodes nouvelles de Probleme_FT_Disc_gen
  //
  void associate_triple_line_model(Triple_Line_Model_FT_Disc& tcl_1);
  const Triple_Line_Model_FT_Disc& tcl() const
  {
    return tcl_ ;
  };
  Triple_Line_Model_FT_Disc& tcl()
  {
    return tcl_ ;
  };
  virtual void associer_equation(Equation_base& eq);
  // Raccourcis pour le Front_Tracking
  virtual const Navier_Stokes_FT_Disc&         equation_hydraulique(const Motcle& nom) const;
  virtual const Transport_Interfaces_FT_Disc& equation_interfaces(const Motcle& nom) const;

  // methodes appelees par resoudre_ft_disc -> par encore codees
  virtual void preparer_mise_a_jour(void);

  // methodes appelees dans le readOn -> on ne sait pas encore s'il faut les modifier
  void discretiser(Discretisation_base&) override;
protected:

private:
  // Par convention : dans le vecteur, N.S. en premier, puis Transport_Interfaces,
  //   puis ConvDiff.
  VECT(REF(Equation_base)) equations_;
  REF(Chimie)  la_chimie_;
  Triple_Line_Model_FT_Disc tcl_;
};

#endif
