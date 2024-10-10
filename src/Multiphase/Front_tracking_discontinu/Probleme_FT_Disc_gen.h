/****************************************************************************
* Copyright (c) 2024, CEA
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

#include <Triple_Line_Model_FT_Disc.h>
#include <Pb_Fluide_base.h>
#include <Equation_base.h>
#include <TRUST_Vector.h>
#include <TRUST_List.h>
#include <TRUST_Ref.h>

class Chimie;
class Equation_base;
class Milieu_base;
class Navier_Stokes_FT_Disc;
class Transport_Interfaces_FT_Disc;

class Probleme_FT_Disc_gen: public Pb_Fluide_base
{
  Declare_instanciable(Probleme_FT_Disc_gen);
public:
  int nombre_d_equations() const override { return equations_.size(); }
  const Equation_base& equation(int i) const override { return equations_[i].valeur(); }
  Equation_base& equation(int i) override { return equations_[i].valeur(); }

  void typer_lire_milieu(Entree& is) override;
  void lire_solved_equations(Entree& is) override;
  Entree& lire_equations(Entree& is, Motcle& dernier_mot) override;

  const Equation_base& get_equation_by_name(const Nom& le_nom) const override;
  Equation_base& getset_equation_by_name(const Nom& le_nom) override;
  void associer_milieu_base(const Milieu_base& milieu) override;
  int associer_(Objet_U& ob) override;
  void completer() override;
  double calculer_pas_de_temps() const override;
  void mettre_a_jour(double temps) override;
  virtual bool updateGivenFields() override;

  // Raccourcis pour le Front_Tracking
  virtual const Navier_Stokes_FT_Disc& equation_hydraulique(const Motcle& nom) const;
  virtual const Transport_Interfaces_FT_Disc& equation_interfaces(const Motcle& nom) const;

  const Triple_Line_Model_FT_Disc& tcl() const { return tcl_; }
  Triple_Line_Model_FT_Disc& tcl() { return tcl_; }

  inline const LIST(OWN_PTR(Equation_base))& get_list_equations() const { return equations_; }

private:
  void add_FT_equation(const Nom& , const Nom& );
  LIST(OWN_PTR(Equation_base)) equations_; // Par convention : dans le vecteur, N.S. en premier, puis Transport_Interfaces, puis ConvDiff.
  OBS_PTR(Chimie) la_chimie_;
  Triple_Line_Model_FT_Disc tcl_;
};

#endif /* Probleme_FT_Disc_gen_included */
