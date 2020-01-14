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
// File:        Convection_Diffusion_Temperature_FT_Disc.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/15
//
//////////////////////////////////////////////////////////////////////////////
#ifndef Convection_Diffusion_Temperature_FT_Disc_included
#define Convection_Diffusion_Temperature_FT_Disc_included

#include <Convection_Diffusion_Temperature.h>
#include <Champ_Fonc.h>
#include <Champ_Don.h>
#include <Ref_Transport_Interfaces_FT_Disc.h>
#include <Ref_Navier_Stokes_std.h>
#include <Ref_Fluide_Diphasique.h>

////Declare_liste(REF(Champ_base));

class Convection_Diffusion_Temperature_FT_Disc: public Convection_Diffusion_Temperature
{
  Declare_instanciable_sans_constructeur(Convection_Diffusion_Temperature_FT_Disc);
public:

  Convection_Diffusion_Temperature_FT_Disc();
  void set_param(Param& titi);
  int lire_motcle_non_standard(const Motcle&, Entree&);
  virtual void preparer_pas_de_temps(void);
  virtual void corriger_pas_de_temps(double dt);
  DoubleTab&   derivee_en_temps_inco(DoubleTab&);
  void         mettre_a_jour(double temps);

  void         calculer_grad_t();
  void         calculer_mpoint(Champ_base& mpoint);
  virtual void suppression_interfaces(const IntVect& num_compo, const ArrOfInt& flags_compo_a_supprimer, int nouvelle_phase);

  void                associer_milieu_base(const Milieu_base& milieu);
  Milieu_base&        milieu();
  const Milieu_base& milieu() const;
  const Champ_base& vitesse_pour_transport();


  void    discretiser(void);

protected:
  void update_rho_cp(double temps);
  // Quelle phase cette equation concerne-t-elle ? 0 ou 1
  int phase_;
  //GB : Ajout de variables :
  int stencil_width_;
  int correction_courbure_ordre_;
  Nom nom_sous_zone_;
  double temp_moy_ini_;
  bool maintien_temperature_;

  REF(Fluide_Diphasique) fluide_dipha_;

  // Champs compris par le postraitement
  LIST(REF(Champ_base)) liste_champs_compris_;

  // Reference a l'equation de transport de l'interface pour l'indicatrice de phase
  REF(Transport_Interfaces_FT_Disc) ref_eq_interface_;
  // Reference a l'equation de navier_stokes pour le champ de vitesse (convection)
  REF(Navier_Stokes_std) ref_eq_ns_;

  // Gradient normal de temperature a l'interface phase 0
  // (grad T scalaire normale a l'interface, normale dirigee
  //  vers la phase 1)
  Champ_Fonc grad_t_;
  Champ_Inc vitesse_convection_;
  Champ_Fonc rho_cp_elem_,rho_cp_comme_T_;


};
#endif
