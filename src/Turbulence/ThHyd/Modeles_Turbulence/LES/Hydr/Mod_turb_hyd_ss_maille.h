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
// File:        Mod_turb_hyd_ss_maille.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/LES/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Mod_turb_hyd_ss_maille_included
#define Mod_turb_hyd_ss_maille_included

#define CSM1 0.063        // Constante viscosite turbulente modele sous maille
#define CSMS1 0.112        // Constante viscosite turbulente modele sous maille selectif
#define CSM2 0.37       // Constante energie cinetique turbulente modele sous maille

#include <Modele_turbulence_hyd_base.h>

/*! @brief Classe Mod_turb_hyd_ss_maille Classe representant le modele de turbulence sous maille pour les
 *
 *     equations de Navier-Stokes.
 *
 * @sa Modele_turbulence_hyd_base Modele_turbulence_hyd_K_Eps
 */
class Mod_turb_hyd_ss_maille : public Modele_turbulence_hyd_base
{

  Declare_base_sans_constructeur(Mod_turb_hyd_ss_maille);

public:

  Mod_turb_hyd_ss_maille();
  void set_param(Param& param) override;
  void discretiser() override;
  void verifie_loi_paroi_diphasique();
  int preparer_calcul() override;
  void completer() override;
  void mettre_a_jour(double ) override;
  inline virtual Champ_Fonc& energie_cinetique_turbulente() { return energie_cinetique_turb_; }
  inline virtual const Champ_Fonc& energie_cinetique_turbulente() const { return energie_cinetique_turb_; }
  virtual void calculer_longueurs_caracteristiques()=0;
  virtual Champ_Fonc& calculer_viscosite_turbulente()=0;
  virtual void calculer_energie_cinetique_turb();

protected:

  Champ_Fonc energie_cinetique_turb_;
  DoubleVect l_;
  Motcle methode;


};

#endif
