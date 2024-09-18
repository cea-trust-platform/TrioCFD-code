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
// File:        Modele_turbulence_hyd_LES_DSGS_VDF.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Modeles_Turbulence/LES/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Modele_turbulence_hyd_LES_DSGS_VDF_included
#define Modele_turbulence_hyd_LES_DSGS_VDF_included

#include <Modele_turbulence_hyd_LES_Smago_VDF.h>

/*! @brief classe Modele_turbulence_hyd_LES_DSGS_VDF Cette classe derivee de Modele_turbulence_hyd_LES_Smago_VDF
 *
 *  correspond a la mise en oeuvre du modele sous maille dynamique en VDF
 *
 *  .SECTION  voir aussi
 *  Modele_turbulence_hyd_LES_Smago_VDF
 *
 */
class Modele_turbulence_hyd_LES_DSGS_VDF: public Modele_turbulence_hyd_LES_Smago_VDF
{
  Declare_instanciable(Modele_turbulence_hyd_LES_DSGS_VDF);
public:
  void associer(const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis& domaine_Cl_dis) override;

protected:
  Champ_Fonc coeff_field;
  DoubleVect model_coeff;

  void calculer_cell_cent_vel(DoubleTab&);
  void calculer_filter_field(const DoubleTab&, DoubleTab&);
  void calculer_filter_tensor(DoubleTab&);
  void calculer_filter_coeff(DoubleVect&);
  void calculer_Lij(const DoubleTab&, const DoubleTab&, DoubleTab&);
  void calculer_Sij(const DoubleTab&, DoubleTab&);
  void calculer_Mij(const DoubleTab&, const DoubleTab&, DoubleTab&);
  void calculer_model_coefficient(const DoubleTab&, const DoubleTab&);
  Champ_Fonc& calculer_viscosite_turbulente() override;
};

#endif /* Modele_turbulence_hyd_LES_DSGS_VDF_included */
