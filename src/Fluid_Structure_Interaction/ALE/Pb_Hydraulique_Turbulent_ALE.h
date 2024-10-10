/****************************************************************************
* Copyright (c) 2019, CEA
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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Pb_Hydraulique_Turbulent_ALE.h
// Directory : $ALE_ROOT/src/New
//
/////////////////////////////////////////////////////////////////////////////


#ifndef Pb_Hydraulique_Turbulent_ALE_included
#define Pb_Hydraulique_Turbulent_ALE_included

#include <Pb_Fluide_base.h>
#include <Navier_Stokes_Turbulent_ALE.h>



/*! @brief classe Pb_Hydraulique_Turbulent_ALE Cette classe represente un probleme d'hydraulique turbulent dans
 *
 *      lequel on resout les equations de Navier Stokes en regime turbulent
 *      pour un fluide incompressible et en maillage mobile (ALE)
 *      La formulation est de type vitesse pression
 *
 * @sa Pb_Fluide_base Navier_Stokes_Turbulent_ALE Pb_Hydraulique
 */
class Pb_Hydraulique_Turbulent_ALE : public Pb_Fluide_base
{
  Declare_instanciable(Pb_Hydraulique_Turbulent_ALE);

public :

  int nombre_d_equations() const override;
  const Equation_base& equation(int) const override;
  Equation_base& equation(int) override;
  inline const Champ_Fonc_base& viscosite_turbulente() const;
  void associer_milieu_base(const Milieu_base& ) override;
  int verifier() override;
  void associer_sch_tps_base(const Schema_Temps_base&) override;


protected :

  Navier_Stokes_Turbulent_ALE eq_hydraulique;

};



/*! @brief Renvoie le champ representant la viscosite turbulente du probleme.
 *
 * @return (Champ_Fonc_base&) le champ representant la viscosite turbulente
 */
inline const Champ_Fonc_base& Pb_Hydraulique_Turbulent_ALE::viscosite_turbulente() const
{
  return eq_hydraulique.viscosite_turbulente();
}

#endif
