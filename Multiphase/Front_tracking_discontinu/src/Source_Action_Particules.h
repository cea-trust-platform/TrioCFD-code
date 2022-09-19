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
// File:        Source_Action_Particules.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/7
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Source_Action_Particules_included
#define Source_Action_Particules_included

#include <Source_base.h>
#include <Transport_Marqueur_FT.h>

/*! @brief Classe Source_Action_Particules Classe mere des classes designant une force exercee par le fluide sur une particule
 *
 *     (ou inversement)
 *     Ces forces interviennent dans le bilan de qdm des particules (one_way_coupling)
 *     et eventuellement dans l equation de Navier_Stokes (two_way_coupling) sous forme
 *     d une source interpolee (voir Source_Reaction_Particules)
 *  .SECTION
 *     Classe non instanciable.
 *
 *
 */
class Source_Action_Particules : public Source_base
{
  Declare_base(Source_Action_Particules);
public :

  void mettre_a_jour(double temps) override;
  inline const DoubleTab& vitesse_fluide() const;
  inline const DoubleTab& rho_fluide() const;
  inline const DoubleTab& visco_dyn_fluide() const;
  inline const DoubleTab& grad_pression() const;

  inline const DoubleTab& vitesse_particules() const;
  inline const DoubleTab& delta_v() const;
  inline const DoubleTab& rho_particules() const;
  inline const DoubleTab& diametre_particules() const;
  inline const DoubleTab& volumes_particules() const;

protected:

  void associer_zones(const Zone_dis& ,const Zone_Cl_dis& ) override;
  void associer_pb(const Probleme_base& ) override;


};

const DoubleTab& Source_Action_Particules::vitesse_fluide() const
{
  const Transport_Marqueur_FT& eq = ref_cast(Transport_Marqueur_FT,equation());
  const DoubleTab& vitesse = eq.vitesse_fluide();
  return vitesse;
}

const DoubleTab& Source_Action_Particules::rho_fluide() const
{
  const Transport_Marqueur_FT& eq = ref_cast(Transport_Marqueur_FT,equation());
  const DoubleTab& rho = eq.rho_fluide();
  return rho;
}

const DoubleTab& Source_Action_Particules::visco_dyn_fluide() const
{
  const Transport_Marqueur_FT& eq = ref_cast(Transport_Marqueur_FT,equation());
  const DoubleTab& visco_dyn = eq.visco_dyn_fluide();
  return visco_dyn;
}

const DoubleTab& Source_Action_Particules::grad_pression() const
{
  const Transport_Marqueur_FT& eq = ref_cast(Transport_Marqueur_FT,equation());
  const DoubleTab& gradient_P = eq.grad_pression();
  return gradient_P;
}


const DoubleTab& Source_Action_Particules::vitesse_particules() const
{
  const Transport_Marqueur_FT& eq = ref_cast(Transport_Marqueur_FT,equation());
  const DoubleTab& vitesse = eq.proprietes_particules().vitesse_particules();
  return vitesse;
}

const DoubleTab& Source_Action_Particules::delta_v() const
{
  const Transport_Marqueur_FT& eq = ref_cast(Transport_Marqueur_FT,equation());
  const DoubleTab& deltav = eq.proprietes_particules().delta_v();
  return deltav;
}

const DoubleTab& Source_Action_Particules::rho_particules() const
{
  const Transport_Marqueur_FT& eq = ref_cast(Transport_Marqueur_FT,equation());
  const DoubleTab& rho = eq.proprietes_particules().masse_vol_particules();
  return rho;
}

const DoubleTab& Source_Action_Particules::diametre_particules() const
{
  const Transport_Marqueur_FT& eq = ref_cast(Transport_Marqueur_FT,equation());
  const DoubleTab& diametre = eq.proprietes_particules().diametre_particules();
  return diametre;
}

const DoubleTab& Source_Action_Particules::volumes_particules() const
{
  const Transport_Marqueur_FT& eq = ref_cast(Transport_Marqueur_FT,equation());
  const DoubleTab& volume = eq.proprietes_particules().volume_particules();
  return volume;
}

#endif
