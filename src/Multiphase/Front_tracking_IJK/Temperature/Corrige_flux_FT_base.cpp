/****************************************************************************
 * Copyright (c) 2019, CEA
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 *modification, are permitted provided that the following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *this list of conditions and the following disclaimer in the documentation
 *and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its contributors
 *may be used to endorse or promote products derived from this software without
 *specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 *AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 *IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 *DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 *FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 *DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 *SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 *CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 *OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 *OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *****************************************************************************/
/////////////////////////////////////////////////////////////////////////////
//
// File      : Corrige_flux_FT.cpp
// Directory : $IJK_ROOT/src/Temperature
//
/////////////////////////////////////////////////////////////////////////////

#include <Corrige_flux_FT.h>
#include <DebogIJK.h>
#include <IJK_FT.h>
#include <IJK_Navier_Stokes_tools.h>
#include <IJK_Thermique.h>
#include <Intersection_Interface_ijk.h>
#include <Param.h>
#include <stat_counters.h>

Implemente_base(Corrige_flux_FT_base, "Corrige_flux_FT_base", Objet_U);

Sortie& Corrige_flux_FT_base::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  return os;
}

Entree& Corrige_flux_FT_base::readOn( Entree& is )
{
  Objet_U::readOn( is );
  return is;
}

void Corrige_flux_FT_base::initialize(const IJK_Splitting& splitting,
                                      const IJK_Field_double& field,
                                      const IJK_Interfaces& interfaces,
                                      const IJK_FT_double& ijk_ft)
{
  /*
   * TODO: est-ce que les pointeurs restent bien toujours les meme ?
   * Si oui je peux ne les initialiser qu'une fois, sinon il faudra
   * les mettre à jour.
   */
  interfaces_ = &interfaces;
  field_ = &field;
  splitting_ = &splitting;
  ref_ijk_ft_ = ijk_ft;
}

void Corrige_flux_FT_base::set_physical_parameters(const double rhocpl,
                                                   const double rhocpv,
                                                   const double ldal,
                                                   const double ldav)
{
  rhocp_l_ = rhocpl;
  rhocp_v_ = rhocpv;
  lda_l_ = ldal;
  lda_v_ = ldav;
}

bool Corrige_flux_FT_base::test_if_stencil_inclut_bout_interface_liquide() const
{
  const int& dir= parcours_.face();

  const double velocity = ref_ijk_ft_->get_velocity(
                          )[dir](parcours_.i(), parcours_.j(), parcours_.k());
  double stencil_inclut_interface = 1.;
  int decal = 0;
  if ( velocity > 0.)
    decal = -1;

  for (int c = 0; c < 3; c++)
    {
      const auto c_elem = parcours_.elem(c + decal);
      stencil_inclut_interface *= ref_ijk_ft_->itfce().I(
                                    c_elem[0],
                                    c_elem[1],
                                    c_elem[2]);
    }

  // On considère que ce n'est pas grave s'il n'y a qu'un tout petit bout d'interface
  return std::abs(1. - stencil_inclut_interface) > 0.05;
}

bool Corrige_flux_FT_base::test_if_stencil_inclut_bout_interface_vapeur() const
{
  const int& dir= parcours_.face();

  const double velocity = ref_ijk_ft_->get_velocity(
                          )[dir](parcours_.i(), parcours_.j(), parcours_.k());
  double stencil_inclut_interface = 1.;
  int decal = 0;
  if ( velocity > 0.)
    decal = -1;

  for (int c = 0; c < 3; c++)
    {
      const auto c_elem = parcours_.elem(c + decal);
      stencil_inclut_interface *=  (
                                     1. - ref_ijk_ft_->itfce().I(
                                       c_elem[0],
                                       c_elem[1],
                                       c_elem[2]));
    }

  // On considère que ce n'est pas grave s'il n'y a qu'un tout petit bout d'interface
  return std::abs(1. - stencil_inclut_interface) > 0.05;
}

