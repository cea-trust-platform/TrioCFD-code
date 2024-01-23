/****************************************************************************
* Copyright (c) 2021, CEA
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
// File:        Diffusion_supplementaire_echelle_temp_turb_VDF.h
// Directory:   $TRUST_ROOT/src/ThHyd/Multiphase/PolyMAC_P0
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Diffusion_supplementaire_echelle_temp_turb_VDF_included
#define Diffusion_supplementaire_echelle_temp_turb_VDF_included

#include <Source_Diffusion_supplementaire_echelle_temp_turb.h>

class Correlation;

/*! @brief Classe Diffusion_supplementaire_echelle_temp_turb_PolyMAC_P0 Cette classe implemente dans PolyMAC_P0 la diffusion supplementaire venant de l'introduction du temps tau
 *
 *
 *
 * @sa Operateur_PolyMAC_P0_base Operateur_base
 */
class Diffusion_supplementaire_echelle_temp_turb_VDF: public Source_Diffusion_supplementaire_echelle_temp_turb
{
  Declare_instanciable(Diffusion_supplementaire_echelle_temp_turb_VDF);
public :
  void dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl = {}) const override { }; // Explicit for now
  void ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl = {}) const override;

protected :
  void face_to_elem(const Domaine_VF& domaine, const DoubleTab& tab_faces,DoubleTab& tab_elems) const;

  double limiter_ = 5 ;
  double limiter_tau_ = 1.e-6;
};

#endif
