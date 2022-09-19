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
// File:        Marqueur_FT.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////
#ifndef Marqueur_FT_included
#define Marqueur_FT_included

#include <Marqueur_Lagrange_base.h>
#include <Marching_Cubes.h>

/*! @brief classe Marqueur_FT La classe Marqueur_FT est une classe instanciable permettant de suivre des particules
 *
 *     -on attribue a cette classe un objet de type Maillage_FT_Disc pour
 *     beneficier de la structure gerant un suivi Lagrangien
 *     On ne retient que le suivi des sommets qui constituent les points suivis
 *
 * @sa Marqueur_Lagrange_base
 */
class Marqueur_FT : public Marqueur_Lagrange_base
{
  Declare_instanciable(Marqueur_FT);
public:

  void discretiser(const Probleme_base& pb, const  Discretisation_base& dis) override;
  inline Ensemble_Lagrange_base& ensemble_points() override;
  inline const Ensemble_Lagrange_base& ensemble_points() const override;
  void calculer_valeurs_champs() override;

protected:

  Maillage_FT_Disc ensemble_points_; //Ensemble des points suivis

};

inline Ensemble_Lagrange_base& Marqueur_FT::ensemble_points()
{
  return ensemble_points_;
  //temporaire
  exit();
}

inline const Ensemble_Lagrange_base& Marqueur_FT::ensemble_points() const
{
  return ensemble_points_;
  //temporaire
  exit();
}

#endif
