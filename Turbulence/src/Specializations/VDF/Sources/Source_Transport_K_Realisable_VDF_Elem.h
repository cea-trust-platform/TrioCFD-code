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
// File:        Source_Transport_K_Realisable_VDF_Elem.h
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/14
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Source_Transport_K_Realisable_VDF_Elem_included
#define Source_Transport_K_Realisable_VDF_Elem_included


#define C21_DEFAULT_K_EPS_REALISABLE  1.9   // dans l'equation de k_eps
#define C3_DEFAULT_K_EPS_REALISABLE 1.0    // de transport de K et Eps source: Chabard et N3S

//
// .DESCRIPTION class Source_Transport_K_Realisable_VDF_Elem
//
// .SECTION voir aussi

#include <Source_Transport_K_Eps_VDF_Elem.h>
#include <Ref_Zone_Cl_VDF.h>
#include <Ref_Transport_Flux_Chaleur_Turbulente.h>
#include <Modele_Fonc_Realisable_base.h>
#include <Ref_Zone_dis.h>
#include <Zone_Cl_dis.h>
#include <Ref_Transport_K_ou_Eps_Realisable.h>
#include <Ref_Convection_Diffusion_Temperature.h>
#include <Ref_Convection_Diffusion_Concentration.h>
#include <Calcul_Production_K_VDF.h>
#include <Source_base.h>

class Probleme_base;
class Champ_Don_base;
class DoubleVect;
class DoubleTab;
class Champ_Face;
class Zone_dis;
class Zone_Cl_dis;
class Zone_Cl_VDF;
class Champ_Face;

#include <Source_Transport_Realisable_VDF_Elem_base.h>

class Source_Transport_K_Realisable_VDF_Elem : public Source_Transport_Realisable_VDF_Elem_base
{
  Declare_instanciable(Source_Transport_K_Realisable_VDF_Elem);
public :
  DoubleTab& ajouter(DoubleTab& ) const;
  void mettre_a_jour(double temps);
  virtual void associer_pb(const Probleme_base& );

protected :
  REF(Transport_K_ou_Eps_Realisable) eqn_k_Rea, eqn_eps_Rea;

private:
  void fill_coeff_matrice(const int , const DoubleTab& , const DoubleVect& , const DoubleVect& , double& , Matrice_Morse& ) const;
};

#endif /* Source_Transport_K_Realisable_VDF_Elem_included */
