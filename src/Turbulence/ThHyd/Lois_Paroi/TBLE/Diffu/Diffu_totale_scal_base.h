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
// File:        Diffu_totale_scal_base.h
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Lois_Paroi/TBLE/Diffu
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Diffu_totale_scal_base_included
#define Diffu_totale_scal_base_included
#include <Diffu_totale_base.h>
#include <TRUST_Ref.h>

class Eq_couch_lim ;

/*! @brief Classe Diffu_totale_scal_base Classe abstraite calculant la diffusivite totale (somme diffusivite
 *
 *     moleculaire et diffusivite turbulente) dans les equations de
 *     couche limite simplifiees necessaires
 *     a l'utilisation des lois de parois de type Balaras
 *
 *
 */
class Diffu_totale_scal_base : public Diffu_totale_base
{
  Declare_base_sans_constructeur(Diffu_totale_scal_base);

public :


  OBS_PTR(Diffu_totale_base) visco_tot;

  Diffu_totale_scal_base() ;
  void set_visco_tot(Diffu_totale_base& a)
  {
    visco_tot=a;
  };

  double getPrandtl()
  {
    return Prandtl;
  };
  void setPrandtl(double p)
  {
    Prandtl=p;
  };


private:
  double Prandtl;

};


#endif
