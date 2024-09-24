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
// File:        Paroi_DWF_hyd_VDF.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/VDF/Turbulence
// Version:     /main/14
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Paroi_DWF_hyd_VDF_included
#define Paroi_DWF_hyd_VDF_included


#include <Pb_MG.h>
#include <Paroi_std_hyd_VDF.h>
#include <Paroi_loi_WW_scal_VDF.h>
#include <Probleme_Couple.h>
#include <TRUST_Ref.h>




class Champ_Fonc_base;

/*! @brief class Paroi_DWF_hyd_VDF cette classe derive de turbulence_paroi_base
 *
 *      qui est la classe standard des lois de parois
 *      pour la turbulence
 *
 *      Surtout ne pas copier sur cette classe pour creer une loi de paroi
 *      avec un couplage !
 *
 *
 *
 */

///////////////////////////////////////////////////////////////////////////
//On conserve la classe car possible developpement dans le futur.
//Actuellement pas utilisee et ajout de commentaires pour
//compilation suite a integration de developpements Couplage_Multi_Echelles
///////////////////////////////////////////////////////////////////////////

class Paroi_DWF_hyd_VDF : public Paroi_std_hyd_VDF
{

  Declare_instanciable(Paroi_DWF_hyd_VDF);

public:


  int init_lois_paroi() override;

  /* calculer_hyd pour k-epsilon qui renvoie une erreur */
  int calculer_hyd(DoubleTab&) override;

  /* calculer_hyd pour L.E.S */
  int calculer_hyd(DoubleTab& , DoubleTab&) override;

protected :

  // Probleme multigrille pour connecter le probleme thermohydraulique avec le sous probleme fin en paroi
  Pb_MG pbMG;
  REF(Probleme_base) pb_fin;

  int thermo,CHT;
};


#endif
