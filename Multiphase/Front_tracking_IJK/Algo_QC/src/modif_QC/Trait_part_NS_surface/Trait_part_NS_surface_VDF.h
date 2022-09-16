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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Trait_part_NS_surface_VDF.h
// Directory : $NEW_ALGO_QC_ROOT/src/modif_QC/Trait_part_NS_surface
//
/////////////////////////////////////////////////////////////////////////////


#ifndef Traitement_particulier_NS_surface_VDF_inclus
#define Traitement_particulier_NS_surface_VDF_inclus

#include <Trait_part_NS_surface.h>

class Navier_Stokes_Turbulent;

/*! @brief classe Traitement_particulier_NS_surface_VDF Cette classe permet de faire les traitements particuliers
 *
 *      pour le calcul d'un canal surface :
 *          * conservation du debit
 *          * calculs de moyennes
 * 
 *
 * @sa Navier_Stokes_Turbulent, Traitement_particulier_base,, Traitement_particulier_VDF 
 */
class Traitement_particulier_NS_surface_VDF : public Traitement_particulier_NS_surface
{
  Declare_instanciable(Traitement_particulier_NS_surface_VDF);

public :

  Entree& lire(const Motcle& , Entree& );
  Entree& lire(Entree& );

protected :

  void remplir_XYZ(DoubleVect& ,DoubleVect& ,DoubleVect&, int& , int&, int&,IntTab&) const;

  void recuperation_grandeurs(DoubleTab&) const;

};

#endif
