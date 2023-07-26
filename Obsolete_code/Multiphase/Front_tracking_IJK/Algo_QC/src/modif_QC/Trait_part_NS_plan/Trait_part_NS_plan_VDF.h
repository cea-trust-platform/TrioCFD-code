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
// File      : Trait_part_NS_plan_VDF.h
// Directory : $NEW_ALGO_QC_ROOT/src/modif_QC/Trait_part_NS_plan
//
/////////////////////////////////////////////////////////////////////////////


#ifndef Traitement_particulier_NS_plan_VDF_inclus
#define Traitement_particulier_NS_plan_VDF_inclus

#include <Trait_part_NS_plan.h>

class Navier_Stokes_Turbulent;

/*! @brief classe Traitement_particulier_NS_plan_VDF Cette classe permet de faire les traitements particuliers
 *
 *      pour le calcul d'un canal plan :
 *          * conservation du debit
 *          * calculs de moyennes
 *
 *
 * @sa Navier_Stokes_Turbulent, Traitement_particulier_base,, Traitement_particulier_VDF
 */
class Traitement_particulier_NS_plan_VDF : public Traitement_particulier_NS_plan
{
  Declare_instanciable(Traitement_particulier_NS_plan_VDF);

public :

  Entree& lire(const Motcle& , Entree& );
  Entree& lire(Entree& );

protected :


  void calculer_valeur_spatiale_vitesse_rho_mu(DoubleTab&,DoubleTab&,double) ;

  double distance_elem(int ,int ) const; // methode pour le calcul de la distance entre deux voisins

  double derivee_un_elem(const DoubleTab& ,int ,int) const; // methode d'ordre deux pour les derivees de grandeurs locialisees au centres des elements
  double derivee_un_elem( DoubleTab& ,int ,int) const; // methode d'ordre deux pour les derivees de grandeurs locialisees au centres des elements
  double derivee_un_elem( DoubleTab& ,int,int ,int) const; // methode d'ordre deux pour les derivees de vitesse localisees au centres des elements

  double derivee_deux_vitesse_coli(const DoubleTab& ,int ,int) const; // methode qui calcule la derivee d'ordre deux lorsque la vitesse et le sens de derivation sont paralleles.
  double derivee_deux_vitesse_nonco(DoubleTab&,int ,int ,int) const; // methode qui calcule la derivee d'ordre deux selon x ou z  lorsque la vitesse et le sens de derivation sont differents.

};

#endif
