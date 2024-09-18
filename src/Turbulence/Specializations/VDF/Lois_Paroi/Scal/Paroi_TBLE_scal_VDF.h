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
// File:        Paroi_TBLE_scal_VDF.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Lois_Paroi/Scal
//
//////////////////////////////////////////////////////////////////////////////

#ifndef Paroi_TBLE_scal_VDF_included
#define Paroi_TBLE_scal_VDF_included

#include <Paroi_std_scal_hyd_VDF.h>
#include <Domaine_Cl_dis.h>
#include <Eq_couch_lim.h>
#include <TRUST_Vector.h>
#include <Milieu_base.h>
#include <Paroi_TBLE_QDM_Scal.h>
#include <TRUST_Ref.h>

#include <Domaine_Cl_dis.h>
class Convection_Diffusion_std;
class Champ_Fonc_base;
class Paroi_TBLE_QDM;

/*! @brief CLASS: Paroi_ODVM_scal_VDF
 *
 * .SECTION  voir aussi
 *  Turbulence_paroi_base
 *
 */
class Paroi_TBLE_scal_VDF : public Paroi_std_scal_hyd_VDF, public Paroi_TBLE_QDM_Scal
{
  Declare_instanciable_sans_constructeur(Paroi_TBLE_scal_VDF);

public:
  int calculer_scal(Champ_Fonc_base& ) override;
  int init_lois_paroi() override;
  void imprimer_nusselt(Sortie&) const override;

  //OC 02/2006 :methodes de reprise/sauvegarde pour TBLE. Pour l'instant les donnees TBLE sont stockes dans des fichiers differents du .sauv => a voir si on met tout cela dans le .sauv a terme
  int sauvegarder(Sortie&) const override ;
  int reprendre(Entree& ) override ;
  const Probleme_base& getPbBase() const override ;
  Paroi_TBLE_QDM& getLoiParoiHydraulique() override;
  MuLambda_TBLE_base& getMuLambda() override;
  Eq_couch_lim& get_eq_couch_lim(int) override;

private:
  int alpha_cv;
  IntVect corresp; //Correspondance compteur-num_face parietale

  int calculer_stats();

  void calculer_convection(int compteur_faces_paroi, int elem, int ndeb, int nfin, int ori, double ts);


};


#endif
