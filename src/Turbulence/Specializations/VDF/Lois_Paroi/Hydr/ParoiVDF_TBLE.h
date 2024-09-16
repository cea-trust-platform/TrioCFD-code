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
// File:        ParoiVDF_TBLE.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Lois_Paroi/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#ifndef ParoiVDF_TBLE_included
#define ParoiVDF_TBLE_included

#include <Paroi_hyd_base_VDF.h>
#include <Domaine_Cl_dis.h>
#include <TRUST_Vector.h>
#include <Eq_couch_lim.h>
#include <Milieu_base.h>
#include <Paroi_TBLE_QDM.h>
#include <Paroi_log_QDM.h>

#include <Domaine_Cl_dis.h>
class Champ_Fonc_base;
class Param;

/*! @brief CLASS: ParoiVDF_TBLE
 *
 * .SECTION  voir aussi
 *  Turbulence_paroi_base
 *
 */
class ParoiVDF_TBLE : public Paroi_hyd_base_VDF, public Paroi_TBLE_QDM, public Paroi_log_QDM
{
  Declare_instanciable_sans_constructeur(ParoiVDF_TBLE);
public:

  ParoiVDF_TBLE();
  void set_param(Param& param) override;
  int lire_motcle_non_standard(const Motcle&, Entree&) override;
  int init_lois_paroi() override;
  int calculer_hyd(DoubleTab& ) override;
  int calculer_hyd_BiK(DoubleTab& , DoubleTab& ) override;
  int calculer_hyd(DoubleTab& , DoubleTab&) override;
  int calculer_hyd(DoubleTab& tab1,int isKeps,DoubleTab& tab2);

  void imprimer_ustar(Sortie&) const override;

  //OC 02/2006 :methodes de reprise/sauvegarde pour TBLE. Pour l'instant les donnees TBLE sont stockes dans des fichiers differents du .sauv => a voir si on met tout cela dans le .sauv a terme
  int sauvegarder(Sortie&) const override ;
  int reprendre(Entree& ) override ;
  const Probleme_base& getPbBase() const override ;

private:

  int alpha_cv;

  int calculer_stats();
  void calculer_convection(int num_face, int face1, int face2, int face3, int face4, int elem,
                           int ndeb, int nfin, int ori, double gradient_de_pression0, double ts0, double gradient_de_pression1, double ts1);

  IntVect corresp; //Correspondance compteur-num_face parietale

};

///////////////////////////////////////////////////////////
//
//  Fonctions inline de la classe ParoiVDF_TBLE
//
///////////////////////////////////////////////////////////

#endif
