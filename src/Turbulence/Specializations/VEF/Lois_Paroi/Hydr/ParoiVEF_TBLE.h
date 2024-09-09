/****************************************************************************
* Copyright (c) 2019, CEA
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
// File:        ParoiVEF_TBLE.h
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Lois_Paroi/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#ifndef ParoiVEF_TBLE_included
#define ParoiVEF_TBLE_included

#include <Paroi_hyd_base_VEF.h>
#include <Domaine_Cl_dis.h>
#include <TRUST_Vector.h>
#include <Eq_couch_lim.h>
#include <Milieu_base.h>
#include <Paroi_TBLE_QDM.h>
#include <Domaine_dis.h>
#include <Domaine_Cl_dis.h>
#include <Paroi_log_QDM.h>

class Champ_Fonc_base;
class Param;

/*! @brief CLASS: ParoiVEF_TBLE
 *
 * .SECTION  voir aussi
 *  Turbulence_paroi_base
 *
 */
class ParoiVEF_TBLE : public Paroi_hyd_base_VEF, public Paroi_TBLE_QDM, public Paroi_log_QDM
{
  Declare_instanciable_sans_constructeur(ParoiVEF_TBLE);
public:

  ParoiVEF_TBLE();
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

protected:


private:

  void traitement_keps(DoubleTab& tab_k_eps, int num_face, const IntTab& face_voisins,
                       const IntTab& elem_faces, int nfac, double dist, double d_plus, double d_visco, double u_star);
  void traitement_keps_BiK(DoubleTab& tab_k, DoubleTab& tab_eps, int num_face, const IntTab& face_voisins,
                           const IntTab& elem_faces, int nfac, double dist, double d_plus, double d_visco, double u_star);
  virtual int calculer_k_eps(double& , double& , double , double , double , double);

  int calculer_stats();
//  int alpha_cv;
  DoubleTab f_tang_moy;
  double alpha; // pour la moyenne glissante du terme source

};

#endif
