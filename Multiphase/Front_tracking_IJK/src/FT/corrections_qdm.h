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
/////////////////////////////////////////////////////////////////////////////
//
// File      : corrections_qdm.h
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#ifndef consigne_initiale_included
#define consigne_initiale_included

#include <FixedVector.h>
#include <IJK_Field.h>
#include <Objet_U.h>
#include <string>
#include <iostream>
#include <math.h>

/*! @brief : class consigne_initiale
 *
 *  <  S'OCCUPE DE CALCULER v_cible_ et qdm_cible_ UNIQUEMENT
 *                 qdm_cible = qdm_du_premier_pas_de_temps
 *                 TODO v_cible_ : est-ce que je le sort ?
 *  >
 *
 *
 *
 */
class consigne_initiale : public Objet_U
{

  Declare_instanciable_sans_constructeur( consigne_initiale ) ;

public :
  consigne_initiale();
  void initialise();
  void set_time(int time_iteration) {time_iteration_ = time_iteration;};
  double get_qdm_cible() {return qdm_cible_;};
  void compute_qdm_cible(double qdm_initiale);

  int get_need_for_vl_vv() {return 0;};
  int get_need_for_rho_l() {return 0;};
protected :
  int time_iteration_;
  double qdm_cible_;
};

#endif /* consigne_initiale_included */

/////////////////////////////////////////////////////////////////////////////////////////
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////


#ifndef cible_donnee_included
#define cible_donnee_included

#include <FixedVector.h>
#include <IJK_Field.h>
#include <Objet_U.h>
#include <string>
#include <iostream>
#include <math.h>

/*! @brief : class cible_donnee
 *
 *  <  S'OCCUPE DE CALCULER v_cible_ et qdm_cible_ UNIQUEMENT
 *                 qdm_cible = alpha_l rho_l v_cible_
 *                 v_cible_ : constante utilisateur ou .sauv
 *  >
 *
 *
 *
 */
class cible_donnee : public Objet_U
{

  Declare_instanciable_sans_constructeur( cible_donnee ) ;

public :
  cible_donnee();
  void set_time() {/* does nothing */;};
  void set_rho_l(double rho_liq) {rho_liq_ = rho_liq;};
  double get_v_cible() {return v_cible_;};
  double get_qdm_cible() {return qdm_cible_;};
  void compute_qdm_cible(double alpha_l);

  int get_need_for_vl_vv() {return 0;};
  int get_need_for_rho_l() {return 1;};
protected :
  double rho_liq_;
  double v_cible_;
  double qdm_cible_;
};

#endif /* cible_donnee_included */

/////////////////////////////////////////////////////////////////////////////////////////
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

#ifndef moyenne_par_morceaux_included
#define moyenne_par_morceaux_included

#include <FixedVector.h>
#include <IJK_Field.h>
#include <Objet_U.h>
#include <string>
#include <iostream>
#include <math.h>

/*! @brief : class moyenne_par_morceaux
 *
 *  <S'OCCUPE DE CALCULER v_cible_ et qdm_cible_ UNIQUEMENT
 *   qdm_cible = alpha_l rho_l v_cible_
 *   Pour (n-1)*duree_morceaux < t < (n)*duree_morceaux :
 * 	  v_cible_ = 1/duree_morceau sum_{morceau_precedent} v_liq-v_vap
 *     qdm_cible_ = 1/duree_morceau sum_{morceau_precedent} a_l rho_l * (v_liq-v_vap)
 *   Pour t = n*duree_morceaux
 *     n <-- n+1 ...
 *   ATTENTION AUX DIMENSION : v_cible_ vaut alpha_l*v_liq (03.03.22)
 *   >
 *
 *
 *
 */
class moyenne_par_morceaux : public Objet_U
{

  Declare_instanciable_sans_constructeur( moyenne_par_morceaux ) ;

public :
  moyenne_par_morceaux();
  void initialise();
  void set_time(double time, double time_step, int time_iteration);
  void set_rho_l(double rho_liq) {rho_liq_ = rho_liq;};
  void compute_qdm_cible(double alpha_l, double vitesse_relative);
  double get_qdm_cible() {return qdm_cible_;};

  int get_need_for_vl_vv() {return 1;};
  int get_need_for_rho_l() {return 1;};

protected :
  // Fetch from IJK_FT
  double timestep_;
  double time_;
  int time_iteration_;
  double rho_liq_;
  // Only here
  //double duree_morceau_effective_;
  double duree_morceau_demandee_;
  double last_update_time_;
  double v_init_guess_;
  double qdm_cible_;
  double v_cible_;
  double qdm_integ_partielle_;
};

#endif /* moyenne_par_morceaux_included */

/////////////////////////////////////////////////////////////////////////////////////////
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

#ifndef moyenne_glissante_included
#define moyenne_glissante_included

#include <FixedVector.h>
#include <IJK_Field.h>
#include <Objet_U.h>
#include <string>
#include <iostream>
#include <math.h>

/*! @brief : class moyenne_glissante
 *
 *  < S'OCCUPE DE CALCULER v_cible_ et qdm_cible_ UNIQUEMENT
 *    qdm_cible = alpha_l rho_l v_cible_
 *    Pour t < duree_demandee_morceau_glissant_ :
 * 		  v_cible_ = 1/t sum_{t} v_liq-v_vap
 * 	 Pour t >= duree_demandee_morceau_glissant_ :
 *         v_cible_ = 1/(t-tm) sum_{tm;t} v_liq-v_vap
 * 	      tm = duree_effective_morceau_glissant_
 *   ATTENTION AUX DIMENSIONS : v_cible_ = \ol{u}^liq - \ol{u}^vap
 *   >
 *
 *
 */
class moyenne_glissante : public Objet_U
{

  Declare_instanciable_sans_constructeur( moyenne_glissante ) ;

public :

  moyenne_glissante();
  void initialise();
  void set_time(double time, double time_step, int time_iteration);
  void set_rho_l(double rho_liq) {rho_liq_ = rho_liq;};
  double compute_v_cible();
  double get_qdm_cible() {return qdm_cible_;};
  double get_v_cible() {return v_cible_;};
  void compute_qdm_cible(double alpha_liq);
  void compute_v_cible(double vitesse_relative);

  int get_need_for_vl_vv() {return 1;};
  int get_need_for_rho_l() {return 1;};

protected :
  int list_index_;
  int offset_list_index_;
  ArrOfDouble liste_instants_;
  ArrOfDouble liste_v_cible_dl_;
  ArrOfDouble liste_v_cible_;
  double rho_liq_;

  int tstep_;
  double time_;
  double timestep_;
  double v_cible_;
  double qdm_cible_;
  double duree_demandee_morceau_glissant_;
  double duree_effective_morceau_glissant_;
};

#endif /* moyenne_glissante_included */

/////////////////////////////////////////////////////////////////////////////////////////
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

#ifndef correction_one_direction_included
#define correction_one_direction_included

#include <FixedVector.h>
#include <IJK_Field.h>
#include <Objet_U.h>
#include <string>
#include <iostream>
#include <math.h>

/*! @brief : class correction_one_direction
 *
 *  < Correction : "vitesse_corrigee_ = vitesse - vitesse_correction_"
 *
 *    A partir de qdm_cible, calcule vitesse_correction_ et vitesse_corrigee_ pour une direction
 *
 *                 u_i * = u_i - {u_correction}_i
 *                 correct_velocity_ = vel_ijk_t - value_correction_
 * 				  value_correction_= ({rho u}_i _moyen - qdm_cible)  /  rho_moyen
 *
 *    qdm_cible : calculee par une des sous-classes.
 *   >
 *
 *
 */
class correction_one_direction : public Objet_U
{

  Declare_instanciable_sans_constructeur( correction_one_direction ) ;

public :
  correction_one_direction();
  void set_correction(correction_one_direction& correction_in);
//  initialise_perp_g();
//  Entree& interpreter(Entree&);
  enum type_correction { CIBLE_CONSTANTE, MOYENNE_PAR_MORCEAUX, MOYENNE_GLISSANTE, REGIME_ETABLI, CONSIGNE_INITIALE };
  type_correction get_type_corr() const;
  double get_value() {return value_correction_;};
  void set_time_for_correction(double time, double time_step, int time_iteration);
  void set_rho_moyen_alpha_l_for_correction(double rho_moyen, double alpha_l);
  void set_rho_vel_moyen_for_correction(double rho_vel_moyen);
  void set_mean_values_for_correction(double rho_vel_moyen, double rho_moyen, double alpha_l);
  void compute_correction_value();
  void compute_correct_velocity(double vel_ijk_t);
  double get_correct_velocity() {return correct_velocity_;};
  double get_velocity_correction() {return value_correction_;};
  double get_correction_value() {return qdm_cible_;};

  int get_need_for_rho_liq();
  int get_need_for_vit_rel();
  int get_need_to_compute_correction_value();

  void set_rho_liq(double rho_liquide);
  void set_vl_vv(double vl_vv);

  cible_donnee get_cible_constante() const;
  moyenne_par_morceaux get_moyenne_par_morceaux() const;
  moyenne_glissante get_moyenne_glissante() const;
  consigne_initiale get_consigne_initiale() const;

protected :
  double value_correction_;
  double correct_velocity_;
  int type_corr_;
  int i_need_vitesse_relative_;
  /* Compute from IJK */
  double rho_vel_moyen_;
  double rho_moyen_;
  double alpha_l_;
  /* Compute here */
  double vitesse_relative_;
  double qdm_cible_;

  cible_donnee parametres_cible_constante_;
  moyenne_par_morceaux parametres_moyenne_par_morceaux_;
  moyenne_glissante parametres_moyenne_glissante_;
  consigne_initiale parametres_consigne_initiale_;
}
;

#endif /* correction_one_direction_included */

/////////////////////////////////////////////////////////////////////////////////////////
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

#ifndef corrections_qdm_included
#define corrections_qdm_included

#include <FixedVector.h>
#include <IJK_Field.h>
#include <Objet_U.h>
#include <string>
#include <iostream>
#include <math.h>

/*! @brief : class corrections_qdm
 *
 *  <Corrections to comply with the mean momentum budget in the three directions.>
 *
 *
 *
 */
class corrections_qdm : public Objet_U
{

  Declare_instanciable_sans_constructeur( corrections_qdm ) ;

public :
  corrections_qdm();
  void set_corrections(correction_one_direction& in_corr_x_ ,
                       correction_one_direction& in_corr_y_ ,
                       correction_one_direction& in_corr_z_);
  void set_time(double time, double time_step, int time_iteration);
  void set_vitesse_relative(int direction, double vitesse_relative);
  void set_rho_liquide(double rho_liquide);
  void set_rho_moyen_alpha_l(double rho_moyen, double alpha_l);
  void set_rho_vel_moyen(int direction, double rho_vel_moyen);
  void set_mean_values_for_corrections(double rho_vel_moyen, double rho_moyen, double alpha_l);
  void compute_correction_value_one_direction(int dir);
  void compute_correct_velocity_one_direction(int dir, double vel_ijk_t);
  Vecteur3 get_correct_velocities();    // u - v_corr
  Vecteur3 get_velocity_corrections();  // mean(rho u) - qdm_cible / mean(rho)
  Vecteur3 get_correction_values();     // qdm_cible
  double get_correct_velocitiy_one_direction(int dir);
  int get_need_for_vitesse_relative(int direction);
  int get_need_to_compute_correction_value_one_direction(int direction);

  enum type_dict_ { GB, GR, NONE_IJK };
  type_dict_ get_type_() const;

  int is_type_gb() const;
  int is_type_gr() const;
  int is_type_none() const;
protected :
  int type_;
  correction_one_direction correction_x_;
  correction_one_direction correction_y_;
  correction_one_direction correction_z_;
};

#endif /* corrections_qdm_included */

/////////////////////////////////////////////////////////////////////////////////////////
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
