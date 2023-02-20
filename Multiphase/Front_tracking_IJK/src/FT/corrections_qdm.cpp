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
// File      : corrections_qdm.cpp
// Directory : $IJK_ROOT/src/FT
//
/////////////////////////////////////////////////////////////////////////////

#include <corrections_qdm.h>
//#include <Interprete_bloc.h>
#include <Param.h>
#include <string>
#include <iostream>
#include <math.h>

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////

Implemente_instanciable_sans_constructeur( consigne_initiale, "consigne_initiale", Objet_U ) ;

consigne_initiale::consigne_initiale()
{
  time_iteration_=0;
  qdm_cible_=0;
}

Sortie& consigne_initiale::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  os  << "{\n"
      << " qdm_cible " << qdm_cible_ << "\n";
  os << " }\n" ;
  return os;
}

Entree& consigne_initiale::readOn( Entree& is )
{
  Objet_U::readOn( is );
  Param param(que_suis_je());
  param.ajouter("qdm_cible",&qdm_cible_ );
  param.lire_avec_accolades(is);
  return is;
}

void consigne_initiale::initialise()
{
  time_iteration_=0;
  qdm_cible_=0;
}

// TODO
void consigne_initiale::compute_qdm_cible(double qdm_init)
{
  if (time_iteration_ == 0)
    qdm_cible_ = qdm_init;
}

/////////////////////////////////////////////////////////////////////////////////////////
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

Implemente_instanciable_sans_constructeur( cible_donnee, "cible_donnee", Objet_U ) ;

cible_donnee::cible_donnee():
  rho_liq_(0.),
  v_cible_(0.),
  qdm_cible_(0.)
{
}

Sortie& cible_donnee::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  os  << "{\n"
      << " vitesse_cible " << v_cible_ << "\n";
  os << " }\n" ;
  return os;
}

Entree& cible_donnee::readOn( Entree& is )
{
  Objet_U::readOn( is );
  Param param(que_suis_je());
  param.ajouter("vitesse_cible",&v_cible_ );
  param.lire_avec_accolades(is);
  return is;
}

void cible_donnee::compute_qdm_cible(double alpha_l)
{
  qdm_cible_ = alpha_l*rho_liq_*v_cible_;
}

/////////////////////////////////////////////////////////////////////////////////////////
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

// IL FAUT ENCORE BIEN INITIALISER LA CLASSE EN GERANT LES SAUVEGARDES, REPRISES //
Implemente_instanciable_sans_constructeur( moyenne_par_morceaux, "moyenne_par_morceaux", Objet_U ) ;

moyenne_par_morceaux::moyenne_par_morceaux() // @suppress("Class members should be properly initialized")
{
  timestep_=0;
  time_=0;
  rho_liq_=0;

  duree_morceau_demandee_=0;
  last_update_time_=0;
  v_init_guess_=0;
  qdm_cible_=0;
  v_cible_=0;
  qdm_integ_partielle_=0;
}


Sortie& moyenne_par_morceaux::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  os  << "{\n";
  /* Sorties utilisateur */
  os  << " duree_morceau " << duree_morceau_demandee_ << "\n"
      << " v_init_guess " << v_init_guess_ << "\n";
  /* Sorties reprises */
  os << " qdm_integ_partielle " << qdm_integ_partielle_ << "\n"
     << " last_update_time " << last_update_time_ << "\n"
     << " qdm_reprise " << qdm_cible_ << "\n"
     << " v_reprise " << v_cible_ << "\n";
  os << " }\n" ;
  return os;
  return os;
}

Entree& moyenne_par_morceaux::readOn( Entree& is )
{
  Objet_U::readOn( is );
  Param param(que_suis_je());
  /* Entrees utilisateur */
  param.ajouter("duree_morceau",&duree_morceau_demandee_);
  param.ajouter("v_init_guess",&v_init_guess_);
  /* Entrees reprise */
  param.ajouter("qdm_integ_partielle",&qdm_integ_partielle_);
  param.ajouter("last_update_time",&last_update_time_);
  param.ajouter("qdm_reprise",&qdm_cible_);
  param.ajouter("v_reprise",&v_cible_);
  param.lire_avec_accolades(is);
  return is;
}

void moyenne_par_morceaux::initialise()
{
  time_iteration_ = 0;
  timestep_=0;
  time_=0;
  rho_liq_=0;

  duree_morceau_demandee_=0;
  last_update_time_=0;
  v_init_guess_=0;
  qdm_cible_=0;
  v_cible_=0;
  qdm_integ_partielle_=0;
}
void moyenne_par_morceaux::set_time(double time, double time_step, int time_iteration)
{
  time_ = time;
  timestep_ = time_step;
  time_iteration_ = time_iteration;
}

void moyenne_par_morceaux::compute_qdm_cible(double alpha_l, double vitesse_relative)
{
  Cout << "time_" << time_ << finl;
  Cout << "duree_morceau_demandee_" << duree_morceau_demandee_ << finl;
  Cout << "last_update_time_" << last_update_time_ << finl;
  Cout << "qdm_integ_partielle_" << qdm_integ_partielle_ << finl;
  Cout << "alpha_l " << alpha_l << finl;
  Cout << "vitesse_relative " << vitesse_relative << finl;
  Cout << "timestep_ " << timestep_ << finl;
  if (time_ < duree_morceau_demandee_ + last_update_time_)
    {
      /* *****************************************************
       * Construit l'integrale pour la cible de qdm
       * *****************************************************/
      Cout << "In the if " << time_<< "/ "<< duree_morceau_demandee_ + last_update_time_ << finl;
      qdm_integ_partielle_ += alpha_l*vitesse_relative*timestep_;
      //duree_morceau_effective_ += timestep_;
    }
  else
    {
      /* ******************************************************************************
       * Passage d'un morceau a l'autre :
       * Mise a jour de la cible de qdm et remise a zero des grandeurs de construction
       * ******************************************************************************/

      /* Finalize integral and update target value */
      Cout << "time_ - last_update_time_" << time_ - last_update_time_ << finl;
      Cout << "rho_liq_" << rho_liq_ << finl;
      Cout << "qdm_integ_partielle_" << qdm_integ_partielle_ << finl;
      Cout << "last_update_time_" << last_update_time_ << finl;
      qdm_integ_partielle_ *= rho_liq_;
      qdm_integ_partielle_ /= (time_ - last_update_time_);
      v_cible_ = qdm_cible_/rho_liq_;
      qdm_cible_ = qdm_integ_partielle_;
      /* Refresh for the next mean */
      qdm_integ_partielle_ = 0;
      last_update_time_ = time_;
    }

}

/////////////////////////////////////////////////////////////////////////////////////////
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

// IL FAUT ENCORE BIEN INITIALISER LA CLASSE EN GERANT LES SAUVEGARDES, REPRISES //
Implemente_instanciable_sans_constructeur( moyenne_glissante, "moyenne_glissante", Objet_U ) ;

moyenne_glissante::moyenne_glissante() // @suppress("Class members should be properly initialized")
{
  list_index_=0;
  offset_list_index_=0; /*è MOLTO IMPORTANTE lasciare il valore 0*/
  liste_instants_.resize_array(0);
  liste_instants_.set_smart_resize(1);
  liste_v_cible_dl_.resize_array(0);
  liste_v_cible_dl_.set_smart_resize(1);
  liste_v_cible_.resize_array(0);
  liste_v_cible_.set_smart_resize(1);
  rho_liq_=0;

  tstep_=0;
  time_=0;
  timestep_=0;
  v_cible_=0;
  qdm_cible_=0;
  duree_demandee_morceau_glissant_=0;
  duree_effective_morceau_glissant_=0;
}

Sortie& moyenne_glissante::printOn( Sortie& os ) const
{
  Objet_U::printOn( os );
  os  << "{\n";
  /* Sorties utilisateurs */
  os  << " duree_intervalle " << duree_demandee_morceau_glissant_ << "\n"
      << " v_init_guess " << v_cible_ << "\n";
  /* Sorties reprise*/
  os << " list_index " << list_index_ << "\n"
     << " offset_list_index " << list_index_+1 << "\n"
     << " liste_instants " << liste_instants_ << "\n"
     << " liste_v_cible_dl " << liste_v_cible_dl_ << "\n"
     << " liste_v_cible " << liste_v_cible_ << "\n"
     << " duree_effective_morceau_glissant " << duree_effective_morceau_glissant_ << "\n"
     << " v_cible_en_cours " << v_cible_ << "\n";
  os << " }\n" ;
  return os;
}

Entree& moyenne_glissante::readOn( Entree& is )
{
  Objet_U::readOn( is );
  Param param(que_suis_je());
  /* Entrees utilisateur */
  param.ajouter("duree_intervalle", &duree_demandee_morceau_glissant_);
  param.ajouter("v_init_guess", &v_cible_);
  /* Entrees reprises */
  param.ajouter("list_index", &list_index_);
  param.ajouter("offset_list_index", &offset_list_index_);
  param.ajouter("liste_instants", &liste_instants_);
  param.ajouter("liste_v_cible_dl", &liste_v_cible_dl_);
  param.ajouter("liste_v_cible", &liste_v_cible_);
  param.ajouter("duree_effective_morceau_glissant", &duree_effective_morceau_glissant_);
  param.ajouter("v_cible_en_cours", &v_cible_);
  param.lire_avec_accolades(is);
  return is;
}

void moyenne_glissante::initialise()
{
  list_index_=0;
  offset_list_index_=0; /*è MOLTO IMPORTANTE lasciare il valore 0*/
  liste_instants_.resize_array(0);
  liste_v_cible_dl_.resize_array(0);
  liste_v_cible_.resize_array(0);
  rho_liq_=0;

  tstep_=0;
  time_=0;
  timestep_=0;
  v_cible_=0;
  qdm_cible_=0;
  duree_demandee_morceau_glissant_=0;
  duree_effective_morceau_glissant_=0;
}

void moyenne_glissante::set_time(double time, double time_step, int time_iteration)
{
  time_ = time;
  tstep_ = time_iteration;
  timestep_ = time_step;
}
void moyenne_glissante::compute_v_cible(double vitesse_relative)
{
  /* ******************************************
   * time_ <= duree_morceau_glissant_ : la moyenne ne glisse pas encore
   * *****************************************/
  Cout << "time_" << time_ << finl;
  Cout << "duree_demandee_morceau_glissant_" << duree_demandee_morceau_glissant_ << finl;
  Cout << "list_index_ " << list_index_ << finl;
  Cout << "tstep_ " << tstep_ << finl;
  Cout << "offset_list_index_ " << offset_list_index_ << finl;
  Cout << "liste_instants_.size_array()" << liste_instants_.size_array() << finl ;
  Cout << "duree_effective_morceau_glissant_" << duree_effective_morceau_glissant_ << finl;
  Cout << "v_cible_" << v_cible_ << finl;

  if (time_ <= duree_demandee_morceau_glissant_)
    {
      Cout << "DEB : time_ <= duree_demandee_morceau_glissant_" << finl;
      /* Construction des listes pour la moyenne */
      liste_instants_.append_array(time_);
      double v_cible_dl = vitesse_relative*timestep_;
      liste_v_cible_dl_.append_array(v_cible_dl);
      duree_effective_morceau_glissant_ = time_;


      if (time_ != 0)
        {
          Cout << "time_ != 0 " << time_ << finl;
          /* Construction de la moyenne au fur et a mesure */
          const int i_last = liste_instants_.size_array()-1;
          v_cible_ = v_cible_*liste_instants_[i_last-1] + liste_v_cible_dl_[i_last]*timestep_;
          v_cible_ /= duree_effective_morceau_glissant_;
          Cout << "liste_instants_.size_array()-1" << liste_instants_.size_array()-1 << finl;
        }
      else
        {
          Cout << "ELSE : time_ != 0 " << time_ << finl;
          /* Instant initial, ou instant de reprise */
          Cout << "v_cible_" << v_cible_ << finl;
          v_cible_ += 0.;
        }
      Cout << "duree_effective_morceau_glissant_" << duree_demandee_morceau_glissant_ << finl;
      Cout << "FIN : time_ <= duree_demandee_morceau_glissant_" << finl;
    }
  /* ******************************************
   * time_ > duree_morceau_glissant_ : la moyenne glisse
   *  si on a rempli N+1 instants :
   *  v_c^n = { v_c^(n-1) T^(n-1) -       vdt^(n-N)          +  vdt^n } / (t^n  -      t^(n-N)         )
   *  OU =>     v_cible_^(mem.)     liste_v_cible_dl_^(prem.)   calc.      mem.  liste_instants_^(prem.)
   * *****************************************/
  if (time_ > duree_demandee_morceau_glissant_)
    {
      const int len_list = liste_instants_.size_array();
      Cout << time_ <<" > " << duree_demandee_morceau_glissant_ << finl;
      Cout << "vitesse_relative" << vitesse_relative << finl;
      Cout << "timestep_" << timestep_ << finl;
      Cout << "liste_v_cible_dl_[list_index_]"<< liste_v_cible_dl_[list_index_] << finl;
      Cout << "liste_instants_[list_index_]" << liste_instants_[list_index_] << finl;
      list_index_ = (tstep_-1+offset_list_index_)%len_list;
      Cout << "list_index_ " << list_index_ << finl;
      Cout << "(list_index_+1)%len_list" << (list_index_+1)%len_list << finl;
      Cout << "liste_v_cible_dl_[list_index_]"<< liste_v_cible_dl_[list_index_] << finl;
      Cout << "liste_instants_[(list_index_+1)%len_list]" << liste_instants_[(list_index_+1)%len_list] << finl;
      Cout << "(vitesse_relative)*timestep_" << (vitesse_relative)*timestep_<< finl;
      /* With old time_interval_ */
      v_cible_ = v_cible_*duree_effective_morceau_glissant_ -liste_v_cible_dl_[(list_index_+1)%len_list] + (vitesse_relative)*timestep_;
      duree_effective_morceau_glissant_ = time_ - liste_instants_[(list_index_+1)%len_list];
      v_cible_ /= duree_effective_morceau_glissant_;

      /* Au lieu de decaler TOUT LE MONDE a gauche, un remplace la case d'indice j=(N_max+i)%N_max = i%N_max
      *  et apres on vient piocher la valeur tstep%N_max dans mes listes */
      liste_instants_[list_index_] = time_;
      liste_v_cible_dl_[list_index_] = (vitesse_relative)*timestep_;
    }
}

void moyenne_glissante::compute_qdm_cible(double alpha_liq)
{
  qdm_cible_ = alpha_liq*rho_liq_*v_cible_;
}

/////////////////////////////////////////////////////////////////////////////////////////
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

Implemente_instanciable_sans_constructeur( correction_one_direction, "correction_one_direction", Objet_U ) ;

correction_one_direction::correction_one_direction() :
  value_correction_(0),
  correct_velocity_(0),
  type_corr_(CONSIGNE_INITIALE), // plutot le mettre en consigne_initiale_ pour un usage par defaut
  i_need_vitesse_relative_(0),
  rho_vel_moyen_(0.),
  rho_moyen_(0.),
  alpha_l_(0.),
  vitesse_relative_(0.),
  qdm_cible_(0.)
{
  // Ces membres sont initialises par le constructeur par defaut:
  //  parametres_cible_constante_.initialise();
  // parametres_moyenne_par_morceaux_.initialise();
  // parametres_moyenne_glissante_.initialise();
}


Sortie& correction_one_direction::printOn( Sortie& os ) const
{
  //Objet_U::printOn( os );
  // TODO : Help, je ne sais pas faire l'operation inverse de get_type_corr_()
  os  << "{\n";
  if ( get_type_corr()==CIBLE_CONSTANTE)    os << " type_correction " << "cible_constante" << "\n";
  if ( get_type_corr()==MOYENNE_PAR_MORCEAUX)    os << " type_correction " << "moyenne_par_morceaux" << "\n";
  if ( get_type_corr()==MOYENNE_GLISSANTE)    os << " type_correction " << "moyenne_glissante" << "\n";
  if ( get_type_corr()==REGIME_ETABLI)    os << " type_correction " << "regime_etabli" << "\n";
  if ( get_type_corr()==CONSIGNE_INITIALE)    os << " type_correction " << "consigne_initiale" << "\n";
  os << " parametres_cible_constante " << parametres_cible_constante_ << "\n"
     << " parametres_moyenne_par_morceaux " << parametres_moyenne_par_morceaux_ << "\n"
     << " parametres_moyenne_glissante " << parametres_moyenne_glissante_ << "\n"
     << " parametres_consigne_initiale " << parametres_consigne_initiale_ << "\n";
  os << " }\n" ;
  return os;

}

Entree& correction_one_direction::readOn( Entree& is )
{
  //Objet_U::readOn( is );
  Param param(que_suis_je());
  param.ajouter("type_correction",&type_corr_);
  param.dictionnaire("cible_constante", CIBLE_CONSTANTE);
  param.dictionnaire("moyenne_par_morceaux", MOYENNE_PAR_MORCEAUX);
  param.dictionnaire("moyenne_glissante", MOYENNE_GLISSANTE);
  param.dictionnaire("regime_etabli", REGIME_ETABLI); // unused (29/06/2022)
  param.dictionnaire("consigne_initiale", CONSIGNE_INITIALE);
  param.ajouter("parametres_cible_constante",&parametres_cible_constante_);
  param.ajouter("parametres_moyenne_par_morceaux",&parametres_moyenne_par_morceaux_);
  param.ajouter("parametres_moyenne_glissante",&parametres_moyenne_glissante_);
  param.ajouter("parametres_consigne_initiale",&parametres_consigne_initiale_);
  param.lire_avec_accolades(is);
  return is;
}

void correction_one_direction::set_correction(correction_one_direction& correction_in)
{
  type_corr_ = correction_in.get_type_corr();
  parametres_cible_constante_ = correction_in.get_cible_constante();
  parametres_moyenne_par_morceaux_ = correction_in.get_moyenne_par_morceaux();
  parametres_moyenne_glissante_ = correction_in.get_moyenne_glissante();
  parametres_consigne_initiale_ = correction_in.get_consigne_initiale();
}

correction_one_direction::type_correction correction_one_direction::get_type_corr() const
{
  return (type_correction) type_corr_;
}

cible_donnee correction_one_direction::get_cible_constante() const
{return parametres_cible_constante_;}
moyenne_par_morceaux correction_one_direction::get_moyenne_par_morceaux() const
{return parametres_moyenne_par_morceaux_;}
moyenne_glissante correction_one_direction::get_moyenne_glissante() const
{return parametres_moyenne_glissante_;}
consigne_initiale correction_one_direction::get_consigne_initiale() const
{return parametres_consigne_initiale_;}

void correction_one_direction::set_time_for_correction(double time, double time_step, int time_iteration)
{
  // switch VOLONTAIREMENT NON UTILISE CAR DIMINUE LA LISIBILITE
  if ( get_type_corr()==CIBLE_CONSTANTE)
    {
      parametres_cible_constante_.set_time();
    }
  else if ( get_type_corr()==MOYENNE_PAR_MORCEAUX)
    {
      parametres_moyenne_par_morceaux_.set_time(time,time_step,time_iteration);
    }
  else if ( get_type_corr()==MOYENNE_GLISSANTE)
    {
      parametres_moyenne_glissante_.set_time(time,time_step,time_iteration);
    }
  else if ( get_type_corr()==CONSIGNE_INITIALE)
    {
      parametres_consigne_initiale_.set_time(time_iteration);
    }
  else
    Cerr << "Error, unknown type_correction_ "<< (int)get_type_corr() <<" for correction_perp_g" << finl;
}
void correction_one_direction::set_rho_moyen_alpha_l_for_correction(double rho_moyen, double alpha_l)
{
  rho_moyen_ = rho_moyen;
  alpha_l_ = alpha_l;
}

void correction_one_direction::set_rho_vel_moyen_for_correction(double rho_vel_moyen)
{
  rho_vel_moyen_ = rho_vel_moyen;
}

void correction_one_direction::set_mean_values_for_correction(double rho_vel_moyen, double rho_moyen, double alpha_l)
{
  rho_vel_moyen_ = rho_vel_moyen;
  rho_moyen_ = rho_moyen;
  alpha_l_ = alpha_l;
}

void correction_one_direction::compute_correction_value()
{

  // switch VOLONTAIREMENT NON UTILISE CAR DIMINUE LA LISIBILITE
  if ( get_type_corr()==CIBLE_CONSTANTE)
    {
      parametres_cible_constante_.compute_qdm_cible(alpha_l_);
      qdm_cible_ = parametres_cible_constante_.get_qdm_cible();
    }
  else if ( get_type_corr()==MOYENNE_PAR_MORCEAUX)
    {
      parametres_moyenne_par_morceaux_.compute_qdm_cible(alpha_l_, vitesse_relative_);
      qdm_cible_ = parametres_moyenne_par_morceaux_.get_qdm_cible();
    }
  else if ( get_type_corr()==MOYENNE_GLISSANTE)
    {
      parametres_moyenne_glissante_.compute_v_cible(vitesse_relative_);
      parametres_moyenne_glissante_.compute_qdm_cible(vitesse_relative_);
      qdm_cible_ = parametres_moyenne_glissante_.get_qdm_cible();
    }
  else if ( get_type_corr()==CONSIGNE_INITIALE)
    {
      parametres_consigne_initiale_.compute_qdm_cible(rho_vel_moyen_);
      qdm_cible_ = parametres_consigne_initiale_.get_qdm_cible();
    }
  else
    Cerr << "Error, unknown type_correction_ "<< (int)get_type_corr() <<" for correction_perp_g" << finl;

  value_correction_= (rho_vel_moyen_ - qdm_cible_) / rho_moyen_;
}

void correction_one_direction::compute_correct_velocity(double vel_ijk_t)
{
  correct_velocity_ = vel_ijk_t - value_correction_;
}

void correction_one_direction::set_rho_liq(double rho_l)
{
  // switch VOLONTAIREMENT NON UTILISE CAR DIMINUE LA LISIBILITE
  if ( get_type_corr()==CIBLE_CONSTANTE)
    {
      parametres_cible_constante_.set_rho_l(rho_l);
    }
  else if ( get_type_corr()==MOYENNE_PAR_MORCEAUX)
    {
      parametres_moyenne_par_morceaux_.set_rho_l(rho_l);
    }
  else if ( get_type_corr()==MOYENNE_GLISSANTE)
    {
      parametres_moyenne_glissante_.set_rho_l(rho_l);
    }
  else if ( get_type_corr()==CONSIGNE_INITIALE)
    {
      // no need for rho_l ?;
    }
  else
    Cerr << "Error, unknown type_correction_ "<< (int)get_type_corr() <<" for correction_perp_g" << finl;
}

int correction_one_direction::get_need_for_vit_rel()
{
  int need_to_compute_relative_velocity = 0;
  // switch VOLONTAIREMENT NON UTILISE CAR FONCTIONNE SUR DES INT UNIQUEMENT, CE QUI DIMINUE LA LISIBILITE
  if ( get_type_corr()==CIBLE_CONSTANTE)
    {
      need_to_compute_relative_velocity = parametres_cible_constante_.get_need_for_vl_vv();
    }
  else if ( get_type_corr()==MOYENNE_PAR_MORCEAUX)
    {
      need_to_compute_relative_velocity = parametres_moyenne_par_morceaux_.get_need_for_vl_vv();
    }
  else if ( get_type_corr()==MOYENNE_GLISSANTE)
    {
      need_to_compute_relative_velocity = parametres_moyenne_glissante_.get_need_for_vl_vv();
    }
  else if ( get_type_corr()==CONSIGNE_INITIALE)
    {
      need_to_compute_relative_velocity = parametres_consigne_initiale_.get_need_for_vl_vv();
    }
  else
    Cerr << "Error, unknown type_correction_ "<< (int)get_type_corr() <<" for correction_perp_g" << finl;
  return need_to_compute_relative_velocity;
}

void correction_one_direction::set_vl_vv(double vl_vv)
{
  vitesse_relative_ = vl_vv;
}

int correction_one_direction::get_need_for_rho_liq()
{
  int need_to_compute_relative_velocity = 0;
  // switch VOLONTAIREMENT NON UTILISE CAR FONCTIONNE SUR DES INT UNIQUEMENT, CE QUI DIMINUE LA LISIBILITE
  if ( get_type_corr()==CIBLE_CONSTANTE)
    {
      need_to_compute_relative_velocity = parametres_cible_constante_.get_need_for_rho_l();
    }
  else if ( get_type_corr()==MOYENNE_PAR_MORCEAUX)
    {
      need_to_compute_relative_velocity = parametres_moyenne_par_morceaux_.get_need_for_rho_l();
    }
  else if ( get_type_corr()==MOYENNE_GLISSANTE)
    {
      need_to_compute_relative_velocity = parametres_moyenne_glissante_.get_need_for_rho_l();
    }
  else if ( get_type_corr()==CONSIGNE_INITIALE)
    {
      need_to_compute_relative_velocity = parametres_consigne_initiale_.get_need_for_rho_l();
    }
  else
    Cerr << "Error, unknown type_correction_ "<< (int)get_type_corr() <<" for correction_perp_g" << finl;
  return need_to_compute_relative_velocity;
}

int correction_one_direction::get_need_to_compute_correction_value()
{
  // TODO : Si on  ne veux corriger que dans une direction de l'espace,
  // pas besoin de definir les correction des autres directions
  int need_to_compute_correction_value = 1;
//	  // switch VOLONTAIREMENT NON UTILISE CAR FONCTIONNE SUR DES INT UNIQUEMENT, CE QUI DIMINUE LA LISIBILITE
//	  if ( get_type_corr()==CONSIGNE_INITIALE)
//	    {
//		  need_to_compute_correction_value = parametres_consigne_initiale_.get_need_to_compute_correction_value();
//	    }
//
  return need_to_compute_correction_value;
}

/////////////////////////////////////////////////////////////////////////////////////////
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////

Implemente_instanciable_sans_constructeur( corrections_qdm, "corrections_qdm", Objet_U ) ;

corrections_qdm::corrections_qdm() :
  type_(2),    // type_ mis a none (-> NONE) par defaut
  write_me_(0) // par defaut on n'ecrit pas les informations
{
}

Sortie& corrections_qdm::printOn( Sortie& os ) const
{
  //Objet_U::printOn( os );
  // TODO : Help : je ne sais pas faire l'operation inverse de get_type
  os  << "{\n";
  if (get_type_()==GR) os << "type " << " gr" << "\n";
  if (get_type_()==GB) os << "type " << " gb" << "\n";
  if (get_type_()==NONE_IJK) os << "type " << " none" << "\n";
  os   << " write_info " << write_me_ << "\n";
  os   << " correction_x " << correction_x_ << "\n"
       << " correction_y " << correction_y_ << "\n"
       << " correction_z " << correction_z_ << "\n";
  os << " }\n" ;
  return os;
  return os;
}

Entree& corrections_qdm::readOn( Entree& is )
{
  //Objet_U::readOn( is );
  Param param(que_suis_je());
  param.ajouter("type",&type_);
  param.dictionnaire("gb",GB);
  param.dictionnaire("gr",GR);
  param.dictionnaire("none",NONE_IJK);
  param.ajouter_flag("write_infos",&write_me_);
  param.ajouter("correction_x",&correction_x_ );
  param.ajouter("correction_y",&correction_y_ );
  param.ajouter("correction_z",&correction_z_ );
  param.lire_avec_accolades(is);
  return is;
}

corrections_qdm::type_dict_ corrections_qdm::get_type_() const
{
  return (type_dict_) type_;
}

void corrections_qdm::set_corrections(correction_one_direction& in_corr_x_,
                                      correction_one_direction& in_corr_y_,
                                      correction_one_direction& in_corr_z_)
{
  correction_x_.set_correction(in_corr_x_);
  correction_y_.set_correction(in_corr_y_);
  correction_z_.set_correction(in_corr_z_);
}

void corrections_qdm::set_time(double time_, double time_step, int time_iteration)
{
  correction_x_.set_time_for_correction(time_, time_step, time_iteration);
  correction_y_.set_time_for_correction(time_, time_step, time_iteration);
  correction_z_.set_time_for_correction(time_, time_step, time_iteration);
}

void corrections_qdm::set_rho_moyen_alpha_l(double rho_moyen, double alpha_l)
{
  correction_x_.set_rho_moyen_alpha_l_for_correction(rho_moyen, alpha_l);
  correction_y_.set_rho_moyen_alpha_l_for_correction(rho_moyen, alpha_l);
  correction_z_.set_rho_moyen_alpha_l_for_correction(rho_moyen, alpha_l);
}

void corrections_qdm::set_rho_vel_moyen(int dir, double rho_vel_moyen)
{
  switch(dir)
    {
    case 0:
      correction_x_.set_rho_vel_moyen_for_correction(rho_vel_moyen);
      break;
    case 1:
      correction_y_.set_rho_vel_moyen_for_correction(rho_vel_moyen);
      break;
    case 2:
      correction_z_.set_rho_vel_moyen_for_correction(rho_vel_moyen);
      break;
    default:
      Cerr << "In corrections_qdm::compute_correction_one_direction, dir " << dir << " Error." << finl;
    }
}

void corrections_qdm::set_mean_values_for_corrections(double rho_vel_moyen, double rho_moyen, double alpha_l)
{
  correction_x_.set_mean_values_for_correction(rho_vel_moyen, rho_moyen, alpha_l);
  correction_y_.set_mean_values_for_correction(rho_vel_moyen, rho_moyen, alpha_l);
  correction_z_.set_mean_values_for_correction(rho_vel_moyen, rho_moyen, alpha_l);
}

void corrections_qdm::compute_correction_value_one_direction(int dir)
{
  switch(dir)
    {
    case 0:
      correction_x_.compute_correction_value();
      break;
    case 1:
      correction_y_.compute_correction_value();
      break;
    case 2:
      correction_z_.compute_correction_value();
      break;
    default:
      Cerr << "In corrections_qdm::compute_correction_one_direction, dir " << dir << " Error." << finl;
    }
}

void corrections_qdm::compute_correct_velocity_one_direction(int dir, double vel_ijk_t)
{
  switch(dir)
    {
    case 0:
      correction_x_.compute_correct_velocity(vel_ijk_t);
      break;
    case 1:
      correction_y_.compute_correct_velocity(vel_ijk_t);
      break;
    case 2:
      correction_z_.compute_correct_velocity(vel_ijk_t);
      break;
    default:
      Cerr << "In corrections_qdm::compute_correction_one_direction, dir " << dir << " Error." << finl;
    }
}

Vecteur3 corrections_qdm::get_correct_velocities()
{
  /* Donne acces a correct_velocity_,
   * ->  correct_velocity_ = v_ijk_t - value_correction */
  Vecteur3 correct_velocities;
  correct_velocities[0] = correction_x_.get_correct_velocity();
  correct_velocities[1] = correction_y_.get_correct_velocity();
  correct_velocities[2] = correction_z_.get_correct_velocity();
  return correct_velocities;
}

Vecteur3 corrections_qdm::get_velocity_corrections()
{
  /* Donne acces a value_correction_,
   * -> value_correction_ =  (rho_vel_moyen_ - qdm_cible_) / rho_moyen_*/
  Vecteur3 velocity_corrections;
  velocity_corrections[0] = correction_x_.get_velocity_correction();
  velocity_corrections[1] = correction_y_.get_velocity_correction();
  velocity_corrections[2] = correction_z_.get_velocity_correction();
  return velocity_corrections;
}

Vecteur3 corrections_qdm::get_correction_values()
{
  /* Donne acces a qdm_cible_,
   * dont l'expression depend du type de correction choisi
   * -> pour VITESSE_CIBLE :  qdm_cible_ = alpha_l rho_l u_cible */
  Vecteur3 correction_values;
  correction_values[0] = correction_x_.get_correction_value();
  correction_values[1] = correction_y_.get_correction_value();
  correction_values[2] = correction_z_.get_correction_value();
  return correction_values;
}


double corrections_qdm::get_correct_velocitiy_one_direction(int dir)
{
  switch(dir)
    {
    case 0:
      return correction_x_.get_correct_velocity();
      break;
    case 1:
      return correction_y_.get_correct_velocity();
      break;
    case 2:
      return correction_z_.get_correct_velocity();
      break;
    default:
      return 0.;
      Cerr << "In corrections_qdm::compute_correction_one_direction, dir " << dir << " Error." << finl;
    }
}

int corrections_qdm::get_need_for_vitesse_relative(int direction)
{
  switch(direction)
    {
    case 0:
      return correction_x_.get_need_for_vit_rel();
      break;
    case 1:
      return correction_y_.get_need_for_vit_rel();
      break;
    case 2:
      return correction_z_.get_need_for_vit_rel();
      break;
    default:
      return 0.;
      Cerr << "In corrections_qdm::compute_correction_one_direction, dir " << direction << " Error." << finl;
    };
}

int corrections_qdm::get_need_to_compute_correction_value_one_direction(int direction)
{
  switch(direction)
    {
    case 0:
      return correction_x_.get_need_to_compute_correction_value();
      break;
    case 1:
      return correction_y_.get_need_to_compute_correction_value();
      break;
    case 2:
      return correction_z_.get_need_to_compute_correction_value();
      break;
    default:
      return 0.;
      Cerr << "In corrections_qdm::compute_correction_one_direction, dir " << direction << " Error." << finl;
    };
}

void corrections_qdm::set_vitesse_relative(int direction, double vl_vv)
{
  if (correction_x_.get_need_for_vit_rel() && direction==0)
    {
      correction_x_.set_vl_vv(vl_vv);
    }
  if (correction_y_.get_need_for_vit_rel() && direction==1)
    {
      correction_y_.set_vl_vv(vl_vv);
    }
  if (correction_z_.get_need_for_vit_rel() && direction==2)
    {
      correction_z_.set_vl_vv(vl_vv);
    }
}

void corrections_qdm::set_rho_liquide(double rho_liquide)
{
  if (correction_x_.get_need_for_rho_liq())
    {
      correction_x_.set_rho_liq(rho_liquide);
    }
  if (correction_y_.get_need_for_rho_liq())
    {
      correction_y_.set_rho_liq(rho_liquide);
    }
  if (correction_z_.get_need_for_rho_liq())
    {
      correction_z_.set_rho_liq(rho_liquide);
    }
}

int corrections_qdm::is_type_gb() const
{
  if (type_ == GB)
    return 1;
  else
    return 0;
}

int corrections_qdm::is_type_gr() const
{
  if (type_ == GR)
    return 1;
  else
    return 0;
}

int corrections_qdm::is_type_none() const
{
  if (type_ == NONE_IJK)
    return 1;
  else
    return 0;
}

int corrections_qdm::write_me() const
{return write_me_;}
/////////////////////////////////////////////////////////////////////////////////////////
// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////
