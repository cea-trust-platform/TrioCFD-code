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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Init_spectral.cpp
// Directory : $IJK_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <init_forcage_THI.h>
#include <Param.h>
#include <string>
#include <iostream>
#include <math.h>
#include <fftw3.h>

#include <Force_ph.h>
#include <Force_sp.h>
#include <Random_process.h>

Implemente_instanciable_sans_constructeur( init_forcage_THI, "init_forcage_THI", Objet_U ) ;

// TODO : Il faudrait aussi initialiser
//Force_sp f_sp_THI;
//Force_ph f_ph_THI;
//Random_process random;
init_forcage_THI::init_forcage_THI():
  type_forcage(0),
  facteur_forcage_(0),
  forced_advection_(0),
  time_to_be_del_(0),
  forcage_ts_stop(-1),
  forcage_t_stop(-1.),
  mode_min(0),
  mode_max(0),
  amplitude(1.0),
  eps_etoile(0.1),
  tL(0.02),
  i_offset(0),
  j_offset(0),
  k_offset(0)
{
  advection_length_.resize_array(3);
  advection_length_ = 0.;
}

Sortie& init_forcage_THI::printOn( Sortie& os ) const
{
  //Objet_U::printOn(os);
  os  << "{\n"
      << "   type " << type_forcage << "\n"
      << "   facteur " << facteur_forcage_ << "\n"
      << "   forced_advection " << forced_advection_ << "\n"
      << "   advection_velocity " << advection_velocity_ << "\n"
      << "   advection_length " << advection_length_ << "\n"
      << "   stops_at_time_step " << forcage_ts_stop << "\n"
      << "   stops_at_time " << forcage_t_stop << "\n";
  os  << "   minimal_forced_mode " << mode_min << "\n"
      << "   maximal_forced_mode " << mode_max << "\n";
  os  << "   amplitude " << amplitude << "\n";
  os  << "   dissipation " << eps_etoile << "\n";
  os  << "   temps_grande_echelle " << tL << "\n";
//  if (random_fixed_)
//    os << "   random_fixed \n" ;
  os  << "   random_process " << random_ << "\n";
  os << " }\n" ;
  return os;
}

Entree& init_forcage_THI::readOn( Entree& is )
{
  // Objet_U::readOn( is );
  Param param(que_suis_je());
  param.ajouter("type",&type_forcage);
  param.ajouter("facteur",&facteur_forcage_);
  param.ajouter("forced_advection",&forced_advection_);
  if (forced_advection_==1)
    {param.ajouter("advection_velocity", &advection_velocity_, Param::REQUIRED);}
  else
    {param.ajouter("advection_velocity", &advection_velocity_);}
  param.ajouter("advection_length", &advection_length_);
  param.ajouter("stops_at_time_step",&forcage_ts_stop);
  param.ajouter("stops_at_time",&forcage_t_stop);
  param.ajouter("minimal_forced_mode",&mode_min);
  param.ajouter("maximal_forced_mode",&mode_max);
  param.ajouter("amplitude",&amplitude);
  param.ajouter("dissipation",&eps_etoile);
  param.ajouter("temps_grande_echelle",&tL);
  if (type_forcage == 3)
    {param.ajouter("random_process",&random_, Param::REQUIRED);}
  else
    {param.ajouter("random_process",&random_);}
  param.lire_avec_accolades(is);
  return is;
}

void init_forcage_THI::compute_initial_chouippe(int my_nproc_tot,
                                                const IJK_Grid_Geometry& my_geom,
                                                int my_ni, int my_nj, int my_nk,
                                                const IJK_Splitting& splitting_,
                                                Nom nom_sauvegarde)
// FixedVector<IJK_Field_double, 3> v_)
{
  double Lx(my_geom.get_domain_length(DIRECTION_I));
  double Ly(my_geom.get_domain_length(DIRECTION_J));
  double Lz(my_geom.get_domain_length(DIRECTION_K));

  const double Ox = my_geom.get_origin(DIRECTION_I) ;
  const double Oy = my_geom.get_origin(DIRECTION_J) ;
  const double Oz = my_geom.get_origin(DIRECTION_K) ;

  // Test parallelisation
  i_offset = splitting_.get_offset_local(DIRECTION_I);
  j_offset = splitting_.get_offset_local(DIRECTION_J);
  k_offset = splitting_.get_offset_local(DIRECTION_K);

  int ni(my_ni), nj(my_nj), nk(my_nk);

  int mmin(mode_min), mmax(mode_max);
  // double kmin(2.*M_PI*mmin/Lx), kmax(2.*M_PI*mmax/Lx);
  double kmin(2.*M_PI*mmin/Lx), kmax(2.*M_PI*mmax/Lx);
  int nl(mmax), nm(mmax), nn(mmax);

  std::string nom_fichier_random("./random_gen.out");
  std::string nom_fichier_spectral("./spectral.out");
  std::string nom_fichier_physique("./physique.out");

  f_sp_THI.initialise(nl,nm,nn, mmin,mmax,kmin,kmax,amplitude, nom_fichier_spectral);
  // f_ph_THI.initialise(my_nproc_tot,ni,nj,nk,nl,nm,nn,Lx,Ly,Lz,Ox,Oy,Oz,mmin,mmax,kmin,kmax, nom_fichier_physique, splitting_);//,i_offset,j_offset,k_offset);
  f_ph_THI.initialise(my_nproc_tot,ni,nj,nk,nl,nm,nn,Lx,Ly,Lz,Ox,Oy,Oz,mmin,mmax,kmin,kmax, nom_fichier_physique, splitting_,i_offset,j_offset,k_offset);
  if (type_forcage == 3)
    {
      random_.initialise(eps_etoile,tL, nl,nm,nn, nom_fichier_random, nom_sauvegarde);//, random_fixed_);//,i_offset,j_offset,k_offset);
      // random.initialise(eps_etoile,tL, nl,nm,nn, nom_fichier_random, random_fixed_,i_offset,j_offset,k_offset);
    }
}


void init_forcage_THI::compute_THI_force(const int time_iteration,
                                         const double tstep,
                                         const double current_time,
                                         const IJK_Splitting& my_splitting
                                         // const int rk_step,
                                        )
{
  /*
   * En fonction de l'entree utilisateur pour type dans le jdd, cette methode va appeler
   * La fonction de forcage adequate.
   * -> type; facteur; forced_advection>advection_velocity; stops_at_time.step; ...
   *   o forced_advection_==-1 : on advecte avec \ol{u}^l
   *   o forced_advection_== 1 : on advecte avec la valeur lue dans le jdd
   * Remarque : dans le cas d'une fonction de forcage deterministe, appelee lors d'une resolution
   *            RK3, est ce qu'on pourrai forcer uniquement un des 3 sous pas de temps plutot que tous les forcer ?
   * Remarque : Pour un calcul avec reprise et pour forced_advection_==-1, il s'ssurer que la vitesse d'advection
   *            correspond bien a \ol{u}^l DES LE PREMIER PAS DE TEMPS.
  */
  int type_forcage_active(activate_forcage(time_iteration,current_time));

  if (type_forcage_active==0)
    {
      Cout << "Normalement on ne rentre pas ici. Que si on a force de la thi, puisqu'on arrete.";
      f_sp_THI.set_zero();
      f_ph_THI.set_zero();
    }
  if (type_forcage_active==100)
    {
      Cout << "On force un dirac en spectral, uniX" << finl;
      f_sp_THI.compute_dirac_point_uniX_alongX();
    }
  if (type_forcage_active==200)
    {
      Cout << "On force un dirac en spectral, uniX" << finl;
      f_sp_THI.compute_dirac_point_uniX_alongY();
    }
  if (type_forcage_active==010)
    {
      Cout << "On force un dirac en spectral, uniY" << finl;
      f_sp_THI.compute_dirac_point_uniY();
    }
  if (type_forcage_active==001 || type_forcage_active==1)
    {
      Cout << "On force un dirac en spectral, uniZ" << finl;
      f_sp_THI.compute_dirac_point_uniZ();
    }
  if (type_forcage_active==101)
    {
      Cout << "On force un dirac en spectral, uniXZ" << finl;
      f_sp_THI.compute_dirac_point_div_nulle();
    }
  if (type_forcage_active==2)
    {
      Cout << "On force une porte en spectral, soit un cube" << finl;
      f_sp_THI.compute_door_cube();
    }
  if (type_forcage_active==3)
    {
//      Cout << "type_forcage_active : " << type_forcage_active << finl;
//      Cout << "On applique le processus Ornstein-Uhlenbeck" << finl;

      random_.next_step2(tstep, time_iteration);
      f_sp_THI.compute_step2(random_.get_b_flt());
    }
  if (type_forcage_active==20)
    {
      f_ph_THI.cheat_function();
    }

  /* PASSAGE DU DOMAINE SPECTRAL AU DOMAINE PHYSIQUE */
  if (type_forcage_active!=20 && type_forcage_active!=0)
    {
      /*type 20 ecrit directement le champ physique sans passer par la TF
       * Remarque : Si on ajoute d'autre "type"s qui ne doivent pas passer par la TF il faut penser
       *            a changer le if () juste au dessus.
      */

      /* ADVECTION DU CHAMP FORCE : forced_advection 1 ou forced_advection -1 dans jdd
       * forced_advection  1 : advection selon advection_velocity (jdd)
       * forced_advection -1 : advection selon moy{u_l}^l         (calcul),
       *                       /!\ Dans le cas d'une reprise, doit on lire la valeur inscrite dans le .sauv ?
       * forced_advection  0 : advcetion nulle
       * */
      if (forced_advection_==-1 || forced_advection_==1)
        {
          f_ph_THI.from_spect_to_phys_opti2_advection(f_sp_THI.get_coeff_flt(),advection_length_);
        }
      else
        f_ph_THI.from_spect_to_phys_opti2(f_sp_THI.get_coeff_flt());
    }

}

int init_forcage_THI::get_semi_gen()
{
  return random_.get_semi_gen();
}

ArrOfDouble init_forcage_THI::get_b_flt()
{
  return random_.get_b_flt();

}

FixedVector<IJK_Field_double, 3> init_forcage_THI::get_force_ph()
{
  return f_ph_THI.get_force_attribute();
}


FixedVector<IJK_Field_double, 3>& init_forcage_THI::get_force_ph2()
{
  return f_ph_THI.get_force_attribute2();
}


int init_forcage_THI::get_type_forcage()
{
  return type_forcage;
}

int init_forcage_THI::get_facteur_forcage()
{
  return facteur_forcage_;
}

void init_forcage_THI::update_advection_velocity(ArrOfDouble& value)
{
  advection_velocity_.resize_array(3);
  for (int dir=0; dir<3; dir++)
    {
      advection_velocity_[dir]= value[dir];
    }
}

void init_forcage_THI::update_advection_length(double dt)
{
  Cout <<"update_advection_length : dt " << dt << finl;
  time_to_be_del_+=dt;
  Cout << "update_advection_length : time_to_be_del_ " << time_to_be_del_ << finl;
  for (int dir=0; dir<3; dir++)
    {
      Cout << "update_advection_length : advection_length_["<<dir<< "] : " << advection_length_[dir] << finl;
      advection_length_[dir]+= (advection_velocity_[dir]*dt);
      Cout << "update_advection_length : advection_length_["<<dir<< "] : " << advection_length_[dir] << finl;
    }
}

int init_forcage_THI::get_forced_advection()
{
  return forced_advection_;
}

int init_forcage_THI::activate_forcage(const int current_time_step, const double current_time)
{
  int stop(0);
  int no_stop(get_type_forcage());

  // J'ai l'impression qu'on peut faire plus malin. Mais bon je ne vois pas trop
  if ((forcage_ts_stop < 0 ) && (forcage_t_stop < 0))
    return no_stop;
  else if ((forcage_ts_stop < 0) && (forcage_t_stop > current_time))
    {
      return no_stop;
    }
  else if ((forcage_ts_stop > current_time_step) && (forcage_t_stop < 0))
    {
      return no_stop;
    }
  else if ((forcage_ts_stop > current_time_step) && (forcage_t_stop > current_time))
    {
      return no_stop;
    }
  else
    {
      return stop;
    }
}
