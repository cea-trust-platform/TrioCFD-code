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
// File      : Tube_base.cpp
// Directory : $IJK_ROOT/src/IBC
//
/////////////////////////////////////////////////////////////////////////////
#include <Parser.h>
#include <Tube_base.h>
#include <Param.h>
#include <IJK_Field.h>

Implemente_base(Tube_base, "Tube_base", Objet_U);
Implemente_instanciable(Faisceau_Tubes, "Faisceau_Tubes", VECT(DERIV(Tube_base)));

static void evaluer_f_et_df_dt_et_df2_dt2(const Noms& expressions, const double t, const double dt, Vecteur3& f, Vecteur3& df_dt, Vecteur3& df2_dt2)
{
  for (int dir = 0; dir < 3; dir++)
    {
      // Pour chaque direction, evaluer l'expression mathematique
      Parser parser;
      parser.setString(std::string(expressions[dir]));
      parser.setNbVar(1); // il n'y a qu'une variable dans l'expression : le temps t; si on avait eu deux variables ca aurait ete :  parser.setNbVar(2);
      parser.addVar("t"); // la premiere variable, d'indice zero sera le temps
      parser.parseString(); // Lance la conversion de la chaine en une succession d'operations mathematiques

      // Calcul de la position au temps courant:
      parser.setVar((int)0, t); // variable numero zero = "t"
      f[dir] = parser.eval(); // permet d'evaluer l'expression donc la position

      // Calcul de la vitesse au temps courant:
      double delta_t = dt; // on fait une difference finie avec ce pas de temps, on peut aussi prendre dt/10 par ex...
      parser.setVar((int)0, t + delta_t); // l'expression est prise au temps t + delta_t
      double f2 = parser.eval(); // on evalue l'expression au temps t + delat_t
      df_dt[dir] = (f2 - f[dir]) / delta_t;

      // Calcul de l'acceleration au temps courant par difference finie centree:
      parser.setVar((int)0, t + delta_t); // l'expression est prise au temps t + delta_t
      double x_p1 = parser.eval(); // position au temps n+1
      parser.setVar((int)0, t - delta_t); // l'expression est prise au temps t - delta_t
      double x_m1 = parser.eval(); // position au temps n-1
      df2_dt2[dir] = (x_p1 - 2 * f[dir] + x_m1) / ( delta_t * delta_t); // a = (xn+1 -2*xn + xn-1) / (dt*dt)
    }
}

void Tube_base::initialize(const IJK_Splitting& x)
{
  ref_splitting_ = x;
}

Sortie& Tube_base::printOn(Sortie& os) const
{
  return os;
}

Entree& Tube_base::readOn(Entree& is)
{
  return is;
}

Implemente_instanciable(Tube_impose, "Tube_impose", Tube_base);

Sortie& Tube_impose::printOn(Sortie& os) const
{
  // On ecrit ce qu'il faut pour relire l'objet lors d'une reprise:
  os << " { ";
  if (nom_ != "??")
    os << " nom " << nom_;
  if (omega_ != 0)
    os << " omega " << omega_;
  os << " pos_x " << expressions_position_[0] ;
  os  << " pos_y " << expressions_position_[1] ;
  os << " pos_z " << expressions_position_[2] ;
  os << " rayon " << rayon_ ;
  os  << " }";
  // On ecrit position, vitesse et acceleration actuelle en commentaires
  os << " # pos=[ " << pos_[0] << " " << pos_[1] << " " << pos_[2]  ;
  os  << " ] v=[ " << v_[0] << " " << v_[1] << " " << v_[2] ;
  os  << " ] a=[ " << acceleration_[0] << " " << acceleration_[1] << " " << acceleration_[2] << " ] # " << finl;
  return os;
}

Entree& Tube_impose::readOn(Entree& is)
{
  omega_ = 0.;
  expressions_position_.dimensionner(3);
  Param param(que_suis_je());
  param.ajouter("pos_x", &expressions_position_[0], Param::REQUIRED);
  param.ajouter("pos_y", &expressions_position_[1], Param::REQUIRED);
  param.ajouter("pos_z", &expressions_position_[2], Param::REQUIRED);
  param.ajouter("rayon", &rayon_, Param::REQUIRED);
  param.ajouter("nom",   &nom_);
  param.ajouter("omega", &omega_);
  param.lire_avec_accolades(is);

  // Initialisation a t=0
  update_vitesse_position(0., 1e-6, Vecteur3(0,0,0));

  return is;
}

void Tube_impose::update_vitesse_position(double current_time, double dt, const Vecteur3& force_appliquee)
{
  evaluer_f_et_df_dt_et_df2_dt2(expressions_position_, current_time, dt, pos_, v_, acceleration_);
}


Implemente_instanciable(Tube_libre, "Tube_libre", Tube_base);

Sortie& Tube_libre::printOn(Sortie& os) const
{
  // On ecrit ce qu'il faut pour relire l'objet lors d'une reprise:
  os << " { ";
  if (nom_ != "??")
    os << " nom " << nom_;
  if (omega_ != 0)
    os << " omega " << omega_;
  os << " pos_x " << pos_[0] ;
  os   << " pos_y " << pos_[1] ;
  os  << " pos_z " << pos_[2] ;
  os  << " v_x " << v_[0] ; ;
  os  << " v_y " << v_[1] ;
  os     << " v_z " << v_[2] ;
  os    << " a_x " << acceleration_[0] ;
  os   << " a_y " << acceleration_[1] ;
  os  << " a_z " << acceleration_[2] ;
  os << " pos_eq_x " << position_equilibre_[0] ;
  os  << " pos_eq_y " << position_equilibre_[1] ;
  os << " pos_eq_z " << position_equilibre_[2]; ;
  if (blocage_[DIRECTION_I]) os << " blocage_i";
  if (blocage_[DIRECTION_J]) os << " blocage_j";
  if (blocage_[DIRECTION_K]) os << " blocage_k";
  os << " rho_s " << rho_cylindre_;
  os  << " c " << amortissement_;
  os  << " k " << raideur_;
  os << " rayon " << rayon_;
  os << " }";
  return os;
}

Entree& Tube_libre::readOn(Entree& is)
{
  omega_ = 0.;
  pos_ = Vecteur3(1e38, 1e38, 1e38);
  // Par defaut, vitesse initiale nulle
  v_ = Vecteur3(0,0,0);
  // Par defaut, acceleration initiale nulle
  acceleration_ = Vecteur3(0,0,0);

  Param param(que_suis_je());
  // Position initiale
  param.ajouter("pos_x", &pos_[0]);
  param.ajouter("pos_y", &pos_[1]);
  param.ajouter("pos_z", &pos_[2]);
  // Vitesse initiale
  param.ajouter("v_x", &v_[0]);
  param.ajouter("v_y", &v_[1]);
  param.ajouter("v_z", &v_[2]);
  // Acceleration initiale
  param.ajouter("a_x", &acceleration_[0]);
  param.ajouter("a_y", &acceleration_[1]);
  param.ajouter("a_z", &acceleration_[2]);

  // Position d'equilibre
  param.ajouter("pos_eq_x", &position_equilibre_[0], Param::REQUIRED);
  param.ajouter("pos_eq_y", &position_equilibre_[1], Param::REQUIRED);
  param.ajouter("pos_eq_z", &position_equilibre_[2], Param::REQUIRED);
  // Degres de liberte du tube
  blocage_.resize_array(3);
  blocage_ = 0; // par defaut, non bloque.
  param.ajouter_flag("blocage_i", &blocage_[DIRECTION_I]);
  param.ajouter_flag("blocage_j", &blocage_[DIRECTION_J]);
  param.ajouter_flag("blocage_k", &blocage_[DIRECTION_K]);
  // Parametres mecaniques
  param.ajouter("rho_s", &rho_cylindre_, Param::REQUIRED);
  param.ajouter("c", &amortissement_, Param::REQUIRED);
  param.ajouter("k", &raideur_, Param::REQUIRED);

  param.ajouter("rayon", &rayon_, Param::REQUIRED);
  param.ajouter("nom",   &nom_);

  param.lire_avec_accolades(is);

  // Si on n'a pas donne de position initiale, prendre la position d'equilibre
  for (int i = 0; i < 3; i++)
    if (pos_[i] > 1e37)
      pos_[i] = position_equilibre_[i];

  return is;
}

void Tube_libre::update_vitesse_position(double current_time, double dt, const Vecteur3& force_appliquee)
{
  // constante du Newmark
  const double beta = 0.25;
  const double gamma = 0.5;
  const double timestep = dt;
  ArrOfDouble a(8);
  a[0] = 1/( beta * timestep * timestep);
  a[1] = gamma /( beta * timestep);
  a[2] =  1 /( beta * timestep);
  a[3] = 1/(2 *beta)-1;
  a[4] = gamma/ beta -1;
  a[5] = 0.5* timestep * (gamma/ beta -2);
  a[6] = timestep * (1- gamma);
  a[7] = gamma * timestep;
  Cout << "Coefficients newmark: " << a;
  // Calcul de la masse du tube:
  const double surface = M_PI * rayon_ * rayon_; // surface analytique
  // Calcul de la hauteur du tube
  const IJK_Splitting& splitting = ref_splitting_.valeur();
  const IJK_Grid_Geometry& geom = splitting.get_grid_geometry();
  ArrOfDouble noeuds_y;
  noeuds_y = geom.get_node_coordinates(1); // rentre les coordonnees des noeuds selon y dans le tableau noeuds_y
  const int n_element_y = noeuds_y.size_array()-1; // nombre d'elements selon y = nombre de noeuds -1
  const double delta_y = geom.get_constant_delta(1);
  const double hauteur_cylindre = n_element_y * delta_y;
  // Calcul du volume du tube
  double volume_cylindre = surface * hauteur_cylindre; // Vol analytique : vol = pi*rayon*rayon*h
  Cout << "Volume du tube: " << volume_cylindre << finl;;
  // Faire le calcul du volume avec la surface discretisee
  //Calcul de la masse du cylindre
  double m_s = rho_cylindre_ * volume_cylindre;
  double c= amortissement_;
  double k= raideur_;
  Vecteur3 pos_newmark;

  Vecteur3 F_droite, K_droite;
  // ATTENTION les matrices de M, C et K sont des matrices identites *m, c et k!!!
  // Donc on ne resout pas une equation matricielle
  double K_gauche = a[0] * m_s + a[1] * c + k; // cas particulier!! Dans le cas general K_gauche est une matrice

  for (int dir = 0; dir < 3; dir++)
    {
      if (blocage_[dir])
        {
          v_[dir] = acceleration_[dir] = 0.;
          acceleration_[dir] = 0.;
          continue;
        }
      // pos_newmark = pos_tube - pos_eq; donc v_newmark=v_ et a_newmark=acceleration_
      pos_newmark[dir] =  pos_[dir] - position_equilibre_[dir];
      // ce qu'on resout : K_gauche * position(n+1) = Force_appliquee - K_droite = D_droite
      // K_gauche est ici la matrice identite * (a[0] * m_s + a[1] * c + k), mais on pourrait mettre une raideur, masse, amortissement different selon la direction
      // K_droite est un vecteur car c'est : M * (a[0] * z_newmark[dir] + a[2] * v_newmark[dir] + a[3] * a_newmark[dir]) + C *(a[1] * z_newmark[dir] + a[4] * v_newmark[dir] + a[5] * a_newmark[dir])
      //ici comme M=Id * m_s et C=Id*c alors ca donne ce cas particulier :
      K_droite[dir] = m_s * (a[0] * pos_newmark[dir] + a[2] * v_[dir] + a[3] * acceleration_[dir]) + c * (a[1] * pos_newmark[dir] + a[4] * v_[dir] + a[5] * acceleration_[dir]);
      F_droite[dir] = force_appliquee[dir] + K_droite[dir];
      double tmp_pos = pos_newmark[dir]; // copie de pos_newmark au temps n
      // Calcul de la position au temps n+1 //
      pos_newmark[dir] = F_droite[dir] / K_gauche; // pos_newmark au temps n+1
      pos_[dir] = pos_newmark[dir] + position_equilibre_[dir]; // pos du tube au temps n+1
      // Calcul de l'acceleration au temps n+1 //
      double tmp_acceleration = acceleration_[dir]; // copie de l'acceleration au temps n
      acceleration_[dir] = a[0] * (pos_newmark[dir] - tmp_pos) - a[2] * v_[dir] - a[3] * acceleration_[dir]; // a(n+1) = a_0 * (pos_nwemark(n+1) - pos_newmark(n)) - a_2 v(n) - a_3 a(n)
      // Calcul de la vitesse au temps n+1 //
      v_[dir] += a[6] * tmp_acceleration + a[7] * acceleration_[dir]; // vz(n+1) = vz(n) + a_6 * az(n) + a_7 * az (n+1)

    }

}


// Les formats d'entree et de sortie sont ceux de VECT(xxx)
// c'est a dire
//  N  ( nombre d'objets dans le vecteur)
//  objet[0]
//  objet[1]
//  ...
Entree& Faisceau_Tubes::readOn(Entree& is)
{
  return VECT(DERIV(Tube_base))::readOn(is);
}

Sortie&   Faisceau_Tubes::printOn(Sortie& is) const
{
  return VECT(DERIV(Tube_base))::printOn(is);
}
