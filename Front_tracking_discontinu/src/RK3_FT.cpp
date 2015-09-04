/****************************************************************************
* Copyright (c) 2015, CEA
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
// File:        RK3_FT.cpp
// Directory:   $TRUST_ROOT/src/Front_tracking_discontinu
// Version:     /main/19
//
//////////////////////////////////////////////////////////////////////////////

#include <RK3_FT.h>
#include <Equation.h>
#include <Probleme_base.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <Ref_Champ_Inc_base.h>
#include <DoubleTabs.h>

Implemente_instanciable(RK3_FT,"RK3_FT",RK3);


// Description:
//    Simple appel a: Schema_Temps_base::printOn(Sortie& )
//    Ecrit le schema en temps sur un flot de sortie.
// Precondition:
// Parametre: Sortie& s
//    Signification: un flot de sortie
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Sortie&
//    Signification: le flot de sortie modifie
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
Sortie& RK3_FT::printOn(Sortie& s) const
{
  return  Schema_Temps_base::printOn(s);
}


// Description:
//    Lit le schema en temps a partir d'un flot d'entree.
//    Simple appel a: Schema_Temps_base::readOn(Entree& )
// Precondition:
// Parametre: Entree& s
//    Signification: un flot d'entree
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Entree&
//    Signification: le flot d'entree modifie
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
Entree& RK3_FT::readOn(Entree& s)
{
  return Schema_Temps_base::readOn(s) ;
}

////////////////////////////////
//                            //
// Caracteristiques du schema //
//                            //
////////////////////////////////


// Description:
//    Renvoie le nombre de valeurs temporelles a conserver.
//    Ici : n et n+1, donc 2.
int RK3_FT::nb_valeurs_temporelles() const
{
  return 2 ;
}

// Description:
//    Renvoie le nombre de valeurs temporelles futures.
//    Ici : n+1, donc 1.
int RK3_FT::nb_valeurs_futures() const
{
  return 1 ;
}

// Description:
//    Renvoie le le temps a la i-eme valeur future.
//    Ici : t(n+1)
double RK3_FT::temps_futur(int i) const
{
  assert(i==1);
  return temps_courant()+pas_de_temps();
}

// Description:
//    Renvoie le temps que doivent rendre les champs a
//    l'appel de valeurs()
//    Ici : t(n+1)
double RK3_FT::temps_defaut() const
{
  return temps_courant()+pas_de_temps();
}

/////////////////////////////////////////
//                                     //
// Fin des caracteristiques du schema  //
//                                     //
/////////////////////////////////////////

int RK3_FT::faire_un_pas_de_temps_eqn_base(Equation_base& eqn)
{
  Cerr << "on ne devrait pas passer ds faire_un_pas_de_temps_eqn_base de RK3_FT" << finl;
  return 0;
}

// Description:
//    Effectue un pas de temps de Runge Kutta d'ordre 3,
//    sur l'equation passee en parametre.
//    Le schema de Runge Kutta  d'ordre 3
//     (cas 7 de Williamson) s'ecrit :
//     q1=h f(x0)
//     x1=x0+b1 q1
//     q2=h f(x1) +a2 q1
//     x2=x1+b2 q2
//     q3=h f(x2)+a3 q2
//     x3=x2+b3 q3
//      avec a1=0, a2=-5/9, a3=-153/128
//                              b1=1/3, b2=15/16, b3=8/15
// Precondition:
// Parametre: Equation_base& eqn
//    Signification: l'equation que l'on veut faire avancer d'un
//                   pas de temps
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: int
//    Signification: renvoie toujours 1
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:

bool RK3_FT::iterateTimeStep(bool& converged)
{

  Probleme_base& prob=pb_base();

  int i;
  const double save_temps_courant=temps_courant_;
  int nb_eqn=prob.nombre_d_equations();
  for(i=0; i<nb_eqn; i++)
    if (prob.equation(i).equation_non_resolue())
      {
        Cerr << "Option equation_non_resolue is not yet supported for the " << que_suis_je() << " scheme." << finl;
        exit();
      }
  double dt_temp;

  double a2=-5./9.;
  double a3=-153./128.;
  double b1=1./3.;
  double b2=15./16.;
  double b3=8./15.;

  const double epsilon_dt = dt_ * 0.001;

  Champ_Inc_base& inconnue_ = prob.equation(0).inconnue().valeur();

  DoubleTab  qNSi(inconnue_.valeurs());
  DoubleTab  qNSj(inconnue_.valeurs());
  DoubleTabFT  qIi;
  DoubleTabFT  qIj;

  //  ufinal =VI.valeurs();
  // exemple
  //  equ_int.calculer_vitesse_transport_interpolee(
  //   champ_eulerien,
  //   equ_int.maillage_interface(),
  //   resu_vitesse_lagrange,
  //   1 /* reprojeter la vitesse L2 a partir de N.S. */ );
  //
  //  equ_int.transporter_sans_changement_topologie(vitesse, coeff);

  //ss pas de temps 1
  //mise a jour des conditions aux limites
  ////////////////////////////////////////////////////////////////////////////////

  for(i=0; i<nb_eqn; i++)
    {
      // Astuce pour calculer la derivee gpoint de la vitesse imposee au bord:
      //  gpoint est calcule par difference finie entre les deux derniers
      //  "mettre_a_jour". On fait une difference finie avec un delta_t tout petit:
      prob.equation(i).zone_Cl_dis().mettre_a_jour(temps_courant_ - epsilon_dt);
      prob.equation(i).zone_Cl_dis().mettre_a_jour(temps_courant_);
      prob.equation(i).zone_Cl_dis()->Gpoint(temps_courant_ - epsilon_dt,temps_courant_);
    }

  temps_courant_ += dt_*1./3.;

  dt_temp=dt_;
  dt_=0.;
  prob.equation(0).initTimeStep(dt_);
  dt_=dt_temp;

  prob.equation(0).derivee_en_temps_inco(qNSi);
  inconnue_.futur() = inconnue_.valeurs();
  // B.M.: calcul identique a l'implementation precedente:
  inconnue_.futur().ajoute_sans_ech_esp_virt(b1 * dt_, qNSi, VECT_ALL_ITEMS);

  if (prob.que_suis_je() == "Probleme_FT_Disc_gen")
    {
      Transport_Interfaces_FT_Disc& equ_int= ref_cast(Transport_Interfaces_FT_Disc,prob.equation(1));
      const Transport_Interfaces_FT_Disc& equ_int_const = equ_int;
      equ_int.calculer_vitesse_transport_interpolee(inconnue_,
                                                    equ_int_const.maillage_interface(),
                                                    qIi,
                                                    1);
      equ_int.transporter_sans_changement_topologie(qIi,b1*dt_,temps_courant_);
    }


  //ss pas de temps 2

  prob.equation(0).mettre_a_jour(temps_courant_);

  //mise a jour des conditions aux limites
  ////////////////////////////////////////////////////////////////////////////////////////////

  for(i=0; i<nb_eqn; i++)
    {
      prob.equation(i).zone_Cl_dis().mettre_a_jour(temps_courant_ - epsilon_dt);
      prob.equation(i).zone_Cl_dis().mettre_a_jour(temps_courant_);
      prob.equation(i).zone_Cl_dis()->Gpoint(temps_courant_ - epsilon_dt,temps_courant_);
    }

  temps_courant_ += dt_*5./12.;

  dt_temp=dt_;
  dt_=0.;
  prob.equation(0).initTimeStep(dt_);
  dt_=dt_temp;


  prob.equation(0).derivee_en_temps_inco(qNSj);
  qNSj.ajoute_sans_ech_esp_virt(a2, qNSi, VECT_ALL_ITEMS);

  inconnue_.futur() = inconnue_.valeurs();
  inconnue_.futur().ajoute_sans_ech_esp_virt(b2 * dt_, qNSj, VECT_ALL_ITEMS);

  if (prob.que_suis_je() == "Probleme_FT_Disc_gen")
    {
      Transport_Interfaces_FT_Disc& equ_int= ref_cast(Transport_Interfaces_FT_Disc,prob.equation(1));
      const Transport_Interfaces_FT_Disc& equ_int_const = equ_int;
      equ_int.calculer_vitesse_transport_interpolee(inconnue_,
                                                    equ_int_const.maillage_interface(),
                                                    qIj,
                                                    1);
      qIi *= -1.*a2;
      qIj -= qIi;
      equ_int.transporter_sans_changement_topologie(qIj,b2*dt_,temps_courant_);
    }


  //ss pas de temps 3

  prob.equation(0).mettre_a_jour(temps_courant_);

  //mise a jour des conditions aux limites
  //////////////////////////////////////////////////////////////////////////////////////////////

  for(i=0; i<nb_eqn; i++)
    {
      prob.equation(i).zone_Cl_dis().mettre_a_jour(temps_courant_ - epsilon_dt);
      prob.equation(i).zone_Cl_dis().mettre_a_jour(temps_courant_);
      prob.equation(i).zone_Cl_dis()->Gpoint(temps_courant_ - epsilon_dt,temps_courant_);
    }

  temps_courant_ += dt_*1./4.;

  dt_temp=dt_;
  dt_=0.;
  prob.equation(0).initTimeStep(dt_);
  dt_=dt_temp;

  prob.equation(0).derivee_en_temps_inco(qNSi);
  qNSi.ajoute_sans_ech_esp_virt(a3, qNSj);

  inconnue_.futur() = qNSi;
  inconnue_.futur() *= b3 * dt_;
  update_critere_statio(inconnue_.futur(), prob.equation(0));
  inconnue_.futur() += inconnue_.valeurs();


  if (prob.que_suis_je() == "Probleme_FT_Disc_gen")
    {
      Transport_Interfaces_FT_Disc& equ_int= ref_cast(Transport_Interfaces_FT_Disc,prob.equation(1));
      const Transport_Interfaces_FT_Disc& equ_int_const = equ_int;
      equ_int.calculer_vitesse_transport_interpolee(inconnue_,
                                                    equ_int_const.maillage_interface(),qIi,1);
      qIj *= -1.*a3;
      qIi -= qIj;
      equ_int.transporter_sans_changement_topologie(qIi,b3*dt_,temps_courant_);
    }

  //fin du pas de temps
  if (prob.que_suis_je() == "Probleme_FT_Disc_gen")
    {
      // pour ne pas deplacer les interfaces:
      inconnue_.valeurs() = 0;
      // on tourne la roue juste apres donc on s'en fout de perdre les valeurs.
      prob.equation(1).mettre_a_jour(temps_courant_);
    }

  prob.equation(0).mettre_a_jour(temps_courant_);
  ////assert(nb_eqn < 3);

  //  pb.postraiter();

  temps_courant_=save_temps_courant;

  // On applique le RK3 sur les equations restantes:
  int premiere_equation;
  if (prob.que_suis_je() == "Probleme_FT_Disc_gen")
    premiere_equation=2;
  else
    premiere_equation=1;
  for(i=premiere_equation; i<nb_eqn; i++)
    RK3::faire_un_pas_de_temps_eqn_base(prob.equation(i));

  converged=true;
  return true;
}


int RK3_FT::faire_un_pas_de_temps_pb_couple(Probleme_Couple& pbc)
{
  // Il s'agit d'un prototype fait pour des problemes couples dont chacun a au plus deux equations :
  // NS et une equation pr le transport de l'interface.
  int i;
  int n_pb= pbc.nb_problemes();
  for(i=0; i<n_pb; i++)
    {
      Probleme_base& pb = ref_cast(Probleme_base,pbc.probleme(i));
      const int nb_eqn=pb.nombre_d_equations();
      if (nb_eqn > 2)
        {
          Cerr << "RK3_FT gere au plus deux equations par pb, votre nb d'equations est : "
               << nb_eqn << finl;
          assert(0);
        }
    }

  double save_temps_courant=temps_courant_;
  double dt_temp;

  double a2=-5./9.;
  double a3=-153./128.;
  double b1=1./3.;
  double b2=15./16.;
  double b3=8./15.;

  REF(Champ_Inc_base) * inconnues = new REF(Champ_Inc_base)[n_pb];
  VECT(DoubleTab) qNSi(n_pb);
  VECT(DoubleTab) qNSj(n_pb);
  DoubleTabFT  qIi;
  DoubleTabFT  qIj;
  for (i=0; i<n_pb; i++)
    {
      // <REF(Champ_Inc_base)> = <Champ_Inc_base>
      Probleme_base& pb = ref_cast(Probleme_base,pbc.probleme(i));
      inconnues[i] = pb.equation(0).inconnue().valeur();
      qNSi[i] = inconnues[i].valeur().valeurs();
      qNSj[i] = inconnues[i].valeur().valeurs();
    }

  //ss pas de temps 1
  //mise a jour des conditions aux limites
  for(i=0; i<n_pb; i++)
    {
      ////Probleme_base& pb = pbc.probleme(i);
      Probleme_base& pb = ref_cast(Probleme_base,pbc.probleme(i));
      const int nb_eqn=pb.nombre_d_equations();
      for(int j=0; j<nb_eqn; j++)
        {
          Equation_base& eqn = pb.equation(j);
          eqn.zone_Cl_dis().mettre_a_jour(temps_courant_);
        }
    }
  temps_courant_ += dt_*1./3.;

  dt_temp=dt_;
  dt_=0.;

  for(i=0; i<n_pb; i++)
    {
      Probleme_base& pb = ref_cast(Probleme_base,pbc.probleme(i));
      pb.equation(0).initTimeStep(dt_);
    }

  dt_=dt_temp;

  for(i=0; i<n_pb; i++)
    {
      Probleme_base& pb = ref_cast(Probleme_base,pbc.probleme(i));
      pb.equation(0).derivee_en_temps_inco(qNSi[i]);

      //on teste si l'etat stationnaire est atteind
      double accroissement_max_abs=qNSi[i].mp_max_abs_vect();
      set_stationnaire_atteint() *= ( accroissement_max_abs < seuil_statio_ );

      inconnues[i].valeur().futur() = inconnues[i].valeur().valeurs();
      inconnues[i].valeur().futur().ajoute_sans_ech_esp_virt(b1 * dt_, qNSi[i], VECT_ALL_ITEMS);

      if (pbc.probleme(i).que_suis_je() == "Probleme_FT_Disc_gen")
        {
          Transport_Interfaces_FT_Disc& equ_int=
            ref_cast(Transport_Interfaces_FT_Disc,pb.equation(1));
          const Transport_Interfaces_FT_Disc& equ_int_const = equ_int;
          equ_int.calculer_vitesse_transport_interpolee(inconnues[i].valeur(),
                                                        equ_int_const.maillage_interface(),
                                                        qIi,
                                                        1);
          equ_int.transporter_sans_changement_topologie(qIi,b1*dt_,temps_courant_);
        }
    }

  //ss pas de temps 2
  //mise a jour des conditions aux limites
  for(i=0; i<n_pb; i++)
    {
      Probleme_base& pb = ref_cast(Probleme_base,pbc.probleme(i));
      pb.equation(0).mettre_a_jour(temps_courant_);
      const int nb_eqn=pb.nombre_d_equations();
      for(int j=0; j<nb_eqn; j++)
        {
          Equation_base& eqn = pb.equation(j);
          eqn.zone_Cl_dis().mettre_a_jour(temps_courant_);
        }
    }
  temps_courant_ += dt_*5./12.;

  dt_temp=dt_;
  dt_=0.;
  for(i=0; i<n_pb; i++)
    {
      Probleme_base& pb = ref_cast(Probleme_base,pbc.probleme(i));
      pb.equation(0).initTimeStep(dt_);
    }
  dt_=dt_temp;

  for(i=0; i<n_pb; i++)
    {
      Probleme_base& pb = ref_cast(Probleme_base,pbc.probleme(i));
      pb.equation(0).derivee_en_temps_inco(qNSj[i]);
      qNSj[i].ajoute_sans_ech_esp_virt(a2, qNSi[i]);
      inconnues[i].valeur().futur() = inconnues[i].valeur().valeurs();
      inconnues[i].valeur().futur().ajoute(b2 * dt_, qNSj[i], VECT_ALL_ITEMS);
      if (pbc.probleme(i).que_suis_je() == "Probleme_FT_Disc_gen")
        {
          Transport_Interfaces_FT_Disc& equ_int=
            ref_cast(Transport_Interfaces_FT_Disc,pb.equation(1));
          const Transport_Interfaces_FT_Disc& equ_int_const = equ_int;
          equ_int.calculer_vitesse_transport_interpolee(inconnues[i].valeur(),
                                                        equ_int_const.maillage_interface(),
                                                        qIj,
                                                        1);
          qIi *= -1.*a2;
          qIj -= qIi;
          equ_int.transporter_sans_changement_topologie(qIj,b2*dt_,temps_courant_);
        }
    }

  //ss pas de temps 3
  //mise a jour des conditions aux limites
  for(i=0; i<n_pb; i++)
    {
      Probleme_base& pb = ref_cast(Probleme_base,pbc.probleme(i));
      pb.equation(0).mettre_a_jour(temps_courant_);
      const int nb_eqn=pb.nombre_d_equations();
      for(int j=0; j<nb_eqn; j++)
        {
          Equation_base& eqn = pb.equation(j);
          eqn.zone_Cl_dis().mettre_a_jour(temps_courant_);
        }
    }
  temps_courant_ += dt_*1./4.;

  dt_temp=dt_;
  dt_=0.;
  for(i=0; i<n_pb; i++)
    {
      Probleme_base& pb = ref_cast(Probleme_base,pbc.probleme(i));
      pb.equation(0).initTimeStep(dt_);
    }
  dt_=dt_temp;

  for(i=0; i<n_pb; i++)
    {
      Probleme_base& pb = ref_cast(Probleme_base,pbc.probleme(i));
      pb.equation(0).derivee_en_temps_inco(qNSi[i]);
      qNSi[i].ajoute_sans_ech_esp_virt(a3, qNSj[i], VECT_ALL_ITEMS);
      inconnues[i].valeur().futur() = inconnues[i].valeur().valeurs();
      inconnues[i].valeur().futur().ajoute_sans_ech_esp_virt(b3 * dt_, qNSi[i], VECT_ALL_ITEMS);
      if (pbc.probleme(i).que_suis_je() == "Probleme_FT_Disc_gen")
        {
          Transport_Interfaces_FT_Disc& equ_int=
            ref_cast(Transport_Interfaces_FT_Disc,pb.equation(1));
          const Transport_Interfaces_FT_Disc& equ_int_const = equ_int;
          equ_int.calculer_vitesse_transport_interpolee(inconnues[i].valeur(),
                                                        equ_int_const.maillage_interface(),
                                                        qIi,
                                                        1);
          qIj *= -1.*a3;
          qIi -= qIi;
          equ_int.transporter_sans_changement_topologie(qIi,b3*dt_,temps_courant_);
        }
    }

  //fin du pas de temps
  for(i=0; i<n_pb; i++)
    {
      Probleme_base& pb = ref_cast(Probleme_base,pbc.probleme(i));

      if ( pb.que_suis_je() == "Probleme_FT_Disc_gen")
        {
          // pour ne pas deplacer les interfaces:
          inconnues[i].valeur().valeurs() = 0;
          // on tourne la roue juste apres donc on s'en fout de perdre les valeurs.
          pb.equation(1).mettre_a_jour(temps_courant_);
        }
      pb.equation(0).mettre_a_jour(temps_courant_);
      // On applique le RK3 sur les equations restantes:
      int premiere_equation;
      if (pb.que_suis_je() == "Probleme_FT_Disc_gen")
        premiere_equation=2;
      else
        premiere_equation=1;
      for(i=premiere_equation; i<pb.nombre_d_equations(); i++)
        RK3::faire_un_pas_de_temps_eqn_base(pb.equation(i));
    }
  temps_courant_=save_temps_courant;
  delete[] inconnues;

  return 1;
}

