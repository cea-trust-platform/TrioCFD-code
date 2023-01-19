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
// File:        Pb_QC_base.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Quasi_Compressible
// Version:     /main/19
//
//////////////////////////////////////////////////////////////////////////////

#include <Pb_QC_base.h>
#include <Equation_base.h>
#include <Fluide_Quasi_Compressible.h>
#include <Schema_Euler_explicite.h>
#include <Pred_Cor.h>
#include <RRK2.h>
#include <RK3.h>
#include <RK3b.h>
#include <Debog.h>
#include <Schema_Euler_Implicite.h>
#include <Zone.h>
#include <Loi_Etat_GP.h>
#include <Navier_Stokes_std.h>
#include <Navier_Stokes_QC_impl.h>
#include <Loi_Fermeture_base.h>
#include <Probleme_Couple.h>

Implemente_base(Pb_QC_base,"Pb_QC_base",Pb_qdm_fluide);

Sortie& Pb_QC_base::printOn(Sortie& os) const
{
  return Pb_qdm_fluide::printOn(os);
}


Entree& Pb_QC_base::readOn(Entree& is)
{
  return Pb_qdm_fluide::readOn(is);
}


void Pb_QC_base::associer_milieu_base(const Milieu_base& mil)
{
  if (sub_type(Fluide_Quasi_Compressible,mil))
    {
      Pb_qdm_fluide::associer_milieu_base(mil);
    }
  else
    {
      Cerr << "Un milieu de type " << mil.que_suis_je()
           << " ne peut etre associe a "<< finl
           << "un probleme quasi-compressible" << finl;
      exit();
    }
}


void Pb_QC_base::associer_sch_tps_base(const Schema_Temps_base& sch)
{
  if (!sub_type(Schema_Euler_explicite,sch)
      && !sub_type(Schema_Euler_Implicite,sch)
      && !sub_type(RRK2,sch)
      && !sub_type(RK3,sch)
      && !sub_type(RK3b,sch) )
    {
      Cerr << finl;
      Cerr << "TRUST can't solve a " << que_suis_je() << finl;
      Cerr << "with a " << sch.que_suis_je() << " time scheme." << finl;
      exit();
    }
  if ( sub_type(Schema_Euler_Implicite,sch) && coupled_==1 )
    {
      Cerr << finl;
      Cerr << "Coupled problem with unique Euler implicit time scheme with quasi-compressible fluid are not allowed!" << finl;
      Cerr << "Choose another time scheme or set a time scheme for each of your problems." << finl;
      exit();
    }
  Pb_qdm_fluide::associer_sch_tps_base(sch);
}


void Pb_QC_base::preparer_calcul()
{
  // WEC : Attention le double appel a le_fluide.preparer_calcul() est necessaire !
  // (je ne sais pas pourquoi, vu sur le cas test Champs_fonc_QC_rayo_semi_transp)
  Fluide_Quasi_Compressible& le_fluide = ref_cast(Fluide_Quasi_Compressible,milieu());

  //Le nom sera a nouveau modifie dans Loi_Etat_GR_rhoT::initialiser()
  if (le_fluide.type_fluide()=="Gaz_Reel")
    {
      equation(1).inconnue()->nommer("temperature");
    }

  le_fluide.completer(*this);
  le_fluide.preparer_calcul();
  Pb_qdm_fluide::preparer_calcul();
  Loi_Etat_GP& la_loi = ref_cast(Loi_Etat_GP,le_fluide.loi_etat().valeur());		//modif YB 28/08/09
  // if (reprise_effectuee()*0) //do not work with gcc 7.1.1
# if 0
  {
    // on a stocke rho dans inco chaleur a la suvegarde , je reviens en T°
    la_loi.remplir_T();
    exit();
  }
#endif
  le_fluide.calculer_masse_volumique();
  la_loi.remplir_T();			//modif YB 28/08/09
  le_fluide.preparer_calcul();
}


bool Pb_QC_base::initTimeStep(double dt)
{
  bool ok=Pb_qdm_fluide::initTimeStep(dt);
  Fluide_Quasi_Compressible& le_fluide = ref_cast(Fluide_Quasi_Compressible,milieu());
  le_fluide.preparer_pas_temps();
  return ok;
}


void Pb_QC_base::mettre_a_jour(double temps)
{
  Debog::set_nom_pb_actuel(le_nom());
  equation(1).mettre_a_jour(temps);
  equation(0).mettre_a_jour(temps);

  for(int i=2; i<nombre_d_equations(); i++) // FACTORISABLE ????
    equation(i).mettre_a_jour(temps);


  // WEC : si ca pouvait venir ici, ce serait factorise avec Pb_base
  //   for(int i=0; i<nombre_d_equations(); i++)
  //     equation(i).mettre_a_jour(temps);
  les_postraitements.mettre_a_jour(temps);
  domaine().mettre_a_jour(temps,domaine_dis(),*this);
  for  (auto& itr : liste_loi_fermeture_)
    {
      Loi_Fermeture_base& loi=itr.valeur();
      loi.mettre_a_jour(temps);
    }
}

// methode simulatant equation_base_mettre_a_jour
// on cherche a en avoir le minimum pour faire comme en incompressible
void equation_base_mettre_a_jour(double temps, Equation_base& mon_eq)
{
  mon_eq.milieu().mettre_a_jour(temps);
  mon_eq.inconnue().mettre_a_jour(temps);
  return;
  mon_eq.zone_Cl_dis()->avancer(temps);
}


bool Pb_QC_base::iterateTimeStep(bool& converged)
{
  Debog::set_nom_pb_actuel(le_nom());

  // WEC 05/01/2005 : reimporte la partie specifique de Resoudre_Qcomp

  Fluide_Quasi_Compressible& le_fluide = ref_cast(Fluide_Quasi_Compressible,milieu());
  Schema_Temps_base& sch=schema_temps();

  if (schema_temps().lsauv())
    {
      sauver();
    }
  int test_a_ecrire=(sch.que_suis_je()!="Runge_Kutta_ordre_3_QC");
  if (test_a_ecrire) // Schemas normaux en temps
    {
      //REMPLACEMENT PAR LE SCHEMA DE RESOLUTION PARTICULIER AU QUASI COMPRESSIBLE
      double temps_present=sch.temps_courant();
      double temps_futur=temps_present+sch.pas_de_temps();


//1. Resoudre la conservation de la masse
      for (int i=1; i<nombre_d_equations(); i++)
        {
          sch.faire_un_pas_de_temps_eqn_base(equation(i));
          equation_base_mettre_a_jour( temps_futur, equation(i));
        }

//1. bis la masse volumique du fluide recupere l'inconnue de inco-ch
      le_fluide.calculer_masse_volumique();

//2. Resoudre l'EDO sur la pression
      le_fluide.Resoudre_EDO_PT();

//3. Calcule T a partir de la loi d'etat sans les proprietes du gaz (lambda, mu, alpha, ...)
      le_fluide.loi_etat().valeur().remplir_T();

//4. Resoudre l'equation de quantite de mouvement
      sch.faire_un_pas_de_temps_eqn_base(equation(0));

//5. Calcule T a partir de la loi d'etat puis calcule les proprietes du gaz (lambda, mu, alpha, ...)
      le_fluide.calculer_coeff_T();

      equation_base_mettre_a_jour( temps_futur, equation(0));

    }
  else // if (test_a_ecrire)
    {
      // Schemas RK3 specifique au QC
      ArrOfDouble a(4),b(4),dt_intermediaire(4),t_i(4),t_i2(4);
      a(1)=0;
      b(1)=1./3.;
      a(2)=-5./9.;
      b(2)=15./16.;
      a(3)=-153./128.;
      b(3)=8./15.;

      double temps_present=sch.temps_courant();
      double dt=sch.pas_de_temps();

      // A CORRIGER quand on corrigera evaluation de Pthn+k
      //      double dtprecedent=dt;
      dt_intermediaire(1)=1./3.*dt;
      dt_intermediaire(2)=5./12.*dt;
      dt_intermediaire(3)=1./4.*dt;
      t_i(1)=1./3.;

      t_i(2)=t_i(1)+5./12.; //3./4.; // 5./12.;
      t_i(3)=1.;

      t_i2(1)=1./3.;
      t_i2(2)=5./12.*3./2.;
      t_i2(3)=1.;
      double temps_futur=temps_present+dt;

      DoubleTab Uref(equation(0).inconnue().valeurs());

      int nb_ss_pas=3;
      DoubleTab& tab_rho_elem = le_fluide.masse_volumique();
      DoubleTab Frho(equation(1).inconnue().valeurs());
      Frho=0;
      DoubleTab Fu(equation(0).inconnue().valeurs());
      Fu=0;
      //double dt_passe=dtprecedent/4.;
      double Fk =0;
      for (int dti=1; dti<nb_ss_pas+1; dti++)
        {
          // RK3
          //1. Resoudre la conservation de la masse
          for (int i=1; i<nombre_d_equations(); i++)
            {
              if (0)
                {
                  // test
                  const Navier_Stokes_std& eqns=ref_cast(Navier_Stokes_std,equation(0));
                  DoubleTab uu(get_champ("pression").valeurs());
                  eqns.operateur_divergence().calculer(equation(0).inconnue().valeurs(),uu);
                  Cerr<<" iii "<<mp_max_abs_vect(uu)<<finl;

                }
              Equation_base& eqn=equation(i);

              DoubleTab& futur   = eqn.inconnue().futur();
              DoubleTab& present   = eqn.inconnue().valeurs();
              futur=present;
              DoubleTab dIdt(futur);

              eqn.derivee_en_temps_inco(dIdt);

              Frho*=a[dti];
              Frho+=dIdt;
              futur.ajoute(dt*b[dti],Frho);
              eqn.zone_Cl_dis().imposer_cond_lim(eqn.inconnue(),temps_futur);

              schema_temps().update_critere_statio(dIdt, eqn);

              assert(i==1);

              present=futur;
              equation_base_mettre_a_jour( temps_futur, equation(i));
            }

          le_fluide.calculer_masse_volumique();
          tab_rho_elem = le_fluide.masse_volumique();

          sch.set_dt()=dt_intermediaire(dti);
          //2. Resoudre l'EDO sur la pression
          // A CORRIGER <---- encore d'actualitee? F.A 24/03/11
          le_fluide.Resoudre_EDO_PT(); // fonction remise en fonctionnement le 3/05/11 pourquoi l'a ton enlever?
          // dt_passe=dt_intermediaire(dti);

          //3. Calcule T a partir de la loi d'etat sans les proprietes du gaz (lambda, mu, alpha, ...)
          le_fluide.loi_etat().valeur().remplir_T();

          //4. Resoudre l'equation de quantite de mouvement
          // important pour la projection
          Equation_base& eqn=equation(0);

          // FAUX <---- encore d'actualitee? F.A 24/03/11
          DoubleTab& futur   = eqn.inconnue().futur();
          DoubleTab& present   = eqn.inconnue().valeurs();
          futur=present;

          eqn.zone_Cl_dis().imposer_cond_lim(eqn.inconnue(),temps_futur);
          //	 { int i=0; Cerr<<i <<" la"<<dti<<" "<< futur(i)<<" "<< present(i)<<" "<<Uref(i)<<endl; }
          int nbval=futur.dimension_tot(0);
          // pour l instant on n essaye pas les val interpolees pour u
          nbval=0;

          for (int i=0; i<nbval; i++)
            {
              if (futur(i)!=present(i))
                {
                  //	 Cerr<<i <<" ici "<< futur(i)<<" "<< present(i)<<" "<<Uref(i)<<endl;
                  //present(i)+=(futur(i)-present(i))*t_i(dti);
                  //present(i)=Uref(i)+(futur(i)-Uref(i))*t_i(dti);
                  present(i)=present(i)+(futur(i)-present(i))*t_i2(dti);
                  Cerr<<i <<" ici2 "<< futur(i)<<" "<< present(i)<<" "<<Uref(i)<<endl;
                  //present(i)=futur(i);
                  futur(i)=present(i);
                }
            }
          present=futur;

          DoubleTab dudt(futur);

          Navier_Stokes_QC_impl& nsimp=dynamic_cast<Navier_Stokes_QC_impl&>(eqn);
          Navier_Stokes_std& eqns=ref_cast(Navier_Stokes_std,eqn);

          nsimp.derivee_en_temps_impl_p1(eqns,dudt,eqns.fluide(),eqns.matrice_pression(),eqns.assembleur_pression(),0);
          // ne peut pas etre commenter sans provoquer une erreur : "flux_bords_ not dimensioned for the operator Op_Conv_Centre_VDF_Face"


          Fu*=a[dti];
          Fk*=a[dti];
          Fu+=dudt;
          Fk+= 1;
          dudt=Fu;
          dudt/=Fk;

          nsimp.derivee_en_temps_impl_projection(eqns,dudt,eqns.fluide(),eqns.matrice_pression(),eqns.assembleur_pression(),0);

          futur.ajoute(dt_intermediaire[dti],dudt);

          sch.set_dt()=dt;
          //eqn.zone_Cl_dis().imposer_cond_lim(eqn.inconnue(),temps_futur);
          schema_temps().update_critere_statio(dudt, eqn);
          present=futur;

          //5. Calcule T a partir de la loi d'etat puis calcule les proprietes du gaz (lambda, mu, alpha, ...)
          le_fluide.calculer_coeff_T();


          equation_base_mettre_a_jour( temps_futur, equation(0));
          equation(0).inconnue().reculer();
        }
      assert(est_egal(temps_futur,temps_present+sch.pas_de_temps()));
    }
  //Resolution des eqations supplmentaires dans le cas d un probleme quasi_compressible
  //avec un liste d quations de scalaire
  //La classe mere (Pb_Thermohydraulique_QC ou Pb_Thermohydraulique_Turbulent_QC) possede deux equations

  for (int i=0; i<nombre_d_equations(); i++)
    {
      // on recule les inconnues ( le pb mettra a jour les equations)
      equation(i).inconnue().reculer();
    }
  //FIN REMPLACEMENT PAR LE SCHEMA DE RESOLUTION PARTICULIER AU QUASI COMPRESSIBLE
  // Calculs coeffs echange sur l'instant sur lequel doivent agir les operateurs.
  double tps=schema_temps().temps_defaut();
  for(int i=0; i<nombre_d_equations(); i++)
    {
      equation(i).zone_Cl_dis()->calculer_coeffs_echange(tps);
    }

  converged=true;
  return true;
}
