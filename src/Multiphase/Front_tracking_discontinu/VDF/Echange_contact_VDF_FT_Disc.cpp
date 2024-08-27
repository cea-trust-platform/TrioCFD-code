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
// File      : Echange_contact_VDF_FT_Disc.cpp
// Directory : $FRONT_TRACKING_DISCONTINU_ROOT/src/VDF
//
/////////////////////////////////////////////////////////////////////////////

#include <Echange_contact_VDF_FT_Disc.h>

#include <Champ_front_calc.h>
#include <Probleme_base.h>
#include <Champ_Uniforme.h>
#include <Schema_Temps_base.h>
#include <Milieu_base.h>
#include <Modele_turbulence_scal_base.h>
#include <Domaine_VDF.h>
#include <Equation_base.h>
#include <Conduction.h>
#include <Param.h>
#include <Probleme_FT_Disc_gen.h>
#include <Triple_Line_Model_FT_Disc.h>
#include <Domaine_VF.h>
#include <Front_VF.h>


Implemente_instanciable( Echange_contact_VDF_FT_Disc, "Echange_contact_VDF_FT_Disc", Echange_contact_VDF ) ;
// XD echange_contact_vdf_ft_disc condlim_base echange_contact_vdf_ft_disc 1 echange_conatct_vdf en prescisant la phase

Sortie& Echange_contact_VDF_FT_Disc::printOn( Sortie& os ) const
{
  Echange_contact_VDF::printOn( os );
  return os;
}

Entree& Echange_contact_VDF_FT_Disc::readOn( Entree& s )
{
  if (app_domains.size() == 0) app_domains = { Motcle("Thermique") };

  //  Echange_contact_VDF::readOn( is );
  Cerr<<"Lecture des parametres du contact (Echange_contact_VDF_FT_Disc::readOn)"<<finl;
  Param param("Echange_contact_VDF_FT_Disc::readOn");
  param.ajouter("autre_probleme",&nom_autre_pb_,Param::REQUIRED); // XD_ADD_P chaine name of other problem
  param.ajouter("autre_bord",&nom_bord,Param::REQUIRED); // XD_ADD_P chaine name of other boundary
  param.ajouter("autre_champ_temperature",&nom_champ,Param::REQUIRED); // XD_ADD_P chaine name of other field
  param.ajouter("nom_mon_indicatrice",&nom_champ_indicatrice_,Param::REQUIRED);  // XD_ADD_P chaine name of indicatrice
  int phase;
  param.ajouter("phase",&phase,Param::REQUIRED); // XD_ADD_P int phase
  param.lire_avec_accolades(s);
  indicatrice_ref_ = double(phase);
  nom_bord_oppose_=nom_bord;

  h_paroi=1e10; // why not git 1/h_paroi = 0....?
  // T_autre_pb(): temp type front from other calculation/Champ dans le domaine
  T_autre_pb().typer("Champ_front_calc");
  // T_ext(): mettre_a_jour utilise des donnees externes,
  // Peut aussi initialized by a champ dans le domaine.
  T_ext().typer("Ch_front_var_instationnaire_dep");
  T_ext()->fixer_nb_comp(1);

  return s;
}
void Echange_contact_VDF_FT_Disc::mettre_a_jour(double temps)
{
  // Champ_front_calc:: mettre_a_jour()
  //  par default distant_= 1
  //  trace the value corresponding from champ inconnu
  //   Champ_Inc_P0_base::trace(), trace the element from distant
  T_autre_pb()->mettre_a_jour(temps);
  indicatrice_->mettre_a_jour(temps);

  const DoubleTab& I = indicatrice_->valeurs_au_temps(temps);

  Champ_front_calc& ch=ref_cast(Champ_front_calc, T_autre_pb().valeur());
  // le_milieu =  SOLID
  const Milieu_base& le_milieu=ch.milieu();
  int nb_comp = le_milieu.conductivite()->nb_comp();
  assert(nb_comp==1);


  int is_pb_fluide=0;


  DoubleTab& mon_h= h_imp_->valeurs();
  int opt=0;
  calculer_h_autre_pb( autre_h, 0., opt);
  // Here, compute h_diff in the fluid side
  calculer_h_mon_pb(mon_h,0.,opt);

  // juste acceder la valeur..., et les remplir
  // pas forcement des chose dedans et valable.
  DoubleTab& mon_Tex= T_ext()-> valeurs();
  calculer_Teta_equiv(mon_Tex,mon_h,autre_h,is_pb_fluide,temps);

  // Attention: Ti_wall_ should be updated after TCL model
  // AS it will be used in TCL model, in particular,
  // influence the heat flux in MESO zone
  // OK here: TCL model is called in UpdateField
  DoubleTab& mon_Ti= Ti_wall_-> valeurs();
  DoubleTab Twalltmp(mon_Ti);
  Twalltmp.detach_vect();
  calculer_Teta_paroi(Twalltmp,mon_h,autre_h,is_pb_fluide,temps);

  int taille=mon_h.dimension(0);

  for (int ii=0; ii<taille; ii++)
    for (int jj=0; jj<nb_comp; jj++)
      {
        if (est_egal(I(ii,0),indicatrice_ref_))
          {
            mon_h(ii,jj)=1./(1./autre_h(ii,jj)+1./mon_h(ii,jj));
            mon_Ti(ii, jj) = Twalltmp(ii, jj);
          }
        else
          {
            // Using Ghost Fluid Method, mon_h and T in pure-L/pure-G cells are NO-PHYSICAL
            // the BC here are not important
            // Setting mon_h = 0 => isothermal BC
            mon_h(ii,jj) = 0.;
            // Update mon_Ti in the same way (to have a good convergence)
            mon_Ti(ii, jj) = Twalltmp(ii, jj);
          }
      }


  Nom nom_pb=mon_dom_cl_dis->equation().probleme().le_nom();
  Probleme_base& pb_gen=ref_cast(Probleme_base, Interprete::objet(nom_pb));
  Probleme_FT_Disc_gen *pbft = dynamic_cast<Probleme_FT_Disc_gen*>(&pb_gen);

  const Domaine_VF& le_dom=ref_cast(Domaine_VF, mon_dom_cl_dis -> domaine_dis().valeur());
  const DoubleVect& surface= le_dom.face_surfaces();


  if (pbft-> tcl().is_activated())
    {

      // phi_ext_ used in *Eval_Diff_VDF_Elem_Gen.tpp* L97
      DoubleTab& mon_phi = phi_ext_-> valeurs();
      mon_phi = 0;
// **************************************To be implemented*******************
      // 2 - phase cells at pb-Boundary when solving T-eq at Liquid side
      //mixed mesh => Text, Twall, mon_h
      if (indicatrice_ref_ == 1)
        {
          Triple_Line_Model_FT_Disc *tcl = pbft ? &pbft->tcl() : nullptr;
          const ArrOfDouble& Q_from_CL = tcl->Q();
          const ArrOfInt& faces_with_CL_contrib = tcl-> boundary_faces();
          const IntTab& face_voisins = le_dom.face_voisins();


          const int nb_contact_line_contribution = faces_with_CL_contrib.size_array();


          for (int jj = 0; jj < nb_comp; jj++)
            {
              // In case of parallelization:
              // raccord in liquid domain and faces_with_CL_contrib are in the same processeur
              // => have the same face Number in the same position
              // use face number to check correspondence
              // fill mon_phi
              for (int idx = 0; idx < nb_contact_line_contribution; idx++)
                {
                  const int facei = faces_with_CL_contrib[idx];
                  bool Not_find_ = true;

                  int ii;
                  for (ii=0; ii<taille; ii++)
                    {
                      const int face = ii + frontiere_dis ().frontiere ().num_premiere_face ();
                      if (facei == face)
                        {
                          Not_find_ = false;
                          break;
                        }
                    }

                  const double sign = (face_voisins (facei, 0) == -1) ? -1. : 1.;
                  const double TCL_wall_flux = Q_from_CL[idx] / surface (facei);
                  const double val = -sign * TCL_wall_flux;

                  if (!Not_find_)
                    mon_phi(ii, jj)  += val;
                  else
                    Process::exit(Nom("Echange_contact_VDF_FT_Disc : missing element corresponding") + Nom(facei) + " ! Check all faces number in TCL are at BC?" );
                }

              // replace mon_h and mon_Ti;
              for (int ii=0; ii<taille; ii++)
                {
                  if (!est_egal(mon_phi(ii, jj), 0.))
                    {
                      mon_Ti(ii, jj) = T_ext()->valeurs()(ii, jj) - mon_phi(ii, jj) /autre_h(ii) ;
                      mon_h(ii,jj) = 0.;
                    }
                }
            }
        }
    }

  // put in the end: to make sure to update the *modified* h_imp_, phi_ext_, and T_ext
  mon_h.echange_espace_virtuel();
  Echange_global_impose::mettre_a_jour(temps);
  mon_Ti.echange_espace_virtuel();
  Ti_wall_->mettre_a_jour(temps);


  // check if to inject a new nuclateion seed
  // only check in the liquid - equation

  if ((pbft->tcl ().reinjection_tcl ()) && (indicatrice_ref_ == 1) )
    {
      bool& ready_inject = pbft->tcl ().ready_inject_tcl ();
      ready_inject = false;

      int BC_has_tcl = 0;

      for (int ii = 0; ii < taille; ii++)
        if (!est_egal (I (ii, 0), indicatrice_ref_))
          {
            BC_has_tcl = 1;
            break;
          }
      BC_has_tcl = Process::mp_max (BC_has_tcl);

      if (BC_has_tcl == 0)
        {
          double sum_T = 0;
          double sum_surface = 0;

          for (int ii = 0; ii < taille; ii++)
            {
              const int face = ii
                               + frontiere_dis ().frontiere ().num_premiere_face ();
              if (le_dom.xv (face, 0) <= pbft->tcl ().Rc_inject ())
                {
                  sum_surface += surface (face);
                  sum_T += surface (face) * mon_Ti (ii, 0);
                }
            }

          sum_T = Process::mp_sum (sum_T);
          sum_surface = Process::mp_sum (sum_surface);

          sum_T = (sum_surface > DMINFLOAT) ? sum_T / sum_surface : 0;

          ready_inject = (sum_T >= pbft->tcl ().tempC_tcl ()) ? true : false;
        }

      ready_inject = Process::mp_max ((int)ready_inject);
    }

}

void Echange_contact_VDF_FT_Disc::completer()
{
  Echange_contact_VDF::completer ();

  indicatrice_.typer("Champ_front_calc");
  Champ_front_calc& ch=ref_cast(Champ_front_calc, indicatrice_.valeur());

  Nom nom_bord_=frontiere_dis().frontiere().le_nom();
  Nom nom_pb = mon_dom_cl_dis->equation ().probleme ().le_nom ();

  int distant=0;
  // when solving pure condution pb for solid
  if (sub_type(Conduction,domaine_Cl_dis().equation()))
    {
      nom_pb=nom_autre_pb_;
      nom_bord_=nom_bord_oppose_;
      distant=1;
    }
  // check the coherance and fixer nb of component
  ch.creer(nom_pb, nom_bord_, nom_champ_indicatrice_);
  // changer distant = 0; for indicatrice_...
  // par default, 1;
  ch.set_distant(distant);
  ch.associer_fr_dis_base(T_ext()->frontiere_dis());
  ch.completer();
  int nb_cases=domaine_Cl_dis().equation().schema_temps().nb_valeurs_temporelles();
  ch.fixer_nb_valeurs_temporelles(nb_cases);

  Probleme_base& pb_gen = ref_cast(Probleme_base, Interprete::objet (nom_pb));

  //ONCE phi_ext_lu_ true, will be completer in father classes

  const Probleme_FT_Disc_gen *pbft =
    dynamic_cast<const Probleme_FT_Disc_gen*> (&pb_gen);

  int nb_faces_raccord1 = frontiere_dis().frontiere().nb_faces ();

  if (pbft->tcl ().is_activated () && phi_ext_lu_ == false)
    {
      phi_ext_lu_ = true;

      derivee_phi_ext_.typer ("Champ_front_fonc");
      derivee_phi_ext_->fixer_nb_comp (1);
      derivee_phi_ext_->associer_fr_dis_base (frontiere_dis ());
      derivee_phi_ext_->valeurs ().resize (nb_faces_raccord1, 1);

      phi_ext_.typer ("Champ_front_fonc");
      phi_ext_->fixer_nb_comp (1);
      phi_ext_->associer_fr_dis_base (frontiere_dis ());
      phi_ext_->valeurs ().resize (nb_faces_raccord1, 1);
    }



  Ti_wall_.typer ("Champ_front_fonc");
  Ti_wall_->fixer_nb_comp (1);
  Ti_wall_->associer_fr_dis_base (frontiere_dis ());
  Ti_wall_->valeurs().resize (nb_faces_raccord1, 1);
  Ti_wall_-> fixer_nb_valeurs_temporelles(nb_cases);

}



/*! @brief Change le i-eme temps futur de la CL.
 *
 */
void Echange_contact_VDF_FT_Disc::changer_temps_futur(double temps,int i)
{
  Echange_contact_VDF::changer_temps_futur(temps,i);
  indicatrice_->changer_temps_futur(temps,i);
  Ti_wall_ -> changer_temps_futur(temps,i);
}

/*! @brief Tourne la roue de la CL
 *
 */
int Echange_contact_VDF_FT_Disc::avancer(double temps)
{
  int ok=Echange_contact_VDF::avancer(temps);
  ok = ok && indicatrice_->avancer(temps);
  ok = ok && Ti_wall_ -> avancer(temps);
  return ok;
}

/*! @brief Tourne la roue de la CL
 *
 */
int Echange_contact_VDF_FT_Disc::reculer(double temps)
{
  int ok=Echange_contact_VDF::reculer(temps);
  ok = ok && indicatrice_->reculer(temps);
  ok = ok && Ti_wall_ -> reculer(temps);
  return ok;
}

int Echange_contact_VDF_FT_Disc::initialiser(double temps)
{

  // T_autre_pb is ALSO created and initialised in the following line
  if (!Echange_contact_VDF::initialiser (temps))
    return 0;

  DoubleTab& mon_Ti = Ti_wall_->valeurs ();

  int is_pb_fluide = 0;
  DoubleTab mon_h (mon_Ti);
  DoubleTab mautre_h (mon_Ti);
  int opt = 0;
  calculer_h_autre_pb (mautre_h, 0., opt);
  // Here, compute h_diffusion in the fluid side
  calculer_h_mon_pb (mon_h, 0., opt);

  // Use Twalltmp to avoid resize distributed array mon_Ti
  // when calling calculer_Teta_paroi
  DoubleTab Twalltmp (mon_Ti);
  Twalltmp.detach_vect ();

  calculer_Teta_paroi (Twalltmp, mon_h, mautre_h, is_pb_fluide, temps);

  int taille = mon_Ti.dimension (0);
  for (int ii = 0; ii < taille; ii++)
    mon_Ti (ii, 0) = Twalltmp (ii);

  Champ_front_calc& chbis=ref_cast(Champ_front_calc, indicatrice_.valeur());
  return chbis.initialiser(temps,domaine_Cl_dis().equation().inconnue());
}

void Echange_contact_VDF_FT_Disc::set_temps_defaut(double temps)
{
  if (Ti_wall_.non_nul())
    Ti_wall_->set_temps_defaut(temps);
  if (indicatrice_.non_nul())
    indicatrice_->set_temps_defaut(temps);
  Echange_global_impose::set_temps_defaut(temps);
}
double Echange_contact_VDF_FT_Disc::Ti_wall(int i) const
{

  if (Ti_wall_->valeurs().size() == 1)
    return Ti_wall_->valeurs()(0, 0);
  else if (Ti_wall_->valeurs().dimension(1) == 1)
    return Ti_wall_->valeurs()(i, 0);
  else
    Cerr << "Echange_contact_VDF_FT_Disc::Ti_wall_ erreur" << finl;
  exit();
  return 0.;
}
