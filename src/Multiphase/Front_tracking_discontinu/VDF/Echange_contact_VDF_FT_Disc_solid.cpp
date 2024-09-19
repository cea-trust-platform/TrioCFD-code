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
// File      : Echange_contact_VDF_FT_Disc_solid.cpp
// Directory : $FRONT_TRACKING_DISCONTINU_ROOT/src/VDF
//
/////////////////////////////////////////////////////////////////////////////

#include <Echange_contact_VDF_FT_Disc_solid.h>

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
#include <Interprete.h>
#include <Probleme_FT_Disc_gen.h>
#include <Triple_Line_Model_FT_Disc.h>
#include <Domaine_Cl_VDF.h>
#include <EcrFicPartage.h>


Implemente_instanciable( Echange_contact_VDF_FT_Disc_solid, "Echange_contact_VDF_FT_Disc_solid", Echange_contact_VDF_FT_Disc ) ;
// XD echange_contact_vdf_ft_disc_solid condlim_base echange_contact_vdf_ft_disc_solid 1 echange_conatct_vdf en prescisant la phase

int meme_point2(const DoubleVect& a,const DoubleVect& b);

Sortie& Echange_contact_VDF_FT_Disc_solid::printOn( Sortie& os ) const
{
  Echange_contact_VDF_FT_Disc::printOn( os );
  return os;
}

Entree& Echange_contact_VDF_FT_Disc_solid::readOn( Entree& s )
{
  if (app_domains.size() == 0) app_domains = { Motcle("Thermique") };

  //  Echange_contact_VDF_FT_Disc::readOn( is );
  Cerr<<"Lecture des parametres du contact (Echange_contact_VDF_FT_Disc_solid::readOn)"<<finl;
  Param param("Echange_contact_VDF_FT_Disc_solid::readOn");
  param.ajouter("autre_probleme",&nom_autre_pb_,Param::REQUIRED); // XD_ADD_P chaine name of other problem
  param.ajouter("autre_bord",&nom_bord,Param::REQUIRED);  // XD_ADD_P chaine name of other boundary
  param.ajouter("autre_champ_temperature_indic1",&nom_champ,Param::REQUIRED); // XD_ADD_P chaine name of temperature indic 1
  param.ajouter("autre_champ_temperature_indic0",&nom_champ_T2_autre_pb_,Param::REQUIRED); // XD_ADD_P chaine name of temperature indic 0
  param.ajouter("autre_champ_indicatrice",&nom_champ_indicatrice_,Param::REQUIRED); // XD_ADD_P chaine name of indicatrice
  param.ajouter("dt_impr_Tw",&dt_impr_Tw_);
  param.lire_avec_accolades(s);

  nom_bord_oppose_=nom_bord;
  h_paroi=1e10;
  numero_T_=0;
  T_autre_pb().typer("Champ_front_calc");
  T_ext().typer("Ch_front_var_instationnaire_dep");
  T_ext()->fixer_nb_comp(1);
  return s;
}
void Echange_contact_VDF_FT_Disc_solid::mettre_a_jour(double temps)
{
  // update T_autre pb1/2
  T2_autre_pb_->mettre_a_jour(temps);
  T_autre_pb_->mettre_a_jour(temps);
  indicatrice_->mettre_a_jour(temps);
  int nb_comp;
  {
    Champ_front_calc& ch=ref_cast(Champ_front_calc, T_autre_pb().valeur());
    const Milieu_base& le_milieu=ch.milieu();
    nb_comp = le_milieu.conductivite()->nb_comp();
    assert(nb_comp==1);
  }
  const DoubleTab& I = indicatrice_->valeurs_au_temps(temps);

  int is_pb_fluide=0;

  DoubleTab& hh_imp= h_imp_->valeurs();
  hh_imp=0;
  DoubleTab mon_h(hh_imp);
  DoubleTab& Text=T_ext()-> valeurs();
  DoubleTab& mon_Ti= Ti_wall_-> valeurs();

  DoubleTab Texttmp(Text);
  DoubleTab Twalltmp(Text);
  Twalltmp.detach_vect();

  int opt=0;
  // h of solid
  calculer_h_mon_pb(mon_h,0.,opt);
  // need to overwrite mon_h by h_micro in 2-phase cells
  // 1, compoute Q_micro

  // numero_T = 0 T_autre_pb_
  // numero_T = 1 T2_autre_pb_
  for( int n=0; n<2; n++)
    {
      numero_T_=n;
      // h of fluid
      calculer_h_autre_pb( autre_h, 0., opt);

      calculer_Teta_paroi(Twalltmp,mon_h,autre_h,is_pb_fluide,temps);
      // calculer_Teta_equiv(Text,mon_h,autre_h,is_pb_fluide,temps);
      calculer_Teta_equiv(Texttmp,mon_h,autre_h,is_pb_fluide,temps);
      // on a calculer Teta paroi, on peut calculer htot dans himp (= mon_h)
      int taille=mon_h.dimension(0);
      double I_ref_=1.;
      if (n==1)
        I_ref_=0;
      for (int ii=0; ii<taille; ii++)
        {
          if (est_egal(I(ii,0),I_ref_))
            for (int jj=0; jj<nb_comp; jj++)
              {
                hh_imp(ii,jj)=1./(1./autre_h(ii,jj)+1./mon_h(ii,jj));

                Text(ii,jj)=Texttmp(ii,jj);
                mon_Ti(ii) = Twalltmp(ii, jj);
              }
        }
    }
  const Domaine_VF& le_dom = ref_cast(
                               Domaine_VF, mon_dom_cl_dis->domaine_dis());


  Probleme_base& pb_gen = ref_cast(Probleme_base, Interprete::objet (nom_autre_pb_));
  const Probleme_FT_Disc_gen *pbft = dynamic_cast<const Probleme_FT_Disc_gen*> (&pb_gen);
  if (pbft->tcl ().is_activated ())
    {
      DoubleTab& mon_phi = phi_ext_->valeurs ();
      mon_phi = 0;

      // Numero_T = 0 corresponding I_ref_=1. Liquid side
      // use T_autre_pb_ for ref_cast
      numero_T_=0;
      const Champ_front_calc& ch = ref_cast(Champ_front_calc, T_autre_pb ().valeur ());
      const Domaine_Cl_dis_base& zcldis = ch.domaine_Cl_dis();
      const Domaine_Cl_VDF& zclvdf = ref_cast(Domaine_Cl_VDF, zcldis);
      const Front_VF& front_vf = ref_cast(Front_VF, ch.front_dis ());


      const Cond_lim_base& la_cl = zclvdf.condition_limite_de_la_frontiere(front_vf.frontiere().le_nom());

      if (sub_type(Echange_contact_VDF_FT_Disc, la_cl))
        {
          const Echange_contact_VDF_FT_Disc& la_cl_typee = ref_cast(
                                                             Echange_contact_VDF_FT_Disc, la_cl);

          const DoubleTab& autre_phi  = la_cl_typee.phi_ext()->valeurs ();



          // trace_elem_xxx(x, y) is elelmet-based function, X should be is indiced with element number
          // when several paves are used (progressive mesh), the elem-number is not continued at BC
          // create an elelmet-based Field phi_filed, filled the value of autre_phi at the boundary


          const Equation_base& autre_eqn = ref_cast(Equation_base, ch.domaine_Cl_dis().equation());
          const DoubleTab& autre_inco =autre_eqn.inconnue().valeurs();

          DoubleTab phi_filed(autre_inco);
          phi_filed = 0.;

          const Domaine_dis_base& domainedis = ref_cast(Domaine_dis_base, ch.domaine_dis());
          const IntTab& face_voisins_loc = domainedis.face_voisins();

          for (int ii = 0; ii < autre_phi.dimension (0); ii++)
            {
              const int face_loc = ii + front_vf.frontiere().num_premiere_face();
              int elem = face_voisins_loc(face_loc, 0);
              if (elem == -1)
                elem = face_voisins_loc(face_loc, 1);
              phi_filed(elem) = autre_phi(ii);
            }

          Nom nom_racc1=frontiere_dis().frontiere().le_nom();
          if (mon_dom_cl_dis -> domaine().raccord(nom_racc1)->que_suis_je() =="Raccord_distant_homogene")
            {
              // front_vf.frontiere ().trace_elem_distant (autre_phi, mon_phi);
              front_vf.frontiere ().trace_elem_distant (phi_filed, mon_phi);
            }
          else // Raccord_local_homogene
            {
              // front_vf.frontiere ().trace_elem_local (autre_phi, mon_phi);
              front_vf.frontiere ().trace_elem_local (phi_filed, mon_phi);
            }
        }

      const Equation_base& mon_eqn = domaine_Cl_dis().equation();
      const DoubleTab& mon_inco = mon_eqn.inconnue().valeurs();
      const IntTab& face_voisins = le_dom.face_voisins();

      // replace mon_h and mon_Ti;
      int taille = mon_h.dimension(0);
      for (int jj = 0; jj < nb_comp; jj++)
        {
          for (int ii = 0; ii < taille; ii++)
            {
              if (!est_egal(mon_phi(ii, jj), 0.))
                {
                  mon_phi(ii, jj) = -mon_phi(ii, jj);

                  hh_imp(ii, jj) = 0.;

                  const int face = ii + frontiere_dis().frontiere().num_premiere_face();

                  const int elemi = face_voisins(face, 0) + face_voisins(face, 1) + 1;
                  mon_Ti(ii, jj) = mon_inco(elemi, 0) + mon_phi(ii, jj) / mon_h(ii);
                }
            }
        }
    }

  numero_T_=0;
  // put in the end: to make sure to update the *modified* h_imp_, phi_ext_, and Text
  Echange_global_impose::mettre_a_jour(temps);
  Ti_wall_->mettre_a_jour(temps);


  // print Twall


  const Schema_Temps_base& sch = mon_dom_cl_dis->equation().schema_temps();
  double temps_courant = sch.temps_courant();
  double temps_prec = sch.temps_precedent();
  double dt= sch.pas_de_temps() ;


  if (dt_impr_Tw_ != DMAXFLOAT)
    {
      bool is_imp = sch.temps_final_atteint () || sch.nb_pas_dt_max_atteint ();
      is_imp = is_imp || (dt_impr_Tw_ <= dt);

      if (!is_imp)
        {
          // Voir Schema_Temps_base::limpr pour information sur epsilon et modf
          double i, j, epsilon = 1.e-8;
          modf (temps_courant / dt_impr_Tw_ + epsilon, &i);
          modf (temps_prec / dt_impr_Tw_ + epsilon, &j);
          is_imp = ( i>j );
        }

      if (is_imp)
        {

          int ndeb = frontiere_dis ().frontiere ().num_premiere_face ();
          int nfin = ndeb + frontiere_dis ().frontiere ().nb_faces ();


          EcrFicPartage filTwall;
          Nom nom_pb=mon_dom_cl_dis->equation().probleme().le_nom();
          Nom fichier=Objet_U::nom_du_cas()+"_"+nom_pb+"_"+frontiere_dis ().frontiere ().le_nom()+"_"+"twall.face";

          // On cree le fichier au premier pas de temps si il n'y a pas reprise
          if ( est_egal(temps_prec, 0) && !pb_gen.reprise_effectuee())
            {
              filTwall.ouvrir(fichier);
            }
          // Sinon on l'ouvre
          else
            {
              filTwall.ouvrir(fichier,ios::app);
            }



          if(je_suis_maitre())
            {
              filTwall << finl;
              if (dimension == 2)
                {
                  filTwall << "--------------------------------------------------------------------------------------------" << finl;
                  filTwall << "Time\t\t| X\t\t\t| Y\t\t\t| Twall" << finl;
                  filTwall << "--------------------------------------------------------------------------------------------" << finl;
                }
            }


          for (int face = ndeb; face < nfin; face++)
            {
              filTwall << temps << "\t| " << le_dom.xv (face, 0) << "\t| " << le_dom.xv (face, 1) << "\t| " << Ti_wall (face - ndeb) << finl;
            }
          filTwall.syncfile ();
        }
    }

}


void Echange_contact_VDF_FT_Disc_solid::completer()
{
  Echange_contact_VDF_FT_Disc::completer();

  // configure T2_autre pb_
  T2_autre_pb_.typer("Champ_front_calc");
  Champ_front_calc& ch=ref_cast(Champ_front_calc, T2_autre_pb_.valeur());

  ch.creer(nom_autre_pb_, nom_bord_oppose_, nom_champ_T2_autre_pb_);

  ch.associer_fr_dis_base(T_ext()->frontiere_dis());
  ch.completer();
  int nb_cases=domaine_Cl_dis().equation().schema_temps().nb_valeurs_temporelles();
  ch.fixer_nb_valeurs_temporelles(nb_cases);
}



/*! @brief Change le i-eme temps futur de la CL.
 *
 */
void Echange_contact_VDF_FT_Disc_solid::changer_temps_futur(double temps,int i)
{
  Echange_contact_VDF_FT_Disc::changer_temps_futur(temps,i);
  T2_autre_pb_->changer_temps_futur(temps,i);
}

/*! @brief Tourne la roue de la CL
 *
 */
int Echange_contact_VDF_FT_Disc_solid::avancer(double temps)
{
  int ok=Echange_contact_VDF_FT_Disc::avancer(temps);
  ok = ok && T2_autre_pb_->avancer(temps);
  return ok;
}

/*! @brief Tourne la roue de la CL
 *
 */
int Echange_contact_VDF_FT_Disc_solid::reculer(double temps)
{
  int ok=Echange_contact_VDF_FT_Disc::reculer(temps);
  ok = ok && T2_autre_pb_->reculer(temps);
  return ok;
}

int Echange_contact_VDF_FT_Disc_solid::initialiser(double temps)
{

  if (!Echange_contact_VDF_FT_Disc::initialiser(temps))
    return 0;

  Champ_front_calc& ch=ref_cast(Champ_front_calc, T2_autre_pb_.valeur());
  return ch.initialiser(temps,domaine_Cl_dis().equation().inconnue());
}

