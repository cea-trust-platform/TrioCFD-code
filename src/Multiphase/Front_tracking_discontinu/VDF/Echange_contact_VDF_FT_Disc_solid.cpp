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


Implemente_instanciable( Echange_contact_VDF_FT_Disc_solid, "Echange_contact_VDF_FT_Disc_solid", Echange_contact_VDF_FT_Disc ) ;
// XD echange_contact_vdf_ft_disc_solid condlim_base echange_contact_vdf_ft_disc_solid 1 echange_conatct_vdf en prescisant la phase

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
  T2_autre_pb_.mettre_a_jour(temps);
  T_autre_pb_.mettre_a_jour(temps);
  indicatrice_.mettre_a_jour(temps);
  int nb_comp;
  {
    Champ_front_calc& ch=ref_cast(Champ_front_calc, T_autre_pb().valeur());
    const Milieu_base& le_milieu=ch.milieu();
    nb_comp = le_milieu.conductivite()->nb_comp();
    assert(nb_comp==1);
  }
  const DoubleTab& I = indicatrice_.valeur().valeurs_au_temps(temps);

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



  Probleme_base& pb_gen = ref_cast(Probleme_base, Interprete::objet (nom_autre_pb_));
  if (sub_type(Probleme_FT_Disc_gen, pb_gen))
    {
      const Probleme_FT_Disc_gen *pbft = dynamic_cast<const Probleme_FT_Disc_gen*> (&pb_gen);
      if (pbft->tcl ().is_activated ())
        {
          DoubleTab& mon_phi = phi_ext_->valeurs ();
          mon_phi = 0;
          for (int n = 0; n < 2; n++)
            {
              numero_T_ = n;
              int taille = mon_h.dimension (0);
              double I_ref_ = 1.;
              if (n == 1)
                I_ref_ = 0;

              for (int ii = 0; ii < taille; ii++)
                {
                  if (I (ii, 0) > 0 && I (ii, 0) < 1 && (I_ref_ == 1))
                    {
                      for (int jj = 0; jj < nb_comp; jj++)
                        {

                          const Triple_Line_Model_FT_Disc *tcl = pbft ? &pbft->tcl () : nullptr;

                          const ArrOfInt& faces_with_CL_contrib = tcl->boundary_faces ();
                          const ArrOfDouble& Q_from_CL = tcl->Q ();

                          const Domaine_VF& le_dom = ref_cast( Domaine_VF, mon_dom_cl_dis->domaine_dis ().valeur ());
                          const IntTab& face_voisins = le_dom.face_voisins ();
                          const DoubleVect& surface = le_dom.face_surfaces ();

                          // const int face = ii+ frontiere_dis().frontiere().num_premiere_face();
                          // const int face = ii+ T_ext().frontiere_dis().frontiere().num_premiere_face();

                          const Champ_front_calc& ch = ref_cast( Champ_front_calc, T_autre_pb ().valeur ());
                          const Front_VF& front_vf = ref_cast(Front_VF, ch.front_dis ());
                          // num_face in Liquid-Domaine
                          const int face = ii + front_vf.num_premiere_face ();

                          const int nb_contact_line_contribution = faces_with_CL_contrib.size_array ();
                          int nb_contrib = 0;
                          double flux_local = 0.;
                          hh_imp (ii, jj) = 0.;

                          for (int idx = 0; idx < nb_contact_line_contribution;
                               idx++)
                            {
                              // element i

                              const int facei = faces_with_CL_contrib[idx];

                              if (facei == face)
                                {
                                  nb_contrib++;
                                  const double sign = (face_voisins (face, 0) == -1) ? -1. : 1.;

                                  // const int elemi = face_voisins(faces, 0)+face_voisins(faces, 1)+1;
                                  const int faces = ii + frontiere_dis ().frontiere ().num_premiere_face ();
                                  const double TCL_wall_flux = Q_from_CL[idx]
                                                               / surface (faces);
                                  // val should be : -rho*Cp * flux(W)
                                  // probably because the whole energy equation is written with rhoCp somewhere...
                                  // and the sign should be negative for incoming flux (towards the fluid) by convention.
                                  const double val = -sign * TCL_wall_flux;
                                  if (nb_contrib == 1)
                                    flux_local = val;
                                  // hh_imp(ii,jj) = val/( mon_inco(elemi, 0) - Text(ii));
                                  else
                                    flux_local += val;
                                  // hh_imp(ii,jj) += val/( mon_inco(elemi, 0) - Text(ii));
                                }
                            }
                          mon_phi (ii, jj) += flux_local;

                          const Equation_base& mon_eqn =
                            domaine_Cl_dis ().equation ();
                          const DoubleTab& mon_inco =
                            mon_eqn.inconnue ().valeurs ();

                          const int faces = ii + frontiere_dis ().frontiere ().num_premiere_face ();
                          const int elemi = face_voisins (faces, 0) + face_voisins (faces, 1) + 1;

                          mon_Ti (ii, jj) = mon_inco (elemi, 0) + flux_local / mon_h (ii);
                        }

                    }
                }
            }
        }
    }

  numero_T_=0;
  // put in the end: to make sure to update the *modified* h_imp_, phi_ext_, and Text
  Echange_global_impose::mettre_a_jour(temps);
  Ti_wall_.mettre_a_jour(temps);

}


void Echange_contact_VDF_FT_Disc_solid::completer()
{
  Echange_contact_VDF_FT_Disc::completer();
  T2_autre_pb_.typer("Champ_front_calc");
  Champ_front_calc& ch=ref_cast(Champ_front_calc, T2_autre_pb_.valeur());


  Nom nom_bord_=frontiere_dis().frontiere().le_nom();
  Nom nom_pb=domaine_Cl_dis().equation().probleme().le_nom();
  int distant=0;
  if (sub_type(Conduction,domaine_Cl_dis().equation()))
    {
      nom_pb=nom_autre_pb_;
      nom_bord_=nom_bord_oppose_;
      distant=1;
    }
  else
    {
      abort() ;
    }
  ch.creer(nom_pb, nom_bord_, nom_champ_T2_autre_pb_);
  ch.set_distant(distant);

  ch.associer_fr_dis_base(T_ext().frontiere_dis());

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





  Champ_front_calc& cha=ref_cast(Champ_front_calc, T_autre_pb().valeur());
  cha.creer(nom_autre_pb_, nom_bord, nom_champ);

  const Milieu_base& le_milieu = cha.milieu ();
  int nb_comp = le_milieu.conductivite ()->nb_comp ();

  Nom nom_racc1 = frontiere_dis ().frontiere ().le_nom ();
  Domaine_dis_base& domaine_dis1 = domaine_Cl_dis ().domaine_dis ().valeur ();
  int nb_faces_raccord1 =
    domaine_dis1.domaine ().raccord (nom_racc1).valeur ().nb_faces ();

  Probleme_base& pb_gen = ref_cast(Probleme_base, Interprete::objet (nom_autre_pb_));

  if (sub_type(Probleme_FT_Disc_gen, pb_gen))
    {
      const Probleme_FT_Disc_gen *pbft = dynamic_cast<const Probleme_FT_Disc_gen*> (&pb_gen);

      if (pbft->tcl ().is_activated ())
        {
          phi_ext_lu_ = true;

          derivee_phi_ext_.typer ("Champ_front_fonc");
          derivee_phi_ext_->fixer_nb_comp (nb_comp);
          derivee_phi_ext_->associer_fr_dis_base (frontiere_dis ());
          derivee_phi_ext_.valeurs ().resize (nb_faces_raccord1, nb_comp);

          phi_ext_.typer ("Champ_front_fonc");
          phi_ext_->fixer_nb_comp (nb_comp);
          phi_ext_->associer_fr_dis_base (frontiere_dis ());
          phi_ext_.valeurs ().resize (nb_faces_raccord1, nb_comp);

        }
    }

  Champ_front_calc& ch=ref_cast(Champ_front_calc, T2_autre_pb_.valeur());
  return ch.initialiser(temps,domaine_Cl_dis().equation().inconnue());
}

