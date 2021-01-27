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
#include <Zone_VDF.h>
#include <Equation_base.h>
#include <Conduction.h>
#include <Param.h>


Implemente_instanciable( Echange_contact_VDF_FT_Disc, "Echange_contact_VDF_FT_Disc", Echange_contact_VDF ) ;
// XD echange_contact_vdf_ft_disc condlim_base echange_contact_vdf_ft_disc 1 echange_conatct_vdf en prescisant la phase

Sortie& Echange_contact_VDF_FT_Disc::printOn( Sortie& os ) const
{
  Echange_contact_VDF::printOn( os );
  return os;
}

Entree& Echange_contact_VDF_FT_Disc::readOn( Entree& s )
{
  //  Echange_contact_VDF::readOn( is );
  Cerr<<"Lecture des parametres du contact (Echange_contact_VDF_FT_Disc::readOn)"<<finl;
  Nom nom_pb, nom_bord;
  Motcle nom_champ;
  Param param("Echange_contact_VDF_FT_Disc::readOn");
  param.ajouter("autre_probleme",&nom_pb,Param::REQUIRED); // XD_ADD_P chaine name of other problem
  param.ajouter("autre_bord",&nom_bord,Param::REQUIRED); // XD_ADD_P chaine name of other boundary
  param.ajouter("autre_champ_temperature",&nom_champ,Param::REQUIRED); // XD_ADD_P chaine name of other field
  param.ajouter("nom_mon_indicatrice",&nom_champ_indicatrice_,Param::REQUIRED);  // XD_ADD_P chaine name of indicatrice
  int phase;
  param.ajouter("phase",&phase,Param::REQUIRED); // XD_ADD_P int phase
  param.lire_avec_accolades(s);
  indicatrice_ref_ = double(phase);
  nom_autre_pb_=nom_pb;
  nom_bord_oppose_=nom_bord;

  h_paroi=1e10;
  T_autre_pb().typer("Champ_front_calc");
  Champ_front_calc& ch=ref_cast(Champ_front_calc, T_autre_pb().valeur());
  ch.creer(nom_pb, nom_bord, nom_champ);
  T_ext().typer("Ch_front_var_instationnaire_dep");
  T_ext()->fixer_nb_comp(1);

  return s;
}
void Echange_contact_VDF_FT_Disc::mettre_a_jour(double temps)
{

  Echange_contact_VDF::mettre_a_jour(temps);

  DoubleTab& mon_h= h_imp_->valeurs();

  indicatrice_.mettre_a_jour(temps);
  const DoubleTab& I = indicatrice_.valeur().valeurs_au_temps(temps);
  //Cerr<<zone_Cl_dis().equation().probleme().le_nom()<<" "<<I<<finl;
  int taille=mon_h.dimension(0);
  int nb_comp=mon_h.dimension(1);
  for (int ii=0; ii<taille; ii++)
    for (int jj=0; jj<nb_comp; jj++ )
      if (!est_egal(I(ii,0),indicatrice_ref_))
        {
          mon_h(ii,jj)=0;                                     ;
        }
}


void Echange_contact_VDF_FT_Disc::completer()
{
  Echange_contact_VDF::completer();
  indicatrice_.typer("Champ_front_calc");
  Champ_front_calc& ch=ref_cast(Champ_front_calc, indicatrice_.valeur());


  Nom nom_bord=frontiere_dis().frontiere().le_nom();
  Nom nom_pb=zone_Cl_dis().equation().probleme().le_nom();
  int distant=0;
  if (sub_type(Conduction,zone_Cl_dis().equation()))
    {
      nom_pb=nom_autre_pb_;
      nom_bord=nom_bord_oppose_;
      distant=1;
    }
  ch.creer(nom_pb, nom_bord, nom_champ_indicatrice_);
  ch.set_distant(distant);

  ch.associer_fr_dis_base(T_ext().frontiere_dis());

  ch.completer();

  int nb_cases=zone_Cl_dis().equation().schema_temps().nb_valeurs_temporelles();
  ch.fixer_nb_valeurs_temporelles(nb_cases);

}



// Description:
//    Change le i-eme temps futur de la CL.
void Echange_contact_VDF_FT_Disc::changer_temps_futur(double temps,int i)
{
  Echange_contact_VDF::changer_temps_futur(temps,i);
  indicatrice_->changer_temps_futur(temps,i);
}

// Description:
//    Tourne la roue de la CL
int Echange_contact_VDF_FT_Disc::avancer(double temps)
{
  int ok=Echange_contact_VDF::avancer(temps);
  ok = ok && indicatrice_->avancer(temps);
  return ok;
}

// Description:
//    Tourne la roue de la CL
int Echange_contact_VDF_FT_Disc::reculer(double temps)
{
  int ok=Echange_contact_VDF::reculer(temps);
  ok = ok && indicatrice_->reculer(temps);
  return ok;
}

int Echange_contact_VDF_FT_Disc::initialiser(double temps)
{
  if (!Echange_contact_VDF::initialiser(temps))
    return 0;
  Champ_front_calc& ch=ref_cast(Champ_front_calc, indicatrice_.valeur());
  return ch.initialiser(temps,zone_Cl_dis().equation().inconnue());
}
