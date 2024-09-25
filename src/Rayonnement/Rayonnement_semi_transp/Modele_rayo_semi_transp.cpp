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
// File:        Modele_rayo_semi_transp.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src
// Version:     /main/19
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_rayo_semi_transp.h>
#include <Flux_radiatif_base.h>
#include <Probleme_base.h>
#include <Schema_Temps_base.h>
#include <Symetrie.h>
#include <Frontiere_dis_base.h>
#include <Fluide_base.h>
#include <Frontiere_ouverte_rayo_semi_transp.h>
#include <Frontiere_ouverte_temperature_imposee_rayo_semi_transp.h>
#include <Neumann_paroi_rayo_semi_transp_VEF.h>
#include <Neumann_paroi_rayo_semi_transp_VDF.h>
#include <Echange_externe_impose_rayo_semi_transp.h>
#include <Echange_contact_rayo_semi_transp_VDF.h>
#include <Echange_global_impose_rayo_semi_transp.h>
#include <Temperature_imposee_paroi_rayo_semi_transp.h>
#include <Source_rayo_semi_transp_base.h>
#include <verif_cast.h>
#include <Champ_Uniforme.h>
#include <Discretisation_base.h>
#include <Domaine.h>

Implemente_instanciable(Modele_rayo_semi_transp,"Modele_rayo_semi_transp",Probleme_base);


double const Modele_rayo_semi_transp::sigma = 5.67e-8;

////////////////////////////////////////////////////////////////////
// Description :
//           Lecture du Modele_rayo_semi_transp
//           On lit les parametres de l'equation de rayonnements
// Precondition :
// Parametre :
//     Signification :
//     Valeurs par defaut :
//     Contraintes :
//     Entree :
//     Entree/Sortie :
//     Sortie :
// Retour :
//     Signification :
//     Contraintes :
// Exception :
// Effets de bord :
// Postcondition :
//
Entree& Modele_rayo_semi_transp::readOn(Entree& is)
{
  return Probleme_base::readOn(is);
}

////////////////////////////////////////////////////////////////////
// Description :
//
// Precondition :
// Parametre :
//     Signification :
//     Valeurs par defaut :
//     Contraintes :
//     Entree :
//     Entree/Sortie :
//     Sortie :
// Retour :
//     Signification :
//     Contraintes :
// Exception :
// Effets de bord :
// Postcondition :
//
Sortie& Modele_rayo_semi_transp::printOn(Sortie& os) const
{
  return os;
}

bool Modele_rayo_semi_transp::initTimeStep(double dt)
{
  return eq_rayo().initTimeStep(dt);
}

bool Modele_rayo_semi_transp::iterateTimeStep(bool& converged)
{
  converged=true;
  return eq_rayo().solve();
}

void Modele_rayo_semi_transp::validateTimeStep()
{
  double temps=probleme().presentTime();
  eq_rayo().mettre_a_jour(temps);
  calculer_flux_radiatif();
  les_postraitements_.mettre_a_jour(temps);
  schema_temps().mettre_a_jour();
}


void Modele_rayo_semi_transp::associer_sch_tps_base(const Schema_Temps_base& un_schema_en_temps)
{
  le_schema_en_temps_ = un_schema_en_temps;
  le_schema_en_temps_->associer_pb(*this);
}

////////////////////////////////////////////////////////////////////
// Description :
//       L'appel de cette methode permet d'obtenir le champ
//       de l'irradiance a l'instant courant
// Precondition :
// Parametre :
//     Signification :
//     Valeurs par defaut :
//     Contraintes :
//     Entree :
//     Entree/Sortie :
//     Sortie :
// Retour :
//     Signification :
//     Contraintes :
// Exception :
// Effets de bord :
// Postcondition :
//
Champ_Inc_base& Modele_rayo_semi_transp::put_irradience()
{
  Champ_Inc_base& irradiance = eq_rayo().inconnue();
  return irradiance;
}


int is_la_cl_rayo(const Cond_lim_base& la_cl, Cond_Lim_rayo_semi_transp*& la_cl_rayo_semi_transp)
{
  if (sub_type(Frontiere_ouverte_rayo_semi_transp,la_cl))
    {
      la_cl_rayo_semi_transp =&( verif_cast(Cond_Lim_rayo_semi_transp&, ref_cast(Frontiere_ouverte_rayo_semi_transp,la_cl)));
      return 1;
    }
  if (sub_type(Frontiere_ouverte_temperature_imposee_rayo_semi_transp,la_cl))
    {
      la_cl_rayo_semi_transp =&( verif_cast(Cond_Lim_rayo_semi_transp&, ref_cast(Frontiere_ouverte_temperature_imposee_rayo_semi_transp,la_cl)));
      return 1;
    }
  if (sub_type(Neumann_paroi_rayo_semi_transp_VEF,la_cl))
    {
      la_cl_rayo_semi_transp =&( verif_cast(Cond_Lim_rayo_semi_transp&, ref_cast(Neumann_paroi_rayo_semi_transp_VEF,la_cl)));
      return 1;
    }
  if (sub_type(Echange_externe_impose_rayo_semi_transp,la_cl))
    {
      la_cl_rayo_semi_transp =&( verif_cast(Cond_Lim_rayo_semi_transp&, ref_cast(Echange_externe_impose_rayo_semi_transp,la_cl)));
      return 1;
    }
  if (sub_type(Neumann_paroi_rayo_semi_transp_VDF,la_cl))
    {
      la_cl_rayo_semi_transp =&( verif_cast(Cond_Lim_rayo_semi_transp&, ref_cast(Neumann_paroi_rayo_semi_transp_VDF,la_cl)));
      return 1;
    }
  if (sub_type(Echange_contact_rayo_semi_transp_VDF,la_cl))
    {
      la_cl_rayo_semi_transp =&( verif_cast(Cond_Lim_rayo_semi_transp&, ref_cast(Echange_contact_rayo_semi_transp_VDF,la_cl)));
      return 1;
    }
  if (sub_type(Temperature_imposee_paroi_rayo_semi_transp,la_cl))
    {
      la_cl_rayo_semi_transp =&( verif_cast(Cond_Lim_rayo_semi_transp&, ref_cast(Temperature_imposee_paroi_rayo_semi_transp,la_cl)));
      return 1;
    }
  if (sub_type(Echange_global_impose_rayo_semi_transp,la_cl))
    {
      la_cl_rayo_semi_transp =&( verif_cast(Cond_Lim_rayo_semi_transp&, ref_cast(Echange_global_impose_rayo_semi_transp,la_cl)));
      return 1;
    }

  return 0;
}

void Modele_rayo_semi_transp::preparer_calcul()
{

  int contient_source_rayo_semi_transp=0;

  for (int j=0; j<probleme().nombre_d_equations(); j++)
    {

      // Associer le modele au CL rayonnantes.
      Domaine_Cl_dis_base& la_zcl = probleme().equation(j).domaine_Cl_dis();
      for (int num_cl=0; num_cl<la_zcl.nb_cond_lim(); num_cl++)
        {
          Cond_lim_base& la_cl = la_zcl.les_conditions_limites(num_cl).valeur();
          Cond_Lim_rayo_semi_transp* la_cl_rayo_semi_transp;

          if (is_la_cl_rayo(la_cl,la_cl_rayo_semi_transp))
            {
              la_cl_rayo_semi_transp->associer_modele(*this);
              la_cl_rayo_semi_transp->recherche_emissivite_et_A();
              // Dans le cas d'un echange contact, il faut aussi completer la CL opposee
              if (sub_type(Echange_contact_rayo_semi_transp_VDF,la_cl))
                {
                  Echange_contact_rayo_semi_transp_VDF& la_cl_ech=ref_cast(Echange_contact_rayo_semi_transp_VDF,la_cl);
                  Echange_contact_rayo_semi_transp_VDF& la_cl_opp=la_cl_ech.la_Cl_opposee();
                  la_cl_opp.associer_modele(*this);
                  la_cl_opp.recherche_emissivite_et_A();
                }
            }
        }

      // Associer le modele au terme source de rayonnement de l'equation de temperature
      Sources& les_sources = probleme().equation(j).sources();
      for(int num_source = 0; num_source < les_sources.size(); num_source++)
        {
          if (  (sub_type(Source_rayo_semi_transp_base,les_sources[num_source].valeur()))
                || (les_sources[num_source]->que_suis_je() == "Source_rayo_semi_transp_QC_VDF_P0_VDF")
                || (les_sources[num_source]->que_suis_je() == "Source_rayo_semi_transp_QC_VEF_P1NC") )
            {
              contient_source_rayo_semi_transp = 1;
            }
        }
    }

  if (contient_source_rayo_semi_transp==0)
    {
      Cerr<<"Attention, vous n'avez pas defini de terme source de rayonnement semi transparent"<<finl;
      Cerr<<"pensez a ajourter le terme source Source_rayo_semi_transp dans la liste des termes "<<finl;
      Cerr<<"sources de l'equation de l'energie"<<finl;
      exit();
    }
  eq_rayo().completer();
}


void Modele_rayo_semi_transp::discretiser(Discretisation_base& dis)
{

  // Typage de l'equation de rayonnement
  Cerr<<"typage de l'equation de rayonnement " ;
  Equation_base& eq_base = probleme().equation(1);
  Nom disc = eq_base.discretisation().que_suis_je(), type = "Eq_rayo_semi_transp_";
  if(disc=="VEFPreP1B")
    disc="VEF";

  type+=disc;
  Cerr << type << finl;

  Eq_rayo_.typer(type);
  //
  // Creation des associations pour l'equation de rayonnement.
  // Necessairement ici car, contrairement aux cas classiques, l'equation
  // a besoin d'etre typee !
  //
  Eq_rayo_->associer_modele_rayonnement(*this);
  Eq_rayo_->associer_sch_tps_base(schema_temps());

  //  Probleme_base::discretiser(dis);
  associer();
  la_discretisation_ = dis;
  Cerr << "Discretisation du domaine associe au probleme " << le_nom() << finl;
  if (!le_domaine_.non_nul())
    Process::exit("ERROR: Discretize - You're trying to discretize a problem without having associated a Domain to it!!! Fix your dataset.");
  // Initialisation du tableau renum_som_perio
  le_domaine_->init_renum_perio();

  dis.associer_domaine(le_domaine_.valeur());
  le_domaine_dis_ = dis.discretiser();
  // Can not do this before, since the Domaine_dis is not typed yet:
  le_domaine_dis_->associer_domaine(le_domaine_);

  Cerr << "Discretisation des equations" << finl;
  for(int i=0; i<nombre_d_equations(); i++)
    {
      equation(i).associer_domaine_dis(domaine_dis());
      equation(i).discretiser();
    }

  // Association du fluide + diverses operations

  if (sub_type(Fluide_base,probleme().milieu()))
    {

      Fluide_base& fluide = ref_cast(Fluide_base,probleme().milieu());
      Eq_rayo_->associer_fluide(fluide);

      if (fluide.is_rayo_semi_transp())
        {
          Champ_Don& l_rayo = fluide.longueur_rayo();
          Champ_Don& coeff_abs = fluide.kappa();

          if (sub_type(Champ_Uniforme,coeff_abs.valeur()))
            {

              // Typage du Champ_Don longueur_rayo_ comme un champ_uniforme
              l_rayo.typer("Champ_Uniforme");
              Champ_Uniforme& ch_l_rayo=ref_cast(Champ_Uniforme,l_rayo.valeur());
              ch_l_rayo.nommer("longueur_de_rayonnement");
              ch_l_rayo.fixer_nb_comp(1);
              // Le nombre de valeurs nodales est fixe a 1 ici car il s'agit d'un
              // Champ_Uniforme, dans les autres cas, il faut le fixer egale au nombre
              // d'elements ou de faces en fonctions de la localisation du champ
              ch_l_rayo.fixer_nb_valeurs_nodales(1);
              ch_l_rayo.fixer_unite("m");
              ch_l_rayo.changer_temps(0);
            }
          else
            {
              Cerr <<"le coefficient d'absorption n'est pas un OWN_PTR(Champ_base) Uniforme "<<finl;
              Cerr <<"mais un "<<coeff_abs.que_suis_je()<<". modifier la methode "<<finl;
              Cerr <<"Pb_Couple_rayo_semi_transp::discretiser pour pouvoir prendre"<<finl;
              Cerr <<"en compte ce type de Champ_Don"<<finl;
            }
          const double temps = probleme().schema_temps().temps_courant();
          fluide.initialiser(temps);
        }
      else
        {
          Cerr<<"Erreur 0 dans Pb_Couple_rayo_semi_transp::discretiser"<<finl;
          Cerr<<"vous n'avez probablement pas renseigne tous les "<<finl;
          Cerr<<"parametres physique de votre fluide incompressible "<<finl;
          Cerr<<"pour pouvoir traiter un probleme de rayonnement semi transparent"<<finl;
          exit();
        }
    }

  else
    {
      Cerr<<"Erreur dans Modele_rayo_semi_transp::readOn "<<finl;
      Cerr<<"Le modele de rayonnement semi transparent ne peut etre utilise"<<finl;
      Cerr<<"qu'avec un Fluide_base et non "<< probleme().milieu().que_suis_je()<<finl;
      exit();
    }

  Cerr << "Modele_rayo_semi_transp::discretiser fin" << finl;
}

// On met a jour le flux radiatif pour tous les bords du probleme
void Modele_rayo_semi_transp::calculer_flux_radiatif()
{
  Conds_lim& les_cl_rayo = eq_rayo().domaine_Cl_dis().les_conditions_limites();

  int num_cl_rayo=0;
  for(num_cl_rayo=0; num_cl_rayo<les_cl_rayo.size(); num_cl_rayo++)
    {
      Cond_lim& la_cl_rayo = eq_rayo().domaine_Cl_dis().les_conditions_limites(num_cl_rayo);
      if(sub_type(Flux_radiatif_base,la_cl_rayo.valeur()))
        {
          Flux_radiatif_base& la_cl_rayon = ref_cast(Flux_radiatif_base,la_cl_rayo.valeur());
          Equation_base& eq_temp = probleme().equation(1);
          la_cl_rayon.calculer_flux_radiatif(eq_temp);
        }
      else if (sub_type(Symetrie,la_cl_rayo.valeur()))
        {
          ;
        }
      else
        {
          Cerr<<"Erreur : les conditions aux limites de l'equation de rayonnement"<<finl;
          Cerr<<"doivent forcement etre du type rayonnantes"<<finl;
          exit();
        }
    }
}


////////////////////////////////////////////////////////////////////
// Description :
//
// Precondition :
// Parametre :
//     Signification :
//     Valeurs par defaut :
//     Contraintes :
//     Entree :
//     Entree/Sortie :
//     Sortie :
// Retour :
//     Signification :
//     Contraintes :
// Exception :
// Effets de bord :
// Postcondition :
//
const Champ_front_base& Modele_rayo_semi_transp::flux_radiatif(const Nom& nom_bord) const
{
  //  Cerr<<"Modele_rayo_semi_transp::flux_radiatif const : Debut"<<finl;
  // On fait une boucle sur les bords pour trouver celui dont le nom est nom_bord
  const Conds_lim& les_cl_rayo = eq_rayo().domaine_Cl_dis().les_conditions_limites();


  int num_cl_rayo=0;
  for(num_cl_rayo=0; num_cl_rayo<les_cl_rayo.size(); num_cl_rayo++)
    {
      const Cond_lim& la_cl_rayo = eq_rayo().domaine_Cl_dis().les_conditions_limites(num_cl_rayo);

      if(la_cl_rayo->frontiere_dis().le_nom()==nom_bord)
        {
          if(sub_type(Flux_radiatif_base,la_cl_rayo.valeur()))
            {
              Flux_radiatif_base& la_cl_rayon = ref_cast_non_const(Flux_radiatif_base,la_cl_rayo.valeur());
              return la_cl_rayon.flux_radiatif();
            }
          else
            {
              Cerr<<"Erreur : les conditions aux limites de l'equation de rayonnement"<<finl;
              Cerr<<"doivent forcement etre du type rayonnantes"<<finl;
              exit();

            }
        }
    }
  Cerr<<"Erreur : Modele_rayo_semi_transp::flux_radiatif"<<finl;
  Cerr<<"il n'y a pas de condition a la limite portant le nom "<<nom_bord<<finl;
  exit();
  //pour les compilos
  const Cond_lim& la_cl_rayo = eq_rayo().domaine_Cl_dis().les_conditions_limites(0);
  Flux_radiatif_base& la_cl_rayon = ref_cast_non_const(Flux_radiatif_base,la_cl_rayo.valeur());
  return la_cl_rayon.flux_radiatif();
  //  Cerr<<"Modele_rayo_semi_transp::flux_radiatif const : Fin"<<finl;
}


