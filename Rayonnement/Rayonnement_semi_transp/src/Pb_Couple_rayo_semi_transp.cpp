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
// File:        Pb_Couple_rayo_semi_transp.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src
// Version:     /main/24
//
//////////////////////////////////////////////////////////////////////////////

#include <Pb_Couple_rayo_semi_transp.h>
#include <Schema_Temps_base.h>
#include <Interprete.h>
#include <Modele_rayo_semi_transp.h>
#include <Source_rayo_semi_transp_base.h>
#include <Fluide_base.h>
#include <Champ_Uniforme.h>
#include <verif_cast.h>
#include <Interprete_bloc.h>

Implemente_instanciable(Pb_Couple_rayo_semi_transp,"Pb_Couple_rayo_semi_transp",Probleme_Couple);


Entree& Pb_Couple_rayo_semi_transp::readOn(Entree& is)
{
  return is;
}

Sortie& Pb_Couple_rayo_semi_transp::printOn(Sortie& os) const
{
  return Probleme_Couple::printOn(os);
}

void Pb_Couple_rayo_semi_transp::initialize( )
{
  Probleme_Couple::initialize();
  Probleme_base& le_pb=modele().probleme();
  // Associer le modele aux sources de rayonnement
  for (int i=0; i<le_pb.nombre_d_equations(); i++)
    {
      Sources& les_sources=le_pb.equation(i).sources();
      for (int j=0; j<les_sources.size(); j++)
        {
          Source& la_source=les_sources[j];
          if (sub_type(Source_rayo_semi_transp_base,la_source.valeur()))   // premier cas
            {
              Source_rayo_semi_transp_base& source_rayo=ref_cast(Source_rayo_semi_transp_base,la_source.valeur());
              Cerr << "Association MODELE a SOURCE" << finl;
              source_rayo.associer_modele_rayo(modele());
            }
        }
    }

  modele().eq_rayo().resoudre(presentTime());
  modele().calculer_flux_radiatif();

  for (int i=0; i<nb_problemes(); i++)
    {
      Probleme_base& pb=ref_cast(Probleme_base,probleme(i));
      for (int j=0; j<pb.nombre_d_equations(); j++)
        {
          pb.equation(j).zone_Cl_dis()->calculer_coeffs_echange(presentTime());
        }
    }
}

// void Pb_Couple_rayo_semi_transp::terminate() {
//   Probleme_Couple::terminate();
//   modele().terminate();
// }

void Pb_Couple_rayo_semi_transp::associer_sch_tps_base(Schema_Temps_base& sch)
{
  Probleme_Couple::associer_sch_tps_base(sch);
  sch_clone=sch;
  if (!le_modele_.non_nul())
    {
      Cerr << "Attention, le modele de rayonnement semi transparent n'est pas encore defini." << finl;
      Cerr << "La definition des problemes a change depuis la 1.5.2. En particulier, le modele" << finl;
      Cerr << "de rayonnement (qui est devenu un probleme avec son propre postraitement) doit etre" << finl;
      Cerr << "associe au probleme couple avant l'association du schema en temps avec le probleme couple." << finl;
      Cerr << "Voir la documentation des problemes thermohydrauliques avec modele de rayonnement" << finl;
      Cerr << "ou contacter le support TRUST pour plus de precisions." << finl;
      exit();
    }
  modele().associer_sch_tps_base(sch_clone); // association
}

int Pb_Couple_rayo_semi_transp::associer_(Objet_U& ob)
{
  Cerr << "Appel a associer_ "  << ob.que_suis_je() << finl;
  if(sub_type(Modele_rayo_semi_transp, ob))
    {
      Cerr << "association du modele au pbc" << finl;
      le_modele_rayo_associe(ref_cast(Modele_rayo_semi_transp, ob));

      // Je rajoute le modele de rayonnement dans la liste des problemes
      int nb_groupes = groupes.size_array();
      if (nb_groupes == 0)
        {
          // on est passe par associer dans le readOn
          // donc les groupes ne sont pas encore crees
          groupes.resize_array(2);
          groupes[0] = nb_problemes();
          Probleme_Couple::ajouter(modele());
          groupes[1] = 1;
        }
      else
        {
          Probleme_Couple::ajouter(modele());

          ArrOfInt new_groupes(nb_groupes+1);
          for (int i=0; i<nb_groupes; ++i)
            {
              new_groupes[i] = groupes[i];
            }
          new_groupes[nb_groupes] = 1;
          groupes.ref_array(new_groupes);
        }

      return 1;
    }
  else if (Probleme_Couple::associer_(ob))
    {
      if (sub_type(Probleme_base,ob))
        {
          Probleme_base& pb=ref_cast(Probleme_base,ob);
          if (sub_type(Fluide_base,pb.milieu()))
            {
              Fluide_base& fluide=ref_cast(Fluide_base,pb.milieu());
              fluide.fixer_type_rayo();
            }
        }
      return 1;
    }
  else
    return 0;
}

void Pb_Couple_rayo_semi_transp::le_modele_rayo_associe(const Modele_rayo_semi_transp& un_modele_de_rayonnement)
{
  if (le_modele_.non_nul())
    {
      Cerr << "Attention : on ne peut associer qu'un modele de rayonnement a un Pb_Couple_rayo_semi_transp" << finl;
      exit();
    }
  le_modele_ = un_modele_de_rayonnement;
  int le_pb_a_associer=-1;
  // On associe au modele de rayonnement
  // le dernier probleme fluide du probleme couple.
  for(int l=0; l<nb_problemes(); l++)
    {
      Probleme_base& pb=ref_cast(Probleme_base,probleme(l));
      if (pb.milieu().is_rayo_semi_transp())
        le_pb_a_associer=l;
    }

  if (le_pb_a_associer==-1)
    {
      Cerr << "Attention : il n'y a aucun probleme fluide auquel associer le modele de rayonnement." << finl;
      exit();
    }
  Probleme_base& le_pb=ref_cast(Probleme_base,probleme(le_pb_a_associer));
  // Le probleme a associer est maintenant reference dans le_pb.

  // Associer le probleme au modele
  modele().associer_probleme(le_pb);

  // Clonage et association du domaine (WEC)
  // Deviendra inutile avec la version de gomtrie de B. Mathieu
  der_domaine_clone.typer("Zone");
  Zone& dom_clone=ref_cast(Zone,der_domaine_clone.valeur());
  dom_clone=le_pb.domaine();
  Nom new_name=dom_clone.le_nom()+"_copy";
  dom_clone.nommer(new_name); // nommage
  modele().associer_domaine(dom_clone); // association
  // Ici ajouter d'ventuelles autres associations ...
  Interprete_bloc::interprete_courant().ajouter(new_name,der_domaine_clone);

}

// bool Pb_Couple_rayo_semi_transp::initTimeStep(double dt) {
//   bool ok = Probleme_Couple::initTimeStep(dt);
//   ok = ok && modele().initTimeStep(dt);
//   return ok;
// }

// bool Pb_Couple_rayo_semi_transp::iterateTimeStep(bool& converged) {
//   bool cv;
//   bool ok = Probleme_Couple::iterateTimeStep(converged);
//   ok = ok && modele().iterateTimeStep(cv);
//   converged = converged && cv;
//   return ok;
// }

// void Pb_Couple_rayo_semi_transp::validateTimeStep() {
//   Probleme_Couple::validateTimeStep();
//   modele().validateTimeStep();
// }

