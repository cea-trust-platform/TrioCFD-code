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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Modele_turbulence_scal_Fluctuation_Temperature.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/Common/Scal
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_scal_Fluctuation_Temperature.h>
#include <Probleme_base.h>
#include <Mod_turb_hyd_base.h>
#include <DoubleTrav.h>
#include <Param.h>

Implemente_instanciable(Modele_turbulence_scal_Fluctuation_Temperature,"Modele_turbulence_scal_Fluctuation_Temperature",Modele_turbulence_scal_base);

//// printOn
//

Sortie& Modele_turbulence_scal_Fluctuation_Temperature::printOn(Sortie& s ) const
{
  return s << que_suis_je() << " " << le_nom();
}


//// readOn
//

Entree& Modele_turbulence_scal_Fluctuation_Temperature::readOn(Entree& is )
{

  return Modele_turbulence_scal_base::readOn(is);
}

void Modele_turbulence_scal_Fluctuation_Temperature::set_param(Param& param)
{
  param.ajouter_non_std("Transport_Fluctuation_Temperature",this);
  param.ajouter_non_std("Transport_Flux_Chaleur_Turbulente",this);
  // on ajoute pas les params de la classe mere car ce n'etait pas fait avt
}
int Modele_turbulence_scal_Fluctuation_Temperature::lire_motcle_non_standard(const Motcle& mot, Entree& s)
{
  Cerr << "Lecture des parametres du modele de fluctuation thermique. Il doit y avoir deux types d'equation." << finl;
  Motcle accouverte = "{" , accfermee = "}" ;
  Motcles les_mots(2);
  {
    les_mots[0] = "Transport_Fluctuation_Temperature";
    les_mots[1] = "Transport_Flux_Chaleur_Turbulente";
  }
  int rang=les_mots.search(mot);
  switch(rang)
    {
    case 0:
      {
        Cerr << "Lecture de l'equation Transport_Fluctuation_Temperature" << finl;
        eqn_transport_Fluctu_Temp.associer_modele_turbulence(*this);
        s >> eqn_transport_Fluctu_Temp;
        break;
      }
    case 1:
      {
        Cerr << "Lecture de l'equation Transport_Flux_Chaleur_Turbulente" << finl;
        eqn_transport_Flux_Chaleur_Turb.associer_modele_turbulence(*this);
        s >> eqn_transport_Flux_Chaleur_Turb;
        break;
      }
    default :
      {
        return Modele_turbulence_scal_base::lire_motcle_non_standard(mot,s);
        break;
      }
    }
  return 1;
}

void Modele_turbulence_scal_Fluctuation_Temperature::associer_viscosite_turbulente(const Champ_Fonc& visc_turb)
{
  la_viscosite_turbulente = visc_turb;
}



int Modele_turbulence_scal_Fluctuation_Temperature::preparer_calcul()
{
  eqn_transport_Fluctu_Temp.preparer_calcul();
  eqn_transport_Flux_Chaleur_Turb.preparer_calcul();
  Modele_turbulence_scal_base::preparer_calcul();
  mettre_a_jour(eqn_transport_Fluctu_Temp.schema_temps().temps_courant());
  return 1;
}

bool Modele_turbulence_scal_Fluctuation_Temperature::initTimeStep(double dt)
{
  bool ok=eqn_transport_Fluctu_Temp.initTimeStep(dt);
  ok = ok && eqn_transport_Flux_Chaleur_Turb.initTimeStep(dt);
  return ok;
}

Champ_Fonc& Modele_turbulence_scal_Fluctuation_Temperature::calculer_diffusivite_turbulente()
{
  DoubleTab& alpha_t = diffusivite_turbulente_.valeurs();
  const DoubleTab& nu_t = la_viscosite_turbulente->valeurs();
  double temps = la_viscosite_turbulente->temps();

  if (temps != diffusivite_turbulente_.temps())
    {
      static const double Prdt_turbulent = 0.9;

      int n= alpha_t.size();
      if (nu_t.size() != n)
        {
          Cerr << "Les DoubleTab des champs diffusivite_turbulente et viscosite_turbulente" << finl;
          Cerr << "doivent avoir le meme nombre de valeurs nodales" << finl;
          exit();
        }

      for (int i=0; i<n; i++)
        alpha_t[i] = nu_t[i]/Prdt_turbulent;

      diffusivite_turbulente_.changer_temps(temps);
    }
  return diffusivite_turbulente_;
}

void Modele_turbulence_scal_Fluctuation_Temperature::mettre_a_jour(double temps)
{
  //  Champ_Inc& ch_Fluctu_Temp = Fluctu_Temperature();
  Schema_Temps_base& sch1 = eqn_transport_Fluctu_Temp.schema_temps();
  // Voir Schema_Temps_base::faire_un_pas_de_temps_pb_base
  eqn_transport_Fluctu_Temp.zone_Cl_dis().mettre_a_jour(temps);
  sch1.faire_un_pas_de_temps_eqn_base(eqn_transport_Fluctu_Temp);
  //eqn_transport_Fluctu_Temp.inconnue().mettre_a_jour(temps);
  eqn_transport_Fluctu_Temp.mettre_a_jour(temps);
  eqn_transport_Fluctu_Temp.controler_grandeur();
  //  Champ_Inc& ch_Flux_Chaleur_Turb = Flux_Chaleur_Turb();
  Schema_Temps_base& sch2 = eqn_transport_Flux_Chaleur_Turb.schema_temps();
  // Voir Schema_Temps_base::faire_un_pas_de_temps_pb_base
  eqn_transport_Flux_Chaleur_Turb.zone_Cl_dis().mettre_a_jour(temps);
  sch2.faire_un_pas_de_temps_eqn_base(eqn_transport_Flux_Chaleur_Turb);
  //eqn_transport_Flux_Chaleur_Turb.inconnue().mettre_a_jour(temps);
  eqn_transport_Flux_Chaleur_Turb.mettre_a_jour(temps);
  eqn_transport_Flux_Chaleur_Turb.controler_grandeur();
  calculer_diffusivite_turbulente();
}

void Modele_turbulence_scal_Fluctuation_Temperature::completer()
{
  eqn_transport_Fluctu_Temp.completer();
  eqn_transport_Flux_Chaleur_Turb.completer();
  const Probleme_base& mon_pb = mon_equation->probleme();
  const RefObjU& modele_turbulence = mon_pb.equation(0).get_modele(TURBULENCE);
  const Mod_turb_hyd_base& mod_turb_hydr = ref_cast(Mod_turb_hyd_base,modele_turbulence.valeur());
  const Champ_Fonc& visc_turb = mod_turb_hydr.viscosite_turbulente();
  associer_viscosite_turbulente(visc_turb);
}

int Modele_turbulence_scal_Fluctuation_Temperature::sauvegarder(Sortie& os) const
{
  Modele_turbulence_scal_base::a_faire(os);
  int bytes=0;
  bytes += eqn_transport_Fluctu_Temp.sauvegarder(os);
  bytes += eqn_transport_Flux_Chaleur_Turb.sauvegarder(os);
  return bytes;

}

int Modele_turbulence_scal_Fluctuation_Temperature::reprendre(Entree& is)
{
  Modele_turbulence_scal_base::reprendre(is);
  if (mon_equation.non_nul())
    {
      eqn_transport_Fluctu_Temp.reprendre(is);
      eqn_transport_Flux_Chaleur_Turb.reprendre(is);
      return 1;
    }
  else
    {
      double dbidon;
      Nom bidon;
      DoubleTrav tab_bidon;
      is >> bidon >> bidon;
      is >> dbidon;
      tab_bidon.jump(is);
      return 1;
    }
}
void Modele_turbulence_scal_Fluctuation_Temperature::imprimer(Sortie&) const
{
}

const Champ_base& Modele_turbulence_scal_Fluctuation_Temperature::get_champ(const Motcle& nom) const
{
  REF(Champ_base) ref_champ;
  try
    {
      return Modele_turbulence_scal_base::get_champ(nom);
    }
  catch (Champs_compris_erreur)
    {
    }

  try
    {
      return eqn_transport_Fluctu_Temp.get_champ(nom);
    }
  catch (Champs_compris_erreur)
    {
    }
  try
    {
      return eqn_transport_Flux_Chaleur_Turb.get_champ(nom);
    }
  catch (Champs_compris_erreur)
    {
    }

  throw Champs_compris_erreur();
  return ref_champ;
}

void Modele_turbulence_scal_Fluctuation_Temperature::get_noms_champs_postraitables(Noms& nom,Option opt) const
{
  Modele_turbulence_scal_base::get_noms_champs_postraitables(nom,opt);

  eqn_transport_Fluctu_Temp.get_noms_champs_postraitables(nom,opt);

  eqn_transport_Flux_Chaleur_Turb.get_noms_champs_postraitables(nom,opt);

}
