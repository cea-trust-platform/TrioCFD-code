/****************************************************************************
* Copyright (c) 2022, CEA
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

#include <Estimateur_Aposteriori_P0_VEF.h>
#include <Navier_Stokes_Aposteriori.h>
#include <Champ_P1_isoP1Bulle.h>
#include <Discretisation_base.h>
#include <Schema_Temps_base.h>
#include <Postraitements.h>
#include <Probleme_base.h>
#include <Champ_P1NC.h>
#include <EChaine.h>

Implemente_instanciable(Navier_Stokes_Aposteriori,"Navier_Stokes_Aposteriori",Navier_Stokes_std);
// XD Navier_Stokes_Aposteriori navier_stokes_standard Navier_Stokes_Aposteriori -1 Modification of the Navier_Stokes_standard class in order to accept the estimateur_aposteriori post-processing. To post-process estimateur_aposteriori, add this keyword into the list of fields to be post-processed. This estimator whill generate a map of aposteriori error estimators; it is defined on each mesh cell and is a measure of the local discretisation error. This will serve for adaptive mesh refinement

Sortie& Navier_Stokes_Aposteriori::printOn(Sortie& is) const { return Navier_Stokes_std::printOn(is); }
Entree& Navier_Stokes_Aposteriori::readOn(Entree& is) { return Navier_Stokes_std::readOn(is); }

void Navier_Stokes_Aposteriori::discretiser()
{
  if (discretisation().que_suis_je() != "VEFPreP1B")
    {
      Cerr << "Error : Pb_Hydraulique_Aposteriori is currently only available for the discretization VEFPreP1B !" << finl;
      Cerr << "Please update your data file ..." << finl;
      Process::exit();
    }
  Navier_Stokes_std::discretiser();
}

void Navier_Stokes_Aposteriori::creer_champ(const Motcle& motlu)
{
  Navier_Stokes_std::creer_champ(motlu);

  if (motlu == "estimateur_aposteriori")
    {
      if (!estimateur_aposteriori_.non_nul())
        {
          estimateur_aposteriori();
          champs_compris_.ajoute_champ(estimateur_aposteriori_);

          Nom chaine = "{ numero_source 0 sources { refChamp { pb_champ ";
          chaine += probleme().le_nom();
          chaine += "  vitesse } } }";
          EChaine echaine(chaine);
          echaine >> champ_src_;
          champ_src_.completer(probleme().postraitements().front());
        }
    }
}


void Navier_Stokes_Aposteriori::get_noms_champs_postraitables(Noms& nom,Option opt) const
{
  Navier_Stokes_std::get_noms_champs_postraitables(nom,opt);
  Noms noms_compris = champs_compris_.liste_noms_compris();
  noms_compris.add("estimateur_aposteriori");
  if (opt==DESCRIPTION)
    Cerr<<" Navier_Stokes_Aposteriori : "<< noms_compris <<finl;
  else
    nom.add(noms_compris);
}

const Champ_base& Navier_Stokes_Aposteriori::get_champ(const Motcle& nom) const
{
  const double temps_init = schema_temps().temps_init();

  if (nom == "estimateur_aposteriori")
    {
      Champ_Fonc_base& ch = ref_cast_non_const(Champ_Fonc_base, estimateur_aposteriori_.valeur());
      if (((ch.temps() != la_vitesse->temps()) || (ch.temps() == temps_init)) && (la_vitesse->mon_equation_non_nul()))
        ch.mettre_a_jour(la_vitesse->temps());

      return champs_compris_.get_champ(nom);
    }
  else
    return Navier_Stokes_std::get_champ(nom);
}

void Navier_Stokes_Aposteriori::estimateur_aposteriori()
{
  const Domaine_VEF& domaine_vef = ref_cast(Domaine_VEF, domaine_dis());
  const Domaine_Cl_VEF& domaine_cl_vef = ref_cast(Domaine_Cl_VEF, domaine_Cl_dis());
  estimateur_aposteriori_.typer("Estimateur_Aposteriori_P0_VEF");
  Estimateur_Aposteriori_P0_VEF& ch = ref_cast(Estimateur_Aposteriori_P0_VEF, estimateur_aposteriori_.valeur());
  ch.associer_domaine_dis_base(domaine_vef);
  const Champ_P1NC& vit = ref_cast(Champ_P1NC, la_vitesse.valeur());
  const Champ_P1_isoP1Bulle& pres = ref_cast(Champ_P1_isoP1Bulle, la_pression.valeur());
  ch.associer_champ(vit, pres, diffusivite_pour_transport() /* viscosite_cinematique */, domaine_cl_vef);
  ch.nommer("estimateur_aposteriori");
  ch.fixer_nb_comp(1);
  ch.fixer_nb_valeurs_nodales(domaine_vef.nb_elem());
  ch.fixer_unite("sans");
  ch.changer_temps(la_vitesse->temps());
}
