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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Navier_Stokes_Turbulent_ALE.cpp
// Directory : $ALE_ROOT/src/New
//
/////////////////////////////////////////////////////////////////////////////

#include <Navier_Stokes_Turbulent_ALE.h>
#include <Op_Diff_Turbulent_base.h>
#include <Probleme_base.h>
#include <Domaine.h>
#include <Champ_Uniforme.h>
#include <Fluide_Incompressible.h>
#include <Avanc.h>
#include <Modele_turbulence_hyd_null.h>
#include <Discretisation_base.h>
#include <Schema_Temps_base.h>
#include <Discret_Thyd.h>
#include <Modele_turbulence_hyd_K_Eps.h>
#include <Modele_turbulence_hyd_K_Eps_Realisable.h>
#include <Param.h>

Implemente_instanciable(Navier_Stokes_Turbulent_ALE,"Navier_Stokes_Turbulent_ALE",Navier_Stokes_std_ALE);
// XD Navier_Stokes_Turbulent_ALE Navier_Stokes_std_ALE Navier_Stokes_Turbulent_ALE -1 Resolution of hydraulic turbulent Navier-Stokes eq. on mobile domain (ALE)
// XD  attr modele_turbulence modele_turbulence_hyd_deriv modele_turbulence 1 Turbulence model for Navier-Stokes equations.


/*! @brief Impression de l'equation sur un flot de sortie.
 *
 * Simple appel a Equation_base::printOn(Sortie&).
 *
 * @param (Sortie& is) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Navier_Stokes_Turbulent_ALE::printOn(Sortie& is) const
{
  return Equation_base::printOn(is);
}


/*! @brief Lit les specifications de l'equation de Navier Stokes a partir d'un flot d'entree.
 *
 *     Simple appel a Navier_Stokes_std_ALE::readOn(Entree&)
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 * @throws pas de modele de turbulence speficie
 */
Entree& Navier_Stokes_Turbulent_ALE::readOn(Entree& is)
{
  Navier_Stokes_std_ALE::readOn(is);
  return is;
}

void Navier_Stokes_Turbulent_ALE::set_param(Param& param)
{
  Navier_Stokes_std_ALE::set_param(param);
  //param.ajouter_non_std("solveur_pression",(this),Param::REQUIRED);
  param.ajouter_non_std("modele_turbulence",(this),Param::REQUIRED);
}

int Navier_Stokes_Turbulent_ALE::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="diffusion")
    {
      Cerr << "Reading and typing of the diffusion operator : ";
      terme_diffusif.associer_diffusivite(diffusivite_pour_transport());
      // Si on a lu le modele de turbulence et qu'il est nul,
      // alors on utilise l'operateur de diffusion standard.
      if (le_modele_turbulence.non_nul() // L'operateur a ete type (donc lu)
          && sub_type(Modele_turbulence_hyd_null, le_modele_turbulence.valeur()))
        is >> terme_diffusif;
      else
        lire_op_diff_turbulent(is);
      terme_diffusif.associer_diffusivite_pour_pas_de_temps(diffusivite_pour_pas_de_temps());
      return 1;
    }
  else if (mot=="modele_turbulence")
    {
      Cerr << "Reading and typing of the turbulence model :" ;
      Motcle typ;
      is >> typ;
      Motcle nom1("Modele_turbulence_hyd_");
      nom1 += typ;
      Nom discr = discretisation().que_suis_je();

      if (typ.debute_par("SOUS_MAILLE") || discr == "VDF_Hyper" || typ.debute_par("LONGUEUR_MELANGE") || (typ == "K_Epsilon_V2"))
        {
          if (dimension == 2 && discr != "VDF_Hyper")
            {
              Cerr << "Vous traitez un cas turbulent en dimension 2 avec un modele sous maille" << finl;
              Cerr << "Attention a l'interpretation des resultats !!" << finl;
            }

          nom1 += "_";
          // les operateurs de diffusion sont communs aux discretisations VEF et VEFP1B
          if (discr == "VEFPreP1B") discr = "VEF";
          nom1 += discr;
        }
      if (nom1 == "MODELE_TURBULENCE_HYD_SOUS_MAILLE_LM_VEF")
        {
          Cerr << "Le mot cle Sous_maille_LM s'appelle desormais Longueur_Melange pour etre coherent en VDF et VEF." << finl;
          Cerr << "Changer votre jeu de donnees." << finl;
          Process::exit();
        }

      Cerr << nom1 << finl;

      le_modele_turbulence.typer(nom1);
      le_modele_turbulence->associer_eqn(*this);
      le_modele_turbulence->associer(domaine_dis(), domaine_Cl_dis());
      is >> le_modele_turbulence.valeur(); // on lit !

      le_modele_turbulence->discretiser();
      RefObjU le_modele;
      le_modele = le_modele_turbulence.valeur();
      liste_modeles_.add_if_not(le_modele);

      return 1;
    }
  else
    return Navier_Stokes_std_ALE::lire_motcle_non_standard(mot,is);
}

const Champ_Don_base& Navier_Stokes_Turbulent_ALE::diffusivite_pour_transport() const
{
  return fluide().viscosite_cinematique();
}

const Champ_base& Navier_Stokes_Turbulent_ALE::diffusivite_pour_pas_de_temps() const
{
  return fluide().viscosite_cinematique();
}

// Lecture et typage de l'operateur diffusion turbulente.
// Attention : il faut avoir fait "terme_diffusif.associer_diffusivite" avant
//             d'enter ici.
Entree& Navier_Stokes_Turbulent_ALE::lire_op_diff_turbulent(Entree& is)
{
  Motcle accouverte = "{" , accfermee = "}" ;
  Nom type="Op_Dift_";

  Nom discr=discretisation().que_suis_je();
  // les operateurs de diffusion sont communs aux discretisations VEF et VEFP1B
  if(discr=="VEFPreP1B") discr="VEF";
  type+=discr;

  Nom nb_inc;
  if (terme_diffusif.diffusivite().nb_comp() == 1)
    nb_inc = "_";
  else
    nb_inc = "_Multi_inco_";
  type+= nb_inc ;

  Nom type_inco=inconnue().que_suis_je();
  type+=(type_inco.suffix("Champ_"));

  if (axi)
    type+="_Axi";

  Motcle motbidon;
  is >>  motbidon;
  if (motbidon!=accouverte)
    {
      Cerr << "A { was expected while reading" << finl;
      Cerr << "the turbulent diffusive term" << finl;
      exit();
    }
  is >>  motbidon;
  if (motbidon=="negligeable")
    {
      type="Op_Dift_negligeable";
      terme_diffusif.typer(type);
      terme_diffusif.l_op_base().associer_eqn(*this);
      terme_diffusif->associer_diffusivite(terme_diffusif.diffusivite());
      is >> motbidon;
      //on lit la fin de diffusion { }
      if ( motbidon != accfermee)
        Cerr << " On ne peut plus entrer d option apres negligeable "<< finl;
    }
  else if (motbidon=="standard")
    {
      type+="_standard";
      terme_diffusif.typer(type);
      terme_diffusif.l_op_base().associer_eqn(*this);
      terme_diffusif->associer_diffusivite(terme_diffusif.diffusivite());
      is>>terme_diffusif.valeur();
    }
  else if (motbidon==accfermee)
    {
      terme_diffusif.typer(type);
      terme_diffusif.l_op_base().associer_eqn(*this);
      Cerr << terme_diffusif->que_suis_je() << finl;
      terme_diffusif->associer_diffusivite(terme_diffusif.diffusivite());
    }
  else if (motbidon=="stab")
    {
      type+="_stab";
      terme_diffusif.typer(type);
      terme_diffusif.l_op_base().associer_eqn(*this);
      terme_diffusif->associer_diffusivite(terme_diffusif.diffusivite());
      is>>terme_diffusif.valeur();
    }
  else
    {
      type += motbidon;
      is >> motbidon;
      if ( motbidon != accfermee)
        {
          Cerr << " No option are now readable "<< finl;
          Cerr << " for the turbulent diffusive term" << finl;

        }
      if ( discr == "VEF" )
        {
          type += "_const_P1NC" ;
        }
      terme_diffusif.typer(type);
      terme_diffusif.l_op_base().associer_eqn(*this);
      Cerr << terme_diffusif->que_suis_je() << finl;
      terme_diffusif->associer_diffusivite(terme_diffusif.diffusivite());
    }
  return is;
}

/*! @brief Prepare le calcul.
 *
 * Simple appe a Modele_turbulence_hyd_base::preparer_caclul() sur
 *     le membre reprresentant la turbulence.
 *
 * @return (int) renvoie toujours 1
 */
int Navier_Stokes_Turbulent_ALE::preparer_calcul()
{

  Turbulence_paroi_base& loipar=le_modele_turbulence->loi_paroi();
  if (le_modele_turbulence->has_loi_paroi_hyd())
    loipar.init_lois_paroi();

  Navier_Stokes_std_ALE::preparer_calcul();
  le_modele_turbulence->preparer_calcul();
  return 1;
}

bool Navier_Stokes_Turbulent_ALE::initTimeStep(double dt)
{
  bool ok=Navier_Stokes_std_ALE::initTimeStep(dt);
  ok = ok && le_modele_turbulence->initTimeStep(dt);
  return ok;
}



/*! @brief Sauvegarde l'equation (et son modele de turbulence) sur un flot de sortie.
 *
 * @param (Sortie& os) un flot de sortie
 * @return (int) renvoie toujours 1
 */
int Navier_Stokes_Turbulent_ALE::sauvegarder(Sortie& os) const
{
  int bytes=0;
  bytes += Navier_Stokes_std_ALE::sauvegarder(os);
  assert(bytes % 4 == 0);
  bytes += le_modele_turbulence->sauvegarder(os);
  assert(bytes % 4 == 0);
  return bytes;
}


/*! @brief Reprise de l'equation et de son modele de turbulence a partir d'un flot d'entree.
 *
 * @param (Entree& is) un flot d'entree
 * @return (int) renvoie toujours 1
 * @throws fin de fichier rencontre pendant la reprise
 */
int Navier_Stokes_Turbulent_ALE::reprendre(Entree& is)
{
  Navier_Stokes_std_ALE::reprendre(is);
  double temps = schema_temps().temps_courant();
  Nom ident_modele(le_modele_turbulence->que_suis_je());
  ident_modele += probleme().domaine().le_nom();
  ident_modele += Nom(temps,probleme().reprise_format_temps());

  avancer_fichier(is, ident_modele);
  le_modele_turbulence->reprendre(is);

  return 1;
}


/*! @brief Appels successifs a: Navier_Stokes_std_ALE::completer()
 *
 *       Mod_Turb_Hyd::completer() [sur le membre concerne]
 *
 */
void Navier_Stokes_Turbulent_ALE::completer()
{
  Navier_Stokes_std_ALE::completer();
  le_modele_turbulence->completer();
}


/*! @brief Effecttue une mise a jour en temps de l'equation.
 *
 * @param (double temps) le temps de mise a jour
 */
void Navier_Stokes_Turbulent_ALE::mettre_a_jour(double temps)
{
  Navier_Stokes_std_ALE::mettre_a_jour(temps);
  le_modele_turbulence->mettre_a_jour(temps);
}

void Navier_Stokes_Turbulent_ALE::creer_champ(const Motcle& motlu)
{
  Navier_Stokes_std_ALE::creer_champ(motlu);

  if (le_modele_turbulence.non_nul())
    le_modele_turbulence->creer_champ(motlu);
}

bool Navier_Stokes_Turbulent_ALE::has_champ(const Motcle& nom, OBS_PTR(Champ_base)& ref_champ) const
{
  if (Navier_Stokes_std_ALE::has_champ(nom))
    return Navier_Stokes_std_ALE::has_champ(nom, ref_champ);

  if (le_modele_turbulence.non_nul())
    if (le_modele_turbulence->has_champ(nom))
      return le_modele_turbulence->has_champ(nom, ref_champ);

  return false; /* rien trouve */
}

bool Navier_Stokes_Turbulent_ALE::has_champ(const Motcle& nom) const
{
  if (Navier_Stokes_std_ALE::has_champ(nom))
    return true;

  if (le_modele_turbulence.non_nul())
    if (le_modele_turbulence->has_champ(nom))
      return true;

  return false; /* rien trouve */
}

const Champ_base& Navier_Stokes_Turbulent_ALE::get_champ(const Motcle& nom) const
{
  if (Navier_Stokes_std_ALE::has_champ(nom))
    return Navier_Stokes_std_ALE::get_champ(nom);

  if (le_modele_turbulence.non_nul())
    if (le_modele_turbulence->has_champ(nom))
      return le_modele_turbulence->get_champ(nom);

  throw std::runtime_error("Field not found !");
}

void Navier_Stokes_Turbulent_ALE::get_noms_champs_postraitables(Noms& nom,Option opt) const
{
  Navier_Stokes_std_ALE::get_noms_champs_postraitables(nom,opt);
  if (le_modele_turbulence.non_nul())
    le_modele_turbulence->get_noms_champs_postraitables(nom,opt);
}

void Navier_Stokes_Turbulent_ALE::imprimer(Sortie& os) const
{
  Navier_Stokes_std_ALE::imprimer(os);
  le_modele_turbulence->imprimer(os);
}

const RefObjU& Navier_Stokes_Turbulent_ALE::get_modele(Type_modele type) const
{
  for(const auto& itr : liste_modeles_)
    {
      const RefObjU&  mod = itr;
      if (mod.non_nul())
        if ((sub_type(Modele_turbulence_hyd_base,mod.valeur())) && (type==TURBULENCE))
          return mod;
    }
  return Equation_base::get_modele(type);
}

// Impression du residu dans fic (generalement dt_ev)
// Cette methode peut etre surchargee par des equations
// imprimant des residus particuliers (K-Eps, Concentrations,...)
void Navier_Stokes_Turbulent_ALE::imprime_residu(SFichier& fic)
{
  Equation_base::imprime_residu(fic);
  // Si K-Eps, on imprime
  if (sub_type(Modele_turbulence_hyd_K_Eps,le_modele_turbulence.valeur()))
    ref_cast(Modele_turbulence_hyd_K_Eps,le_modele_turbulence.valeur()).eqn_transp_K_Eps().imprime_residu(fic);
  else if (sub_type(Modele_turbulence_hyd_K_Eps_Realisable,le_modele_turbulence.valeur()))
    ref_cast(Modele_turbulence_hyd_K_Eps_Realisable,le_modele_turbulence.valeur()).eqn_transp_K_Eps().imprime_residu(fic);
}

// Retourne l'expression du residu (de meme peut etre surcharge)
Nom Navier_Stokes_Turbulent_ALE::expression_residu()
{
  Nom tmp=Equation_base::expression_residu();
  if (sub_type(Modele_turbulence_hyd_K_Eps,le_modele_turbulence.valeur()))
    tmp+=ref_cast(Modele_turbulence_hyd_K_Eps,le_modele_turbulence.valeur()).eqn_transp_K_Eps().expression_residu();
  else if (sub_type(Modele_turbulence_hyd_K_Eps_Realisable,le_modele_turbulence.valeur()))
    tmp+=ref_cast(Modele_turbulence_hyd_K_Eps_Realisable,le_modele_turbulence.valeur()).eqn_transp_K_Eps().expression_residu();
  return tmp;
}

