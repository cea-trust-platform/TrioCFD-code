/****************************************************************************
* Copyright (c) 2024, CEA
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
// File:        Navier_Stokes_FT_Disc.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/55
//
//////////////////////////////////////////////////////////////////////////////

#include <Parametre_implicite.h>
#include <Navier_Stokes_FT_Disc.h>
#include <Probleme_FT_Disc_gen.h>
#include <Modele_turbulence_hyd_null.h>
#include <Discret_Thyd.h>
#include <Operateur_Diff_base.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <Fluide_Diphasique.h>
#include <Fluide_Incompressible.h>
#include <Assembleur_base.h>
#include <TRUST_Vector.h>
#include <Schema_Temps_base.h>
#include <Domaine_VDF.h>
#include <Debog.h>
#include <Domaine_VF.h>
#include <Terme_Source_Constituant_Vortex_VEF_Face.h>
#include <TRUSTTrav.h>
#include <Matrice_Morse_Sym.h>
#include <Matrice_Bloc.h>
#include <Param.h>
#include <Vecteur3.h>
#include <Domaine_Cl_VDF.h>

#define NS_VERBOSE 0 // To activate verbose mode on err ...

Implemente_instanciable(Navier_Stokes_FT_Disc, "Navier_Stokes_FT_Disc", Navier_Stokes_Turbulent);

/*! @brief Calcul de champ_rho_faces_ et champ_rho_elem_ en fonction de l'indicatrice: rho_elem_ = indicatrice * ( rho(phase_1) - rho(phase_0) ) + rho(phase_0)
 *
 *    rho_faces_ = 0.5 * (rho_elem(voisin_0) + rho_elem(voisin_1))
 *   Le calcul des viscosites cinematique et dynamique est pour l'instant le suivant:
 *    nu_elem = indicatrice * ( nu(phase_1) - nu(phase_0) ) + nu(phase_0)
 *    mu_elem = nu_elem * rho_elem
 *
 */
static void FT_disc_calculer_champs_rho_mu_nu_dipha(const Domaine_dis_base& domaine_dis_base, const Fluide_Diphasique& fluide, const DoubleVect& indicatrice_elem, DoubleVect& rho_elem,
                                                    DoubleVect& nu_elem, DoubleVect& mu_elem, DoubleVect& rho_faces)
{
  const Fluide_Incompressible& phase_0 = fluide.fluide_phase(0);
  const Fluide_Incompressible& phase_1 = fluide.fluide_phase(1);
  const DoubleTab& tab_rho_phase_0 = phase_0.masse_volumique()->valeurs();
  const DoubleTab& tab_rho_phase_1 = phase_1.masse_volumique()->valeurs();
  const double rho_phase_0 = tab_rho_phase_0(0, 0);
  const double rho_phase_1 = tab_rho_phase_1(0, 0);
  const double delta_rho = rho_phase_1 - rho_phase_0;
  const DoubleTab& tab_nu_phase_0 = phase_0.viscosite_cinematique()->valeurs();
  const DoubleTab& tab_nu_phase_1 = phase_1.viscosite_cinematique()->valeurs();
  const double nu_phase_0 = tab_nu_phase_0(0, 0);
  const double nu_phase_1 = tab_nu_phase_1(0, 0);
  const double delta_nu = nu_phase_1 - nu_phase_0;
  const double mu_phase_0 = nu_phase_0 * rho_phase_0;
  const double mu_phase_1 = nu_phase_1 * rho_phase_1;
  const double delta_mu = mu_phase_1 - mu_phase_0;
  double mu = 0.;
  const int formule_mu = fluide.formule_mu();

  // Calcul de rho, nu, mu aux elements
  {
    const int n = indicatrice_elem.size();
    assert(n == rho_elem.size());
    assert(n == nu_elem.size());
    assert(n == mu_elem.size());
    for (int i = 0; i < n; i++)
      {
        const double indic = indicatrice_elem[i];
        const double rho = indic * delta_rho + rho_phase_0;
        rho_elem[i] = rho;
        double nu = indic * delta_nu + nu_phase_0;
        switch(formule_mu)
          {
          case 0: // standard default method (will be removed)
            {
              mu = nu * rho;
            }
            break;
          case 1: // Arithmetic average
            {
              mu = indic * delta_mu + mu_phase_0;
            }
            break;
          case 2: // Harmonic average
            {
              mu = (mu_phase_0 * mu_phase_1) / (mu_phase_1 - indic * delta_mu);
            }
            break;
          default:
            {
              Cerr << "The method specified for formule_mu in not recognized. \n" << finl;
              Cerr << "you can choose : standard, arithmetic or harmonic. \n" << finl;
              Process::exit();
            }
          }
        mu_elem[i] = mu;
        nu = mu / rho;
        nu_elem[i] = nu;
      }
    rho_elem.echange_espace_virtuel();
    mu_elem.echange_espace_virtuel();
    nu_elem.echange_espace_virtuel();
  }

  // Calcul de rho aux faces (on suppose que la vitesse est aux faces)
  {
    assert(rho_elem.size() == domaine_dis_base.nb_elem());
    const IntTab& face_voisins = domaine_dis_base.face_voisins();
    const int nfaces = face_voisins.dimension(0);
    assert(rho_faces.size() == nfaces);
    for (int i = 0; i < nfaces; i++)
      {
        const int elem0 = face_voisins(i, 0);
        const int elem1 = face_voisins(i, 1);
        double rho = 0.;
        if (elem0 >= 0)
          rho = rho_elem[elem0];
        if (elem1 >= 0)
          rho += rho_elem[elem1];
        if (elem0 >= 0 && elem1 >= 0)
          rho *= 0.5;
        rho_faces[i] = rho;
      }
    rho_faces.echange_espace_virtuel();
  }
}

static void FT_disc_calculer_champs_rho_mu_nu_mono(const Domaine_dis_base& zdis, const Fluide_Incompressible& fluide, Champ_Fonc& champ_rho_elem_, Champ_Don& champ_mu_, Champ_Don& champ_nu_,
                                                   Champ_Fonc& champ_rho_faces_)
{

  if (sub_type(Champ_Uniforme,champ_rho_elem_.valeur()) && (sub_type(Champ_Uniforme, champ_nu_.valeur())))
    {
      const DoubleTab& tab_rho_phase_0 = fluide.masse_volumique()->valeurs();
      const double rho = tab_rho_phase_0(0, 0);
      const DoubleTab& tab_nu_phase_0 = fluide.viscosite_cinematique()->valeurs();
      const double nu = tab_nu_phase_0(0, 0);
      const double mu = nu * rho;
      champ_rho_elem_->valeurs() = rho;
      champ_nu_->valeurs() = nu;
      champ_mu_->valeurs() = mu;
      champ_rho_faces_->valeurs() = rho;
    }
  else
    {
      const int nb_el = zdis.domaine().nb_elem();
      const Domaine_VF& zvf = ref_cast(Domaine_VF, zdis);
      const IntTab& face_vois = zvf.face_voisins();
      const DoubleVect& vol = zvf.volumes();
      const int nb_faces = zvf.nb_faces();
      int elem1, elem2;
      double volume;
      //

      IntVect les_polys(nb_el);
      for (int i = 0; i < nb_el; i++)
        les_polys(i) = i;

      const DoubleTab& cg = zvf.xp();
      DoubleTab& val_rho = champ_rho_elem_->valeurs();
      DoubleTab& val_nu = champ_nu_->valeurs();
      DoubleTab& val_mu = champ_mu_->valeurs();
      DoubleTab& val_rho_faces = champ_rho_faces_->valeurs();

      fluide.masse_volumique()->valeur_aux_elems(cg, les_polys, val_rho);
      fluide.viscosite_dynamique()->valeur_aux_elems(cg, les_polys, val_mu);
      val_rho.echange_espace_virtuel();
      val_mu.echange_espace_virtuel();

      for (int el = 0; el < nb_el; el++)
        val_nu(el) = val_mu(el) / val_rho(el);
      val_nu.echange_espace_virtuel();

      for (int fac = 0; fac < nb_faces; fac++)
        {
          elem1 = face_vois(fac, 0);
          elem2 = face_vois(fac, 1);
          val_rho_faces(fac) = 0.;
          volume = 0.;
          if (elem1 != -1)
            {
              val_rho_faces(fac) += val_rho(elem1) * vol(elem1);
              volume += vol(elem1);
            }
          if (elem2 != -1)
            {
              val_rho_faces(fac) += val_rho(elem2) * vol(elem2);
              volume += vol(elem2);
            }
          val_rho_faces(fac) /= volume;
        }

      val_rho_faces.echange_espace_virtuel();
    }
}

Sortie& Navier_Stokes_FT_Disc::printOn(Sortie& os) const
{
  return Navier_Stokes_std::printOn(os);
}

Entree& Navier_Stokes_FT_Disc::readOn(Entree& is)
{
  terme_diffusif.associer_diffusivite(champ_mu_);
  return Navier_Stokes_Turbulent::readOn(is);
}

void Navier_Stokes_FT_Disc::set_param(Param& param)
{
  Navier_Stokes_Turbulent::set_param(param);
  param.ajouter_flag("boussinesq_approximation", &variables_internes().is_boussinesq_);
  param.ajouter_flag("new_mass_source", &variables_internes().new_mass_source_);
  param.ajouter_flag("matrice_pression_invariante", &variables_internes().matrice_pression_invariante);
  param.ajouter_non_std("equation_interfaces_vitesse_imposee", (this));
  param.ajouter_non_std("equations_interfaces_vitesse_imposee", (this));
  param.ajouter_non_std("equation_interfaces_proprietes_fluide", (this));
  param.ajouter("clipping_courbure_interface", &variables_internes().clipping_courbure_interface);
  param.ajouter_non_std("equation_temperature_mpoint", (this));
  param.ajouter_non_std("equation_temperature_mpoint_vapeur", (this));
  param.ajouter_non_std("terme_gravite", (this));
  param.ajouter("equations_concentration_source_vortex", &variables_internes().equations_concentration_source_fluide_);
  param.ajouter_non_std("repulsion_aux_bords", (this));
  param.ajouter_non_std("penalisation_forcage", (this));
  param.ajouter_flag("mpoint_inactif_sur_qdm", &variables_internes().mpoint_inactif);
  param.ajouter_flag("mpoint_vapeur_inactif_sur_qdm", &variables_internes().mpointv_inactif);
  param.ajouter("correction_courbure_ordre", &variables_internes().correction_courbure_ordre_);
  param.ajouter_non_std("interpol_indic_pour_dI_dt", (this));
  param.ajouter_non_std("OutletCorrection_pour_dI_dt", (this));
}

int Navier_Stokes_FT_Disc::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot == "diffusion")
    {
      terme_diffusif.associer_diffusivite(diffusivite_pour_transport());
      Cerr << "Reading and typing of the diffusion operator : " << finl;
      // Si on a lu le modele de turbulence et qu'il est nul,
      // alors on utilise l'operateur de diffusion standard.
      if (le_modele_turbulence.non_nul() // L'operateur a ete type (donc lu)
          && sub_type(Modele_turbulence_hyd_null, le_modele_turbulence.valeur()))
        {
          is >> terme_diffusif; // Operateur de diffusion standard (non turbulent)
          if (Process::je_suis_maitre())
            {
              Cerr << " WARNING : standard diffusion operator used for front_tracking\n";
              Cerr << "           the transposed term grad(v) is missing !!!" << finl;
            }
        }
      else
        lire_op_diff_turbulent(is);

      // Le coefficient de diffusion est une viscosite dynamique.
      // Il faut le diviser par rho pour calculer le pas de temps de stabilite.
      terme_diffusif->associer_champ_masse_volumique(champ_rho_elem_.valeur());

      // Pareil avec le modele de turbulence
      if (le_modele_turbulence.non_nul())
        le_modele_turbulence->associer_champ_masse_volumique(champ_rho_elem_.valeur());
      return 1;
    }
  else if ((mot == "equation_interfaces_vitesse_imposee") || (mot == "equation_interfaces_proprietes_fluide"))
    {
      Motcle nom_equation;
      is >> nom_equation;
      Cerr << mot << " equation : " << nom_equation << finl;
      const Transport_Interfaces_FT_Disc& eq = probleme_ft().equation_interfaces(nom_equation);

      if (mot == "equation_interfaces_vitesse_imposee")
        {
          if (variables_internes().ref_eq_interf_vitesse_imposee.size() > 0)
            {
              Cerr << "Error: You have already a equation_interfaces_vitesse_imposee defined." << finl;
              Cerr << "Since 1.6.4 TRUST version, you need to use the following syntax when" << finl;
              Cerr << "using several equation_interfaces_vitesse_imposee equations:" << finl;
              Cerr << "equations_interfaces_vitesse_imposee number_of_equations equation_name_one equation_name_two ..." << finl;
              exit();
            }
          else
            {
              Cerr << "===============================================================================" << finl;
              Cerr << "Warning: You are using a future obsolete syntax for defining a solid interface:" << finl;
              Cerr << "equation_interfaces_vitesse_imposee " << nom_equation << finl;
              Cerr << "Should be written:" << finl;
              Cerr << "equations_interfaces_vitesse_imposee 1 " << nom_equation << finl;
              Cerr << "===============================================================================" << finl;
            }
          variables_internes().ref_eq_interf_vitesse_imposee.add(eq);
        }
      else if (mot == "equation_interfaces_proprietes_fluide")
        variables_internes().ref_eq_interf_proprietes_fluide = eq;
      return 1;
    }
  else if (mot == "equations_interfaces_vitesse_imposee")
    {
      Noms na;
      is >> na;
      variables_internes().ref_eq_interf_vitesse_imposee.reset();
      for (int i = 0; i < na.size(); i++)
        {
          const Transport_Interfaces_FT_Disc& eq = probleme_ft().equation_interfaces(na[i]);
          variables_internes().ref_eq_interf_vitesse_imposee.add(eq);
        }
      return 1;
    }
  else if (mot == "equation_temperature_mpoint")
    {
      Nom nom_equation;
      is >> nom_equation;
      Cerr << " equation : " << nom_equation << finl;
      const Equation_base& eq = probleme_ft().get_equation_by_name(nom_equation);
      if (!sub_type(Convection_Diffusion_Temperature_FT_Disc, eq))
        {
          Cerr << " Error : equation is not of type Convection_Diffusion_Temperature_FT_Disc" << finl;
          exit();
        }
      variables_internes().ref_equation_mpoint_ = ref_cast(Convection_Diffusion_Temperature_FT_Disc, eq);
      return 1;
    }
  else if (mot == "equation_temperature_mpoint_vapeur")
    {
      Nom nom_equation;
      is >> nom_equation;
      Cerr << " equation : " << nom_equation << finl;
      const Equation_base& eq = probleme_ft().get_equation_by_name(nom_equation);
      if (!sub_type(Convection_Diffusion_Temperature_FT_Disc, eq))
        {
          Cerr << " Error : equation is not of type Convection_Diffusion_Temperature_FT_Disc" << finl;
          Process::exit();
        }
      variables_internes().ref_equation_mpoint_vap_ = ref_cast(Convection_Diffusion_Temperature_FT_Disc, eq);
      return 1;
    }
  else if (mot == "terme_gravite")
    {
      Motcles mots;
      mots.add("rho_g"); // Terme source volumique traditionnel avec courants parasites
      mots.add("grad_i"); // Terme source "front-tracking" sans courants parasites mais pbs en VEF.
      Motcle motbis;
      is >> motbis;
      Cerr << "Reading terme_gravite : " << motbis << finl;
      const int r = mots.search(motbis);
      switch(r)
        {
        case 0:
          variables_internes().terme_gravite_ = Navier_Stokes_FT_Disc_interne::GRAVITE_RHO_G;
          break;
        case 1:
          variables_internes().terme_gravite_ = Navier_Stokes_FT_Disc_interne::GRAVITE_GRAD_I;
          break;
        default:
          Cerr << "Error " << mots << "was expected whereas " << motbis << " has been found." << finl;
          barrier();
          exit();
        }
      return 1;
    }
  else if (mot == "repulsion_aux_bords")
    {
      is_repulsion = 1;
      is >> minx;
      is >> maxx;
      is >> pente;
      Cerr << "Interfaces repulsion on the boundaries for : " << minx << " " << maxx << finl;
      return 1;
    }
  else if (mot == "penalisation_forcage")
    {
      variables_internes().is_penalized = 1;
      variables_internes().is_explicite = 0;
      variables_internes().eta = 1.e-12;
      Cerr << "Navier_Stokes_FT_Disc : option penalisation_forcage" << finl;
      Motcles mots;
      mots.add("pression_reference"); // Valeur reference pour penalisation L2 pression
      mots.add("domaine_flottant_fluide"); // Traitement local Dirichlet pression si les CL pression sont toutes en Neumann homogene
      Motcle motbis;
      is >> motbis;
      Motcle accouverte = "{", accfermee = "}";
      if (motbis == accouverte)
        {
          is >> motbis;
          while (motbis != accfermee)
            {
              int rang = mots.search(motbis);
              switch(rang)
                {
                case 0:
                  {
                    is >> variables_internes().p_ref_pena;
                    Cerr << "Lecture de la pression de reference : " << variables_internes().p_ref_pena << finl;
                    break;
                  }
                case 1:
                  {
                    variables_internes().is_pfl_flottant = 1;
                    is >> variables_internes().x_pfl_imp;
                    is >> variables_internes().y_pfl_imp;
                    if (Objet_U::dimension == 3)
                      is >> variables_internes().z_pfl_imp;
                    Cerr << "Domaine flottant fluide" << finl;
                    Cerr << "Lecture de la position du point de reference pression fluide : " << variables_internes().x_pfl_imp << " " << variables_internes().y_pfl_imp << " "
                         << variables_internes().z_pfl_imp << finl;
                    break;
                  }
                default:
                  Cerr << "Erreur, on attendait " << mots << "On a trouve : " << motbis << finl;
                  barrier();
                  exit();
                }
              is >> motbis;
            }
        }
      else
        {
          Cerr << "Erreur a la lecture des parametres de la penalisation du forcage " << finl;
          Cerr << "On attendait : " << accouverte << finl;
          Process::exit();
        }
    }
  else if (mot == "interpol_indic_pour_dI_dt")
    {
      Motcles motcles2(9);
      motcles2[0] = "interp_standard";
      motcles2[1] = "interp_modifiee";
      motcles2[2] = "interp_ai_based";
      motcles2[3] = "interp_standard_uvext";
      motcles2[4] = "interp_modifiee_uvext";
      motcles2[5] = "interp_ai_based_uvext";
      motcles2[6] = "interp_standard_uiext";
      motcles2[7] = "interp_modifiee_uiext";
      motcles2[8] = "interp_ai_based_uiext";
      Motcle motlu;
      is >> motlu;
      Cerr << mot << " " << motlu << finl;
      Cout << "Setting the type of interpolation for calculer_dI_dt to " << motlu << "." << finl;
      int rang = motcles2.search(motlu);
      switch(rang)
        {
        case Navier_Stokes_FT_Disc_interne::INTERP_STANDARD:
          {
            variables_internes_.type_interpol_indic_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::INTERP_STANDARD;
            if (Process::je_suis_maitre())
              Cerr << " The interpolation of indicatrice to faces in calculer_dI_dt is based on the historical way" << " where a mean value + upwind is used." << finl;
            return 1;
          }
        case Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE:
          {
            variables_internes_.type_interpol_indic_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE;
            if (Process::je_suis_maitre())
              Cerr << " The interpolation of indicatrice to faces in calculer_dI_dt is based on the field indicatrice_faces" << " as defined by the interfacial transport option." << finl;
            return 1;
          }
        case Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED:
          {
            variables_internes_.type_interpol_indic_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED;
            if (Process::je_suis_maitre())
              Cerr << " The interpolation of indicatrice to faces in calculer_dI_dt is based on the interfacial area" << " and on the normal to the interface." << finl;
            return 1;
          }
        case Navier_Stokes_FT_Disc_interne::INTERP_STANDARD_UVEXT:
          {
            variables_internes_.type_interpol_indic_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::INTERP_STANDARD_UVEXT;
            if (Process::je_suis_maitre())
              Cerr << " The interpolation of indicatrice to faces in calculer_dI_dt is based on the historical way" << " where a mean value + upwind is used." << " Additionally, uv_ext is used."
                   << finl;
            return 1;
          }
        case Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE_UVEXT:
          {
            variables_internes_.type_interpol_indic_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE_UVEXT;
            if (Process::je_suis_maitre())
              Cerr << " The interpolation of indicatrice to faces in calculer_dI_dt is based on the field indicatrice_faces" << " as defined by the interfacial transport option."
                   << " Additionally, uv_ext is used." << finl;
            return 1;
          }
        case Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED_UVEXT:
          {
            variables_internes_.type_interpol_indic_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED_UVEXT;
            if (Process::je_suis_maitre())
              Cerr << " The interpolation of indicatrice to faces in calculer_dI_dt is based on the interfacial area" << " and on the normal to the interface." << " Additionally, uv_ext is used."
                   << finl;
            return 1;
          }
        case Navier_Stokes_FT_Disc_interne::INTERP_STANDARD_UIEXT:
          {
            variables_internes_.type_interpol_indic_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::INTERP_STANDARD_UIEXT;
            if (Process::je_suis_maitre())
              Cerr << " The interpolation of indicatrice to faces in calculer_dI_dt is based on the historical way" << " where a mean value + upwind is used." << " Additionally, ui_ext is used."
                   << finl;
            return 1;
          }
        case Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE_UIEXT:
          {
            variables_internes_.type_interpol_indic_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE_UIEXT;
            if (Process::je_suis_maitre())
              Cerr << " The interpolation of indicatrice to faces in calculer_dI_dt is based on the field indicatrice_faces" << " as defined by the interfacial transport option."
                   << " Additionally, ui_ext is used." << finl;
            return 1;
          }
        case Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED_UIEXT:
          {
            variables_internes_.type_interpol_indic_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED_UIEXT;
            if (Process::je_suis_maitre())
              Cerr << " The interpolation of indicatrice to faces in calculer_dI_dt is based on the interfacial area" << " and on the normal to the interface." << " Additionally, ui_ext is used."
                   << finl;
            return 1;
          }
        default:
          Cerr << "Transport_Interfaces_FT_Disc::lire\n" << "The options for methode_transport are :\n" << motcles2;
          Process::exit();
        }
    }
  else if (mot == "OutletCorrection_pour_dI_dt")
    {
      Motcles motcles2(4);
      motcles2[0] = "no_correction";
      motcles2[1] = "CORRECTION_GHOST_INDIC";
      motcles2[2] = "zero_net_flux_on_mixed_cells";
      motcles2[3] = "zero_out_flux_on_mixed_cells";
      Motcle motlu;
      is >> motlu;
      Cerr << mot << " " << motlu << finl;
      Cout << "Setting the type of correction at outlet BC for calculer_dI_dt to " << motlu << "." << finl;
      int rang = motcles2.search(motlu);
      switch(rang)
        {
        case Navier_Stokes_FT_Disc_interne::NO_CORRECTION:
          {
            variables_internes_.OutletCorrection_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::NO_CORRECTION;
            if (Process::je_suis_maitre())
              Cerr << " No correction of div(chi u) at exit (historical way)" << finl;
            return 1;
          }
        case Navier_Stokes_FT_Disc_interne::CORRECTION_GHOST_INDIC:
          {
            variables_internes_.OutletCorrection_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::CORRECTION_GHOST_INDIC;
            if (Process::je_suis_maitre())
              Cerr << " Correction of chi in ghost cells (virtually)" << finl;
            return 1;
          }
        case Navier_Stokes_FT_Disc_interne::ZERO_NET_FLUX_ON_MIXED_CELLS:
          {
            variables_internes_.OutletCorrection_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::ZERO_NET_FLUX_ON_MIXED_CELLS;
            if (Process::je_suis_maitre())
              Cerr << " correction of div(chi u) at exit : zero divergence on cells touching outlet." << finl;
            return 1;
          }
        case Navier_Stokes_FT_Disc_interne::ZERO_OUT_FLUX_ON_MIXED_CELLS:
          {
            variables_internes_.OutletCorrection_pour_dI_dt_ = Navier_Stokes_FT_Disc_interne::ZERO_OUT_FLUX_ON_MIXED_CELLS;
            if (Process::je_suis_maitre())
              {
                Cerr << " correction of div(chi u) at exit : zero vapour mass flux on cells touching outlet." << finl;
                Cerr << " This is a bad option because it does not let vapour get out explicitly (prevents interface contact)" << finl;
                Cerr << " Should not be used or with great care." << finl;
                // Process::exit();
              }
            return 1;
          }
        default:
          Cerr << "Transport_Interfaces_FT_Disc::lire\n" << "The options for methode_transport are :\n" << motcles2;
          Process::exit();
        }
    }
  else
    return Navier_Stokes_Turbulent::lire_motcle_non_standard(mot, is);
  return 1;
}

const Champ_Don& Navier_Stokes_FT_Disc::diffusivite_pour_transport() const
{
  return champ_mu_;
}

void Navier_Stokes_FT_Disc::associer_pb_base(const Probleme_base& un_probleme)
{
  if (!sub_type(Probleme_FT_Disc_gen, un_probleme))
    {
      Cerr << "Error for the method Navier_Stokes_FT_Disc::associer_pb_base\n";
      Cerr << " Navier_Stokes_FT_Disc equation must be associated to\n";
      Cerr << " a Probleme_FT_Disc_gen problem type\n";
      exit();
    }
  probleme_ft_ = ref_cast(Probleme_FT_Disc_gen, un_probleme);
  Navier_Stokes_std::associer_pb_base(un_probleme);
}

const Fluide_Diphasique& Navier_Stokes_FT_Disc::fluide_diphasique() const
{
  const REF(Fluide_Diphasique) &fluide_dipha = variables_internes().ref_fluide_diphasique;
  if (!fluide_dipha.non_nul())
    {
      Cerr << "Error for the method Navier_Stokes_FT_Disc::fluide_diphasique()\n";
      Cerr << " The fluid has not been associated\n";
      Cerr << " (check that the fluid is of Fluide_Diphasique type)" << finl;
      assert(0);
      exit();
    }
  return fluide_dipha.valeur();
}

void Navier_Stokes_FT_Disc::associer_milieu_base(const Milieu_base& un_fluide)
{
  if (sub_type(Fluide_Diphasique, un_fluide))
    variables_internes().ref_fluide_diphasique = ref_cast(Fluide_Diphasique, un_fluide);
  else
    Navier_Stokes_Turbulent::associer_milieu_base(un_fluide);
}

const Milieu_base& Navier_Stokes_FT_Disc::milieu() const
{
  if (variables_internes().ref_fluide_diphasique.non_nul())
    return variables_internes().ref_fluide_diphasique.valeur();
  else
    return Navier_Stokes_Turbulent::milieu();
}

Milieu_base& Navier_Stokes_FT_Disc::milieu()
{
  if (variables_internes().ref_fluide_diphasique.non_nul())
    return variables_internes().ref_fluide_diphasique.valeur();
  else
    return Navier_Stokes_Turbulent::milieu();
}

/*! @brief Discretisation des champs utilises dans cette equation.
 *
 * Fonction appelee par Probleme_base::discretiser.
 *   B. Mathieu : a titre experimental, au lieu de dupliquer les noms
 *                des champs ici et dans "a_pour_champ_Fonc", on stocke
 *                les champs dans une liste. (voir a_pour_champ_fonc).
 *
 */
void Navier_Stokes_FT_Disc::discretiser()
{
  Navier_Stokes_Turbulent::discretiser();
  const Discretisation_base& dis = discretisation();
  const double temps = schema_temps().temps_courant();
  const Domaine_dis_base& mon_dom_dis = domaine_dis();
  LIST(REF(Champ_base)) &champs_compris = variables_internes().liste_champs_compris;

  dis.discretiser_champ("champ_elem", mon_dom_dis, "diffusivite", "m^2/s", 1 /* composantes */, temps, champ_nu_);
  champs_compris.add(champ_nu_.valeur());
  //Nouvelle formulation
  champs_compris_.ajoute_champ(champ_nu_);

  dis.discretiser_champ("champ_elem", mon_dom_dis, "viscosite_dynamique", "kg/m.s", 1 /* composantes */, temps, champ_mu_);
  champs_compris.add(champ_mu_.valeur());
  champs_compris_.ajoute_champ(champ_mu_);

  dis.discretiser_champ("champ_elem", mon_dom_dis, "masse_volumique", "kg/m3", 1 /* composantes */, temps, champ_rho_elem_);
  champs_compris.add(champ_rho_elem_.valeur());
  champs_compris_.ajoute_champ(champ_rho_elem_);

  // La masse volumique discretisee sur les volumes de controle de la vitesse
  dis.discretiser_champ("vitesse", mon_dom_dis, "masse_volumique_vnodes", "kg/m3", 1 /* composantes */, temps, champ_rho_faces_);
  champs_compris.add(champ_rho_faces_.valeur());
  champs_compris_.ajoute_champ(champ_rho_faces_);

  // Variables internes
  dis.discretiser_champ("pression", mon_dom_dis, "second_membre_projection", "", 1 /* composantes */, temps, variables_internes().second_membre_projection);
  champs_compris.add(variables_internes().second_membre_projection.valeur());
  champs_compris_.ajoute_champ(variables_internes().second_membre_projection);
  dis.discretiser_champ("champ_elem", mon_dom_dis, "second_membre_projection_jump", "", 1 /* composantes */, temps, variables_internes().second_membre_projection_jump_);
  champs_compris.add(variables_internes().second_membre_projection_jump_.valeur());
  champs_compris_.ajoute_champ(variables_internes().second_membre_projection_jump_);
  dis.discretiser_champ("vitesse", mon_dom_dis, "gradient_pression_interne", "m/s2", -1 /* nb composantes par defaut */, temps, variables_internes().gradient_pression);
  champs_compris.add(variables_internes().gradient_pression.valeur());
  champs_compris_.ajoute_champ(variables_internes().gradient_pression);
  dis.discretiser_champ("vitesse", mon_dom_dis, "derivee_u_etoile", "", -1 /* nb composantes par defaut */, temps, variables_internes().derivee_u_etoile);
  champs_compris.add(variables_internes().derivee_u_etoile.valeur());
  champs_compris_.ajoute_champ(variables_internes().derivee_u_etoile);
  dis.discretiser_champ("vitesse", mon_dom_dis, "terme_diffusion_vitesse", "", -1 /* nb composantes par defaut */, temps, variables_internes().terme_diffusion);
  champs_compris.add(variables_internes().terme_diffusion.valeur());
  champs_compris_.ajoute_champ(variables_internes().terme_diffusion);
  dis.discretiser_champ("vitesse", mon_dom_dis, "terme_convection_vitesse", "", -1 /* nb composantes par defaut */, temps, variables_internes().terme_convection);
  champs_compris.add(variables_internes().terme_convection.valeur());
  champs_compris_.ajoute_champ(variables_internes().terme_convection);
  dis.discretiser_champ("vitesse", mon_dom_dis, "terme_source_vitesse", "", -1 /* nb composantes par defaut */, temps, variables_internes().terme_source);
  champs_compris.add(variables_internes().terme_source.valeur());
  champs_compris_.ajoute_champ(variables_internes().terme_source);
  dis.discretiser_champ("vitesse", mon_dom_dis, "terme_source_interfaces", "", -1 /* nb composantes par defaut */, temps, variables_internes().terme_source_interfaces);
  champs_compris.add(variables_internes().terme_source_interfaces.valeur());
  champs_compris_.ajoute_champ(variables_internes().terme_source_interfaces);
  if (dis.que_suis_je() == "VEFPreP1B")
    {
      dis.discretiser_champ("pression", mon_dom_dis, "indicatrice_p1b", "", 1 /* composantes */, temps, variables_internes().indicatrice_p1b);
      champs_compris.add(variables_internes().indicatrice_p1b.valeur());
      champs_compris_.ajoute_champ(variables_internes().indicatrice_p1b);
    }

  dis.discretiser_champ("vitesse", mon_dom_dis, "gradient_indicatrice", "", -1 /* nb composantes par defaut */, temps, variables_internes().gradient_indicatrice);
  champs_compris.add(variables_internes().gradient_indicatrice.valeur());
  champs_compris_.ajoute_champ(variables_internes().gradient_indicatrice);
  dis.discretiser_champ("vitesse", mon_dom_dis, "potentiel_faces", "", 1 /* composante */, temps, variables_internes().potentiel_faces);
  champs_compris.add(variables_internes().potentiel_faces.valeur());
  champs_compris_.ajoute_champ(variables_internes().potentiel_faces);
  dis.discretiser_champ("champ_elem", mon_dom_dis, "potentiel_elements", "", 1 /* composante */, temps, variables_internes().potentiel_elements);
  champs_compris.add(variables_internes().potentiel_elements.valeur());
  champs_compris_.ajoute_champ(variables_internes().potentiel_elements);
  dis.discretiser_champ("vitesse", mon_dom_dis, "vitesse_delta_interface", "m/s", -1 /* nb composantes par defaut */, 1, temps, variables_internes().delta_u_interface);
  // Pour pouvoir faire filtrer_L2:
  variables_internes().delta_u_interface->associer_eqn(*this);
  champs_compris.add(variables_internes().delta_u_interface.valeur());
  champs_compris_.ajoute_champ(variables_internes().delta_u_interface);
  dis.discretiser_champ("pression", mon_dom_dis, "pression_laplacien_d", "", 1 /* composante */, temps, variables_internes().laplacien_d);
  champs_compris.add(variables_internes().laplacien_d.valeur());
  champs_compris_.ajoute_champ(variables_internes().laplacien_d);
  dis.discretiser_champ("temperature", mon_dom_dis, "temperature_mpoint", "", 1 /* composante */, temps, variables_internes().mpoint);
  champs_compris.add(variables_internes().mpoint.valeur());
  champs_compris_.ajoute_champ(variables_internes().mpoint);
  dis.discretiser_champ("temperature", mon_dom_dis, "temperature_mpointv", "", 1 /* composante */, temps, variables_internes().mpoint_vap);
  champs_compris.add(variables_internes().mpoint_vap.valeur());
  champs_compris_.ajoute_champ(variables_internes().mpoint_vap);
  // Pour variation temporelle dI_dt
  dis.discretiser_champ("pression", mon_dom_dis, "derivee_temporelle_indicatrice", "", 1 /* composante */, temps, variables_internes().derivee_temporelle_indicatrice);
  champs_compris.add(variables_internes().derivee_temporelle_indicatrice.valeur());
  champs_compris_.ajoute_champ(variables_internes().derivee_temporelle_indicatrice);

  // Eulerian Interfacial area :
  dis.discretiser_champ("champ_elem", mon_dom_dis, "interfacial_area", "m2", 1 /* composante */, temps, variables_internes().ai);
  champs_compris.add(variables_internes().ai.valeur());
  champs_compris_.ajoute_champ(variables_internes().ai);

  // Velocity jump "u0" computed for phase 0 :
  Nom nom = Nom("vitesse_jump0_") + le_nom();
  dis.discretiser_champ("vitesse", mon_dom_dis, nom, "m/s", -1 /* nb composantes par defaut */, 1, temps, variables_internes().vitesse_jump0_);
  variables_internes().vitesse_jump0_->associer_eqn(*this);
  champs_compris.add(variables_internes().vitesse_jump0_.valeur());
  champs_compris_.ajoute_champ(variables_internes().vitesse_jump0_);
}

/*! @brief Methode surchargee de Navier_Stokes_std, appelee par Navier_Stokes_std::discretiser().
 *
 *   L'assembleur pression est particulier pour le front-tracking
 *   en VEF (en attendant qu'on factorise tous ces assembleurs pression)
 *
 */
void Navier_Stokes_FT_Disc::discretiser_assembleur_pression()
{
  const Discretisation_base& dis = discretisation();
  if (dis.que_suis_je() == "VDF")
    {
      // Assembleur pression generique (prevu pour rho_variable)
      assembleur_pression_.typer("Assembleur_P_VDF");
    }
  else if (dis.que_suis_je() == "VEFPreP1B")
    {
      // Assembleur pression generique (prevu pour rho_variable)
      assembleur_pression_.typer("Assembleur_P_VEFPreP1B");
    }

  Assembleur_base& assembleur = assembleur_pression_.valeur();
  assembleur.associer_domaine_dis_base(domaine_dis());
  // B. Mathieu, 07_2004 : premier jet de la methode, on resout en pression.
  // Version actuelle : pas d'increment de pression
  assembleur.set_resoudre_increment_pression(0);
  assembleur.set_resoudre_en_u(1);
}

/*! @brief methode appelee par Navier_Stokes_std::preparer_calcul.
 *
 * ..
 *
 */
void Navier_Stokes_FT_Disc::projeter()
{
  if (Process::je_suis_maitre() && limpr())
    Cerr << "Navier_Stokes_FT_Disc::projeter does nothing" << finl;
}

/*! @brief methode appelee par Probleme_base::preparer_calcul()
 *
 */
int Navier_Stokes_FT_Disc::preparer_calcul()
{
  Cerr << "Navier_Stokes_FT_Disc::preparer_calcul()" << finl;
  Equation_base::preparer_calcul();

  le_modele_turbulence->preparer_calcul();

  {
    // Si le mot cle equation_interfaces_proprietes_fluide a ete specifie,
    // l'indicatrice de phase correspondant a l'equation sert a calculer
    // les proprietes physiques du milieu.
    // Sinon, le milieu doit etre de type Fluide_Incompressible.

    REF(Transport_Interfaces_FT_Disc) &ref_equation = variables_internes().ref_eq_interf_proprietes_fluide;

    if (ref_equation.non_nul())
      {

        // Couplage Navier-Stokes / Fluide : les interfaces determinent la position des phases.
        // On utilise les proprietes du fluide diphasique

        if (Process::je_suis_maitre())
          {
            Cerr << "Initialisation of the physical properties (rho, mu, ...)\n" << " based on the indicatrice field of the equation " << ref_equation->le_nom() << finl;

          }
        FT_disc_calculer_champs_rho_mu_nu_dipha(domaine_dis(), fluide_diphasique(), ref_equation->get_update_indicatrice().valeurs(), // indicatrice
                                                champ_rho_elem_->valeurs(), champ_nu_->valeurs(), champ_mu_->valeurs(), champ_rho_faces_->valeurs());
      }
    else
      {

        // Pas de couplage Navier-Stokes / Fluide : le fluide est monophasique.

        const Fluide_Incompressible& phase_0 = ref_cast(Fluide_Incompressible, milieu());
        const Domaine_dis_base& zdis = domaine_dis();
        FT_disc_calculer_champs_rho_mu_nu_mono(zdis, phase_0, champ_rho_elem_, champ_mu_, champ_nu_, champ_rho_faces_);
      }

    // Premiere version :les termes sources sont homogenes a rho*v
    //   sources().associer_champ_rho(champ_rho_elem_.valeur());
    // Nouvelle version: les termes "sources()" sont homogenes a v.
    //   on ne fait rien.
  }

  DoubleTab& tab_vitesse = inconnue()->valeurs();

  // On assemble la matrice de pression pour la premiere fois.
  assembleur_pression()->assembler_rho_variable(matrice_pression_, champ_rho_faces_.valeur());
  // Informe le solveur que la matrice a change :
  solveur_pression()->reinit();
  // Calcul du second membre :
  //  div(vpoint) a l'interieur du domaine,
  //  prise en compte des conditions aux limites de pression/vitesse
  DoubleTab& secmem = variables_internes().second_membre_projection->valeurs();
  divergence.calculer(tab_vitesse, secmem);
  secmem *= -1;
  // Il faut faire ceci car on ne resout pas en "increment de pression":
  assembleur_pression_->modifier_secmem(secmem);

  // Ajout pour la sauvegarde au premier pas de temps si reprise
  la_pression->changer_temps(schema_temps().temps_courant());

  // Resolution du systeme en pression : calcul de la_pression
  solveur_pression_.resoudre_systeme(matrice_pression_.valeur(), secmem, la_pression->valeurs());
  assembleur_pression_->modifier_solution(la_pression->valeurs());
  // Calcul d(u)/dt = vpoint + 1/rho*grad(P)
  DoubleTab& gradP = variables_internes().gradient_pression->valeurs();
  gradient.calculer(la_pression->valeurs(), gradP);
  solveur_masse->appliquer(gradP);
  // Correction de la vitesse :
  if (projection_a_faire()) // Temporaire pour permettre de ne pas resoudre NS avec mettant operateurs nuls et projection_initiale 0
    {
      int i, j;
      const DoubleTab& rho_faces = champ_rho_faces_->valeurs();
      const int n = tab_vitesse.dimension(0);
      const int m = tab_vitesse.line_size();

      for (i = 0; i < n; i++)
        for (j = 0; j < m; j++)
          tab_vitesse(i, j) -= gradP(i, j) / rho_faces(i);

      tab_vitesse.echange_espace_virtuel();
    }

  if (le_traitement_particulier.non_nul())
    {
      if (le_traitement_particulier->support_ok())
        le_traitement_particulier->associer_champ_masse_volumique(champ_rho_faces_.valeur());
      le_traitement_particulier->preparer_calcul_particulier();
    }

  return 1;
}

const int& Navier_Stokes_FT_Disc::get_is_penalized() const
{
  return variables_internes().is_penalized;
}

const int& Navier_Stokes_FT_Disc::get_new_mass_source() const
{
  return variables_internes().new_mass_source_;
}

const DoubleTab& Navier_Stokes_FT_Disc::get_interfacial_area() const
{
  return variables_internes().ai->valeurs();
}

DoubleTab& Navier_Stokes_FT_Disc::get_set_interfacial_area()
{
  return variables_internes().ai->valeurs();
}

const DoubleTab& Navier_Stokes_FT_Disc::get_mpoint() const
{
  return variables_internes().mpoint->valeurs();
}

DoubleTab& Navier_Stokes_FT_Disc::get_set_mpoint()
{
  return variables_internes().mpoint->valeurs();
}

void Navier_Stokes_FT_Disc::preparer_pas_de_temps()
{
}

void Navier_Stokes_FT_Disc::mettre_a_jour(double temps)
{
  Navier_Stokes_Turbulent::mettre_a_jour(temps);

  champ_rho_elem_->mettre_a_jour(temps);
  champ_rho_faces_->mettre_a_jour(temps);
  champ_mu_->mettre_a_jour(temps);
  champ_nu_->mettre_a_jour(temps);
}

const SolveurSys& Navier_Stokes_FT_Disc::get_solveur_pression() const
{
  return solveur_pression_;
}

/*! @brief Calcul des forces de tension superficielles (F_sigma) et de la partie a rotationnel non nul de la gravite (G_rot) (si GRAVITE_GRAD_I) :
 *
 *   F_sigma = INTEGRALE sur le volume de controle (
 *             sigma_aux_faces * courbure_aux_faces * gradient(indicatrice)
 *             + gradient_sigma )
 *   G_rot   = INTEGRALE sur le volume de controle (
 *                phi * gradient(rho) )   (avec phi = potentiel de pesanteur)
 *
 * @param (gradient_indicatrice) le gradient de l'indicatrice issu de l'operateur "gradient", donc homogene a l'integrale du gradient sur les volumes de controle de la vitesse.
 * @param (potentiel_faces) un champ aux faces a une composante, ou on stocke le "potentiel aux faces"
 * @param (champ) le champ aux faces (meme discretisation que la vitesse) ou on stocke le terme source des forces superficielles.
 */
//
//  The method is no longer const because it changes a member of variables_internes() to store the Eulerian interfacial Area.
void Navier_Stokes_FT_Disc::calculer_champ_forces_superficielles(const Maillage_FT_Disc& maillage, const Champ_base& gradient_indicatrice, Champ_base& potentiel_elements, Champ_base& potentiel_faces,
                                                                 Champ_base& champ)
{
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis());
  // Nombre de faces
  const int nb_faces = domaine_vf.nb_faces();
  {
    // Le champ et le gradient de l'indicatrice doivent etre aux faces
    assert(champ.valeurs().dimension(0) == nb_faces);
    assert(potentiel_faces.valeurs().dimension(0) == nb_faces);
    assert(gradient_indicatrice.valeurs().dimension(0) == nb_faces);
  }
  const int dim = Objet_U::dimension;

  // Calcul du "potentiel aux sommets du maillage"
  ArrOfDouble potentiel_sommets;
  potentiel_sommets.resize_array(maillage.nb_sommets(), RESIZE_OPTIONS::NOCOPY_NOINIT);

  // (ce potentiel est constant sur chaque portion connexe d'interface si
  //  l'interface est a l'equilibre).
  //  courbure * sigma + potentiel_gravite_phase_1 - potentiel_gravite_phase_0
  {
    // Pour l'instant, la tension de surface est constante
    const Fluide_Diphasique& fluide_dipha = fluide_diphasique();
    const double sigma = fluide_dipha.sigma();

    const int n = maillage.nb_sommets();
    const ArrOfDouble& courbure_sommets = maillage.get_update_courbure_sommets();
    //Calcul du "potentiel aux sommets"
    int i;
    {
      const double clipping_courbure_max = variables_internes().clipping_courbure_interface;
      int clip_counter = 0;
      for (i = 0; i < n; i++)
        {
          double c = courbure_sommets[i];
          // Clipping de la courbure: si la courbure est superieure a la
          // valeur maxi autorisee, on limite (permet de ne pas plomber le
          // pas de temps s'il y a une singularite geometrique dans le maillage)
          if (std::fabs(c) > clipping_courbure_max)
            {
              clip_counter++;
              c = ((c > 0) ? 1. : -1.) * clipping_courbure_max;
            }
          potentiel_sommets[i] = c * sigma;
        }
      clip_counter = mp_sum(clip_counter);
      if (clip_counter > 0 && je_suis_maitre())
        {
          Cerr << "Navier_Stokes_FT_Disc::calculer_champ_forces_superficielles : clip_count " << clip_counter << finl;
        }
    }

    //Ajout des termes de gravite
    if (variables_internes().terme_gravite_ == Navier_Stokes_FT_Disc_interne::GRAVITE_GRAD_I)
      {
        if (milieu().a_gravite())
          {
            const double rho_0 = fluide_dipha.fluide_phase(0).masse_volumique()->valeurs()(0, 0);
            const double rho_1 = fluide_dipha.fluide_phase(1).masse_volumique()->valeurs()(0, 0);
            const double delta_rho = rho_1 - rho_0;

            // Pour l'instant : gravite uniforme g => phi(s) = - x scalaire g
            const DoubleTab& gravite = milieu().gravite().valeurs();
            if (gravite.nb_dim() != 2 || gravite.line_size() != dim)
              {
                Cerr << "Error for the method calculer_champ_forces_superficielles\n";
                Cerr << "gravite.dimension(1) != Objet_U::dimension" << finl;
                Process::exit();
              }
            const DoubleTab& sommets = maillage.sommets();

            for (i = 0; i < n; i++)
              {
                double p = 0.;
                for (int j = 0; j < dim; j++)
                  p += sommets(i, j) * gravite(0, j);

                if (is_repulsion)
                  {
                    double dx = 0.;
                    double px = sommets(i, 0);
                    if (px > maxx)
                      dx = px - maxx;
                    else if (px < minx)
                      dx = minx - px;

                    double dy = 0.;
                    double py = sommets(i, 1);
                    if (py > maxx)
                      dy = py - maxx;
                    else if (py < minx)
                      dy = minx - py;

                    p += sqrt(dx * dx + dy * dy) * pente;
                  }

                potentiel_sommets[i] -= p * delta_rho;
              }
          }
      }
  }
#if 0
  // Calcul du "potentiel aux faces" : interpolation sur les faces des elements
  // traverses par l'interface.
  {
    // Tableau des facettes du maillage interfaces
    const IntTab& facettes = maillage.facettes();
    // Tableau des faces des elements :
    const IntTab& elem_faces = domaine_vf.elem_faces();
    const int nb_faces_par_element = elem_faces.line_size();
    // Le tableau "potentiel aux faces" a remplir :
    DoubleTab& valeurs_potentiel_faces = potentiel_faces.valeurs();
    valeurs_potentiel_faces = 0.;
    // Somme des poids des contributions ajoutees aux faces
    ArrOfDouble poids(nb_faces);
    poids = 0.;

    // Boucle sur les faces de l'interface
    const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
    const ArrOfInt& index_facette_elem = intersections.index_facette();
    double pot[3] = {0., 0., 0.};

    int fa7;
    const int nb_facettes = maillage.nb_facettes();
    for (fa7 = 0; fa7 < nb_facettes; fa7++)
      {
        // Potentiel aux sommets de la facette:
        pot[0] = potentiel_sommets[facettes(fa7,0)];
        pot[1] = potentiel_sommets[facettes(fa7,1)];
        if (dim == 3)
          pot[2] = potentiel_sommets[facettes(fa7,2)];
        // Boucle sur les elements traverses par la facette
        int index = index_facette_elem[fa7];
        while (index >= 0)
          {
            const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
            // Calcul du potentiel au centre de l'intersection
            double p = 0.;
            int i;
            for (i = 0; i < 3; i++)
              p += data.barycentre_[i] * pot[i];
            // La contribution aux faces est ponderee par la surface d'intersection
            p *= data.surface_intersection_;

            // Ajout de la contribution aux faces de l'element
            const int element = data.numero_element_;
            for (i = 0; i < nb_faces_par_element; i++)
              {
                const int face = elem_faces(element, i);
                poids[face] += data.surface_intersection_;
                valeurs_potentiel_faces(face) += p;
              }
            index = data.index_element_suivant_;
          }
      }

    // Il reste a diviser les valeurs aux face par la somme des poids
    for (int face = 0; face < nb_faces; face++)
      {
        double p = poids[face];
        if (p > 0.)
          valeurs_potentiel_faces(face) /= p;
      }
    valeurs_potentiel_faces.echange_espace_virtuel();
  }
#else
  // Calcul du "potentiel aux elements" :
  DoubleTab poids(potentiel_elements.valeurs());
  {
    // Tableau des facettes du maillage interfaces
    const IntTab& facettes = maillage.facettes();
    // Le tableau "potentiel aux faces" a remplir :
    DoubleTab& valeurs_potentiel_elements = potentiel_elements.valeurs();
    const ArrOfDouble& surface_facettes = maillage.get_update_surface_facettes();
    valeurs_potentiel_elements = 0.;
    // Somme des poids des contributions ajoutees aux faces
    poids = 0.;

    // Boucle sur les faces de l'interface
    const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
    const ArrOfInt& index_facette_elem = intersections.index_facette();
    double pot[3] = { 0., 0., 0. };

    int fa7;
    const int nb_facettes = maillage.nb_facettes();
    for (fa7 = 0; fa7 < nb_facettes; fa7++)
      {
        // Potentiel aux sommets de la facette:
        pot[0] = potentiel_sommets[facettes(fa7, 0)];
        pot[1] = potentiel_sommets[facettes(fa7, 1)];
        if (dim == 3)
          pot[2] = potentiel_sommets[facettes(fa7, 2)];
        // Boucle sur les elements traverses par la facette
        int index = index_facette_elem[fa7];
        const double surface_facette = surface_facettes[fa7];
        while (index >= 0)
          {
            const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
            // Calcul du potentiel au centre de l'intersection
            double p = 0.;
            int i;
            for (i = 0; i < 3; i++)
              p += data.barycentre_[i] * pot[i];
            // La contribution aux faces est ponderee par la surface d'intersection
#ifdef AVEC_BUG_SURFACES
            const double surf = data.surface_intersection_;
#else
            const double surf = data.fraction_surface_intersection_ * surface_facette;
#endif
            p *= surf;

            // Ajout de la contribution a l'element
            const int element = data.numero_element_;
            valeurs_potentiel_elements(element) += p;
            poids(element) += surf;

            index = data.index_element_suivant_;
          }
      }
    valeurs_potentiel_elements.echange_espace_virtuel();
    poids.echange_espace_virtuel();
    if (champ.valeurs().line_size() > 1)   // VEF ?
      {
        // Compute a potential on elements that are neighbours of
        // elements containing an interface :
        //  For VEF, the gradient of the indicator function can be
        //  non zero on the faces of these elements.
        int element;
        const int nb_elements = poids.dimension(0);
        const IntTab& face_voisins = le_dom_dis.valeur().face_voisins();
        const Domaine_VF& domainevf = ref_cast(Domaine_VF, le_dom_dis.valeur());
        const IntTab& elem_faces = domainevf.elem_faces();
        const int nb_faces_par_element = elem_faces.line_size();
        DoubleVect copie_valeurs_potentiel_elements(valeurs_potentiel_elements);
        DoubleVect copie_poids(poids);
        for (element = 0; element < nb_elements; element++)
          {
            double potential = 0.; // sum of potentials of neighbouring elements
            double p = 0.;   // sum of weights of neighbouring elements
            int i_face;
            for (i_face = 0; i_face < nb_faces_par_element; i_face++)
              {
                const int face = elem_faces(element, i_face);
                const int elem_voisin_0 = face_voisins(face, 0);
                const int elem_voisin_1 = face_voisins(face, 1);
                const int elem_voisin = elem_voisin_0 + elem_voisin_1 - element;
                if (elem_voisin >= 0)   // Not a boundary of the domain ?
                  {
                    potential += copie_valeurs_potentiel_elements(elem_voisin);
                    p += copie_poids(elem_voisin);
                  }
              }
            const double old_weight = copie_poids(element);
            // Do not change values for elements that contain an interface:
            if (p > 0. && old_weight == 0.)
              {
                // Decrease weight so that values on faces adjacent to elements
                // containing an interface are almost untouched.
                static double reduction_factor = 0.1;
                potential = potential * reduction_factor;
                p = p * reduction_factor;
                // Assign value
                valeurs_potentiel_elements(element) = potential;
                poids(element) = p;
              }
          }
        valeurs_potentiel_elements.echange_espace_virtuel();
        poids.echange_espace_virtuel();
        Debog::verifier("Navier_Stokes_FT_Disc::calculer_champ_forces_superficielles poids:", poids);
      }
  }

  // I take the opportunity to store the Eulerian Interfacial Area...
  DoubleTab& interfacial_area = get_set_interfacial_area();
  interfacial_area = poids;
  {
    // Boucle sur les faces
    int face;
    DoubleTab& valeurs_potentiel_faces = potentiel_faces.valeurs();
    valeurs_potentiel_faces = 0.;
    const DoubleTab& valeurs_potentiel_elements = potentiel_elements.valeurs();
    const int nb_faces_pot = valeurs_potentiel_faces.dimension(0);
    const IntTab& face_voisins = le_dom_dis.valeur().face_voisins();
    for (face = 0; face < nb_faces_pot; face++)
      {
        double p = 0.; // Somme des poids des deux elements voisins
        double pot = 0.; // Somme des potentiels
        // Boucle sur les deux elements voisins de la face
        for (int i = 0; i < 2; i++)
          {
            int element = face_voisins(face, i);
            if (element >= 0)
              {
                // If neighbour exists (not a boundary face)
                p += poids(element);
                pot += valeurs_potentiel_elements(element);
              }
          }
        if (p > 0.)
          valeurs_potentiel_faces(face) = pot / p;
      }
    valeurs_potentiel_faces.echange_espace_virtuel();
    Debog::verifier("Navier_Stokes_FT_Disc::calculer_champ_forces_superficielles valeurs_potentiel_faces:", valeurs_potentiel_faces);
  }
#endif
  // Derniere operation : calcul du champ
  //   champ = potentiel_faces * gradient_indicatrice
  {
    const DoubleTab& valeurs_potentiel_faces = potentiel_faces.valeurs();
    const DoubleTab& valeurs_gradient_i = gradient_indicatrice.valeurs();
    DoubleTab& valeurs_champ = champ.valeurs();
    const int nb_compo = valeurs_champ.line_size(); // 1 for VDF, 3 for VEF

    for (int face = 0; face < nb_faces; face++)
      {
        const double p = valeurs_potentiel_faces(face);
        for (int i = 0; i < nb_compo; i++)
          valeurs_champ(face, i) = valeurs_gradient_i(face, i) * p;
      }

    valeurs_champ.echange_espace_virtuel();
    Debog::verifier("Navier_Stokes_FT_Disc::calculer_champ_forces_superficielles valeurs_champ:", valeurs_champ);
  }
}

double compute_indic_ghost(const int elem, const int num_face, const double indic, const int normale_sortante_au_domaine, const Domaine_VF& domVF, const Maillage_FT_Disc& maillage)
{
  double indic_ghost = 0.;
  ArrOfDouble normale(3), normale_face(3); // Normale sortante de I=0 vers I=1
  const double face_surface = domVF.face_surfaces(num_face); //  ==  0. ? 1. : domVF.face_surfaces(num_face);
  if (est_egal(face_surface, 0.))
    {
      // On est sur l'axe de rotation du bidim_axi :
      return indic;
    }
  for (int i = 0; i < Objet_U::dimension; i++)
    normale_face[i] = domVF.face_normales(num_face, i) / face_surface; // pour normalisation
  const double norm = maillage.compute_normale_element(elem, true /* NORMALIZE */, normale);
  if (est_egal(norm, 0.))
    // L'element est pure, non traverse.
    return indic;
  // normale_sortante_au_domaine pour regler le signe du prodscal avec la normale_SORTANTE
  const double prodscal = normale_sortante_au_domaine * dotproduct_array(normale, normale_face);
  double val = 1. / sqrt(2.); // Could be improved. Corresponds to cubic cell and Indic=0.5
  if (prodscal < -val)
    indic_ghost = 0.;
  else if (prodscal < 0)
    // Linear interp to indic if 0
    indic_ghost = indic * (1 + prodscal / val); // prodscal is negative
  else if (prodscal < val)
    indic_ghost = indic + (1. - indic) * (prodscal / val);
  else
    indic_ghost = 1.;
  return indic_ghost;
}

/*! @brief Calcul du gradient de l'indicatrice.
 *
 * Ce gradient est utilise pour calculer le second membre de l'equation de qdm,
 *   contenant les termes de tension de surface.
 *   En VEF, on commence par creer un champ P1B a partir du champ P0
 *   et on calcule le gradient.
 *   Design de classe a revoir pour separer VDF et VEF...
 *
 *  Parametre : indicatrice
 *  Signification : un champ aux elements (l'espace virtuel doit etre a jour)
 *  Parametre : gradient_i
 *  Signification : un champ discretise comme la vitesse dans lequel
 *   on met gradient(indicatrice).
 *
 */
#include<Neumann_sortie_libre.h>
#include<Dirichlet.h>
#include<Dirichlet_homogene.h>
#include<Periodique.h>
#include<Symetrie.h>
#include<Sortie_libre_rho_variable.h>
void Navier_Stokes_FT_Disc::calculer_gradient_indicatrice(const Champ_base& indicatrice, const DoubleTab& distance_interface_sommets, Champ_base& gradient_i)
{
  if (gradient_i.que_suis_je() == "Champ_Fonc_Face")
    {
      gradient.calculer(indicatrice.valeurs(), gradient_i.valeurs());
    }
  else
    {
      int i;
      const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis());
      const Domaine& domaine = domaine_vf.domaine();
      const int nb_elem_tot = domaine.nb_elem_tot();
      //const int nb_sommets = domaine.nb_som();
      //const IntTab & les_elems = domaine.les_elems();
      //const int nb_sommets_par_element = les_elems.dimension(1);

      // Calcul d'une indicatrice p1bulle
      DoubleTab& indic_p1b = variables_internes().indicatrice_p1b->valeurs();
      // Verification du support du champ indicatrice_p1b
      if (dimension == 2 && indic_p1b.size_totale() != nb_elem_tot + domaine.nb_som_tot())
        {
          Cerr << "The method Navier_Stokes_FT_Disc::calculer_gradient_indicatrice is developped" << finl;
          Cerr << "only for an indicatrice field discretized like P0+P1 for the 2D dimension case." << finl;
          Cerr << "Please change the discretization." << finl;
          Process::exit();
        }
      if (dimension == 3 && indic_p1b.size_totale() != nb_elem_tot + domaine.nb_som_tot() + domaine.nb_aretes_tot() && indic_p1b.size_totale() != nb_elem_tot + domaine.nb_som_tot())
        {
          Cerr << "The method Navier_Stokes_FT_Disc::calculer_gradient_indicatrice is developped" << finl;
          Cerr << "only for an indicatrice field discretized like P0+P1 or P0+P1+Pa for the 3D dimension case." << finl;
          Cerr << "Please change the discretization." << finl;
          Process::exit();
        }

      // Le champ P1bulle contient
      //   nelements valeurs au centre des elements
      // suivies de
      //   nsommets valeurs aux sommets
      indic_p1b = 0.;
      // Valeurs aux centres des elements = indicatrice a l'element
      // On recopie des valeurs sur les elements virtuels car on en
      // a besoin pour le calcul des valeurs aux sommets.
      for (i = 0; i < nb_elem_tot; i++)
        {
          indic_p1b(i) = indicatrice.valeurs()(i);
        }
#if 0
      // Premiere version : pas terrible parce que le stencil est tres large
      //  et il faut calculer la "courbure de l'interface" sur une grande largeur
      //  autour de l'interface.

      // Valeurs aux sommets = moyenne des valeurs des elements adjacents
      //  ponderee par le volume des elements.
      // Il y a la un choix a faire qui peut etre important...
      // (autre idee : projection L2 du champ discontinu sur l'espace
      // P1bulle)
      const DoubleVect& volume = domaine_vf.volumes();

      // Somme des poids aux sommets:
      ArrOfDouble poids(nb_sommets);
      poids = 0.;
      // Boucle sur les elements reels et virtuels :
      for (i = 0; i < nb_elem_tot; i++)
        {
          // Ajout d'une contribution a chaque sommet reel de l'element
          double p = volume(i);
          double x = indicatrice(i);
          int j;
          for (j = 0; j < nb_sommets_par_element; j++)
            {
              const int sommet = les_elems(i, j);
              if (sommet < nb_sommets)   // sommet reel ?
                {
                  // Le tableau des ddl contient nb_elem_tot valeurs aux elements
                  // suivies de nb_sommets_tot valeurs aux sommets.
                  // L'inconnue associee au "sommet" se trouve a l'indice
                  // nb_elem_tot + sommet.
                  indic_p1b(nb_elem_tot + sommet) += p * x;
                  poids[sommet] += p;
                }
            }
        }
      // Division par le poids
      for (i = 0; i < nb_sommets; i++)
        {
          const double p = poids[i];
          indic_p1b(nb_elem_tot + i) /= p;
        }
#else
#if 0
      // Calcul de l'indicatrice aux sommets. Deuxieme version, pour
      // limiter le stencil de gradient_indicatrice :
      // l'indicatrice vaut 1 si le sommet est adjacent a un element
      // de phase 1, 0 si le sommet est adjacent a un element de phase 0

      // 1) on met une valeur invalide dans les inconnues:
      for (i = 0; i < nb_sommets; i++)
        indic_p1b(nb_elem_tot + i) = -1;
      // 2) On met 1. ou 0. pour les sommets adjacents a un element
      //    dont l'indicatrice vaut 1. ou 0.
      for (i = 0; i < nb_elem_tot; i++)
        {
          const double indic_element = indic_p1b(i);
          if (indic_element == 0. || indic_element == 1.)
            {
              int j;
              for (j = 0; j < nb_sommets_par_element; j++)
                {
                  const int sommet = les_elems(i,j);
                  indic_p1b(nb_elem_tot + sommet) = indic_element;
                }
            }
        }
      // 3) Pour les sommets indecis, on prend 0 ou 1 selon le signe
      //    de la fonction distance pour ce sommet.
      static const double valeur_distance_invalide = -1e30;
      int error_count = 0;
      for (i = 0; i < nb_sommets; i++)
        {
          const double valeur = indic_p1b(nb_elem_tot + i);
          if (valeur < 0.)
            {
              double newval;
              const double d = distance_interface_sommets(i);
              if (d < valeur_distance_invalide)
                {
                  error_count++;
                  newval = 0.;
                }
              else
                {
                  newval = (d >= 0.) ? 1. : 0.;
                }
              indic_p1b(nb_elem_tot + i) = newval;
            }
        }
      if (error_count > 0)
        {
          Process::Journal()
              << "Navier_Stokes_FT_Disc::calculer_gradient_indicatrice\n"
              << error_count << " sommets ont une valeur indeterminee" << finl;
        }
#else
      // On met l'indicatrice sur P0 et 0 sur P1
#endif
#endif
      indic_p1b.echange_espace_virtuel();

      gradient.calculer(indic_p1b, gradient_i.valeurs());
    }

  const bool ghost_correction = (variables_internes_.OutletCorrection_pour_dI_dt_ == Navier_Stokes_FT_Disc_interne::CORRECTION_GHOST_INDIC);
  if (ghost_correction)
    {
      // Correction du gradient a la ligne de contact :
      const DoubleTab& inco = indicatrice.valeurs();
      DoubleTab& resu = gradient_i.valeurs();
      const int nbdimr = resu.line_size();
      //      assert_espace_virtuel_vect(inco);
      const Domaine_VF& zvf = ref_cast(Domaine_VF, domaine_dis());
      const Domaine_Cl_dis_base& zcldis = domaine_Cl_dis();
      const DoubleVect& face_surfaces = zvf.face_surfaces();
      const IntTab& face_voisins = domaine_dis().face_voisins();

      double coef;
      int n0, n1;
      // Boucle sur les bords pour traiter les conditions aux limites
      int ndeb, nfin, num_face;
      for (int n_bord = 0; n_bord < zvf.nb_front_Cl(); n_bord++)
        {
          // pour chaque Condition Limite on regarde son type
          // Si face de Dirichlet ou de Symetrie on ne fait rien
          // Si face de Neumann on calcule la contribution au terme source

          const Cond_lim& la_cl = zcldis.les_conditions_limites(n_bord);
          Cerr << que_suis_je() << "::calculer_gradient_indicatrice() correction du gradient a la CL : " << la_cl.valeur() << finl;
          const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
          ndeb = le_bord.num_premiere_face();
          nfin = ndeb + le_bord.nb_faces();

          // Correction periodicite :
          if (sub_type(Periodique, la_cl.valeur()))
            {
              for (num_face = ndeb; num_face < nfin; num_face++)
                {
                  n0 = face_voisins(num_face, 0);
                  n1 = face_voisins(num_face, 1);
                  if (!est_egal(n0, n1))
                    {
                      Cerr << "Periodic boundary condition with FT is not supported yet." << finl;
                      Process::exit();
                    }
                  coef = face_surfaces(num_face);          //*porosite_surf(num_face);
                  for (int k = 0; k < nbdimr; k++)
                    {
                      const int normale_sortante_au_domaine = (n0 == -1) ? 1 : -1; // Si on a le ghost dans la case 0, la normale sortante pointe vers "x-", on met donc -1 dans normale_sortante_au_domaine
                      const double dSk = normale_sortante_au_domaine * zvf.face_normales(num_face, k); // its magnitude is the surface.
                      // dSk/coef should be either 0 or 1...
                      assert(dSk / coef > -10. * Objet_U::precision_geom * Objet_U::precision_geom);
                      Cerr << "Check if sign of nk is compatible with the expression" << finl;
                      const double nk = zvf.face_normales(num_face, k);
                      Cerr << "Check if sign of nk is compatible with the expression" << finl;
                      Cerr << "nk=" << nk << finl;
                      Cerr << "dSk=" << dSk << finl;
                      Cerr << "num_face=" << num_face << finl;
                      Cerr << "voisin 0/1 =" << zvf.face_voisins(num_face, 0) << " " << zvf.face_voisins(num_face, 1) << finl;
                      Cerr << "code to validate" << finl;
                      Process::exit();
                      resu(num_face, k) = dSk * (inco(n1) - inco(n0));
                    }
                }
            }
          else if (sub_type(Symetrie, la_cl.valeur())) { /* Do nothing */ }
          else if ((sub_type(Dirichlet, la_cl.valeur())) || (sub_type(Neumann_sortie_libre, la_cl.valeur())) || (sub_type(Dirichlet_homogene, la_cl.valeur())))
            {
              REF(Transport_Interfaces_FT_Disc) &refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
              const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
              const Maillage_FT_Disc& maillage = eq_transport.maillage_interface();
              const Intersections_Elem_Facettes& intersections = maillage.intersections_elem_facettes();
              const ArrOfInt& index_elem = intersections.index_elem();
              for (num_face = ndeb; num_face < nfin; num_face++)
                {
                  n0 = face_voisins(num_face, 0);
                  n1 = face_voisins(num_face, 1);
                  int elem = n0 + n1 + 1;
                  //double surface_totale = 0.;
                  for (int k = 0; k < nbdimr; k++)
                    resu(num_face, k) = 0.;
                  if (index_elem[elem] < 0)
                    // element is pure
                    continue;
                  {
                    const double indic = inco(elem);
                    int normale_sortante_au_domaine = (n0 == -1) ? -1 : 1; // Si on a le ghost dans la case 0, la normale sortante pointe vers "x-", on met donc -1 dans normale_sortante_au_domaine
                    const double indic_ghost = compute_indic_ghost(elem, num_face, indic, normale_sortante_au_domaine, zvf, maillage);
                    //const double delta = zvdf.dim_elem(elem, orientation(num_face));
                    //const double volume = zvdf.volumes(elem);
                    coef = face_surfaces(num_face);                    //*porosite_surf(num_face);
                    if (nbdimr == 1)
                      {
                        resu(num_face) = coef * normale_sortante_au_domaine * (indic_ghost - indic);
                        assert(std::abs(coef - zvf.face_normales(num_face, zvf.orientation(num_face))) < Objet_U::precision_geom * Objet_U::precision_geom);
                      }
                    else
                      {
                        for (int k = 0; k < nbdimr; k++)
                          {
                            normale_sortante_au_domaine = 1; // En VEF, la normale dSk est toujours sortante?
                            const double dSk = normale_sortante_au_domaine * zvf.face_normales(num_face, k); // its magnitude is the surface.
                            // dSk/coef should be either 0 or 1...
                            assert(dSk / coef > -10. * Objet_U::precision_geom * Objet_U::precision_geom);
                            Cerr << "Check if sign of nk is compatible with the expression" << finl;
                            Cerr << "nk=" << dSk / coef << finl;
                            Cerr << "num_face=" << num_face << finl;
                            Cerr << "voisin 0/1 =" << zvf.face_voisins(num_face, 0) << " " << zvf.face_voisins(num_face, 1) << finl;
                            //Process::exit();
                            resu(num_face, k) = dSk * (indic_ghost - indic);
                          }
                      }
                  }
                }
            }
          // Fin de la boucle for
        }
    }
}

void Navier_Stokes_FT_Disc::correct_at_exit_bad_gradient(DoubleTab& u0) const
{
  // Correction du gradient a la sortie (car celui-ci ne doit pas se baser sur la CL de pression):
  // On prefere mettre la valeur d'en face qu'une valeur abherante. -> comme cela, la divergence (plus tard) tendra vers 0)
  const Domaine_VF& zvf = ref_cast(Domaine_VF, domaine_dis());
  const IntTab& elem_faces = zvf.elem_faces();
  const IntTab& face_voisins = zvf.face_voisins();
  const Domaine_Cl_dis_base& zcldis = domaine_Cl_dis();
  // Boucle sur les bords pour traiter les conditions aux limites
  for (int n_bord = 0; n_bord < zvf.nb_front_Cl(); n_bord++)
    {
      // pour chaque Condition Limite on regarde son type
      // Si face de Dirichlet ou de Symetrie on ne fait rien
      // Si face de Neumann on calcule la contribution au terme source
      const Cond_lim& la_cl = zcldis.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
      const int ndeb = le_bord.num_premiere_face();
      const int nfin = ndeb + le_bord.nb_faces();
      if ((sub_type(Dirichlet, la_cl.valeur())) || (sub_type(Neumann_sortie_libre, la_cl.valeur())) || (sub_type(Dirichlet_homogene, la_cl.valeur())))
        {
          // const Domaine_VDF& zvdf = ref_cast(Domaine_VDF, zvf);
          // Cerr << que_suis_je() << "::calculer_delta_u_interface() correction de u0 gradient(phi) a la CL : " <<  la_cl.valeur() << finl;
          const int nb_faces = elem_faces.dimension(1);
          for (int num_face = ndeb; num_face < nfin; num_face++)
            {
              const int elem = face_voisins(num_face, 0) + face_voisins(num_face, 1) + 1;
              int idx_face_de_lelem = 0;
              for (idx_face_de_lelem = 0; idx_face_de_lelem < nb_faces; idx_face_de_lelem++)
                {
                  if (elem_faces(elem, idx_face_de_lelem) == num_face)
                    break; // Face found
                }
              if (nb_faces == idx_face_de_lelem)
                {
                  Cerr << "Face is not found!! " << finl;
                  Process::exit();
                }
              const int num_face_den_face = elem_faces(elem, (idx_face_de_lelem + Objet_U::dimension) % nb_faces);
              // const int elem_voisin = face_voisins(num_face_den_face,0)+face_voisins(num_face_den_face,1)-elem;
              // const double d_elem = dist[elem];
              // const double d_voisin = dist[elem_voisin];
              // const double delta = zvdf.dist_face(num_face, num_face_den_face , zvdf.orientation(num_face));
              // Il faudrait lineariser a l'ordre 2 car cela sert aussi au calcul de laplacien_d.
              // On extrapole le gradient mais on ne peut pas faire mieux car on a un petit stencil a 2 points sur dist.
              for (int c = 0; c < u0.line_size(); c++)
                u0(num_face, c) = u0(num_face_den_face, c);              //*(1.+xxxxx);
            }
        }
    }
}

int get_num_face_den_face(const int num_face, const int elem, const IntTab& elem_faces)
{
  const int nb_faces = elem_faces.dimension(1);
  int idx_face_de_lelem = 0;
  for (idx_face_de_lelem = 0; idx_face_de_lelem < nb_faces; idx_face_de_lelem++)
    {
      if (elem_faces(elem, idx_face_de_lelem) == num_face)
        break; // Face found
    }
  if (nb_faces == idx_face_de_lelem)
    {
      Cerr << "Face is not found!! " << finl;
      Process::exit();
    }
  const int num_face_den_face = elem_faces(elem, (idx_face_de_lelem + Objet_U::dimension) % nb_faces);
  return num_face_den_face;
}

/*! @brief Calcul du saut de vitesse a l'interface du au changement de phase
 *
 *   phase_pilote = -1: u-u0 = champ de vitesse de deplacement de l'interface
 *   phase_pilote = 0 : u+u0 = champ de vitesse de la phase 0
 *   phase_pilote = 1 : u+u0 = champ de vitesse de la phase 1
 *   ordre = 0 : pas de prise en compte de la correction en courbure
 *   ordre = 1 : prise en compte de la correction en courbure a l'ordre 1
 *   ordre = 2 : prise en compte de la correction en courbure a l'ordre 2
 *
 */
void Navier_Stokes_FT_Disc::calculer_delta_u_interface(Champ_base& champ_u0, int phase_pilote, int ordre)
{
  //  static const Stat_Counter_Id count = statistiques().new_counter(1, "calculer_delta_u_interface", 0);
  //  statistiques().begin_count(count);
  DoubleTab& u0 = champ_u0.valeurs();
  const Fluide_Diphasique& fluide_dipha = fluide_diphasique();
  const Fluide_Incompressible& phase_0 = fluide_dipha.fluide_phase(0);
  const Fluide_Incompressible& phase_1 = fluide_dipha.fluide_phase(1);
  const DoubleTab& tab_rho_phase_0 = phase_0.masse_volumique()->valeurs();
  const DoubleTab& tab_rho_phase_1 = phase_1.masse_volumique()->valeurs();
  const double rho_0 = tab_rho_phase_0(0, 0);
  const double rho_1 = tab_rho_phase_1(0, 0);
  //const double delta_un_sur_rho = 1. / rho_1 - 1. / rho_0;
  const double un_sur_rho_0 = 1. / rho_0;
  const double un_sur_rho_1 = 1. / rho_1;

  REF(Transport_Interfaces_FT_Disc) &refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
  const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();

  // Pour permettre de calculer mpoint, mais sans l'utiliser pour le deplacement de l'interface,
  // il suffit de ne pas mettre le mot cle " equation_temperature_mpoint temp" dans le jdd.
  const int nn = variables_internes().second_membre_projection->valeurs().dimension(0); // nombre d'elems
  DoubleTab mpoint;
  mpoint.resize(nn);
  mpoint = 0.; //pour initialiser
  if (variables_internes().ref_equation_mpoint_.non_nul())
    {
      const DoubleTab& mp = variables_internes().ref_equation_mpoint_->get_mpoint();
      // Si inactif, on ne prend pas en compte sa contribution dans le calcul de delta_u:
      if (!variables_internes().mpoint_inactif)
        mpoint = mp;
    }
  if (variables_internes().ref_equation_mpoint_vap_.non_nul())
    {
      const DoubleTab& mpv = variables_internes().ref_equation_mpoint_vap_->get_mpoint();
      // Si inactif, on ne prend pas en compte sa contribution dans le calcul de delta_u:
      if (!variables_internes().mpointv_inactif)
        mpoint += mpv;
    }

  // GB2023.10.10 : I don't understand why I did a distinction only on phase_pilote == 1 (should be the same when it's phase_pilote == 0
  //                for instance in the convection of temperature in the vapour...
  //        BESIDES, the "switch (phase_pilote)" makes no sense if only one phase_pilote is used.
  if ((variables_internes().new_mass_source_) && (phase_pilote != 1))
    {
      const DoubleTab& normale_elements = eq_transport.get_update_normale_interface().valeurs();
      const DoubleTab& interfacial_area = variables_internes().ai->valeurs();

      const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis());
      const IntTab& face_voisins = domaine_vf.face_voisins();
      const int nb_faces = face_voisins.dimension(0);
      const int dim = Objet_U::dimension;
      const DoubleTab& xp = domaine_vf.xp(); // centres de gravite des elements
      const DoubleTab& xv = domaine_vf.xv(); // centres de gravite des faces.
      u0 = 0.;
      if (u0.line_size() == dim) // vef
        {
          Cerr << "Using option new_mass_source is not possible yet in VEF. Contact us. " << finl;
          Process::exit();
          for (int face = 0; face < nb_faces; face++)
            for (int j = 0; j < dim; j++)
              {
                double x = 0.; // a coder...
                u0(face, j) *= x;
              }
        }
      else
        {
          const Domaine_VDF& zvdf = ref_cast(Domaine_VDF, domaine_dis());
          const IntVect& orientation = zvdf.orientation();
          for (int face = 0; face < nb_faces; face++)
            {
              const int dir = orientation[face];
              const double surface = domaine_vf.face_surfaces(face);
              const int e1 = face_voisins(face, 0);
              const int e2 = face_voisins(face, 1);
              const double xf = xv(face, dir);
              if (Objet_U::bidim_axi && fabs(xf) <= DMINFLOAT && (dir == 0))
                {
                  // We are on the symmetry axis => surface = 0.;
                  u0(face) = 0.; // there is no normal velocity at a symmetry axis
                  continue;
                }
              //double x = 0.;
              double xx = 0.;
              // Si on n'est pas au bord...
              if (e1 >= 0)
                {
                  const double nx = normale_elements(e1, dir);
                  //x = c*secmem2[e1]*nx;
                  const double ai = interfacial_area(e1);
                  // nx pointe vers le liquide (sortant de phase 0)
                  if ((fabs(ai) > DMINFLOAT) && (fabs(nx) > DMINFLOAT))
                    {
                      // distance positive on the vapour side (chi_0 = 1, ie indicatrice = 0)
                      const double d = (xf - xp(e1, dir)) * nx;
                      switch(phase_pilote)
                        {
                        case -1:
                          {
                            // Champ de vitesse tel que u + u0 soit continu et egal a la
                            // vitesse de deplacement de l'interface
                            const double un_sur_rho = (d > 0.) ? un_sur_rho_0 : un_sur_rho_1;
                            xx = ai / surface * un_sur_rho * mpoint[e1] * nx;
                            //Cerr << "diff " << x << " " << xx << finl;
                            break;
                          }
                        case 0:
                          {
                            // Champ de vitesse tel que u + u0 soit continu et
                            // u+u0 = la vitesse de la phase 0 dans la phase 0
                            const double p = (d > 0.) ? 0. : (un_sur_rho_1 - un_sur_rho_0);
                            xx = ai / surface * p * mpoint[e1] * nx;
                            //Cerr << "face " << face << " " << xx << finl;
                            break;
                          }
                        case 1:
                          {
                            // Champ de vitesse tel que u + u0 soit continu et
                            // u+u0 = la vitesse de la phase 1 dans la phase 1
                            const double p = (d > 0.) ? (un_sur_rho_0 - un_sur_rho_1) : 0.;
                            xx = ai / surface * p * mpoint[e1] * nx;
                            //Cerr << "face " << face << " " << xx << finl;
                            break;
                          }
                        default:
                          Cerr << "Error for the method Navier_Stokes_FT_Disc::calculer_delta_u_interface phase_pilote" << finl;
                          Process::exit();
                        }
                    }
                }
              // We ADD contribution of e2 if not a boundary to xx
              if (e2 >= 0)
                {
                  const double nx = normale_elements(e2, dir);
                  //x += c*secmem2[e2]*normale_elements(e2, dir);
                  const double ai = interfacial_area(e2);
                  if ((fabs(ai) > DMINFLOAT) && (fabs(nx) > DMINFLOAT))
                    {
                      // distance positive on the vapour side (chi_0 = 1, ie indicatrice = 0)
                      const double d = (xf - xp(e2, dir)) * nx;
                      switch(phase_pilote)
                        {
                        case -1:
                          {
                            // Champ de vitesse tel que u + u0 soit continu et egal a la
                            // vitesse de deplacement de l'interface
                            const double un_sur_rho = (d > 0.) ? un_sur_rho_0 : un_sur_rho_1;
                            xx += ai / surface * un_sur_rho * mpoint[e2] * nx;
                            //Cerr << "diff2 " << x << " " << xx << finl;
                            break;
                          }
                        case 0:
                          {
                            // Champ de vitesse tel que u + u0 soit continu et
                            // u+u0 = la vitesse de la phase 0 dans la phase 0
                            const double p = (d > 0.) ? 0. : (un_sur_rho_1 - un_sur_rho_0);
                            xx += ai / surface * p * mpoint[e2] * nx;
                            //Cerr << "face2 " << face << " " << xx << finl;
                            break;
                          }
                        case 1:
                          {
                            // Champ de vitesse tel que u + u0 soit continu et
                            // u+u0 = la vitesse de la phase 1 dans la phase 1
                            const double p = (d > 0.) ? (un_sur_rho_0 - un_sur_rho_1) : 0.;
                            xx += ai / surface * p * mpoint[e2] * nx;
                            //Cerr << "face2 " << face << " " << xx << finl;
                            break;
                          }
                        default:
                          Cerr << "Error for the method Navier_Stokes_FT_Disc::calculer_delta_u_interface phase_pilote" << finl;
                          Process::exit();
                        }
                    }

                }
              u0(face) = xx;
            }

          // GB 2023.10.10. On inclined interfaces, we realised that some configurations (interface topologies)
          //                can lead to using faces where the velocity jump delta_u has not been extended
          //                for the velocity interpolation at the markers. This leads to a misprediction of delta_u_i
          u0.echange_espace_virtuel();

          //  Pour parcourir les elements qui sont coupes par une facette "facette":
          const Maillage_FT_Disc& mesh = eq_transport.maillage_interface();
          const Intersections_Elem_Facettes& intersections = mesh.intersections_elem_facettes();
          const ArrOfInt& index_facette = intersections.index_facette();
          const Domaine_VF& zvf = ref_cast(Domaine_VF, domaine_dis());
          const IntTab& elem_faces = zvf.elem_faces();
          const int nb_faces_per_elem = elem_faces.line_size();
          for (int facette = 0; facette < mesh.nb_facettes(); facette++)
            {
              int index = index_facette[facette];
              if (index >= 0)
                {
                  const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
                  const int elem = data.numero_element_;
                  // elem is mixed (crossed by an interface)
                  // loop on its faces to find a neighbour:
                  for (int idx_face_de_lelem = 0; idx_face_de_lelem < nb_faces_per_elem; idx_face_de_lelem++)
                    {
                      const int face_of_elem = elem_faces(elem, idx_face_de_lelem);
                      // The neighbour :
                      const int e1 = face_voisins(face_of_elem, 0) + face_voisins(face_of_elem, 1) - elem;
                      if ((e1 >= 0) && (fabs(interfacial_area(e1)) < DMINFLOAT))
                        {
                          // e1 exists (ie face is not a boundary) and is pure
                          // face_of_elem is between elem and e1
                          const int num_face_den_face = get_num_face_den_face(face_of_elem, e1, elem_faces);
                          if (fabs(u0(num_face_den_face)) < DMINFLOAT)
                            {
                              // There is no velocity on the face in front of "face_of_e1" (which is adjacent to a mixed element)
                              // we should extrapolate the other one there.
                              // (ie we assume that this other face was zero because the 2nd neighbour (rang 2) is pure too.
                              u0(num_face_den_face) = u0(face_of_elem);
                            }
                        }
                    }
                }
            }
        }
      u0.echange_espace_virtuel();
      return;
    }

  // Distance a l'interface discretisee aux elements:
  const DoubleTab& dist = eq_transport.get_update_distance_interface().valeurs();
  DoubleTab phi = calculer_div_normale_interface().valeurs();
  {
    const int n = phi.dimension(0);
    for (int i = 0; i < n; i++)
      {
        double d = dist(i);
        double p = 0.;
        if (d >= -1e20)
          {
            const double div_n = phi(i);
            // Distance calculee pour cet element ?
            // Calcul de la fonction phi pour cet element :
            const double mp = mpoint[i];
            switch(ordre)
              {
              case 0:
                // Pas de prise en compte de la correction en courbure
                p = d * mp;
                break;
              case 1:
                //  Prise en compte de la correction en courbure a l'ordre 1
                p = d * (1. - 0.5 * div_n * d) * mp;
                break;
              case 2:
                //  Prise en compte de la correction en courbure a l'ordre 2
                p = d * (1. - 0.5 * div_n * d + div_n * div_n * d * d / 6.) * mp;
                break;
              default:
                Cerr << "Error for the method Navier_Stokes_FT_Disc::calculer_delta_u_interface ordre" << finl;
                Process::exit();
              }
            switch(phase_pilote)
              {
              case -1:
                // Champ de vitesse tel que u + u0 soit continu et egal a la
                // vitesse de deplacement de l'interface
                if (d < 0)
                  p *= un_sur_rho_0;
                else
                  p *= un_sur_rho_1;
                break;
              case 0:
                // Champ de vitesse tel que u + u0 soit continu et
                // u+u0 = la vitesse de la phase 0 dans la phase 0
                if (d < 0)
                  p = 0.; // dans la phase 0
                else
                  p *= (un_sur_rho_0 - un_sur_rho_1); // GB BugFix 2020/10/09
                break;
              case 1:
                // Champ de vitesse tel que u + u0 soit continu et
                // u+u0 = la vitesse de la phase 1 dans la phase 1
                if (d < 0)
                  p *= (un_sur_rho_1 - un_sur_rho_0); // GB BugFix 2020/10/09
                else
                  p = 0.; // dans la phase 1
                break;
              default:
                Cerr << "Error for the method Navier_Stokes_FT_Disc::calculer_delta_u_interface phase_pilote" << finl;
                Process::exit();
              }
          }
        phi(i) = p;
      }
    phi.echange_espace_virtuel();
  }

  // Gradient de phi:
  if (champ_u0.que_suis_je() == "Champ_Face")
    {
      gradient.calculer(phi, u0);
      correct_at_exit_bad_gradient(u0);
    }
  else
    {
      Cerr << "Error for the method Navier_Stokes_FT_Disc::calculer_delta_u_interface\n" << " Non code pour " << champ_u0.que_suis_je() << finl;
      Process::exit();
    }

  // On annule la vitesse calculee pour les faces adjacentes a un element
  // invalide.
  {
    const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis());
    const IntTab& face_voisins = domaine_vf.face_voisins();
    const int nb_faces = domaine_vf.nb_faces();
    ;
    for (int i = 0; i < nb_faces; i++)
      {
        for (int j = 0; j < 2; j++)
          {
            const int elem = face_voisins(i, j);
            if (elem >= 0 && dist(elem) < -1e20)
              {
                u0(i) = 0.;
                break;
              }
          }
      }
  }
  u0.echange_espace_virtuel();
  solveur_masse->appliquer(u0);
}

double calculer_indicatrice_face_privilegie_pure(const DoubleTab& indicatrice, const IntTab& face_voisins, const int num_face)
{
  const int elem0 = face_voisins(num_face, 0);
  const int elem1 = face_voisins(num_face, 1);
  double indic_0 = (elem0 >= 0) ? indicatrice[elem0] : indicatrice[elem1];
  double indic_1 = (elem1 >= 0) ? indicatrice[elem1] : indicatrice[elem0];
  // Decentrage : si la face est adjacente a une face monophasique,
  //  prendre la phase pure pour toute la face:
  if (indic_0 == 0. || indic_0 == 1.)
    indic_1 = indic_0;
  if (indic_1 == 0. || indic_1 == 1.)
    indic_0 = indic_1;
  const double indic_face = (indic_0 + indic_1) * 0.5;
  return indic_face;
}

double calculer_indicatrice_face_based_on_ai(const DoubleTab& indicatrice, const DoubleTab& indicatrice_faces, const IntTab& face_voisins, const Domaine_VF& domaine_vf,
                                             const DoubleTab& interfacial_area, const DoubleTab& normale_elements, const int face, const int dim)
{
  double indic_face = 0.;
  int v;
  for (v = 0; v < 2; v++)
    {
      const int elem = face_voisins(face, v);
      if (elem >= 0)
        {
          // If a neighbour is pure, we use that value at the face and stop further calculation.
          const double indic = indicatrice[elem]; // This is the value of chi_1 (ie =1 in phase 1!)
          //if (indic == 0. || indic == 1.)
          if (indic <= 5e-3 || indic >= 1. - 5e-3)
            {
              indic_face = 1 - indic; // We want chi of phase_0
              break;
            }
          else
            {
              const double surface = domaine_vf.face_surfaces(face);
              //const DoubleVect& volumes = domaine_vf.volumes();
              const double ai = interfacial_area(elem); // nx pointe vers le liquide (sortant de phase 0)
              if (fabs(ai) > DMINFLOAT)
                {
                  double x = 0.;
                  if (dim > 1) // dim==1 en VDF; sinon VEF
                    {
                      for (int j = 0; j < dim; j++)
                        {
                          const double dSj = domaine_vf.face_normales(face, j);
                          const double nx = normale_elements(elem, j);
                          // produit scalaire :
                          x += dSj * nx;
                          x *= ai / surface;
                          // Que/comment Choisir?
                          indic_face += 1 - x;
                          // indic_face += x;
                          Cerr << "Never tested. To be verified. It should depend on a scalar product with the vect (xp-xv)" << finl;
                          Process::exit();
                        }
                    }
                  else
                    {
                      // En VDF, l'acces a orientation permet d'eviter le calcul du produit scalaire.
                      const Domaine_VDF& zvdf = ref_cast(Domaine_VDF, domaine_vf);
                      const IntVect& orientation = zvdf.orientation();
                      const int dir = orientation[face];
                      const double nx = normale_elements(elem, dir);
                      // Assumes a cube, nx larger than diag means we can use the method rather safely
                      if (nx > 0.707)
                        {
                          x = ai / surface * nx;
                          // On suppose que v0 est a gauche et v1 a droite!!!
                          if (v == 0)
                            indic_face += 1 - x; // This way, we build chi_0 because normale points towards chi_1
                          else
                            indic_face += x;
                        }
                      else
                        {
                          // L'interface croise probablement la face d'en face et la methode ne marche plus.
                          // on revient a la methode classique :
                          indic_face += (1. - indicatrice_faces[face]);
                        }
                    }
                }
              else
                {
                  Cerr << " WTF, c'est impossible" << finl;
                  Process::exit();
                }
            }
        }
      else
        {
          // The only neighbour to the face :
          const int elem_voisin = face_voisins(face, 1 - v); // The other one is accessed by 1-v
          const double indic = indicatrice[elem_voisin]; // This is the value of chi_1 (ie =1 in phase 1!)
          indic_face = 1 - indic; // We want chi of phase_0
          break; // c'est important pour le if d'apres.
        }
    }
  if (v == 2)
    // On n'a pas touche le break, on est donc passe 2 fois. donc :
    indic_face *= 0.5;

  // assert((indic_face >=0) && (indic_face<=1.));
  // ca arrive des petits derapages..
  if (indic_face < 0)
    indic_face = 0.;
  if (indic_face > 1.)
    indic_face = 1.;
  return indic_face;
}

void correct_indicatrice_face_bord(const int num_face, const Maillage_FT_Disc& maillage, const Domaine_VF& zvf, const IntTab& face_voisins, const DoubleTab& indicatrice, const bool privilegie_pure,
                                   double& indic_face)
{
  // Correction de l'indicatrice face :
  const int n0 = face_voisins(num_face, 0);
  const int n1 = face_voisins(num_face, 1);
  if ((n0 == -1) or (n1 == -1))
    {
      // On a boundary face
      const int outward_normal = (n0 == -1) ? -1 : 1;
      const int elem = n0 + n1 + 1;
      const double indic_ghost = compute_indic_ghost(elem, num_face, indicatrice(elem), outward_normal, zvf, maillage);
      indic_face = indic_ghost;
      if ((privilegie_pure) && (est_egal(indicatrice(elem), 1.) || est_egal(indicatrice(elem), 0.)))
        {
          indic_face = indicatrice(elem);
        }
    }
}

// Calcul de l'integrale de dI_dt sur chaque element du maillage.
// Le tableau dI_dt doit avoir la bonne structure. L'espace virtuel est
// mis a jour. La method n'est plus const a cause des options
// INTERP_MODIFIEE et AI_BASED qui recalculent indicatrice_faces.
void Navier_Stokes_FT_Disc::calculer_dI_dt(DoubleVect& dI_dt) //const
{
  const double rho_0 = fluide_diphasique().fluide_phase(0).masse_volumique()->valeurs()(0, 0);
  const double rho_1 = fluide_diphasique().fluide_phase(1).masse_volumique()->valeurs()(0, 0);
  const double delta_rho = rho_0 - rho_1;

  double rho_0_sur_delta_rho = 0.;
  if (delta_rho != 0)
    rho_0_sur_delta_rho = rho_0 / delta_rho;

  const DoubleTab& tab_vitesse = inconnue()->valeurs();
  const IntTab& face_voisins = domaine_dis().face_voisins();

  REF(Transport_Interfaces_FT_Disc) &refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
  const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
  const DoubleTab& indicatrice = eq_transport.inconnue()->valeurs();
  const Maillage_FT_Disc& maillage = eq_transport.maillage_interface();
  const Domaine_VF& domVF = ref_cast(Domaine_VF, domaine_dis());
  //const IntVect& orientation = ref_cast(Domaine_VF, domaine_dis()).orientation();
  //  const DoubleTab& indicatrice = variables_internes().ref_eq_interf_proprietes_fluide->inconnue()->valeurs();
  DoubleTab tmp(tab_vitesse); // copie du tableau des vitesses de ns
  const int dim = tab_vitesse.line_size();
  const int n = tab_vitesse.dimension(0);

  // On cree un tableau avec la meme structure que la pression
  DoubleTab resu;
  resu.copy(variables_internes().second_membre_projection->valeurs(), RESIZE_OPTIONS::NOCOPY_NOINIT);
  resu = 0.;

  //On utilise un operateur de divergence temporaire et pas celui porte par l equation
  //pour ne pas modifier les flux_bords_ rempli au cours de Navier_Stokes_std::mettre_a_jour
  Operateur_Div div_tmp;
  div_tmp.associer_eqn(*this);
  div_tmp.typer();
  div_tmp.l_op_base().associer_eqn(*this);
  div_tmp->completer();

  if ((variables_internes_.type_interpol_indic_pour_dI_dt_ == Navier_Stokes_FT_Disc_interne::INTERP_STANDARD_UVEXT)
      || (variables_internes_.type_interpol_indic_pour_dI_dt_ == Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE_UVEXT)
      || (variables_internes_.type_interpol_indic_pour_dI_dt_ == Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED_UVEXT))
    {
      // Avec changement de phase, on veut reconstruire u_vap (ie phase 0)
      // Prise en compte du terme source div(u) du changement de phase
      if (variables_internes().ref_equation_mpoint_.non_nul() || variables_internes().ref_equation_mpoint_vap_.non_nul())
        {
          calculer_delta_u_interface(variables_internes().vitesse_jump0_, 0, variables_internes().correction_courbure_ordre_ /* ordre de la correction en courbure */);
          // u+u0 = champ de vitesse de la phase 0
          tmp += variables_internes().vitesse_jump0_->valeurs();
        }
    }
  if ((variables_internes_.type_interpol_indic_pour_dI_dt_ == Navier_Stokes_FT_Disc_interne::INTERP_STANDARD_UIEXT)
      || (variables_internes_.type_interpol_indic_pour_dI_dt_ == Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE_UIEXT)
      || (variables_internes_.type_interpol_indic_pour_dI_dt_ == Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED_UIEXT))
    {
      // Reconstruction d'un champ de vitesse interfaciale (ie phase 1)
      // Prise en compte du terme source div(u) du changement de phase
      if (variables_internes().ref_equation_mpoint_.non_nul() || variables_internes().ref_equation_mpoint_vap_.non_nul())
        {
          calculer_delta_u_interface(variables_internes().delta_u_interface, -1, variables_internes().correction_courbure_ordre_ /* ordre de la correction en courbure */);
          // u-dui = champ de vitesse d'interface
          tmp -= variables_internes().delta_u_interface->valeurs();

          // Question: il y a un assert_espace_virtuel_vect dans divergence.calculer,
          //  mais l'operateur n'a normalement pas besoin de l'espace virtuel !
          //  La ligne suivante devrait pouvoir etre retiree:
          tmp.echange_espace_virtuel();

          div_tmp->calculer(tmp, resu);
          for (int i = 0; i < resu.size_array(); i++)
            resu[i] *= indicatrice[i];
        }
    }

  const bool ghost_correction = (variables_internes_.OutletCorrection_pour_dI_dt_ == Navier_Stokes_FT_Disc_interne::CORRECTION_GHOST_INDIC);
  switch(variables_internes_.type_interpol_indic_pour_dI_dt_)
    {
    case Navier_Stokes_FT_Disc_interne::INTERP_STANDARD:
    case Navier_Stokes_FT_Disc_interne::INTERP_STANDARD_UVEXT:
      {
        for (int i = 0; i < n; i++)
          {
            double indic_face = calculer_indicatrice_face_privilegie_pure(indicatrice, face_voisins, i);
            if (ghost_correction)
              correct_indicatrice_face_bord(i, maillage, domVF, face_voisins, indicatrice, true /* privilegie_pure */, indic_face);
            const double x = rho_0_sur_delta_rho - indic_face;
            for (int j = 0; j < dim; j++)
              tmp(i, j) *= x;
          }
        break;
      }
    case Navier_Stokes_FT_Disc_interne::INTERP_STANDARD_UIEXT:
      {
        for (int i = 0; i < n; i++)
          {
            double indic_face = calculer_indicatrice_face_privilegie_pure(indicatrice, face_voisins, i);
            if (ghost_correction)
              correct_indicatrice_face_bord(i, maillage, domVF, face_voisins, indicatrice, true /* privilegie_pure */, indic_face);
            const double x = -indic_face;
            for (int j = 0; j < dim; j++)
              tmp(i, j) *= x;
          }
        break;
      }
    case Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE:
    case Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE_UVEXT:
      {
        const DoubleTab& indicatrice_faces = refeq_transport->get_compute_indicatrice_faces().valeurs();
        for (int i = 0; i < n; i++)
          {
            double indic_face = indicatrice_faces(i);
            if (ghost_correction)
              correct_indicatrice_face_bord(i, maillage, domVF, face_voisins, indicatrice, false, indic_face);
            const double x = rho_0_sur_delta_rho - indic_face;
            for (int j = 0; j < dim; j++)
              tmp(i, j) *= x;
          }
        break;
      }
    case Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE_UIEXT:
      {
        const DoubleTab& indicatrice_faces = refeq_transport->get_compute_indicatrice_faces().valeurs();
        for (int i = 0; i < n; i++)
          {
            double indic_face = indicatrice_faces(i);
            if (ghost_correction)
              correct_indicatrice_face_bord(i, maillage, domVF, face_voisins, indicatrice, false, indic_face);
            const double x = -indic_face;
            for (int j = 0; j < dim; j++)
              tmp(i, j) *= x;
          }
        break;
      }
    case Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED:
    case Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED_UVEXT:
    case Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED_UIEXT:
      {
        const DoubleTab& indicatrice_faces = refeq_transport->get_compute_indicatrice_faces().valeurs();
        const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis());
        if (Process::je_suis_maitre())
          Cerr << " The interpolation of indicatrice to faces in calculer_dI_dt is based on the interfacial area" << " and on the normal to the interface." << finl;

        const DoubleTab& normale_elements = eq_transport.get_update_normale_interface().valeurs();
        const DoubleTab& interfacial_area = variables_internes().ai->valeurs();

#if NS_VERBOSE
        const DoubleTab& xp = domaine_vf.xp(); // centres de gravite des elements
        const DoubleTab& xv = domaine_vf.xv(); // centres de gravite des faces.
#endif
        // On fait la moyenne des 2 valeurs calculees sur les voisins
        // ATTENTION, ici on veut la valeur de chiv (cad chi_0) a la face.
        for (int face = 0; face < n; face++)
          {
            double indic_face = calculer_indicatrice_face_based_on_ai(indicatrice, indicatrice_faces, face_voisins, domaine_vf, interfacial_area, normale_elements, face, dim);
            if (ghost_correction)
              correct_indicatrice_face_bord(face, maillage, domVF, face_voisins, indicatrice, false, indic_face);

#if NS_VERBOSE
            const double val = 1.-indicatrice_faces[face]; // indicatrice_faces=chi_1  whereas indic_face=chi_0 !!! Hence, the "1.-"
            if (fabs(indic_face-val)>DMINFLOAT)
              {
                const int elem0 = face_voisins(face, 0);
                const int elem1 = face_voisins(face, 1);
                double indic_0 = (elem0 >= 0) ? indicatrice[elem0] : indicatrice[elem1];
                double indic_1 = (elem1 >= 0) ? indicatrice[elem1] : indicatrice[elem0];
                if (elem0>=0)
                  Cerr << "xp0["<< elem0<< "]: " << xp(elem0, 0) << " "  << xp(elem0, 1) << " "  << finl;
                else
                  Cerr << "xp0: bord!" <<finl;
                if (elem1>=0)
                  Cerr << "xp1["<< elem1<< "]: " << xp(elem1, 0) << " "  << xp(elem1, 1) << " "  << finl;
                else
                  Cerr << "xp1: bord!" <<finl;

                Cerr << "xv: " << xv(face, 0) << " "  << xv(face, 1) << " "  << finl;
                Cerr << "voisins (ou ghost): " << indic_0 << " " << indic_1 << finl;
                Cerr << "GB whats up?face="<< face <<" indic / val / diff " << indic_face << " " << val << " " << indic_face-val << finl;
              }
#endif
            // chi_v * u_v_ext
            for (int j = 0; j < dim; j++)
              tmp(face, j) *= indic_face;
          }
        break;
      }
    default:
      Cerr << " Navier_Stokes_FT_Disc::calculer_dI_dt \n" << " unknown case?" << finl;
      Process::exit();
    }

  if (variables_internes_.type_interpol_indic_pour_dI_dt_ == Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED_UIEXT)
    tmp *= -1;

  // Question: il y a un assert_espace_virtuel_vect dans divergence.calculer,
  //  mais l'operateur n'a normalement pas besoin de l'espace virtuel !
  //  La ligne suivante devrait pouvoir etre retiree:
  tmp.echange_espace_virtuel();

  div_tmp->ajouter(tmp, resu);

  // Correction des flux bords :
  // resu = int_V div(tmp) dv    avec:   tmp = (rho_0_sur_delta_rho - indic_face)*vitesse_ns
  // 		(a voir qui est indic_face -> celle de phase1?)
  // Les flux aux bords dans les mailles diphasiques avec sortie libre sont mal calcules si l'indic_phase l'est mal.
  // (car il est difficile de savoir si l'indic_face est en faite encore pure et qu'il ne sort que d'une phase
  //  dans certaines mailles mixtes).
  // Or, les flux_bord ont deja affecte resu par l'operation :
  //   resu += flux // ...
  // Plusieurs options sont testees :
  //   - On fait rien (NO_CORRECTION)
  //   - On fait la correction de l'indic_face avant (avant le calcul de tmp qui est en entree de la div) (CORRECTION_GHOST_INDIC)
  //   - On suppose que globalement, les cellules touchant le bord sortie sont
  // a flux de masse phasique nulle (ce qui rentre sort). (ZERO_NET_FLUX_ON_MIXED_CELLS)
  // Cela signifie que le flux sur la face de bord est l'oppose des flux internes.
  // Concretement, il suffit de faire "resu(elem) =0." dans les mailles concernees.
  //   - On recalcule le flux sur la face de bord des mailles mixtes et on retranche celui qui avait ete ajoute (ZERO_OUT_FLUX_ON_MIXED_CELLS)
  // Cela ne laisse pas sortir la vapeur, ce n'est donc pas bon.
  switch(variables_internes_.OutletCorrection_pour_dI_dt_)
    {
    case Navier_Stokes_FT_Disc_interne::NO_CORRECTION:
    case Navier_Stokes_FT_Disc_interne::CORRECTION_GHOST_INDIC:
      {
        break;
      }
    case Navier_Stokes_FT_Disc_interne::ZERO_NET_FLUX_ON_MIXED_CELLS:
    case Navier_Stokes_FT_Disc_interne::ZERO_OUT_FLUX_ON_MIXED_CELLS:
      {
        const Domaine_VDF& zvdf = ref_cast(Domaine_VDF, domaine_dis());
        const Domaine_Cl_dis_base& zcldis = domaine_Cl_dis();
        const Domaine_Cl_VDF& zclvdf = ref_cast(Domaine_Cl_VDF, zcldis);
        const DoubleVect& face_surfaces = zvdf.face_surfaces();
        // Boucle sur les bords pour traiter les conditions aux limites
        int ndeb, nfin;
        for (int n_bord = 0; n_bord < zvdf.nb_front_Cl(); n_bord++)
          {
            // pour chaque Condition Limite on regarde son type
            // Si face de Dirichlet ou de Symetrie on ne fait rien
            // Si face de Neumann on calcule la contribution au terme source
            const Cond_lim& la_cl = zclvdf.les_conditions_limites(n_bord);
            //Cerr << que_suis_je() << "::calculer_dI_dt() correction du dIdt a la CL : " <<  la_cl.valeur() << finl;
            const Front_VF& le_bord = ref_cast(Front_VF, la_cl->frontiere_dis());
            ndeb = le_bord.num_premiere_face();
            nfin = ndeb + le_bord.nb_faces();

            // Correction sortie libre :
            // On ne sait pas bien calculer la correction de volume liee a dIdt a la sortie libre.
            // On prefere donc l'annuler dans cet element:
            if ((sub_type(Dirichlet, la_cl.valeur())) || (sub_type(Neumann_sortie_libre, la_cl.valeur())) || (sub_type(Dirichlet_homogene, la_cl.valeur())))
              {
                for (int num_face = ndeb; num_face < nfin; num_face++)
                  {
                    const int n0 = face_voisins(num_face, 0);
                    const int n1 = face_voisins(num_face, 1);
                    const int elem = n0 + n1 + 1;
                    const double indic = indicatrice(elem);
                    if (indic * (1 - indic) > 1e-6) // In a mixed cell!
                      {
                        const double coef = face_surfaces(num_face); //*porosite_surf(num_face);
                        if (variables_internes_.OutletCorrection_pour_dI_dt_ == Navier_Stokes_FT_Disc_interne::ZERO_NET_FLUX_ON_MIXED_CELLS)
                          resu(elem) = 0.;
                        else if (variables_internes_.OutletCorrection_pour_dI_dt_ == Navier_Stokes_FT_Disc_interne::ZERO_OUT_FLUX_ON_MIXED_CELLS)
                          {
                            for (int j = 0; j < dim; j++)
                              {
                                const int outward_normal = (n0 == -1) ? -1 : 1;
                                double flux_bord_calcule_par_operateur_div = 0.;
                                flux_bord_calcule_par_operateur_div = outward_normal * coef * tmp(num_face, j);
                                resu(elem) -= flux_bord_calcule_par_operateur_div; // On RETRANCHE le flux qui avait ete pris a l'etape de calcul du div
                              }
                          }
                        else
                          {
                            Process::exit();
                          }

                      }
                  }
              }
          }
        // Fin de la boucle for
        break;
      }
    default:
      Cerr << "unexpected" << finl;
      Process::exit();
    }

  // Extraction des valeurs
  const int nb_elem = domaine_dis().nb_elem();
  assert(nb_elem == dI_dt.size());

  // Simple copie
  dI_dt.inject_array(resu, nb_elem);

  // L'integrale de div sur l'element est dimension * la valeur aux elements
  //  renvoyee par l'operateur.
  // Les valeurs aux elements sont au debut du tableau resu.
  if (tab_vitesse.line_size() > 1) // i.e. VEF
    dI_dt *= dimension;

#if NS_VERBOSE
  {
    Cerr << "[BEFORE-PCH] Locally, the maximum of dI_dt is : " << dI_dt.mp_max_abs_vect() << finl;
    const double temps = schema_temps().temps_courant();
    double sum = 0.;
    for (int i=0; i< nb_elem; i++)
      sum += dI_dt[i];
    Cerr << "[BEFORE-PCH] " << temps << " The sum is : " << sum << " [not valid in //]" << finl;
  }
#endif

  switch(variables_internes_.type_interpol_indic_pour_dI_dt_)
    {
    case Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED:
    case Navier_Stokes_FT_Disc_interne::INTERP_STANDARD_UVEXT:
    case Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE_UVEXT:
    case Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED_UVEXT:
      {
        // dI_dt contains div(chi_v * u_v^ext) * vol_cell (sign to be checked!!)
        // Integrated, this is equal to 0 or to the integral of chi_f entering or leaving the domain through boundaries...
        // It is lacking the phase-change contribution
#define CHECK 0
#if CHECK
        {
          Cerr << "Checking if secmem2 is correctly filled in that case. Check it if you need it..." <<finl;
          if (0) Process::exit();
          DoubleTab& secmem2 = variables_internes().second_membre_projection_jump_.valeurs();

          //if (variables_internes().ref_equation_mpoint_.non_nul())
          //  variables_internes().ref_equation_mpoint_->calculer_mpoint(variables_internes().mpoint.valeur());
          //if (variables_internes().ref_equation_mpoint_vap_.non_nul())
          //  variables_internes().ref_equation_mpoint_vap_->calculer_mpoint(variables_internes().mpoint_vap.valeur());

          const Fluide_Diphasique& fluide = fluide_diphasique();
          const Fluide_Incompressible& phase_0 = fluide.fluide_phase(0);
          const Fluide_Incompressible& phase_1 = fluide.fluide_phase(1);
          const DoubleTab& tab_rho_phase_0 = phase_0.masse_volumique()->valeurs();
          const DoubleTab& tab_rho_phase_1 = phase_1.masse_volumique()->valeurs();
          const double rho_phase_0 = tab_rho_phase_0(0,0);
          const double rho_phase_1 = tab_rho_phase_1(0,0);
          const double jump_inv_rho = 1./rho_phase_1 - 1./rho_phase_0;

          const DoubleTab& interfacial_area = variables_internes().ai->valeurs();
          // Pas une ref, mais un tableau de travail local dans lequel on peut ajouter mointv
          // DoubleTab mpoint = variables_internes().mpoint->valeurs();
          DoubleTab mpoint = variables_internes().ref_equation_mpoint_->get_mpoint();
          if (variables_internes().ref_equation_mpoint_vap_.non_nul())
            {
              //const DoubleTab& mpointv = variables_internes().mpoint_vap->valeurs();
              const DoubleTab& mpointv = variables_internes().ref_equation_mpoint_vap_->get_mpoint();
              for (int elem = 0; elem < nb_elem; elem++)
                mpoint[elem] += mpointv[elem];
            }
          for (int elem = 0; elem < nb_elem; elem++)
            {
              double diff = (secmem2[elem] - jump_inv_rho*interfacial_area[elem]*mpoint[elem]);
              if (fabs(diff)>DMINFLOAT)
                {
                  Cerr << "Problem, secmem2 is not filled properly in that case? "
                  << "elem= " << elem << " diff=" << diff <<finl;
                  Process::exit();
                }
            }
        }
#endif
        if ((variables_internes().ref_equation_mpoint_.non_nul()) && !variables_internes().mpoint_inactif)
          {
            // Is it necessary to recompute them? no
            //variables_internes().ref_equation_mpoint_->calculer_mpoint(variables_internes().mpoint.valeur());
            const double un_sur_rho_0 = 1. / rho_0;

            // We can't use secmem2 because we can't undo jump_inv_rho as it may be 0. !!!
            // DoubleTab& secmem2 = variables_internes().second_membre_projection_jump_.valeurs();
            const DoubleTab& interfacial_area = variables_internes().ai->valeurs();
            // const DoubleVect& volumes = domaine_vf.volumes();
            // Pas une ref, mais un tableau de travail local dans lequel on peut ajouter mointv
            DoubleTab mpoint = variables_internes().ref_equation_mpoint_->get_mpoint();
            if (variables_internes().ref_equation_mpoint_vap_.non_nul())
              {
                // Is it necessary to recompute them? no
                //variables_internes().ref_equation_mpoint_vap_->calculer_mpoint(variables_internes().mpoint_vap.valeur());
                const DoubleTab& mpointv = variables_internes().ref_equation_mpoint_vap_->get_mpoint();
                for (int elem = 0; elem < nb_elem; elem++)
                  mpoint[elem] += mpointv[elem];
              }
#if TCL_MODEL
            // At this point in the algorithm, mpoint has not been augmented by the CL contribution yet, even though
            // the TCL contribution (micro and meso) are already computed and associated to elements.
            // It is done so because just earlier in this method, we call calculer_delta_u_interface
            // to compute a continuous extension of velocity (to displace the interface markers).

            // Adding the TCL contribution to the **local** DoubleTab mpoint
            // (it is thus temporary for the moment and will be done on the real table later at the second call to
            // Triple_Line_Model_FT_Disc::corriger_mpoint)
            // The difference if we correct it here directly is subtile. It is probably fine for the extension of
            // the interface velocity, but the extension for the convective velocity in the temperature has not been done yet.
            // It may cause trouble in the convection then, but the GFM method kind of erase those cells.
            if (probleme_ft().tcl().is_activated())
              {
                Cerr << "[TCL] Contact line model activated in volume correction" << finl;
                probleme_ft().tcl().corriger_mpoint(mpoint);
              }
#endif
            for (int elem = 0; elem < nb_elem; elem++)
              {
                // By convention, mpoint is positive in condensation. Hence, mpoint >0 is responsible for dIv_dt < 0  => a minus sign!
                //                But, \nabla \chi_v = -n_v \delta_i. => another minus sign!
                //                ==> Consequently, it's a "+"
                // Besides, ai = \int_cell delta^i dv => It's homogeneous to the integral, there's no need for an additional "*volumes[elem]"
                //                                       It can be directly summed to the divergence computed before.
                const double x = mpoint[elem] * interfacial_area[elem] * un_sur_rho_0;
                dI_dt[elem] += x;
              }
          }
        break;
      }
    case Navier_Stokes_FT_Disc_interne::INTERP_STANDARD_UIEXT:
    case Navier_Stokes_FT_Disc_interne::INTERP_MODIFIEE_UIEXT:
    case Navier_Stokes_FT_Disc_interne::INTERP_AI_BASED_UIEXT:
      {
        if (probleme_ft().tcl().is_activated())
          {
            const double un_sur_rho_0 = 1. / rho_0;
            const DoubleTab& interfacial_area = variables_internes().ai->valeurs();
            Cerr << "[TCL] Contact line model activated in volume correction" << finl;
            ArrOfInt& tcl_elems = probleme_ft().tcl().elems();
            ArrOfDouble& tcl_mp = probleme_ft().tcl().mp();
            // We accounted for the contribution that is in Ui but not the contribution from TCL yet:
            for (int idx = 0; idx < tcl_elems.size_array(); idx++)
              {
                const int elem = tcl_elems[idx];
                // By convention, mpoint is positive in condensation. Hence, mpoint >0 is responsible for dIv_dt < 0  => a minus sign!
                //                But, \nabla \chi_v = -n_v \delta_i. => another minus sign!
                //                ==> Consequently, it's a "+"
                // Besides, ai = \int_cell delta^i dv => It's homogeneous to the integral, there's no need for an additional "*volumes[elem]"
                //                                       It can be directly summed to the divergence computed before.
                const double x = tcl_mp[idx] * interfacial_area[elem] * un_sur_rho_0;
                dI_dt[elem] += x;
              }
          }
        break;
      }
    default:
      break;
    }

#if NS_VERBOSE
  {
    Cerr << "[AFTER-PCH] Locally, the maximum of dI_dt is : " << dI_dt.mp_max_abs_vect() << finl;
    const double temps = schema_temps().temps_courant();
    double sum = 0.;
    for (int i=0; i< nb_elem; i++)
      sum += dI_dt[i];
    Cerr << "[AFTER-PCH] " << temps << " The sum is : " << sum << " [not valid in //]" << finl;
  }
#endif
  dI_dt.echange_espace_virtuel();
}

// Description:
//  Compute in one Eulerian cell, the average of Front properties in it :
//  normale, bary_facettes_dans_elem and surface_tot are area-weighted averages in the given cell.
//  Warning : normale is not a UNIT vector!
void compute_normale_barycenter_area_in_cell(const int elem, const Maillage_FT_Disc& mesh, Vecteur3& normale, Vecteur3& bary_facettes_dans_elem, double& surface_tot)
{
  const Intersections_Elem_Facettes& intersections = mesh.intersections_elem_facettes();
  const IntTab& facettes = mesh.facettes();
  const DoubleTab& sommets = mesh.sommets();
  const ArrOfDouble& surface_facettes = mesh.get_update_surface_facettes();
  const DoubleTab& normale_facettes = mesh.get_update_normale_facettes();

  surface_tot = 0.;
  normale = 0.;
  bary_facettes_dans_elem = 0.;
  // Get the begining index defining the position in the list index_elem that contains information
  // relative to the current element :
  int index = intersections.index_elem()[elem];
  if (index < 0)
    return; // No facette in this element.

  // Loop over the facettes crossing the element
  while (index >= 0)
    {
      // Accessing the structure containing all the relevant information for facette number fa7
      // Beware, fraction_surface_intersection_ gives the fraction of the facette that is within the current element.
      const Intersections_Elem_Facettes_Data& data = intersections.data_intersection(index);
      const int fa7 = data.numero_facette_;
      const double surface_facette = surface_facettes[fa7];
      const double surf = data.fraction_surface_intersection_ * surface_facette;
      // We compute the (real) coordinates of barycenter of the fraction of facette within the elem
      //    from the Barycentric Coordinates stored in data.barycentre_
      //    (see the web for more information on Barycentric or Areal coordinates)
      //    coord_barycentre_fraction will contain real coordinates (x,y,z) of the barycenter
      Vecteur3 coord_barycentre_fraction(0., 0., 0.);
      for (int dir = 0; dir < 3; dir++)
        {
          const double nx = normale_facettes(fa7, dir);
          normale[dir] += nx * surf;
        }
      for (int isom = 0; isom < 3; isom++)
        {
          const int num_som = facettes(fa7, isom); // numero du sommet dans le tableau sommets
          const double bary_som = data.barycentre_[isom];
          for (int dir = 0; dir < 3; dir++)
            coord_barycentre_fraction[dir] += bary_som * sommets(num_som, dir);

        }
      coord_barycentre_fraction *= surf;
      surface_tot += surf;
      bary_facettes_dans_elem += coord_barycentre_fraction; // This is done for all the 3 components.

      index = data.index_facette_suivante_;
    }

  if (surface_tot > 0.)
    {
      normale *= 1. / surface_tot;
      bary_facettes_dans_elem *= 1. / surface_tot;
    }
  else
    {
      normale = 0.;
      Cerr << " Error in compute_normale_barycenter_area_in_cell (Navier_Stokes_FT_Disc.cpp)." << finl;
      Cerr << "The element " << elem << " only contains facettes of surface=0, so that surface_totale is zero!" << finl;
      Cerr << "What a mess for the barycentre? ..." << finl;
      //    assert(0);
      Process::exit();
      bary_facettes_dans_elem = 0.;
    }
#if NS_VERBOSE
  const double norm =  normale[0]*normale[0] +  normale[1]*normale[1] +  normale[2]*normale[2];
  if (norm<0.9)
    {
      Cerr << " In Navier_Stokes_FT_Disc.cpp compute_normale_barycenter_area_in_cell." << finl;
      Cerr << "Small normal : " << count << " facettes dans l'element " << elem
           << ". surface_tot = "<< surface_tot << "Norm**2 = " << norm << finl;
    }
#endif

}

void Navier_Stokes_FT_Disc::compute_boussinesq_additional_gravity(const Convection_Diffusion_Temperature_FT_Disc& eq, const Fluide_Diphasique& fluide_dipha, const IntTab& face_voisins,
                                                                  const DoubleVect& volumes_entrelaces, const IntVect& orientation, const DoubleTab& indicatrice, const ArrOfDouble& g, // Vect3
                                                                  DoubleTab& gravite_face) const
{
  const int phase_eq = eq.get_phase();
  const DoubleTab& temperature_eq = eq.inconnue()->valeurs();
  const Fluide_Incompressible& fluide_phase_eq = fluide_dipha.fluide_phase(phase_eq);
  const DoubleTab& tab_beta_th_phase_eq = fluide_phase_eq.beta_t()->valeurs();
  const double beta_th_phase_eq = tab_beta_th_phase_eq(0, 0);

  for (int face = 0; face < gravite_face.dimension(0); face++)
    {
      const int elem0 = face_voisins(face, 0);
      const int elem1 = face_voisins(face, 1);
      double coef = 0.;
      // On suppose la ref T0 egale a Tsat
      //
      // Pour les mailles monophasiques, on peut faire simplement l'hypothese que T = chi_k T_k
      // Dans les mailles diphasiques, T = chi_k T_k est une hypothese discutable.
      // Pour les mailles diphasiques, on pourrait envisager une reconstruction plus precise de la
      // temperature monofluide a partir du gradient (cad de mpoint).
      // Neglected in first approximation. we simply compute T = chi_k T_k
      if (elem0 >= 0)
        {
          double chi = (2 * phase_eq - 1) * indicatrice[elem0] + 1 - phase_eq;
          double T_eq = temperature_eq[elem0];
          coef = chi * T_eq;
        }
      if (elem1 >= 0)
        {
          double chi = (2 * phase_eq - 1) * indicatrice[elem1] + 1 - phase_eq;
          double T_eq = temperature_eq[elem1];
          coef += chi * T_eq;
        }
      if (elem0 >= 0 && elem1 >= 0) // Not a boundary of the domain ?
        coef *= 0.5;
      gravite_face(face) -= volumes_entrelaces(face) * g(orientation[face]) * coef * beta_th_phase_eq;
    }
}

/*! @brief Calcul de la derivee en temps de la vitesse.
 *
 */
DoubleTab& Navier_Stokes_FT_Disc::derivee_en_temps_inco(DoubleTab& vpoint)
{
  // Preparation des champs utilises pour le calcul des derivees en temps
  // S'il n'y a pas d'equation de transport des interfaces entre phases fluides,
  // on ne recalcule pas les proprietes.
  {
    REF(Transport_Interfaces_FT_Disc) &refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
    if (refeq_transport.non_nul())
      {
        FT_disc_calculer_champs_rho_mu_nu_dipha(domaine_dis(), fluide_diphasique(), refeq_transport->inconnue()->valeurs(),
                                                // (indicatrice)
                                                champ_rho_elem_->valeurs(), champ_nu_->valeurs(), champ_mu_->valeurs(), champ_rho_faces_->valeurs());
      }
    else
      {
        const Fluide_Incompressible& phase_0 = ref_cast(Fluide_Incompressible, milieu());
        const Domaine_dis_base& zdis = domaine_dis();
        FT_disc_calculer_champs_rho_mu_nu_mono(zdis, phase_0, champ_rho_elem_, champ_mu_, champ_nu_, champ_rho_faces_);
      }
  }

  vpoint = 0.;

  // =====================================================================
  // Methode de projection :
  // Premiere etape : calcul de u_etoile (tous les termes de N.S. sauf la pression)

  // Contribution des operateurs diffusion et convection :
  // Operateur de diffusion : valeurs discretes homogenes a
  //                                             / d             \    //
  //                INTEGRALE                    | -- (rho * v)  |    //
  //                (sur le volume de controle)  \ dt            /    //
  // B.M. 08/2004 : on envoie la vitesse "v" a l'operateur qui calcule
  //                div (mu * (grad(v)+tr(grad(v))))
  //                (on a associe "mu" a la "diffusivite" de l'operateur,
  //                 voir Navier_Stokes_FT_Disc::lire)
  terme_diffusif.calculer(la_vitesse->valeurs(), variables_internes().terme_diffusion->valeurs());
  solveur_masse->appliquer(variables_internes().terme_diffusion->valeurs());

  // Termes sources : gravite et tension de surface,
  // valeurs discretes homogenes a
  //                                             / d             \    //
  //                INTEGRALE                    | -- (rho * v)  |    //
  //                (sur le volume de controle)  \ dt            /    //

  {
    // Si une equation de transport est associee aux proprietes du fluide,
    // on ajoute le terme de tension de surface.
    REF(Transport_Interfaces_FT_Disc) &refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;

    if (refeq_transport.non_nul())
      {
        const Champ_base& indicatrice = refeq_transport->get_update_indicatrice();
        Champ_base& gradient_i = variables_internes().gradient_indicatrice.valeur();
        // Note:
        // On appelle la version const de maillage_interface() (qui est publique) car
        // on passe par const Transport_Interfaces_FT_Disc :
        const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
        const Maillage_FT_Disc& maillage = eq_transport.maillage_interface();
        const DoubleTab& distance_interface_sommets = eq_transport.get_update_distance_interface_sommets();

        calculer_gradient_indicatrice(indicatrice, distance_interface_sommets, gradient_i);
        calculer_champ_forces_superficielles(maillage, gradient_i, variables_internes().potentiel_elements, variables_internes().potentiel_faces, variables_internes().terme_source_interfaces);
      }
    else
      {
        variables_internes().terme_source_interfaces->valeurs() = 0.;
      }
  }
  solveur_masse->appliquer(variables_internes().terme_source_interfaces->valeurs());

  // Autres termes sources (acceleration / repere mobile)
  //  Valeurs homogenes a
  //                                             / d       \          //
  //                INTEGRALE                    | -- (v)  |          //
  //                (sur le volume de controle)  \ dt      /          //
  // (voir "preparer_calcul", commentaire sources().associer_rho...)
  variables_internes().terme_source->valeurs() = 0.;
  les_sources.ajouter(variables_internes().terme_source->valeurs());
  solveur_masse->appliquer(variables_internes().terme_source->valeurs());

  // Operateur de convection : valeurs discretes homogenes a
  //                                             / d       \          //
  //                INTEGRALE                    | -- (v)  |          //
  //                (sur le volume de controle)  \ dt      /          //
  // B.M. 08/2004 : on transporte "v" et non "rho_v"...

  DoubleTab& terme_convection_valeurs = variables_internes().terme_convection->valeurs();
  bool calcul_explicite = false;
  if (parametre_equation_.non_nul() && sub_type(Parametre_implicite, parametre_equation_.valeur()))
    {
      Parametre_implicite& param2 = ref_cast(Parametre_implicite, parametre_equation_.valeur());
      calcul_explicite = param2.calcul_explicite();
    }
  if (schema_temps().diffusion_implicite() && !calcul_explicite)
    {
      terme_convection_valeurs = 0;
      derivee_en_temps_conv(terme_convection_valeurs, la_vitesse->valeurs());
    }
  else
    {
      terme_convectif.calculer(la_vitesse->valeurs(), terme_convection_valeurs);
    }
  solveur_masse->appliquer(variables_internes().terme_convection->valeurs());

  // Ajout des differentes contributions a vpoint :
  const DoubleTab& tab_rho_faces = champ_rho_faces_->valeurs();
  const DoubleVect& volumes_entrelaces = ref_cast(Domaine_VF, domaine_dis()).volumes_entrelaces();
  const DoubleTab& tab_diffusion = variables_internes().terme_diffusion->valeurs();
  const DoubleTab& termes_sources_interf = variables_internes().terme_source_interfaces->valeurs();
  const DoubleTab& termes_sources = variables_internes().terme_source->valeurs();
  const DoubleTab& tab_convection = variables_internes().terme_convection->valeurs();
  const int n = vpoint.dimension(0);
  const int nbdim1 = (vpoint.line_size() == 1);
  const int m = vpoint.line_size();

  DoubleTab gravite_face(inconnue()->valeurs());
  if (milieu().a_gravite())
    {
      ArrOfDouble g(dimension);
      // Pour l'instant : gravite uniforme g => phi(s) = - x scalaire g
      const DoubleTab& gravite = milieu().gravite().valeurs();
      if (gravite.nb_dim() != 2 || gravite.line_size() != dimension)
        {
          Cerr << "Error for calculer_champ_forces_superficielles\n";
          Cerr << " gravite.line_size() != Objet_U::dimension" << finl;
          Process::exit();
        }
      for (int j = 0; j < dimension; j++)
        g[j] = gravite(0, j);

      // On multiplie par les volumes entrelaces et on applique ensuite le solveur masse
      //  (traitement special des CL de Dirichlet)
      if (nbdim1)
        {
          const IntTab& face_voisins = le_dom_dis.valeur().face_voisins();
          const IntVect& orientation = ref_cast(Domaine_VDF, domaine_dis()).orientation();
          if (variables_internes().terme_gravite_ == Navier_Stokes_FT_Disc_interne::GRAVITE_RHO_G)
            {
              for (int face = 0; face < n; face++)
                gravite_face(face, 0) = volumes_entrelaces(face) * g[orientation[face]];
            }
          else
            {
              gravite_face = 0.; // En gradI, on ne prend pas directement la gravite
            }
          // Boussinesq Approximation in use :
          if (variables_internes().is_boussinesq_)
            {
              REF(Transport_Interfaces_FT_Disc) &refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
              if (!refeq_transport.non_nul())
                {
                  Cerr << "Trying to use Boussinesq approximation on a 2phase flow when the transport equation is not specified" << finl;
                  Process::exit();
                }
              const DoubleTab& indicatrice = refeq_transport->get_update_indicatrice().valeurs();
              // First phase with temperature :
              if (variables_internes().ref_equation_mpoint_.non_nul())
                {
                  compute_boussinesq_additional_gravity(variables_internes().ref_equation_mpoint_.valeur(), fluide_diphasique(), face_voisins, volumes_entrelaces, orientation, indicatrice, g,
                                                        gravite_face);
                }

              // Second phase with temperature :
              if (variables_internes().ref_equation_mpoint_vap_.non_nul())
                {
                  compute_boussinesq_additional_gravity(variables_internes().ref_equation_mpoint_vap_.valeur(), fluide_diphasique(), face_voisins, volumes_entrelaces, orientation, indicatrice, g,
                                                        gravite_face);
                }
              // The end of boussinesq force source term.
            }
          // End of if nbdim1 that is VDF.
        }
      else
        {
          // VEF case:
          if (variables_internes().is_boussinesq_)
            {
              Cerr << "Trying to use Boussinesq approximation on a 2phase flow in VEF? Not yet available. Ask TRUST support." << finl;
              Process::exit();
            }
          if (variables_internes().terme_gravite_ == Navier_Stokes_FT_Disc_interne::GRAVITE_RHO_G)
            {
              for (int face = 0; face < n; face++)
                for (int dim = 0; dim < m; dim++)
                  gravite_face(face, dim) = volumes_entrelaces(face) * g[dim];
            }
          else
            {
              gravite_face = 0.; // En gradI, on ne prend pas directement la gravite
            }
        }
      solveur_masse->appliquer(gravite_face);
    }
  else
    {
      // Pas de gravite :
      gravite_face = 0.;
    }

  IntTab flag_gradP(n);
  IntTab coef_TSF(n);
  coef_TSF = 1;
  int flag_diff;
  DoubleTab& gradP = variables_internes().gradient_pression->valeurs();
  if (schema_temps().diffusion_implicite())
    {
      //on calcule gradP (pour le qdm)
      gradient.calculer(la_pression->valeurs(), gradP);
      solveur_masse->appliquer(gradP);
      // si variables_internes().is_penalized=1 alors variables_internes().is_explicite=0
      if (variables_internes().is_penalized)
        flag_gradP = 0; //(grad P raide)
      else
        flag_gradP = 1;
      // on ajoute pas la diffusion cela sera fait par Gradient_conjugue_diff_impl
      // sauf si !is_explicite (terme forcage vitesse implicite; resolution conjointe forcage diffusion)
      // car dans ce cas la vitesse a imposer a besoin d'une bonne approximation de vpoint
      if (!variables_internes().is_explicite && !calcul_explicite)
        flag_diff = 1;
      else
        flag_diff = 0;
    }
  else
    {
      flag_gradP = 0;
      flag_diff = 1;
    }

  bool interf_vitesse_imposee_ok = false;
  int nb_eqs = variables_internes().ref_eq_interf_vitesse_imposee.size();
  int nb_eq_non_nul = 0;
  for (int k = 0; k < nb_eqs; k++)
    {
      REF(Transport_Interfaces_FT_Disc) &refeq_transport = variables_internes().ref_eq_interf_vitesse_imposee[k];

      if (refeq_transport.non_nul())
        nb_eq_non_nul += 1;
    }
  DoubleTab terme_mul;
  if (nb_eq_non_nul == nb_eqs && nb_eqs != 0)
    {
      interf_vitesse_imposee_ok = true;
      terme_mul.copy(champ_rho_faces_->valeurs(), RESIZE_OPTIONS::COPY_INIT);
      terme_mul = 0.;
    }

  REF(Transport_Interfaces_FT_Disc) &refeq_transport_2pha = variables_internes().ref_eq_interf_proprietes_fluide;
  if (refeq_transport_2pha.non_nul() && interf_vitesse_imposee_ok && variables_internes().is_penalized)
    {
      const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis());
      const Domaine_VF& domaine_vdf = ref_cast(Domaine_VDF, domaine_dis());
      const IntTab& face_voisins = domaine_vf.face_voisins();
      const IntTab& elem_faces = domaine_vf.elem_faces();
      const IntVect& orientation = domaine_vdf.orientation();
      const int nb_faces_elem = elem_faces.line_size();
      for (int k = 0; k < nb_eqs; k++)
        {
          REF(Transport_Interfaces_FT_Disc) &refeq_transport = variables_internes().ref_eq_interf_vitesse_imposee[k];
          const DoubleTab& indicatrice_faces = refeq_transport->get_compute_indicatrice_faces().valeurs();
          for (int i = 0; i < face_voisins.dimension(0); i++)
            {
              if (indicatrice_faces(i) > 0.)
                {
                  flag_gradP(i) = 0;
                  coef_TSF(i) = 0; //annulation local terme source interf
                  const int ori = orientation(i);
                  const int voisin0 = face_voisins(i, 0);
                  if (voisin0 >= 0)
                    {
                      int face_visavi = elem_faces(voisin0, ori) + elem_faces(voisin0, ori + Objet_U::dimension) - i;
                      for (int i_face = 0; i_face < nb_faces_elem; i_face++)
                        {
                          const int face = elem_faces(voisin0, i_face);
                          if (indicatrice_faces(face) == 0. && face == face_visavi)
                            {
                              flag_gradP(face) = 0;
                              coef_TSF(face) = 0; //annulation local terme source interf
                            }
                        }
                    }
                  const int voisin1 = face_voisins(i, 1);
                  if (voisin1 >= 0)
                    {
                      int face_visavi = elem_faces(voisin1, ori) + elem_faces(voisin1, ori + Objet_U::dimension) - i;
                      for (int i_face = 0; i_face < nb_faces_elem; i_face++)
                        {
                          const int face = elem_faces(voisin1, i_face);
                          if (indicatrice_faces(face) == 0. && face == face_visavi)
                            {
                              flag_gradP(face) = 0;
                              coef_TSF(face) = 0; //annulation local terme source interf
                            }
                        }
                    }
                }
            }
        }
    }

  // Ajout des differentes contributions a vpoint :
  for (int i = 0; i < n; i++)
    {
      const double rho_face = tab_rho_faces(i);

      for (int j = 0; j < m; j++)
        vpoint(i, j) = (-flag_gradP(i) * gradP(i, j) + flag_diff * tab_diffusion(i, j) + coef_TSF(i) * termes_sources_interf(i, j)) / rho_face + tab_convection(i, j) + termes_sources(i, j)
                       + gravite_face(i, j);
    }
  vpoint.echange_espace_virtuel();

  //  si terme forcage vitesse explicite => J'ai tout, je peux resoudre (Gradient_conjugue_diff_impl)
  if (schema_temps().diffusion_implicite() && !calcul_explicite && variables_internes().is_explicite)
    {
      DoubleTab derivee(la_vitesse->valeurs());
      // on indique au solveur masse de diviser par rho en plus du volume car l'operateur de diffusion renvoit qqqc en rho d u/dt
      solveur_masse->set_name_of_coefficient_temporel(champ_rho_faces_->le_nom());

      DoubleTrav tt(vpoint);
      tt = vpoint;
      derivee = inconnue()->valeurs();
      Equation_base::Gradient_conjugue_diff_impl(tt, derivee);

      solveur_masse->set_name_of_coefficient_temporel("no_coeff");

      vpoint = derivee;
      // on retire le gradient si on ne penalise pas:
      for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++)
          vpoint(i, j) += gradP(i, j) / tab_rho_faces(i);
    }

  const int nfaces = vpoint.dimension_tot(0);

  if (interf_vitesse_imposee_ok)
    {
      int compteur_vimp_regul = 0;
      DoubleTrav vpoint0(vpoint);
      vpoint0 = vpoint;
      terme_mul = 1.0;
      // S'il y a une equation de transport avec vitesse imposee, on impose:
      DoubleTrav forces_tot(vpoint);
      int nb_eqs_bis = variables_internes().ref_eq_interf_vitesse_imposee.size();
      for (int k = 0; k < nb_eqs_bis; k++)
        {
          REF(Transport_Interfaces_FT_Disc) &refeq_transport = variables_internes().ref_eq_interf_vitesse_imposee[k];

          Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();

          const DoubleTab& inco_val = inconnue()->valeurs();
          const DoubleTab& rho_faces = champ_rho_faces_->valeurs();
          DoubleTab& source_val = variables_internes().terme_source->valeurs();
          const double temps = schema_temps().temps_courant();
          const double dt = schema_temps().pas_de_temps();

          //On ajoute un terme source a vpoint pour imposer au fluide la vitesse de l interface
          //source_val est rempli (peut etre postraite)
          eq_transport.modifier_vpoint_pour_imposer_vit(inco_val, vpoint0, vpoint, rho_faces, source_val, temps, dt, variables_internes().is_explicite, variables_internes().eta);
          source_val.echange_espace_virtuel();
          forces_tot += source_val;

          // Afin de savoir s il existe des interfaces IBC/fluide regularisees.
          // Si oui, on ne penalise que pour indic=1.
          if (eq_transport.get_vimp_regul())
            compteur_vimp_regul++;
          //calcul du terme 1 + somme Xs / eta
          if (!variables_internes().is_explicite)
            {
              const DoubleTab& indicatrice_faces = refeq_transport->get_indicatrice_faces().valeurs();
              for (int i = 0; i < nfaces; i++)
                {
                  if ((eq_transport.get_vimp_regul() == 0 && indicatrice_faces(i) > 0.) || (eq_transport.get_vimp_regul() == 1 && indicatrice_faces(i) == 1.))
                    terme_mul(i) += 1. / variables_internes().eta;
                }
            }
        }
      terme_mul.echange_espace_virtuel();
      Debog::verifier("Navier_Stokes_FT_Disc::derivee_en_temps_inco terme_mul:", terme_mul);

      if (schema_temps().diffusion_implicite())
        {
          // si !variables_internes().is_explicite (terme forcage implicite) on retire la diffusion explicite
          // de vpoint (implicitee dans Gradient_conjugue_diff_impl_IBC)
          if (!variables_internes().is_explicite && !calcul_explicite)
            {
              const DoubleTab& rho_faces = champ_rho_faces_->valeurs();
              const DoubleTab& diffusion = variables_internes().terme_diffusion->valeurs();

              for (int i = 0; i < vpoint.dimension(0); i++)
                for (int j = 0; j < vpoint.line_size(); j++)
                  vpoint(i, j) -= diffusion(i, j) / rho_faces(i);

              vpoint.echange_espace_virtuel();
              DoubleTab derivee(la_vitesse->valeurs());
              // on indique au solveur masse de diviser par rho en plus du volume car l'operateur de diffusion renvoit qqqc en rho d u/dt
              {
                solveur_masse->set_name_of_coefficient_temporel(champ_rho_faces_->le_nom());

                DoubleTrav tt(vpoint);
                tt = vpoint;
                Equation_base::Gradient_conjugue_diff_impl(tt, derivee, terme_mul);

                solveur_masse->set_name_of_coefficient_temporel("no_coeff");
              }
              vpoint = derivee;

              // on retire le gradient
              const int nbis = vpoint.dimension(0);
              const int mbis = vpoint.line_size();
              for (int i = 0; i < nbis; i++)
                for (int j = 0; j < mbis; j++)
                  vpoint(i, j) += (flag_gradP(i) * gradP(i, j)) / rho_faces(i);
            }
          vpoint.echange_espace_virtuel();
        }
      else
        {
          // si on implicite le calcul de vpoint
          // On divise vpoint par le terme multiplicatif calcule avant
          if (!variables_internes().is_explicite)
            {
              const int mbis = vpoint.line_size();
              // calcul de vpoint : vpoint / (1 + somme Xs/eta )
              for (int i = 0; i < nfaces; i++)
                for (int j = 0; j < mbis; j++)
                  vpoint(i, j) /= terme_mul(i);

              vpoint.echange_espace_virtuel();
            }

        }
      Debog::verifier("Navier_Stokes_FT_Disc::derivee_en_temps_inco vpoint:", vpoint);

      // Dans le cas penalize + vitesse imposee regularisee,
      // seuls les ddl tels que indic_faces=1 sont penalisees, les
      // autres sont forces avec un DF :
      if (!variables_internes().is_explicite && compteur_vimp_regul)
        {
          vpoint0 = vpoint;
          // S'il y a une equation de transport avec vitesse imposee, on impose:
          DoubleTrav forces_totbis(vpoint);
          for (int k = 0; k < nb_eqs_bis; k++)
            {
              REF(Transport_Interfaces_FT_Disc) &refeq_transport = variables_internes().ref_eq_interf_vitesse_imposee[k];

              Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();

              const DoubleTab& inco_val = inconnue()->valeurs();
              const DoubleTab& rho_faces = champ_rho_faces_->valeurs();
              DoubleTab& source_val = variables_internes().terme_source->valeurs();
              const double temps = schema_temps().temps_courant();
              const double dt = schema_temps().pas_de_temps();

              //On ajoute un terme source a vpoint pour imposer au fluide la vitesse de l interface
              //source_val est rempli (peut etre postraite)
              eq_transport.modifier_vpoint_pour_imposer_vit(inco_val, vpoint0, vpoint, rho_faces, source_val, temps, dt,/* is_explicite */1, /* eta */1.);
              source_val.echange_espace_virtuel();
              forces_totbis += source_val;
            }
        }
      vpoint.echange_espace_virtuel();

      //Calcul des efforts exerces par le fluide sur chaque interface
      //Attention valable si les ibc ne se chevauchent pas
      if (limpr())
        {
          const DoubleTab& tab_rho_facesbis = champ_rho_faces_->valeurs();
          for (int k = 0; k < nb_eqs_bis; k++)
            {
              REF(Transport_Interfaces_FT_Disc) &refeq_transport = variables_internes().ref_eq_interf_vitesse_imposee[k];
              Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
              eq_transport.calcul_effort_fluide_interface(vpoint, tab_rho_facesbis, forces_tot, variables_internes().is_explicite, variables_internes().eta);
              forces_tot.echange_espace_virtuel();
            }
        }
    }

  // Assemblage de la matrice INTEGRALE            ( div(1/rho * grad(P)) ) :
  //                          (volume de controle pression)
  // Si l'option "matrice_pression_invariante" est activee, on ne recalcule
  // pas la matrice :
  if (!interf_vitesse_imposee_ok)
    {
      if (!variables_internes().matrice_pression_invariante)
        {
          assembleur_pression()->assembler_rho_variable(matrice_pression_, champ_rho_faces_.valeur());
          // On a modifie la matrice, il faut reinitialiser le solveur
          //  (recalcul des preconditionnement, factorisation pour Cholesky, etc...)
          solveur_pression()->reinit();
        }
    }

  // ====================================================================
  // Methode de projection :
  // Deuxieme etape : projection du champ de vitesse sur le sous-espace
  //  a divergence nulle.
  // Trouver la pression P telle que
  //   d(u)/dt = vpoint - 1/rho * grad(P)
  //   div( d(u)/dt ) = 0
  // Soit :
  //   div(1/rho * grad(P)) = div(vpoint)

  // Calcul du second membre :
  //  div(vpoint) a l'interieur du domaine,
  //  prise en compte des conditions aux limites de pression/vitesse

  DoubleTab& secmem = variables_internes().second_membre_projection->valeurs();
  DoubleTab& secmem2 = variables_internes().second_membre_projection_jump_->valeurs();
  const int nb_elem = secmem2.dimension(0);
  const double dt = schema_temps().pas_de_temps();
  const DoubleTab& inco = inconnue()->valeurs();
  // secmem = div(U/dt+vpoint) = div(U(n+1)/dt)
  DoubleTab du(inco);
  du /= dt;
  du += vpoint;
  divergence.calculer(du, secmem);
  secmem *= -1;
#if NS_VERBOSE
  double int_sec_mem = 0;
  for (int elem = 0; elem < secmem.dimension(0); elem++)
    int_sec_mem +=secmem(elem);
  Cerr << "Secmem before tcl= " << int_sec_mem << finl;
#endif
  // Prise en compte des sources de constituant (terme source dans une equation concentration)
  // Pour ne pas faire figurer la source deux fois, je la mets uniquement dans l'equation
  // de concentration... elle produit alors un terme source de div_u dans Navier_Stokes:
  const Noms& noms_eq = variables_internes().equations_concentration_source_fluide_;
  for (int i_eq = 0; i_eq < noms_eq.size(); i_eq++)
    {
      const Equation_base& eq = probleme().get_equation_by_name(noms_eq[i_eq]);
      for (int i_source = 0; i_source < eq.sources().size(); i_source++)
        {
          if (!sub_type(Terme_Source_Constituant_Vortex_VEF_Face, eq.sources()[i_source].valeur()))
            continue;

          const Terme_Source_Constituant_Vortex_VEF_Face& src = ref_cast(Terme_Source_Constituant_Vortex_VEF_Face, eq.sources()[i_source].valeur());
          src.ajouter_terme_div_u(secmem, schema_temps().pas_de_temps());
        }

    }
  secmem.echange_espace_virtuel();

  // Prise en compte du terme source div(u) du changement de phase
  if (variables_internes().ref_equation_mpoint_.non_nul() || variables_internes().ref_equation_mpoint_vap_.non_nul())
    {
      // GB2016 : Le calcul de mpoint ci-dessous me semble inutile car il est fait au debut de calculer_delta_u_interface:
      // GB2016 : Mais si je ne le fais pas, j'ai parfois : 'vx.get_md_vector() == md' dans calculer_delta_u_interface

      // GB2022 : Je ne comprend toujours pas bien pourquoi, mais il faut calculer les mpoints
      //          pour avoir un cas FTD_Boiling_bubble avec une extension correcte
      //          (sinon le champ postraite de T est moche dans l'extension)
      if (variables_internes().ref_equation_mpoint_.non_nul())
        variables_internes().ref_equation_mpoint_->calculer_mpoint(variables_internes().mpoint.valeur());
      if (variables_internes().ref_equation_mpoint_vap_.non_nul())
        variables_internes().ref_equation_mpoint_vap_->calculer_mpoint(variables_internes().mpoint_vap.valeur());

      // Pas une ref, mais un tableau de travail local dans lequel on peut ajouter mointv
      DoubleTab mpoint = variables_internes().ref_equation_mpoint_->get_mpoint();
      if (variables_internes().ref_equation_mpoint_vap_.non_nul())
        {
          const DoubleTab& mpointv = variables_internes().ref_equation_mpoint_vap_->get_mpoint();
          for (int elem = 0; elem < nb_elem; elem++)
            mpoint[elem] += mpointv[elem];
        }
      // We can compute delta_u_interface:
      // depending on the option, either historical or new, the calculation may be based on the values filled in secmem2
      calculer_delta_u_interface(variables_internes().delta_u_interface, -1 /* vitesse de l'interface */, variables_internes().correction_courbure_ordre_ /* ordre de la correction en courbure */);

      const Fluide_Diphasique& fluide_diph = fluide_diphasique();
      const Fluide_Incompressible& phase_0 = fluide_diph.fluide_phase(0);
      const Fluide_Incompressible& phase_1 = fluide_diph.fluide_phase(1);
      const DoubleTab& tab_rho_phase_0 = phase_0.masse_volumique()->valeurs();
      const DoubleTab& tab_rho_phase_1 = phase_1.masse_volumique()->valeurs();
      const double rho_phase_0 = tab_rho_phase_0(0, 0);
      const double rho_phase_1 = tab_rho_phase_1(0, 0);
      const double jump_inv_rho = 1. / rho_phase_1 - 1. / rho_phase_0;
      if (variables_internes().new_mass_source_)
        {

          const DoubleTab& interfacial_area = variables_internes().ai->valeurs();
          for (int elem = 0; elem < nb_elem; elem++)
            secmem2[elem] = jump_inv_rho * interfacial_area[elem] * mpoint[elem];
        }
      else
        {
          const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis());
          const IntTab& face_voisins = domaine_vf.face_voisins();
          const IntTab& elem_faces = domaine_vf.elem_faces();
          const int nb_faces_elem = elem_faces.line_size();
          REF(Transport_Interfaces_FT_Disc) &refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
          const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
#if NS_VERBOSE
          const DoubleTab& indicatrice = eq_transport.inconnue()->valeurs();
#endif
          // Distance a l'interface discretisee aux elements:
          const DoubleTab& distance = eq_transport.get_update_distance_interface().valeurs();
          divergence.calculer(variables_internes().delta_u_interface->valeurs(), secmem2);
          // On ne conserve que la divergence des elements traverses par l'interface
          for (int elem = 0; elem < nb_elem; elem++)
            {
              const double dist = distance(elem);
              int i_face = -1;
              if (dist < -1e20)
                {
                  // Distance invalide: on est loin de l'interface
                  i_face = nb_faces_elem;
                }
              else
                {
                  // Y a-t-il un voisin pour lequel la distance est de signe oppose
                  for (i_face = 0; i_face < nb_faces_elem; i_face++)
                    {
                      const int face = elem_faces(elem, i_face);
                      const int voisin = face_voisins(face, 0) + face_voisins(face, 1) - elem;
                      if (voisin >= 0)
                        {
                          const double d = distance(voisin);
                          if (d > -1e20 && d * dist < 0.)
                            {
#if NS_VERBOSE
                              Cerr << "Compa "<< secmem2[elem] << " " << jump_inv_rho*interfacial_area[elem]*mpoint[elem] << finl;
#endif
                              break; // Changement de signe
                            }
                        }
                    }
                }
              if (i_face == nb_faces_elem)
                {
                  // Tous les voisins sont du meme cote de l'interface
                  secmem2(elem) = 0.;
#if NS_VERBOSE
                  if(interfacial_area[elem]>DMINFLOAT)
                    {
                      Cerr << "[WARNING] secmem2 is set to zero in element whereas phase is not pure (indic= "
                           << indicatrice[elem] << "). This is because the choice is based on the signs of distance." << finl;
                      // Indeed, in a diagonal case, for the cell with "x", indic can be close to (but not) pure and all neighbouring
                      // cells can still have the same distance.
                      // __________
                      // |  | /|  |
                      // |__|/_|__|
                      // | /| x|  |
                      // |/_|__|__|
                    }
#endif
                }
            }
        }
#if TCL_MODEL

      /* double int_sec_mem2 = 0.;
       for (int elem = 0; elem < nb_elem; elem++)
       {
       int_sec_mem2 +=secmem2(elem);
       int_sec_mem +=secmem(elem);
       }
       Cerr << "Integral of secmem2 before TCL and /DT : " << int_sec_mem2 << finl; */

      // Now that the correction "corriger_mpoint" is performed directly into Convection_Diffusion_Temperature_FT_Disc,
      // it is no longer required to compute it here. The correction should propagate to calculer_delta_vitesse (in the discreete sense),
      // and subsequently to secmem2. So that in the end, the whole TCL contribution should be accounted for naturally.
      // However, in the discretized version, It is still required to correct secmem2 even though mpoint was corrected itself.
      // It is because of the way secmem is computed (as the div( delta u)) that is using part of cells not crossed by the interface
      // (in the wall-normal direction). It results in an underestimation of the real TCL contribution...
      if (probleme_ft().tcl().is_activated())
        {
          Cerr << "[TCL] Contact line model is activated" << finl;

          const Triple_Line_Model_FT_Disc& tcl = ref_cast(Triple_Line_Model_FT_Disc, probleme_ft().tcl());

          const ArrOfInt& elems_with_CL_contrib = tcl.elems();
          // const ArrOfInt& faces_with_CL_contrib = probleme_ft().tctl().boundary_faces();
          const ArrOfDouble& Q_from_CL = tcl.Q();
          // const ArrOfDouble& mpoint_from_CL = probleme_ft().tcl().mp();

          // ArrOfInt elems_with_CL_contrib;
          // ArrOfInt faces_with_CL_contrib;
          // ArrOfDouble mpoint_from_CL;
          // ArrOfDouble Q_from_CL;
          // GB. 18/12/19. This call is actually the one filling the TCL tables (elems_, mp_ and Q_);
          // probleme_ft().tcl().compute_TCL_fluxes_in_all_boundary_cells(elems_with_CL_contrib,
          //                                                             faces_with_CL_contrib,
          //                                                             mpoint_from_CL,
          //                                                             Q_from_CL);

          // Correct the field mpoint in wall-adjacent cells to account for TCL model:
          // ---> It's not added to mpoint now. Its contribution is added in
          //      Convection_Diffusion_Temperature_FT_Disc::derivee_en_temps_inco
          //      to be after the evaluation of the extended velocities (interfacial and liquid)
          //      and their interpolation. That way, interpolation can still operate on a smooth field.
          // DoubleTab& mpoint = variables_internes().mpoint->valeurs();
          // probleme_ft().tcl().corriger_mpoint(elems_with_CL_contrib,mpoint_from_CL,mpoint);

          const double Lvap = fluide_diph.chaleur_latente();
          const double coef = jump_inv_rho / Lvap;
          // Correct the secmem2 contribution due to TCL :
          if (!variables_internes().mpoint_inactif)
            probleme_ft().tcl().corriger_secmem(coef, secmem2);

          const int check_consistency = 1; // local option to check that secmem2 in near-wall cell is actually well calculated
          if (check_consistency)
            {
              Cerr << "Verifying Contact line model consistency" << finl;
              double error = 0.;
              const int nb_contact_line_contribution = elems_with_CL_contrib.size_array();
              for (int idx = 0; idx < nb_contact_line_contribution; idx++)
                {
                  const int elem = elems_with_CL_contrib[idx];
                  const double sec = secmem2(elem);
                  double Q = 0.;
                  // Go through the list to find all occurences of elem;
                  for (int idx2 = 0; idx2 < nb_contact_line_contribution; idx2++)
                    {
                      if (elem == elems_with_CL_contrib[idx2])
                        {
                          Q += Q_from_CL[idx2];
                        }
                    }
                  const double value = coef * Q;

                  // sec and value should be the same:
                  error += fabs(sec - value);
                  if (fabs(sec - value) > 1.e-12) // changed from 1.e-12 to 1.e-7 ---- for test
                    {
                      Cerr << "local difference sec-value=" << sec << " - " << value << " = " << (sec - value) << finl;
                    }
                }

              if (error > 1.e-8)
                {
                  Cerr << "Final error : " << error << " is fatal!" << finl;
                  Process::exit();
                }

            }
        }
#endif
      secmem2 /= schema_temps().pas_de_temps();
      secmem += secmem2;
      secmem.echange_espace_virtuel();
#if NS_VERBOSE
      double int_sec_mem2 = 0;
      double int_sec_mem = 0;
      for (int elem = 0; elem < nb_elem; elem++)
        {
          int_sec_mem2 +=secmem2(elem);
          int_sec_mem +=secmem(elem);
        }
      Cerr << "Integral of secmem2 after TCL and /DT : " << int_sec_mem2 << finl;
      Cerr << "Integral of secmem after TCL and /DT : " << int_sec_mem << finl;
#endif
    }

  Champ_Fonc champ_rho_faces_modifie(champ_rho_faces_);
  DoubleTab& rho_faces_modifie = champ_rho_faces_modifie->valeurs();

  if (interf_vitesse_imposee_ok)
    {

      // On verifie s il existe des interfaces IBC/fluide regularisees.
      int compteur_vimp_regul = 0;
      int nb_eqs_bis = variables_internes().ref_eq_interf_vitesse_imposee.size();
      for (int k = 0; k < nb_eqs_bis; k++)
        {
          REF(Transport_Interfaces_FT_Disc) &refeq_transport = variables_internes().ref_eq_interf_vitesse_imposee[k];
          Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();

          if (eq_transport.get_vimp_regul() == 1)
            compteur_vimp_regul++;
        }

      // Creation d'un nouveau rho qui prend en compte les vitesses imposees dans
      // le terme de forcage
      // Si des interfaces IBC/fluide sont regularisees, la projection devient classique
      // (a ameliorer potentiellement).
      // Mais dans le cas ou on a une interface diphasique, on bloque toutes les vitesses impossees
      // en PDF (regularise ou non)
      int modif_rho_true = 0;
      if (variables_internes().is_penalized && (compteur_vimp_regul == 0 || refeq_transport_2pha.non_nul()))
        {
          for (int i = 0; i < nfaces; i++)
            {
              if (compteur_vimp_regul != 0)
                {
                  for (int k = 0; k < nb_eqs; k++)
                    {
                      REF(Transport_Interfaces_FT_Disc) &refeq_transport = variables_internes().ref_eq_interf_vitesse_imposee[k];
                      Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
                      const DoubleTab& indicatrice_faces = eq_transport.get_indicatrice_faces().valeurs();
                      if (eq_transport.get_vimp_regul() == 1 && indicatrice_faces(i) > 0.0 && indicatrice_faces(i) != 1.)
                        terme_mul(i) += 1.0 / variables_internes().eta;
                    }
                }
              rho_faces_modifie(i) *= terme_mul(i);
            }
          rho_faces_modifie.echange_espace_virtuel();
          modif_rho_true = 1;
        }
      Debog::verifier("Navier_Stokes_FT_Disc::derivee_en_temps_inco rho_faces_modifie:", rho_faces_modifie);

      // Assemblage de la matrice INTEGRALE            ( div(1/rho * grad(P)) ) :
      //                          (volume de controle pression)

      assembleur_pression()->assembler_rho_variable(matrice_pression_, champ_rho_faces_modifie.valeur());

      //Penalization L2 de la pression si necessaire
      if (modif_rho_true == 1 && variables_internes().p_ref_pena != -1.e40)
        {
          // On se base sur le nombre de composantes par faces pour la discretisation
          Matrice_Morse_Sym& matrice_valeurs = (vpoint.line_size() == 1 ? ref_cast(Matrice_Morse_Sym, (ref_cast(Matrice_Bloc, matrice_pression_.valeur())).get_bloc(0,0).valeur()) // VDF (1)
                                                :
                                                ref_cast(Matrice_Morse_Sym, matrice_pression_.valeur()));                                              // VEF (>1)
          DoubleTab& pressu = la_pression->valeurs();
          assert(nb_elem == champ_rho_elem_->valeurs().dimension(0));
          const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis());
          const Domaine& mon_dom = domaine_dis().domaine();
          const IntTab& elem_faces = domaine_vf.elem_faces();
          const int nb_faces_elem = elem_faces.line_size();
          const int nb_sommet = mon_dom.nb_som_elem();
          int numero_global_som, ligne_mat;
          int point_fluide_dirichlet = -1;
          if (variables_internes().is_pfl_flottant)
            {
              if (Objet_U::dimension == 3)
                {
                  point_fluide_dirichlet = mon_dom.chercher_elements(variables_internes().x_pfl_imp, variables_internes().y_pfl_imp, variables_internes().z_pfl_imp);
                }
              else
                {
                  point_fluide_dirichlet = mon_dom.chercher_elements(variables_internes().x_pfl_imp, variables_internes().y_pfl_imp);
                }
              if (mp_max(point_fluide_dirichlet) == -1)
                {
                  Cerr << "Point de reference pression fluide situe en dehors du domaine !" << finl;
                  Process::exit();
                }
            }
          for (int e = 0; e < nb_elem; e++)
            {
              int nbfglob = 0;
              int nbfpena = 0;
              for (int f = 0; f < nb_faces_elem; f++)
                {
                  int fglob = elem_faces(e, f);
                  if (fglob >= 0)
                    {
                      nbfglob += 1;
                      int dejafait = 0;
                      for (int k = 0; k < nb_eqs_bis && dejafait == 0; k++)
                        {
                          REF(Transport_Interfaces_FT_Disc) &refeq_transport = variables_internes().ref_eq_interf_vitesse_imposee[k];
                          const DoubleTab& indicatrice_faces = refeq_transport->get_indicatrice_faces().valeurs();
                          if (indicatrice_faces(fglob) > 0.0)
                            {
                              dejafait = 1;
                              nbfpena += 1;
                            }
                        }
                    }
                }
              if (nbfpena == nbfglob && nbfpena != 0)
                {
                  matrice_valeurs(e, e) += 1.0 / variables_internes().eta;
                  secmem(e) += variables_internes().p_ref_pena / variables_internes().eta;
                  pressu(e) = variables_internes().p_ref_pena;
                  if (vpoint.line_size() > 1) // VEF
                    {
                      for (int somloc = 0; somloc < nb_sommet; somloc++)
                        {
                          numero_global_som = mon_dom.sommet_elem(e, somloc);
                          ligne_mat = nb_elem + numero_global_som;
                          //                         matrice_valeurs(ligne_mat,ligne_mat) += 1.0 / variables_internes().eta;
                          pressu(ligne_mat) = 0.0;
                        }
                    }
                }
              if (point_fluide_dirichlet == e)
                {
                  if (nbfpena != nbfglob && nbfpena != 0)
                    {
                      // On impose une reference de pression (p_ref)sur un bord
                      secmem(e) += matrice_valeurs(e, e) * variables_internes().p_ref_pena / (float(nbfglob - nbfpena));
                      matrice_valeurs(e, e) *= float(nbfglob - nbfpena + 1) / (float(nbfglob - nbfpena));
                    }
                  else
                    {
                      Cerr << "Point de reference pression fluide non situe dans une cellule fluide voisin d'une IBC !" << finl;
                      Cerr << "Nb faces IBC = " << nbfpena << finl;
                      Process::exit();
                    }
                }
            }
          secmem.echange_espace_virtuel();
          pressu.echange_espace_virtuel();
        }

      // On a modifie la matrice, il faut reinitialiser le solveur
      //  (recalcul des preconditionnement, factorisation pour Cholesky, etc...)
      solveur_pression()->reinit();
    }

  assembleur_pression_->modifier_secmem(secmem);

  // Resolution du systeme en pression : calcul de la_pression
  solveur_pression_->resoudre_systeme(matrice_pression_.valeur(), secmem, la_pression->valeurs());
  assembleur_pression_->modifier_solution(la_pression->valeurs());
  // Calcul d(u)/dt = vpoint + 1/rho*grad(P)
  gradient.calculer(la_pression->valeurs(), gradP);
  solveur_masse->appliquer(gradP);

  // Correction de vpoint :
  if (projection_a_faire()) // Temporaire pour permettre de ne pas resoudre NS avec mettant operateurs nuls et projection_initiale 0
    {
      const int nbis = vpoint.dimension(0);
      const int mbis = vpoint.line_size();
      for (int i = 0; i < nbis; i++)
        for (int j = 0; j < mbis; j++)
          vpoint(i, j) -= gradP(i, j) / rho_faces_modifie(i);

      vpoint.echange_espace_virtuel();
    }

  // Calcul des efforts sur les IBCs et impression dans un fichier
  if (interf_vitesse_imposee_ok && limpr())
    {
      const DoubleTab& rho = champ_rho_faces_->valeurs();

      DoubleTrav forces_tot_2(vpoint);
      DoubleTrav pressure_part(vpoint);
      DoubleTrav diffusion_part(vpoint);

      //Calcul de la vitesse au temps n+1
      DoubleTab vv(vpoint);
      vv *= schema_temps().pas_de_temps();
      vv += inconnue()->valeurs();

      // Terme de diffusion
      terme_diffusif.calculer(vv, variables_internes().terme_diffusion->valeurs());
      variables_internes().terme_diffusion->valeurs().echange_espace_virtuel();
      solveur_masse->appliquer(variables_internes().terme_diffusion->valeurs());
      const DoubleTab& diffusion = variables_internes().terme_diffusion->valeurs();
      // Terme de convection
      DoubleTrav trav(variables_internes().terme_convection->valeurs());
      derivee_en_temps_conv(trav, la_vitesse->valeurs());
      variables_internes().terme_convection->valeurs() = trav;
      solveur_masse->appliquer(variables_internes().terme_convection->valeurs());
      const DoubleTab& convection = variables_internes().terme_convection->valeurs();
      const int nbis = vpoint.dimension(0);
      const int mbis = vpoint.line_size();

      for (int i = 0; i < nbis; i++)
        for (int j = 0; j < mbis; j++)
          {
            pressure_part(i, j) = gradP(i, j);
            diffusion_part(i, j) = -rho(i) * convection(i, j) - diffusion(i, j);
            forces_tot_2(i, j) = pressure_part(i, j) + diffusion_part(i, j);
          }

      int nb_eqs_bis = variables_internes().ref_eq_interf_vitesse_imposee.size();
      for (int k = 0; k < nb_eqs_bis; k++)
        {
          REF(Transport_Interfaces_FT_Disc) &refeq_transport = variables_internes().ref_eq_interf_vitesse_imposee[k];
          Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
          // Impression des efforts
          eq_transport.impr_effort_fluide_interface(forces_tot_2, pressure_part, diffusion_part);
        }
    }

  return vpoint;
}

const Probleme_FT_Disc_gen& Navier_Stokes_FT_Disc::probleme_ft() const
{
  return probleme_ft_.valeur();
}

Probleme_FT_Disc_gen& Navier_Stokes_FT_Disc::probleme_ft()
{
  return probleme_ft_.valeur();
}

/*! @brief In Front Tracking, pression is in Pa and so pression_pa field <=> pression field
 *
 */
void Navier_Stokes_FT_Disc::calculer_la_pression_en_pa()
{
  la_pression_en_pa->valeurs() = la_pression->valeurs();
}

const Navier_Stokes_FT_Disc_interne& Navier_Stokes_FT_Disc::variables_internes() const
{
  return variables_internes_;
}
Navier_Stokes_FT_Disc_interne& Navier_Stokes_FT_Disc::variables_internes()
{
  return variables_internes_;
}

/*! @brief Si le champ de vitesse est discontinu (calcul avec changement de phase), renvoie un pointeur vers le champ delta_v de "discontinuite", tel que
 *
 *   inconnue - delta_v = vitesse de deplacement des interfaces
 *   (voir Transport_Interfaces_FT_Disc::deplacer_maillage_v_fluide())
 *   Si pas de changement de phase, renvoie un pointeur nul.
 *
 */
const Champ_base* Navier_Stokes_FT_Disc::get_delta_vitesse_interface() const
{
  if (variables_internes().ref_equation_mpoint_.non_nul() || variables_internes().ref_equation_mpoint_vap_.non_nul())
    return &(variables_internes().delta_u_interface.valeur());
  else
    return 0;
}
// Description : calcul de div(n) (la courbure discretisee sur le maillage volumique)
// Faudrait deplacer cette methode dans transport interfaces...
const Champ_base& Navier_Stokes_FT_Disc::calculer_div_normale_interface()
{
  REF(Transport_Interfaces_FT_Disc) &refeq_transport = variables_internes().ref_eq_interf_proprietes_fluide;
  const Transport_Interfaces_FT_Disc& eq_transport = refeq_transport.valeur();
  // Distance a l'interface discretisee aux elements:
  const DoubleTab& dist = eq_transport.get_update_distance_interface().valeurs();

  DoubleTab& phi = variables_internes().laplacien_d->valeurs();
  DoubleTab u0(inconnue()->valeurs());

  //  static const Stat_Counter_Id count = statistiques().new_counter(1, "calculer_div_normale", 0);
  //  statistiques().begin_count(count);

  phi = dist;

  const int n = phi.dimension(0);
  for (int i = 0; i < n; i++)
    {
      if (phi(i) < -1e20)
        phi(i) = 0;
    }
  phi.echange_espace_virtuel();
  //  static const Stat_Counter_Id count2 = statistiques().new_counter(1, "calculer_gradient", 0);
  //  statistiques().begin_count(count2);
  gradient.calculer(phi, u0);
  correct_at_exit_bad_gradient(u0);
  //  statistiques().end_count(count2);
  //  static const Stat_Counter_Id count4 = statistiques().new_counter(1, "calculer_solveur_masse", 0);
  //  statistiques().begin_count(count4);
  solveur_masse->appliquer(u0);
  //  statistiques().end_count(count4);

  // Calcul de integrale(div(u0)) sur les mailles:
  //  static const Stat_Counter_Id count3 = statistiques().new_counter(1, "calculer_div", 0);
  //  statistiques().begin_count(count3);
  divergence.calculer(u0, phi);
  //  statistiques().end_count(count3);
  // Division par le volume des mailles:
  const DoubleVect& volumes = ref_cast(Domaine_VF, domaine_dis()).volumes();
  for (int i = 0; i < n; i++)
    {
      const double p = phi(i);
      if (p != 0.)
        {
          const double v = volumes[i];
          phi(i) = p / v;
        }
    }

  //  statistiques().end_count(count);

  return variables_internes().laplacien_d.valeur();
}

const Champ_Fonc& Navier_Stokes_FT_Disc::champ_rho_faces() const
{
  return champ_rho_faces_;
}

//Renvoie 1 si l option GRAVITE_RHO_G est activee 0 sinon
int Navier_Stokes_FT_Disc::is_terme_gravite_rhog() const
{
  if (variables_internes().terme_gravite_ == Navier_Stokes_FT_Disc_interne::GRAVITE_RHO_G)
    return 1;
  else
    return 0;
}
/*
 void Navier_Stokes_FT_Disc::corriger_mpoint()
 {
 if (probleme_ft().tcl().is_activated())
 {
 DoubleTab& mpoint = variables_internes().mpoint->valeurs();
 probleme_ft().tcl().corriger_mpoint(mpoint);
 }
 }
 */
