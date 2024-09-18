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
// File:        Convection_Diffusion_Concentration_FT_Disc.cpp
// Directory:   $TrioCFD_ROOT/Front_tracking_discontinu/src
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////
#include <Convection_Diffusion_Concentration_FT_Disc.h>
#include <Champ_P1NC.h>
#include <Probleme_base.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <Domaine_VEF.h>
#include <Champ_Inc_P0_base.h>
#include <Check_espace_virtuel.h>
#include <Debog.h>
#include <Domaine.h>
#include <Sous_Domaine.h>
#include <SFichier.h>
#include <Param.h>
#include <Constituant.h>
#include <Fluide_Diphasique.h>
#include <Fluide_Incompressible.h>

Implemente_instanciable_sans_constructeur(Convection_Diffusion_Concentration_FT_Disc,"Convection_Diffusion_Concentration_FT_Disc",Convection_Diffusion_Concentration);

Convection_Diffusion_Concentration_FT_Disc::Convection_Diffusion_Concentration_FT_Disc()
{
  option_ = RIEN;
  phase_a_conserver_ = -1;
  constante_cinetique_ = 0.;
  constante_cinetique_nu_t_ = 0.;
  // default kinematic model is 1
  modele_cinetique_ = 1;
}

/*! @brief cf Convection_Diffusion_Concentration::readOn(is)
 *
 */
Entree& Convection_Diffusion_Concentration_FT_Disc::readOn(Entree& is)
{
  Convection_Diffusion_Concentration::readOn(is);
  if (equations_source_chimie_.size() > 0)
    {
      Cerr << "Chemical reaction using these equations:" << equations_source_chimie_ << finl;
      if (constante_cinetique_nu_t_ > 0.)
        {
          if (nom_equation_nu_t_ == "??")
            {
              Cerr << "Missing EQUATION_NU_T" << finl;
              Process::exit();
            }
        }
      Cerr << "The kinematic model used is : modele_cinetique " << modele_cinetique_ << finl;
    }

  if (je_suis_maitre())
    {
      Nom nom_fic(le_nom());
      nom_fic += "_bilan.out";
      SFichier f(nom_fic);
      char s[1000];
      snprintf(s, 1000, "#%15s %16s %16s %16s %16s",
               "temps",
               "integrale0",
               "integrale1",
               "bilan_sortie",
               "debit_sortie");
      f << s << finl;
    }

  return is;
}

Sortie& Convection_Diffusion_Concentration_FT_Disc::printOn(Sortie& os) const
{
  Convection_Diffusion_Concentration::printOn(os);
  return os;
}

void Convection_Diffusion_Concentration_FT_Disc::set_param(Param& param)
{
  Convection_Diffusion_Concentration::set_param(param);
  param.ajouter("constante_cinetique",&constante_cinetique_);
  param.ajouter("equations_source_chimie",&equations_source_chimie_);
  param.ajouter("modele_cinetique",&modele_cinetique_);
  param.ajouter("constante_cinetique_nu_t",&constante_cinetique_nu_t_);
  param.ajouter("equation_nu_t",&nom_equation_nu_t_);
  param.ajouter("domaine_sortie",&nom_domaine_sortie_);
  param.ajouter("phase",&phase_a_conserver_);
  param.ajouter_condition("(value_of_phase_eq_0)_or_(value_of_phase_eq_1)","phase must be set to 0 or 1");
  param.ajouter_non_std("equation_interface",(this),Param::REQUIRED);
  param.ajouter_non_std("option",(this));
}

int Convection_Diffusion_Concentration_FT_Disc::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="equation_interface")
    {
      Nom nom;
      is >> nom;
      ref_eq_transport_ = mon_probleme->get_equation_by_name(nom);
      if (!sub_type(Transport_Interfaces_FT_Disc,ref_eq_transport_.valeur()))
        {
          Cerr << " Error: the equation indicated for keyword equation_interface "<<finl;
          Cerr << " is not of type (or sub type) of Transport_Interfaces_FT_Disc." << finl;
          barrier();
          exit();
        }
      return 1;
    }
  else if (mot=="option")
    {
      Motcle motbis;
      is >> motbis;
      if (motbis == "RAMASSE_MIETTES_SIMPLE")
        option_ = RAMASSE_MIETTES_SIMPLE;
      else if (motbis == "RIEN")
        option_ = RIEN;
      else
        {
          Cerr << "Error: expected RAMASSE_MIETTES_SIMPLE or RIEN after keyword OPTION" << finl;
          barrier();
          exit();
        }
      return 1;
    }
  else
    return Convection_Diffusion_Concentration::lire_motcle_non_standard(mot,is);
}

/*! @brief calcul de l'integrale du champ sur tout le domaine.
 *
 * A factoriser un jour quelque part...
 *
 */
static double integrale(const Champ_base& ch)
{
  double somme = 0.;
  if (sub_type(Champ_P1NC, ch))
    {
      const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, ch.domaine_dis_base());
      const DoubleTab& val = ch.valeurs();
      assert(val.line_size() == 1);
      const int sz = val.dimension(0);
      const DoubleVect& volumes = domaine_vf.volumes_entrelaces();
      const ArrOfInt& faces_doubles = domaine_vf.faces_doubles();
      for (int i = 0; i < sz; i++)
        {
          const double v = volumes[i];
          const double c = val(i);
          const double fact = faces_doubles[i] ? 0.5 : 1.;
          somme += v * c * fact;
        }
    }
  else
    {
      Cerr << "Method integrale is not developped for " << ch.que_suis_je() << finl;
      Process::barrier();
      Process::exit();
    }
  somme = Process::mp_sum(somme);
  return somme;
}

/*! @brief calcul de l'indicatrice aux faces.
 *
 * .. A factoriser un jour...
 *
 */
static void calculer_indicatrice_comme(const Champ_base& src, DoubleTab& dest, const Champ_base& champ_modele)
{
  assert(sub_type(Champ_Inc_P0_base, src));
  if (sub_type(Champ_P1NC, champ_modele))
    {
      const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, champ_modele.domaine_dis_base());
      const IntTab& face_voisins = domaine_vf.face_voisins();
      const int nb_faces = domaine_vf.nb_faces();
      dest.resize(nb_faces);
      const DoubleTab& indic = src.valeurs();
      assert_espace_virtuel_vect(indic);
      for (int i = 0; i < nb_faces; i++)
        {
          const int elem0 = face_voisins(i, 0);
          const int elem1 = face_voisins(i, 1);
          double s = 0.;
          if (elem0 >= 0)
            s += indic(elem0);
          if (elem1 >= 0)
            s += indic(elem1);
          s *= ((elem0>=0) && (elem1>=0)) ? 0.5 : 1.;
          dest[i] = s;
        }
    }
  else
    {
      Cerr << "Error calculer_indicatrice_comme" << finl;
      Process::exit();
    }
}

//  Application d'un traitement particulier pour conserver la masse du constituant
//  tout en l'empechant de "baver" d'une phase a l'autre.
//  Annule la concentration dans la phase qui n'est pas a conserver et
//  rescale celle de la phase a conserver pour que l'integrale soit conservee.
//  La concentration est discontinue a l'interface liquide-vapeur.
void Convection_Diffusion_Concentration_FT_Disc::ramasse_miette_simple(double temps)
{
  const Transport_Interfaces_FT_Disc& eq = ref_cast(Transport_Interfaces_FT_Disc, ref_eq_transport_.valeur());
  DoubleTab tmp;
  Champ_base& champ_inco = inconnue().valeur();
  DoubleTab&   inco       = champ_inco.valeurs();

  // Calcule une indicatrice de phase localisee comme l'inconnue de l'equation
  calculer_indicatrice_comme(eq.inconnue(), tmp, champ_inco);

  // Calcul de l'integrale de la concentration dans tout le volume:
  const double masse_totale = integrale(champ_inco);

  // On met a zero la concentration dans la phase qui n'est pas a conserver
  {
    const int sz = inco.size();
    assert(sz == tmp.dimension(0));
    for (int i = 0; i < sz; i++)
      {
        double facteur = tmp[i];
        if (phase_a_conserver_ == 0)
          facteur = 1. - facteur;
        assert(facteur > -1e-6 && facteur < 1.+1e-6);
        inco[i] *= facteur;
      }
  }
  // Que reste-t-il comme masse ?
  const double masse_restante = integrale(champ_inco);
  if (je_suis_maitre())
    {
      Journal() << "C_D_concentration_ft_disc::ramasse_miette_simple t= " << temps
                << " masse_totale= " << masse_totale
                << " masse_restante= " << masse_restante << finl;
    }
  // On multiplie la concentration par une constante pour obtenir la masse totale initiale
  if (masse_restante > 0.)
    inco *= (masse_totale / masse_restante);

  inco.echange_espace_virtuel();
}

void Convection_Diffusion_Concentration_FT_Disc::mettre_a_jour_chimie()
{
  // On va ecraser les inconnues des equations de chimie:
  if (equations_source_chimie_.size() != 2)
    {
      Cerr << "Error for method Convection_Diffusion_Concentration_FT_Disc::mettre_a_jour_chimie()\n"
           << "It implies two equations_source_chimie" << finl;
      barrier();
      exit();
    }

  Probleme_base& pb = mon_probleme.valeur();
  DoubleTab& champ1 = pb.getset_equation_by_name(equations_source_chimie_[0]).inconnue()->valeurs();
  DoubleTab& champ2 = pb.getset_equation_by_name(equations_source_chimie_[1]).inconnue()->valeurs();
  DoubleTab& champ3 = inconnue()->valeurs();

  const double dt = pb.schema_temps().pas_de_temps();
  const int nb_faces = champ1.dimension(0);

  const IntTab& face_voisins = ref_cast(Domaine_VF, domaine_dis()).face_voisins();
  const DoubleTab *champ_nu_t = 0;
  if (nom_equation_nu_t_ != "??")
    {
      const Equation_base& eq =pb.get_equation_by_name(nom_equation_nu_t_);
      const RefObjU& modele_turbulence_hydr = eq.get_modele(TURBULENCE);
      const Modele_turbulence_hyd_base& le_modele = ref_cast(Modele_turbulence_hyd_base,modele_turbulence_hydr.valeur());
      champ_nu_t = &(le_modele.viscosite_turbulente()->valeurs());
    }

  const DoubleTab& tab_diffusivite = constituant().diffusivite_constituant()->valeurs();
  double coeff_diffusivite = tab_diffusivite(0,0);
  const DoubleVect& volumes_entrelaces = ref_cast(Domaine_VF, domaine_dis()).volumes_entrelaces();
  const int dim = Objet_U::dimension;

  const Fluide_Diphasique& fluide = ref_cast(Fluide_Diphasique, pb.milieu());
  const Fluide_Incompressible& phase_0 = fluide.fluide_phase(0);
  const Fluide_Incompressible& phase_1 = fluide.fluide_phase(1);
  const DoubleTab& tab_nu_phase_0 = phase_0.viscosite_cinematique()->valeurs();
  const DoubleTab& tab_nu_phase_1 = phase_1.viscosite_cinematique()->valeurs();
  const double nu_phase_0 = tab_nu_phase_0(0,0);
  const double nu_phase_1 = tab_nu_phase_1(0,0);
  const double delta_nu = nu_phase_1 - nu_phase_0;

  DoubleTab indic_faces;
  calculer_indicatrice_comme(ref_eq_transport_->inconnue().valeur(),
                             indic_faces, inconnue().valeur());

  for (int face = 0; face < nb_faces; face++)
    {
      double nu_t = 0.;
      if (champ_nu_t)
        {
          int elem1 = face_voisins(face,0);
          int elem2 = face_voisins(face,1);
          if (elem1 < 0) elem1 = elem2;
          if (elem2 < 0) elem2 = elem1;
          nu_t = ((*champ_nu_t)(elem1) + (*champ_nu_t)(elem2)) * 0.5;
        }
      double c1 = std::max(0., champ1(face));
      double c2 = std::max(0., champ2(face));
      double c3 = std::max(0., champ3(face));

      // ES243900
      // Modele 1 classique
      // Finite rate formulation.
      if (modele_cinetique_ == 1)
        {
          double omega = constante_cinetique_ * (1. + constante_cinetique_nu_t_ * nu_t) * c1 * c2;
          double Rdt = std::min(omega * dt, c1);
          Rdt = std::min(Rdt, c2);

          /* Modif des concentrations */
          champ1(face) = c1 - Rdt;
          champ2(face) = c2 - Rdt;
          champ3(face) = c3 + Rdt;
        }

      // Modele 2 : eddy disspertion concept (EDC) Bertrand et al. 2016 Chemical Engineering Journal
      // traitement particulier pour n'avoir pas de c < 0
      else if (modele_cinetique_ == 2)
        {
          // Besoin du modele LES !!!!
          Cerr << "MODELE_CINETIQUE 2 requires a turbulent Concentration problem" << finl;
          Cerr << "Replace Convection_Diffusion_Concentration_FT_Disc by Convection_Diffusion_Concentration_Turbulent_FT_Disc" << finl;
          barrier();
          Process::exit();
        }

      // Modele 3 : eddy disspertion concept (EDC) Bertrand et al. 2016 Chemical Engineering Journal
      // traitement particulier pour n'avoir pas de c < 0
      else if (modele_cinetique_ == 3)
        {
          const double indic = indic_faces[face];
          const double nu = indic * delta_nu + nu_phase_0;

          // On recalcule la longueur caracteristique de maille
          double l_carac = 2.0*pow(6.*volumes_entrelaces(face),1./double(dim)) ;

          const double Ci = 0.07; // from Yoshizawa's model
          double schm = nu/coeff_diffusivite;

          double tm = l_carac*l_carac*Ci*schm/(nu + nu_t);

          double omega = std::min(c1,c2);
          omega = omega/std::max(dt,tm);
          double Rdt = omega * dt;

          /* Modif des concentrations */
          champ1(face) = c1 - Rdt;
          champ2(face) = c2 - Rdt;
          champ3(face) = c3 + Rdt;
        }

      // Modele 4 : eddy disspertion concept (EDC) Bertrand et al. 2016 Chemical Engineering Journal
      // traitement particulier pour n'avoir pas de c < 0
      // conclusion ==> formulation pour que modeles 2 et 3 soient equivalents
      else if (modele_cinetique_ == 4)
        {
          // On recalcule la longueur caracteristique de maille
          double l_carac = 2.0*pow(6.*volumes_entrelaces(face),1./double(dim)) ;
          // merge de modele 2 et 3  ......
          double tm = l_carac*l_carac/nu_t;

          double omega = std::min(c1,c2);
          omega = omega/std::max(dt,tm);
          double Rdt = omega * dt;

          /* Modif des concentrations */
          champ1(face) = c1 - Rdt;
          champ2(face) = c2 - Rdt;
          champ3(face) = c3 + Rdt;
        }

      else
        {
          Cerr << "MODELE_CINETIQUE not implemented ==> Please specify model 1 -- 4" << finl;
          barrier();
          Process::exit();
        }

    }
  champ1.echange_espace_virtuel();
  champ2.echange_espace_virtuel();
  champ3.echange_espace_virtuel();
}

// Multiplie les valeurs des faces marquees par facteur
//  et calcule l'integrale des valeurs avant modifications
void Convection_Diffusion_Concentration_FT_Disc::multiplier_valeurs_faces(const ArrOfBit marqueurs_faces, double facteur, double& integrale_avant, DoubleTab& tab)
{
  const Domaine_VEF&     domaine_vef = ref_cast(Domaine_VEF, domaine_dis());
  const DoubleVect& volumes = domaine_vef.volumes_entrelaces();
  const DoubleVect& volumes_cl = ref_cast(Domaine_Cl_VEF, domaine_Cl_dis().valeur()).volumes_entrelaces_Cl();
  const int nb_faces = domaine_vef.nb_faces();
  const int nb_faces_non_std = domaine_vef.nb_faces_non_std();
  const ArrOfInt& faces_doubles = domaine_vef.faces_doubles();
  int modif = (facteur != 1.);
  double integrale = 0.;
  for (int i = 0; i < nb_faces; i++)
    {
      if (marqueurs_faces[i])
        {
          double vol = (i < nb_faces_non_std) ? volumes_cl[i] : volumes[i];
          double f = faces_doubles[i] ? 0.5 : 1.;
          double x = tab(i);
          integrale += vol * f * x;
          if (modif)
            tab(i) = facteur * x;
        }
    }
  integrale_avant = mp_sum(integrale);
}

void Convection_Diffusion_Concentration_FT_Disc::marquer_faces_sous_domaine(const Nom& nom_sous_domaine,
                                                                            ArrOfBit& marqueur,
                                                                            int sans_faces_non_std_volume_nul) const
{
  const Domaine_VEF& domaine_vef = ref_cast(Domaine_VEF, domaine_dis());
  const IntTab& elem_faces = domaine_vef.elem_faces();
  const int nb_faces_tot = domaine_vef.nb_faces_tot();
  const int nb_faces_non_std = domaine_vef.nb_faces_non_std();
  const int nb_faces_elem = elem_faces.dimension(1);
  marqueur.resize_array(nb_faces_tot);
  marqueur = 0;
  const Sous_Domaine& sous_domaine = domaine_dis().domaine().ss_domaine(nom_sous_domaine);
  const int nb_elem_sous_domaine = sous_domaine.nb_elem_tot();
  const DoubleVect& volumes_cl = ref_cast(Domaine_Cl_VEF, domaine_Cl_dis().valeur()).volumes_entrelaces_Cl();
  int i;
  for (i = 0; i < nb_elem_sous_domaine; i++)
    {
      if (sans_faces_non_std_volume_nul && i < nb_faces_non_std)
        {
          if (volumes_cl[i] == 0.)
            // Ne pas marquer les faces non std de volume nul
            continue;
        }
      const int elem = sous_domaine[i];
      for (int j = 0; j < nb_faces_elem; j++)
        {
          const int face = elem_faces(elem,j);
          marqueur.setbit(face);
        }
    }
}

/*! @brief appel a mettre_a_jour() de l'ancetre et application de l' "option_" choisie pour traiter la presence d'une interface.
 *
 */
void Convection_Diffusion_Concentration_FT_Disc::mettre_a_jour(double temps)
{
  Debog::verifier("Convection_Diffusion_Concentration_FT_Disc::mettre_a_jour c1", inconnue()->valeurs());
  Convection_Diffusion_Concentration::mettre_a_jour(temps);

  switch(option_)
    {
    case RIEN:
      break;
    case RAMASSE_MIETTES_SIMPLE:
      ramasse_miette_simple(temps);
      break;
    default:
      Cerr << "Internal error for method Convection_Diffusion_Concentration_FT_Disc::mettre_a_jour" << finl;
      barrier();
      exit();
    }

  if (equations_source_chimie_.size() > 0)
    {
      mettre_a_jour_chimie();
    }

  double integrale_sortie = 0.;
  double integrale0 = 0.;
  double integrale1 = 0.;

  if (nom_domaine_sortie_ != "??")
    {
      ArrOfBit marqueur;
      marquer_faces_sous_domaine(nom_domaine_sortie_, marqueur,
                                 0 /* annuler les valeurs sur les faces non std */);
      DoubleTab& tab = inconnue()->valeurs();
      multiplier_valeurs_faces(marqueur, 0., integrale_sortie, tab);
    }
  {
    DoubleTab indic_faces;
    calculer_indicatrice_comme(ref_eq_transport_->inconnue().valeur(),
                               indic_faces, inconnue().valeur());
    const Domaine_VEF&     domaine_vef = ref_cast(Domaine_VEF, domaine_dis());
    const DoubleVect& volumes = domaine_vef.volumes_entrelaces();
    const DoubleVect& volumes_cl = ref_cast(Domaine_Cl_VEF, domaine_Cl_dis().valeur()).volumes_entrelaces_Cl();
    const int nb_faces = domaine_vef.nb_faces();
    const int nb_faces_non_std = domaine_vef.nb_faces_non_std();
    const ArrOfInt& faces_doubles = domaine_vef.faces_doubles();
    const DoubleTab& inco = inconnue()->valeurs();
    for (int i = 0; i < nb_faces; i++)
      {
        double vol = (i < nb_faces_non_std) ? volumes_cl[i] : volumes[i];
        double f = faces_doubles[i] ? 0.5 : 1.;
        double indic = indic_faces(i);
        double x = inco(i) * f * vol;
        integrale0 += x * (1.-indic);
        integrale1 += x * indic;
      }
    integrale0 = mp_sum(integrale0);
    integrale1 = mp_sum(integrale1);
  }
  const double dt = probleme().schema_temps().pas_de_temps();
  if (je_suis_maitre())
    {
      Nom nom_fic(le_nom());
      nom_fic += "_bilan.out";
      SFichier f(nom_fic, ios::out | ios::app);
      char s[1000];
      snprintf(s, 1000, "%16.8f %16.8g %16.8g %16.8g %16.8g",
               temps,
               integrale0,
               integrale1,
               integrale_sortie,
               integrale_sortie / dt);
      f << s << finl;
    }
  Debog::verifier("Convection_Diffusion_Concentration_FT_Disc::mettre_a_jour c2", inconnue()->valeurs());
}
