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
// File:        Probleme_FT_Disc_gen.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/15
//
//////////////////////////////////////////////////////////////////////////////

#include <Convection_Diffusion_Temperature_FT_Disc.h>
#include <Convection_Diffusion_Concentration.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <Triple_Line_Model_FT_Disc.h>
#include <Dirichlet_paroi_defilante.h>
#include <Probleme_FT_Disc_gen.h>
#include <Dirichlet_paroi_fixe.h>
#include <Constituant.h>
#include <TRUST_List.h>
#include <Chimie.h>

Implemente_instanciable(Probleme_FT_Disc_gen, "Probleme_FT_Disc_gen|Pb_FrontTracking_Disc", Pb_Fluide_base);

// XD listdeuxmots_acc listobj nul 1 deuxmots 0 List of groups of two words (with curly brackets).

// XD Pb_FrontTracking_Disc problem_read_generic probleme_ft_disc_gen -1 The generic Front-Tracking problem in the discontinuous version. It differs from the rest of the TRUST code : The problem does not state the number of equations that are enclosed in the problem. Two equations are compulsory : a momentum balance equation (alias Navier-Stokes equation) and an interface tracking equation. The list of equations to be solved is declared in the beginning of the data file. Another difference with more classical TRUST data file, lies in the fluids definition. The two-phase fluid (Fluide_Diphasique) is made with two usual single-phase fluids (Fluide_Incompressible). As the list of equations to be solved in the generic Front-Tracking problem is declared in the data file and not pre-defined in the structure of the problem, each equation has to be distinctively associated with the problem with the Associer keyword.
// XD attr solved_equations listdeuxmots_acc solved_equations 0 List of sovled equations in the form 'equation_type' 'equation_alias'
// XD attr fluide_incompressible fluide_incompressible fluide_incompressible 1 The fluid medium associated with the problem.
// XD attr fluide_diphasique fluide_diphasique fluide_diphasique 1 The diphasic fluid medium associated with the problem.
// XD attr constituant constituant constituant 1 Constituent.
// XD attr Triple_Line_Model_FT_Disc Triple_Line_Model_FT_Disc Triple_Line_Model_FT_Disc 1 not_set

Sortie& Probleme_FT_Disc_gen::printOn(Sortie& os) const { return os; }
Entree& Probleme_FT_Disc_gen::readOn(Entree& is) { return Pb_Fluide_base::readOn(is); }

int Probleme_FT_Disc_gen::associer_(Objet_U& ob)
{
  if (sub_type(Chimie, ob))
    {
      la_chimie_ = ref_cast(Chimie, ob);
      return 1;
    }
  else
    return Pb_Fluide_base::associer_(ob);
}

void Probleme_FT_Disc_gen::lire_solved_equations(Entree& is)
{
  /* Step 1 : special FT : read the list of equations to solve ... */
  // Here are all possible equations !!!
  Noms noms_eq, noms_eq_maj;
  Type_info::les_sous_types(Nom("Equation_base"), noms_eq);
  for (auto &itr : noms_eq) noms_eq_maj.add(Motcle(itr)); //ha ha ha

  Motcle read_mc;
  Nom nom_eq;

  is >> read_mc;

  if (read_mc != "SOLVED_EQUATIONS")
    {
      Cerr << "Error in Probleme_FT_Disc_gen::typer_lire_milieu !!! We expected reading the SOLVED_EQUATIONS bloc instead of " << read_mc << " !!!" << finl;
      Cerr << "Fix your data file !!!" << finl;
      Process::exit();
    }

  is >> read_mc;
  if (read_mc != "{")
    {
      Cerr << "Error in Probleme_FT_Disc_gen::typer_lire_milieu !!! We expected { instead of " << read_mc << " !!!" << finl;
      Cerr << "Fix your data file !!!" << finl;
      Process::exit();
    }

  std::vector<Nom> eq_types, eq_name;

  for (is >> read_mc; read_mc != "}"; is >> read_mc)
    {
      if (noms_eq_maj.rang(read_mc) == -1 || (!read_mc.contient("_FT") && !read_mc.debute_par("CONVECTION_DIFFUSION_CONCENTRATION") && !read_mc.debute_par("CONVECTION_DIFFUSION_TEMPERATURE")))
        {
          Cerr << "Error in Probleme_FT_Disc_gen::typer_lire_milieu !!! The equation " << read_mc << " could not be used with a problem of type " << que_suis_je() << " !!!" << finl;
          Cerr << "You can only use the following equations :" << finl;
          for (auto &itr : noms_eq_maj)
            if (itr.contient("_FT") || itr.debute_par("CONVECTION_DIFFUSION_CONCENTRATION") || itr.debute_par("CONVECTION_DIFFUSION_TEMPERATURE"))
              Cerr << "  - " << itr << finl;
          Process::exit();
        }

      eq_types.push_back(read_mc);
      is >> nom_eq;
      eq_name.push_back(nom_eq);
    }

  if (eq_types.size() != eq_name.size())
    {
      Cerr << "Error in Probleme_FT_Disc_gen::typer_lire_milieu !!! The number of strings read in the bloc SOLVED_EQUATIONS is not correct !!!" << finl;
      Cerr << "Fix your data file !!!" << finl;
      Process::exit();
    }

  /* Step 2 : add equations to the list ... */
  /* Add NAVIER_STOKES_FT_DISC at first */
  for (int i = 0; i < static_cast<int>(eq_types.size()); i++)
    if (eq_types[i] == "NAVIER_STOKES_FT_DISC")
      add_FT_equation(eq_name[i], eq_types[i]);

  /* Add Transport_Interfaces at second */
  for (int i = 0; i < static_cast<int>(eq_types.size()); i++)
    if (eq_types[i].debute_par("TRANSPORT_INTERFACES"))
      add_FT_equation(eq_name[i], eq_types[i]);

  /* Add the remaining */
  for (int i = 0; i < static_cast<int>(eq_types.size()); i++)
    if (eq_types[i] != "NAVIER_STOKES_FT_DISC" && !eq_types[i].debute_par("TRANSPORT_INTERFACES"))
      add_FT_equation(eq_name[i], eq_types[i]);
}

void Probleme_FT_Disc_gen::typer_lire_milieu(Entree& is)
{
  bool needs_constituant = false;

  auto& list_stl = equations_.get_stl_list();
  for (auto& itr : list_stl)
    if (Motcle(itr->que_suis_je()).contient("CONCENTRATION"))
      {
        needs_constituant = true; // pb contient concentration !
        break;
      }

  const int nb_milieu = needs_constituant ? 2 : 1;
  le_milieu_.resize(nb_milieu);

  for (int i = 0; i < nb_milieu; i++)
    {
      is >> le_milieu_[i]; // On commence par la lecture du milieu
      associer_milieu_base(le_milieu_[i].valeur()); // On l'associe a chaque equations (methode virtuelle pour chaque pb ...)
    }

  // Milieu(x) lu(s) ... Lets go ! On discretise les equations
  discretiser_equations();

  // remontee de l'inconnue vers le milieu
  for (int i = 0; i < nombre_d_equations(); i++)
    equation(i).associer_milieu_equation();

  /* On discretise le milieu de l'eq 1 (NS) */
  equations_(0)->milieu().discretiser((*this), la_discretisation_.valeur());

  if (needs_constituant) // Et l'eq de concentration !
    for (auto& itr : equations_.get_stl_list())
      if (sub_type(Convection_Diffusion_Concentration, itr.valeur()))
        {
          itr->milieu().discretiser((*this), la_discretisation_.valeur());
          break;
        }
}

Entree& Probleme_FT_Disc_gen::lire_equations(Entree& is, Motcle& dernier_mot)
{
  Pb_Fluide_base::lire_equations(is, dernier_mot);

  if (dernier_mot == "TRIPLE_LINE_MODEL_FT_DISC")
    {
      is >> tcl_; // on lit
      tcl_.associer_pb(*this);
      tcl_.initialize(); // on initialise !!
      is >> dernier_mot; // et on lit le dernier mot
    }
  return is;
}

void Probleme_FT_Disc_gen::add_FT_equation(const Nom& eq_name, const Nom& eq_type)
{
  equations_.add(OWN_PTR(Equation_base)());
  equations_.dernier().typer(eq_type);
  Equation_base& eq = equations_.dernier().valeur();
  eq.associer_pb_base(*this);
  eq.associer_sch_tps_base(le_schema_en_temps_);
  eq.nommer(Nom(eq_name));
  Cerr << "Equation " << eq_type << " added to the list and renamed to : " << eq_name << " ..." << finl;
}

/*! @brief Verifie que le milieu est de type Fluide_Diphasique et associe le milieu aux equations.
 *
 * Precondition: Toutes les equations doivent avoir ete associees.
 */
void Probleme_FT_Disc_gen::associer_milieu_base(const Milieu_base& un_milieu)
{
  // Le milieu est associe aux equations en fonction de son type :
  // Si le milieu est un "Constituant", on l'associe uniquement aux equations de convection_diffusion_concentration,
  // Sinon, on associe le milieu a toutes les autres equations.
  const int n = nombre_d_equations();
  const int is_constituant = sub_type(Constituant, un_milieu);
  for (int i = 0; i < n; i++)
    {
      Equation_base& eq = equation(i);
      int is_conv_diff = sub_type(Convection_Diffusion_Concentration, eq);
      if ((is_conv_diff && is_constituant) || (!is_conv_diff && !is_constituant))
        eq.associer_milieu_base(un_milieu);
    }
}

const Equation_base& Probleme_FT_Disc_gen::get_equation_by_name(const Nom& un_nom) const
{
  const int n = nombre_d_equations();
  const Equation_base *eq = 0;
  int i;
  for (i = 0; i < n; i++)
    {
      eq = &equation(i);
      if (eq->le_nom() == un_nom)
        break;
    }
  if (i == n)
    {
      Cerr << "Erreur dans Probleme_FT_Disc_gen::get_equation_by_name : " << un_nom << " n'est pas le nom d'une equation !!" << finl;
      Cerr << "Les equations du problemes sont les suivantes :" << finl;
      for (i = 0; i < n; i++)
        Cerr << equation(i).le_nom() << finl;
      Process::exit();
    }
  return *eq;
}

Equation_base& Probleme_FT_Disc_gen::getset_equation_by_name(const Nom& un_nom)
{
  const int n = nombre_d_equations();
  Equation_base *eq = 0;
  int i;
  for (i = 0; i < n; i++)
    {
      eq = &equation(i);
      if (eq->le_nom() == un_nom)
        break;
    }
  if (i == n)
    {
      Cerr << "Erreur dans Probleme_FT_Disc_gen::getset_equation_by_name : " << un_nom << " n'est pas le nom d'une equation !!" << finl;
      Cerr << "Les equations du problemes sont les suivantes :" << finl;
      for (i = 0; i < n; i++)
        Cerr << equation(i).le_nom() << finl;
      Process::exit();
    }
  return *eq;
}

const Transport_Interfaces_FT_Disc& Probleme_FT_Disc_gen::equation_interfaces(const Motcle& un_nom) const
{
  const int n = nombre_d_equations();
  const Transport_Interfaces_FT_Disc *eq_ptr = 0;
  for (int i = 0; i < n; i++)
    {
      const Equation_base& eq = equation(i);
      if (un_nom == eq.le_nom() && sub_type(Transport_Interfaces_FT_Disc, eq))
        {
          eq_ptr = &ref_cast(Transport_Interfaces_FT_Disc, eq);
          break;
        }
    }
  if (eq_ptr == 0)
    {
      if (Process::je_suis_maitre())
        {
          Cerr << "Erreur dans Probleme_FT_Disc_gen::equation_interfaces : Le probleme ne contient pas d'equation Transport_Interfaces_FT_Disc de nom " << un_nom << finl;
          Cerr << "Liste des equations du probleme:" << finl;
          for (int i = 0; i < n; i++)
            {
              const Equation_base& eq = equation(i);
              Cerr << eq.que_suis_je() << "  -->  " << eq.le_nom() << finl;
            }
        }
      Process::exit();
    }
  return *eq_ptr;
}

const Navier_Stokes_FT_Disc& Probleme_FT_Disc_gen::equation_hydraulique(const Motcle& un_nom) const
{
  const int n = nombre_d_equations();
  const Navier_Stokes_FT_Disc *eq_ptr = 0;
  for (int i = 0; i < n; i++)
    {
      const Equation_base& eq = equation(i);
      if (un_nom == eq.le_nom() && sub_type(Navier_Stokes_FT_Disc, eq))
        {
          eq_ptr = &ref_cast(Navier_Stokes_FT_Disc, eq);
          break;
        }
    }
  if (eq_ptr == 0)
    {
      if (Process::je_suis_maitre())
        {
          Cerr << "Erreur dans Probleme_FT_Disc_gen::equation_hydraulique : Le probleme ne contient pas d'equation Navier_Stokes_FT_Disc de nom " << un_nom << finl;
          Cerr << "Liste des equations du probleme:" << finl;
          for (int i = 0; i < n; i++)
            {
              const Equation_base& eq = equation(i);
              Cerr << eq.que_suis_je() << "  -->  " << eq.le_nom() << finl;
            }
        }
      Process::exit();
    }
  return *eq_ptr;
}

double Probleme_FT_Disc_gen::calculer_pas_de_temps(void) const
{
  // prend le min des pas de temps de chaque equation
  double dt = Pb_Fluide_base::calculer_pas_de_temps();

  if (la_chimie_.non_nul())
    {
      dt = std::min(dt, la_chimie_->calculer_pas_de_temps());
      dt = mp_min(dt);
    }
  return dt;
}

void Probleme_FT_Disc_gen::mettre_a_jour(double temps)
{
  if (schema_temps().que_suis_je() == "RK3_FT")
    {
      int nb_eqn = nombre_d_equations();
      for (int i = 2; i < nb_eqn; i++)
        equation(i).mettre_a_jour(temps);

      les_postraitements_.mettre_a_jour(temps);
    }
  else
    {
      Pb_Fluide_base::mettre_a_jour(temps);
    }
  if (la_chimie_.non_nul())
    la_chimie_->mettre_a_jour(temps);
}

void Probleme_FT_Disc_gen::completer(void)
{
  Pb_Fluide_base::completer();

  if (la_chimie_.non_nul())
    {
      la_chimie_->discretiser(*this); /* XXX : Elie Saikali : j'ai retarde ca et mis ici ... */
      la_chimie_->completer(*this);
    }

  if (tcl_.is_activated()) tcl_.completer();
}

bool Probleme_FT_Disc_gen::updateGivenFields()
{
  // add: fill the lists of TCL model if actived, Which will be used for updated the BCs
  if (tcl().is_activated())
    {
      ArrOfInt elems_with_CL_contrib, faces_with_CL_contrib;
      ArrOfDouble mpoint_from_CL, Q_from_CL;
      tcl().compute_TCL_fluxes_in_all_boundary_cells(elems_with_CL_contrib, faces_with_CL_contrib, mpoint_from_CL, Q_from_CL);
    }
  return Probleme_base::updateGivenFields();
}
