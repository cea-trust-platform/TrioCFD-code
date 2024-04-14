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

Implemente_instanciable(Probleme_FT_Disc_gen, "Probleme_FT_Disc_gen", Pb_Fluide_base);

Sortie& Probleme_FT_Disc_gen::printOn(Sortie& os) const { return os; }
Entree& Probleme_FT_Disc_gen::readOn(Entree& is) { return Pb_Fluide_base::readOn(is); }

int Probleme_FT_Disc_gen::associer_(Objet_U& ob)
{
  if (sub_type(Chimie, ob))
    {
      la_chimie_ = ref_cast(Chimie, ob);
      return 1;
    }
  if (sub_type(Triple_Line_Model_FT_Disc, ob))
    {
      associate_triple_line_model(ref_cast(Triple_Line_Model_FT_Disc, ob));
      return 1;
    }
  return Pb_Fluide_base::associer_(ob);
}

void Probleme_FT_Disc_gen::typer_lire_milieu(Entree& is)
{
  /* Step 1 : special FT : read the list of equations to solve ... */
  // Here are all possible equations !!!
  Noms noms_eq, noms_eq_maj;
  Type_info::les_sous_types(Nom("Equation_base"), noms_eq);
  for (auto &itr : noms_eq) noms_eq_maj.add(Motcle(itr)); //ha ha ha

  Motcle read_mc;
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

  std::map<std::string, std::string> solved_eqs;
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

      Nom nom_eq;
      is >> nom_eq;
      solved_eqs.insert( { nom_eq.getString(), read_mc.getString() });
    }

  /* Step 2 : add equations to the list ... */
  /* Add NAVIER_STOKES_FT_DISC at first */
  for (auto& itr : solved_eqs)
    if (itr.second == "NAVIER_STOKES_FT_DISC")
      add_FT_equation(itr.first, itr.second);

  /* Add Transport_Interfaces at second */
  for (auto& itr : solved_eqs)
    if (Nom(itr.second).debute_par("TRANSPORT_INTERFACES"))
      add_FT_equation(itr.first, itr.second);

  /* Add the remaining */
  for (auto& itr : solved_eqs)
    if (itr.second != "NAVIER_STOKES_FT_DISC" && !Nom(itr.second).debute_par("TRANSPORT_INTERFACES"))
      add_FT_equation(itr.first, itr.second);

  /* Step 3 : read medium(media) ... */
  const std::string conc = "CONCENTRATION";
  bool needs_constituant = false;

  for (auto& itr : solved_eqs)
    if (itr.second.find(conc) != std::string::npos)
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

void Probleme_FT_Disc_gen::add_FT_equation(const std::string& eq_name, const std::string& eq_type)
{
  equations_.add(Equation());
  equations_.dernier().typer(eq_type.c_str());
  Equation_base& eq = equations_.dernier().valeur();
  eq.associer_pb_base(*this);
  eq.associer_sch_tps_base(le_schema_en_temps_);
  eq.nommer(Nom(eq_name));
  Cerr << "Equation " << eq_type << " added to the list and renamed to : " << eq_name << " ..." << finl;
}

void Probleme_FT_Disc_gen::associate_triple_line_model(Triple_Line_Model_FT_Disc& tcl_1)
{
  // TODO: A verifier. Working fine?
  // tcl.initialize(); // ca ne marche pas, apres le '=' ci-dessous, le tcl_.elems_.smart_resize_ est toujours a 0!!!
  tcl_ = tcl_1;
  tcl_.initialize();
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
      dt = std::min(dt, la_chimie_.valeur().calculer_pas_de_temps());
      dt = mp_min(dt);
    }
  return dt;
}

void Probleme_FT_Disc_gen::discretiser(Discretisation_base& dis)
{
  Pb_Fluide_base::discretiser(dis);
  if (la_chimie_.non_nul())
    la_chimie_.valeur().discretiser(*this);
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
    la_chimie_.valeur().mettre_a_jour(temps);
}

void Probleme_FT_Disc_gen::completer(void)
{
  Pb_Fluide_base::completer();
  if (la_chimie_.non_nul())
    la_chimie_.valeur().completer(*this);

  if (tcl_.is_activated())
    tcl_.completer();
}

bool Probleme_FT_Disc_gen::updateGivenFields()
{
  // add: fill the lists of TCL model if actived, Which will be used for updated the BCs
  if (tcl().is_activated())
    {

      ArrOfInt elems_with_CL_contrib;
      ArrOfInt faces_with_CL_contrib;
      ArrOfDouble mpoint_from_CL;
      ArrOfDouble Q_from_CL;
      tcl().compute_TCL_fluxes_in_all_boundary_cells(elems_with_CL_contrib, faces_with_CL_contrib, mpoint_from_CL, Q_from_CL);
    }
  return Probleme_base::updateGivenFields();
}
