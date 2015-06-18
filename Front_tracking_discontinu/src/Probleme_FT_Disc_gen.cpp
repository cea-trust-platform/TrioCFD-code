/****************************************************************************
* Copyright (c) 2015, CEA
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
// Directory:   $TRUST_ROOT/src/Front_tracking_discontinu
// Version:     /main/15
//
//////////////////////////////////////////////////////////////////////////////

#include <Probleme_FT_Disc_gen.h>
#include <List_Nom.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <Convection_Diffusion_Temperature_FT_Disc.h>
#include <Constituant.h>
#include <Convection_Diffusion_Concentration.h>
#include <Chimie.h>

Implemente_instanciable(Probleme_FT_Disc_gen,"Probleme_FT_Disc_gen",Pb_qdm_fluide);

Implemente_ref(Probleme_FT_Disc_gen);

Entree& Probleme_FT_Disc_gen::readOn(Entree& is)
{
  return Pb_qdm_fluide::readOn(is);
}

Sortie& Probleme_FT_Disc_gen::printOn(Sortie& os) const
{
  return os;
}

// Pour construire le probleme, il faut lui associer les differentes
// equations dans cet ordre :
//
//  Probleme_FT_Disc_gen pb
//  Navier_Stokes_FT_Disc eq_ns
//  Associer pb eq_ns
//  Transport_Interfaces_FT_Disc eq_transport
//  Associer pb eq_transport
//  ...
//  Discretiser pb
//  Lire pb {
//     comme d'habitude
//  }

int Probleme_FT_Disc_gen::associer_(Objet_U& ob)
{
  if (sub_type(Equation_base, ob))
    {
      associer_equation(ref_cast(Equation_base, ob));
      return 1;
    }
  if (sub_type(Chimie,ob))
    {
      la_chimie_=ref_cast(Chimie,ob);
      return 1;
    }
  return Pb_qdm_fluide::associer_(ob);
}

void Probleme_FT_Disc_gen::associer_equation(Equation_base& eq)
{
  // Ajoute une reference a l'equation dans le tableau

  if (equations_.size()==1)
    if (!sub_type(Navier_Stokes_FT_Disc,equation(0)))
      {
        Cerr<<"L'equation d'hydraulique n'a pas ete associee comme premiere equation"<<finl;
        exit();
      }
  equations_.add(eq);
}

// Description: Verifie que le milieu est de type Fluide_Diphasique
//  et associe le milieu aux equations.
// Precondition: Toutes les equations doivent avoir ete associees.
void Probleme_FT_Disc_gen::associer_milieu_base(const Milieu_base& un_milieu)
{
  // Le milieu est associe aux equations en fonction de son type :
  // Si le milieu est un "Constituant", on l'associe uniquement
  //  aux equations de convection_diffusion_concentration,
  // Sinon, on associe le milieu a toutes les autres equations.
  const int n = nombre_d_equations();
  const int is_constituant = sub_type(Constituant, un_milieu);
  int i;
  for (i = 0; i < n; i++)
    {
      Equation_base& eq = equation(i);
      int is_conv_diff = sub_type(Convection_Diffusion_Concentration, eq);
      if ((is_conv_diff && is_constituant) || (!is_conv_diff && !is_constituant))
        eq.associer_milieu_base(un_milieu);
    }
}

int Probleme_FT_Disc_gen::nombre_d_equations(void) const
{
  return equations_.size();
}

const Equation_base& Probleme_FT_Disc_gen::equation(int i) const
{
  return equations_[i].valeur();
}

Equation_base& Probleme_FT_Disc_gen::equation(int i)
{
  return equations_[i].valeur();
}

const Equation_base& Probleme_FT_Disc_gen::get_equation_by_name(const Nom& un_nom) const
{
  const int n = nombre_d_equations();
  const Equation_base * eq = 0;
  int i;
  for (i = 0; i < n; i++)
    {
      eq = &equation(i);
      if (eq->le_nom() == un_nom)
        break;
    }
  if (i == n)
    {
      Cerr << "Erreur dans Probleme_FT_Disc_gen::get_equation_by_name :\n";
      Cerr << un_nom << " n'est pas le nom d'une equation\n";
      Cerr << "Les equations du problemes sont les suivantes :\n";
      for (i = 0; i < n; i++)
        Cerr << equation(i).le_nom() << finl;
      Cerr << "(Attention : le format de lecture de Probleme_FT_Disc_gen est particulier:\n";
      Cerr << " Navier_Stokes_FT_Disc mon_equation_hydraulique\n";
      Cerr << " Associer pb mon_equation_hydraulique\n";
      Cerr << " Lire pb {\n";
      Cerr << "   mon_equation_hydraulique { ... }\n";
      Cerr << " }" << finl;
      assert(0);
      exit();
    }
  return *eq;
}

Equation_base& Probleme_FT_Disc_gen::getset_equation_by_name(const Nom& un_nom)
{
  const int n = nombre_d_equations();
  Equation_base * eq = 0;
  int i;
  for (i = 0; i < n; i++)
    {
      eq = &equation(i);
      if (eq->le_nom() == un_nom)
        break;
    }
  if (i == n)
    {
      Cerr << "Erreur dans Probleme_FT_Disc_gen::get_equation_by_name :\n";
      Cerr << un_nom << " n'est pas le nom d'une equation\n";
      Cerr << "Les equations du problemes sont les suivantes :\n";
      for (i = 0; i < n; i++)
        Cerr << equation(i).le_nom() << finl;
      Cerr << "(Attention : le format de lecture de Probleme_FT_Disc_gen est particulier:\n";
      Cerr << " Navier_Stokes_FT_Disc mon_equation_hydraulique\n";
      Cerr << " Associer pb mon_equation_hydraulique\n";
      Cerr << " Lire pb {\n";
      Cerr << "   mon_equation_hydraulique { ... }\n";
      Cerr << " }" << finl;
      assert(0);
      exit();
    }
  return *eq;
}

const Transport_Interfaces_FT_Disc&
Probleme_FT_Disc_gen::equation_interfaces(const Motcle& un_nom) const
{
  const int n = nombre_d_equations();
  int i;
  const Transport_Interfaces_FT_Disc * eq_ptr = 0;
  for (i = 0; i < n; i++)
    {
      const Equation_base& eq = equation(i);
      if (un_nom == eq.le_nom()
          && sub_type(Transport_Interfaces_FT_Disc, eq))
        {
          eq_ptr = & ref_cast(Transport_Interfaces_FT_Disc, eq);
          break;
        }
    }
  if (eq_ptr == 0)
    {
      if (Process::je_suis_maitre())
        {
          Cerr << "Erreur dans Probleme_FT_Disc_gen::equation_interfaces(const Motcle & nom)\n"
               << " Le probleme ne contient pas d'equation Transport_Interfaces_FT_Disc\n"
               << " de nom " << un_nom << finl;
          Cerr << "Liste des equations du probleme:" << finl;
          for (i = 0; i < n; i++)
            {
              const Equation_base& eq = equation(i);
              Cerr << eq.que_suis_je() << " " << eq.le_nom() << finl;
            }
        }
      barrier();
      assert(0);
      exit();
    }
  return *eq_ptr;
}

const Navier_Stokes_FT_Disc&
Probleme_FT_Disc_gen::equation_hydraulique(const Motcle& un_nom) const
{
  const int n = nombre_d_equations();
  int i;
  const Navier_Stokes_FT_Disc * eq_ptr = 0;
  for (i = 0; i < n; i++)
    {
      const Equation_base& eq = equation(i);
      if (un_nom == eq.le_nom()
          && sub_type(Navier_Stokes_FT_Disc, eq))
        {
          eq_ptr = & ref_cast(Navier_Stokes_FT_Disc, eq);
          break;
        }
    }
  if (eq_ptr == 0)
    {
      if (Process::je_suis_maitre())
        {
          Cerr << "Erreur dans Probleme_FT_Disc_gen::equation_hydraulique(const Motcle & nom)\n"
               << " Le probleme ne contient pas d'equation Navier_Stokes_FT_Disc\n"
               << " de nom " << un_nom << finl;
          Cerr << "Liste des equations du probleme:" << finl;
          for (i = 0; i < n; i++)
            {
              const Equation_base& eq = equation(i);
              Cerr << eq.que_suis_je() << " " << eq.le_nom() << finl;
            }
        }
      barrier();
      assert(0);
      exit();
    }
  return *eq_ptr;
}

double Probleme_FT_Disc_gen::calculer_pas_de_temps(void) const
{


  // remaillage, calcul indicatrice, et autres
  //mon_equation_interfaces->preparer_pas_de_temps();

  // mise a jour de rho et mu
  //mon_equation_hydraulique->preparer_pas_de_temps();

  // mise a jour de rho_cp et T_sat
  //if (mon_equation_thermique.non_nul())
  //  mon_equation_thermique->preparer_pas_de_temps();

  // prend le min des pas de temps de chaque equation
  double dt = Pb_qdm_fluide::calculer_pas_de_temps();

  // calcul delta_h et limite le pas de temps si evaporation
  //if (mon_equation_thermique.non_nul())
  //  mon_equation_thermique->corriger_pas_de_temps(dt);
  if (la_chimie_.non_nul())

    {
      dt=min(dt,la_chimie_.valeur().calculer_pas_de_temps());
      dt=mp_min(dt);
    }
  return dt;
}


void Probleme_FT_Disc_gen::discretiser(const Discretisation_base& dis)
{
  Pb_qdm_fluide::discretiser(dis);
  if (la_chimie_.non_nul())
    la_chimie_.valeur().discretiser(*this);
}
void Probleme_FT_Disc_gen::mettre_a_jour(double temps)
{
  if (schema_temps().que_suis_je() == "RK3_FT")
    {
      int nb_eqn = nombre_d_equations();
      for (int i=2; i<nb_eqn; i++)
        equation(i).mettre_a_jour(temps);

      les_postraitements.mettre_a_jour(temps);
    }
  else
    {
      Pb_qdm_fluide::mettre_a_jour(temps);
    }
  if (la_chimie_.non_nul())
    la_chimie_.valeur().mettre_a_jour(temps);
}

void Probleme_FT_Disc_gen::preparer_mise_a_jour(void)
{
  //mon_equation_interfaces->preparer_mise_a_jour();
  // delegue aux interfaces les calculs de
  //   -> debit massique
  //   -> tension superficielle
  //   -> force interfaciale
  //   -> etalement de la courbure
}


void Probleme_FT_Disc_gen::completer(void)
{
  Pb_qdm_fluide::completer();
  if (la_chimie_.non_nul())
    la_chimie_.valeur().completer(*this);
}

int Probleme_FT_Disc_gen::verifier(void)
{
  return Pb_qdm_fluide::verifier();
}

void Probleme_FT_Disc_gen::preparer_calcul(void)
{
  Cerr<<"Probleme_FT_Disc_gen::preparer_calcul"<<finl;
  Pb_qdm_fluide::preparer_calcul();
}
