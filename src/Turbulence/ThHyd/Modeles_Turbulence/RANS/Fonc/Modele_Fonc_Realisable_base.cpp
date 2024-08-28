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
// File:        Modele_Fonc_Realisable_base.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Fonc
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_Fonc_Realisable_base.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include<Discretisation_base.h>


Implemente_base(Modele_Fonc_Realisable_base,"Modele_Fonc_Realisable_base",Objet_U);

// XD Modele_Fonc_Realisable_base class_generic Modele_Fonc_Realisable_base 1 Base class for Functions necessary to Realizable K-Epsilon Turbulence Model

// printOn et readOn

Sortie& Modele_Fonc_Realisable_base::printOn(Sortie& s ) const
{
  return s << que_suis_je() << " " << le_nom();;
}

Entree& Modele_Fonc_Realisable_base::readOn(Entree& is )
{
  return is;
}

void Modele_Fonc_Realisable_base::completer()
{
  ;
}

void Modele_Fonc_Realisable_base::typer_lire_Modele_Fonc_Realisable(OWN_PTR(Modele_Fonc_Realisable_base) &mod, const Equation_base& eqn, Entree& is)
{
  Motcle typ, nom1("Modele_");
  is >> typ;

  nom1 += typ;
  {
    nom1 += "_";
    Cerr << nom1 << finl;
    Nom discr = eqn.discretisation().que_suis_je();
    if (discr == "VEFPreP1B")
      discr = "VEF";
    nom1 += discr;
  }

  mod.typer(nom1);
  mod->associer_eqn(eqn);
  mod->associer(eqn.domaine_dis(), eqn.domaine_Cl_dis());
  if (mod->has_seconde_equation())
    {
      mod->associer_eqn_2(mod->seconde_equation());
      mod->associer(mod->seconde_equation().domaine_dis(), mod->seconde_equation().domaine_Cl_dis());
    }
  mod->associer_pb(eqn.probleme());
  is >> mod.valeur();
}

void Modele_Fonc_Realisable_base::associer_pb(const Probleme_base& pb )
{
  Cerr << "Modele_Fonc_Realisable_base::associer_pb" << finl;
  eq_hydraulique = pb.equation(0);
}

void Modele_Fonc_Realisable_base::associer_eqn( const Equation_base& eqn)
{
  mon_equation = eqn;
}

void Modele_Fonc_Realisable_base::associer_eqn_2( const Equation_base& eqn)
{
  ma_seconde_equation = eqn;
}

void Modele_Fonc_Realisable_base::discretiser()
{
  const Discretisation_base& dis=mon_equation->discretisation();
  double temps=0;

  dis.discretiser_champ("champ_elem", mon_equation->domaine_dis(),"distance_paroi","m",1,temps,BR_wall_length_);
  champs_compris_.ajoute_champ(BR_wall_length_);

  Cerr << "Discretisation du modele K Espilon Realisable terminee" << finl;
}


int Modele_Fonc_Realisable_base::preparer_calcul()
{
  return 0;
}


int Modele_Fonc_Realisable_base::sauvegarder(Sortie& ) const
{
  return 0;
}

int Modele_Fonc_Realisable_base::reprendre(Entree& )
{
  return 0;
}

void Modele_Fonc_Realisable_base::creer_champ(const Motcle& motlu)
{
}

int Modele_Fonc_Realisable_base::Calcul_is_Cmu_constant() const
{
  return 0;
}

int Modele_Fonc_Realisable_base::Calcul_is_Reynolds_stress_isotrope() const
{
  return 1;
}


const Champ_base& Modele_Fonc_Realisable_base::get_champ(const Motcle& nom) const
{
  return champs_compris_.get_champ(nom);
}

void Modele_Fonc_Realisable_base::get_noms_champs_postraitables(Noms& nom,Option opt) const
{
  if (opt==DESCRIPTION)
    Cerr<<"Modele_Fonc_Realisable_base : "<<champs_compris_.liste_noms_compris()<<finl;
  else
    nom.add(champs_compris_.liste_noms_compris());
}



