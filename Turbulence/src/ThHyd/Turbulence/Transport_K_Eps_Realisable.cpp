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
// File:        Transport_K_Eps_Realisable.cpp
// Directory:   $TRUST_ROOT/src/ThHyd/Turbulence
// Version:     /main/24
//
//////////////////////////////////////////////////////////////////////////////

#include <Transport_K_Eps_Realisable.h>
#include <Modele_turbulence_hyd_K_Eps_Realisable.h>
#include <Les_Pb_Turb.h>
#include <Param.h>
#include <Fluide_base.h>

Implemente_instanciable(Transport_K_Eps_Realisable,"Transport_K_Eps_Realisable",Transport_K_Eps_base);

// XD Transport_K_Eps_Realisable eqn_base Transport_K_Eps_Realisable -1 Realizable K-Epsilon Turbulence Model Transport Equations for K and Epsilon.

// printOn et readOn

Sortie& Transport_K_Eps_Realisable::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

Entree& Transport_K_Eps_Realisable::readOn(Entree& is )
{
  // Lecture des attributs de l'equation
  Transport_K_Eps_base::readOn(is);
  return is;
}


void Transport_K_Eps_Realisable::set_param(Param& param)
{
  Transport_K_Eps_base::set_param(param);
  param.ajouter("with_nu",&with_nu_);
  param.dictionnaire("no",0);
  param.dictionnaire("yes",1);
}

int Transport_K_Eps_Realisable::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="diffusion")
    {
      Cerr << "Reading and typing of the diffusion operator : " << finl;
      terme_diffusif.associer_diffusivite(diffusivite_pour_transport());
      is >> terme_diffusif;
      return 1;
    }
  else if (mot=="convection")
    {
      Cerr << "Reading and typing of the convection operator : " << finl;
      const Champ_base& vit_transp = vitesse_pour_transport();
      associer_vitesse(vit_transp);
      terme_convectif.associer_vitesse(vit_transp);
      is >> terme_convectif;
      return 1;
    }
  else
    return Transport_K_Eps_base::lire_motcle_non_standard(mot,is);
  return 1;
}

// Retour: int
// Signification: le nombre d'operateurs de l'equation
// Contraintes: toujours egal a 2
int Transport_K_Eps_Realisable::nombre_d_operateurs() const
{
  return 2;
}

// Description:
// renvoie terme_diffusif si i=0
// renvoie terme_convectif si i=1
// exit si i>1
const Operateur& Transport_K_Eps_Realisable::operateur(int i) const
{
  switch(i)
    {
    case 0:
      return terme_diffusif;
    case 1:
      return terme_convectif;
    default :
      Cerr << "Error for "<<que_suis_je()<<"::operateur("<<i<<") !! " << finl;
      Cerr << que_suis_je()<<" has " << nombre_d_operateurs() <<" operators "<<finl;
      Cerr << "and you are trying to access the " << i <<" th one."<< finl;
      exit();
    }
  // Pour les compilos!!
  return terme_diffusif;
}

// Description:
// renvoie terme_diffusif si i=0
// renvoie terme_convectif si i=1
// exit si i>1
Operateur& Transport_K_Eps_Realisable::operateur(int i)
{
  switch(i)
    {
    case 0:
      return terme_diffusif;
    case 1:
      return terme_convectif;
    default :
      Cerr << "Error for "<<que_suis_je()<<"::operateur("<<i<<") !! " << finl;
      Cerr << que_suis_je()<<" has " << nombre_d_operateurs() <<" operators "<<finl;
      Cerr << "and you are trying to access the " << i <<" th one."<< finl;
      exit();
    }
  // Pour les compilos!!
  return terme_diffusif;
}

const Champ_Don& Transport_K_Eps_Realisable::diffusivite_pour_transport()
{
  Fluide_base& fluide_inc = ref_cast(Fluide_base,le_fluide.valeur());
  return fluide_inc.viscosite_cinematique();
}

const Champ_base& Transport_K_Eps_Realisable::vitesse_pour_transport()
{
  const Champ_base& vitesse_transportante = probleme().equation(0).inconnue();
  return vitesse_transportante;
}


void Transport_K_Eps_Realisable::completer()
{
  // Ajout automatique du terme source
  if (les_sources.est_vide())
    {
      Source t;
      Source& so=les_sources.add(t);
      const Probleme_base& pb = probleme();
      Cerr << "Construction and typing for the source term of the Transport_K_Eps_Realisable equation." << finl;
      if (sub_type(Pb_Hydraulique_Turbulent,pb) || milieu().que_suis_je()=="Fluide_Quasi_Compressible")
        {
          Nom typ = "Source_Transport_K_Eps_Realisable";
          so.typer(typ,*this);
        }
      else if (sub_type(Pb_Thermohydraulique_Turbulent,pb))
        {
          Nom typ = "Source_Transport_K_Eps_Realisable_anisotherme";
          so.typer(typ,*this);
        }
      else if (sub_type(Pb_Hydraulique_Concentration_Turbulent,pb))
        {
          Nom typ = "Source_Transport_K_Eps_Realisable_aniso_concen";
          so.typer(typ,*this);
        }
      else if (sub_type(Pb_Thermohydraulique_Concentration_Turbulent,pb))
        {
          Nom typ = "Source_Transport_K_Eps_Realisable_aniso_therm_concen";
          so.typer(typ,*this);
        }
      else
        {
          Cerr<<"The equation "<<que_suis_je()<<" cannot be associated to a problem "<<pb.que_suis_je()<<finl;
          abort();
        }
      so->associer_eqn(*this);
    }
  Equation_base::completer();
}

void Transport_K_Eps_Realisable::associer_modele_turbulence(const Mod_turb_hyd_RANS& modele)
{
  const Equation_base& eqn_hydr = modele.equation();
  associer(eqn_hydr);
  associer_milieu_base(eqn_hydr.milieu());
  associer_vitesse(eqn_hydr.inconnue());
  mon_modele = ref_cast(Modele_turbulence_hyd_K_Eps_Realisable,modele);
  discretiser();
}

void Transport_K_Eps_Realisable::associer_milieu_base(const Milieu_base& un_milieu)
{
  le_fluide = ref_cast(Fluide_base, un_milieu);
}

const Motcle& Transport_K_Eps_Realisable::domaine_application() const
{
  static Motcle domaine = "Transport_Keps_Rea";
  return domaine;
}

