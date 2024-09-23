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
// File:        Transport_K_ou_Eps.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Incompressible/Equations/RANS
//
//////////////////////////////////////////////////////////////////////////////

#include <Transport_K_ou_Eps.h>
#include <Modele_turbulence_hyd_K_Eps_Bicephale.h>
#include <Les_Pb_Turb.h>
#include <Param.h>
#include <EChaine.h>
#include <Fluide_Quasi_Compressible.h>


Implemente_instanciable(Transport_K_ou_Eps,"Transport_K_ou_Eps",Transport_K_ou_Eps_base);

/*! @brief Imprime le type de l'equation sur un flot de sortie.
 *
 * @param (Sortie& s) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Transport_K_ou_Eps::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}

void Transport_K_ou_Eps::set_param(Param& param)
{
  Transport_K_ou_Eps_base::set_param(param);
  param.ajouter("with_nu",&with_nu_);
  param.dictionnaire("no",0);
  param.dictionnaire("yes",1);
}

/*! @brief Lit les specifications d'une equation de transport K-epsilon a partir d'un flot d'entree.
 *
 *     Controle dynamique du type du terme source.
 *
 * @param (Entree& s) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Transport_K_ou_Eps::readOn(Entree& s )
{
  with_nu_=0;
  // Lecture des attributs de l'equation
  Transport_K_ou_Eps_base::readOn(s);

  // Ajout automatique du terme source si pas instancie dans le jeu de donnees
  if (les_sources.est_vide())
    {
      if (transporte_K_)
        {
          Source t;
          Source& so=les_sources.add(t);
          const Probleme_base& pb = probleme();
          Cerr << "Construction and typing for the source term of the Transport_K equation." << finl;
          if (sub_type(Pb_Hydraulique_Turbulent,pb) || milieu().que_suis_je()=="Fluide_Quasi_Compressible")
            {
              Nom typ = "Source_Transport_K";
              Cerr << "TYPAGE DES SOURCES : this " << *this << finl;
              so.typer(typ,*this);
            }
          else if (sub_type(Pb_Thermohydraulique_Turbulent,pb))
            {
              Nom typ = "Source_Transport_K_anisotherme";
              so.typer(typ,*this);
            }
          else if (sub_type(Pb_Hydraulique_Concentration_Turbulent,pb))
            {
              Nom typ = "Source_Transport_K_aniso_concen";
              so.typer(typ,*this);
            }
          else if ( (sub_type(Pb_Thermohydraulique_Concentration_Turbulent,pb)) ) //|| (sub_type(Probleme_Combustion,pb)) )
            {
              Nom typ = "Source_Transport_K_aniso_therm_concen";
              so.typer(typ,*this);
            }
          so->associer_eqn(*this);
        }
      else
        {
          Source t;
          Source& so=les_sources.add(t);
          const Probleme_base& pb = probleme();
          Cerr << "Construction and typing for the source term of the Transport_Eps equation." << finl;
          if (sub_type(Pb_Hydraulique_Turbulent,pb) || milieu().que_suis_je()=="Fluide_Quasi_Compressible")
            {
              Nom typ = "Source_Transport_Eps";
              Cerr << "TYPAGE DES SOURCES : this " << *this << finl;
              so.typer(typ,*this);
            }
          else if (sub_type(Pb_Thermohydraulique_Turbulent,pb))
            {
              Nom typ = "Source_Transport_Eps_anisotherme";
              so.typer(typ,*this);
            }
          else if (sub_type(Pb_Hydraulique_Concentration_Turbulent,pb))
            {
              Nom typ = "Source_Transport_Eps_aniso_concen";
              so.typer(typ,*this);
            }
          else if ( (sub_type(Pb_Thermohydraulique_Concentration_Turbulent,pb)) ) //|| (sub_type(Probleme_Combustion,pb)) )
            {
              Nom typ = "Source_Transport_Eps_aniso_therm_concen";
              so.typer(typ,*this);
            }
          so->associer_eqn(*this);
        }
    }
  return s;
}


int Transport_K_ou_Eps::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="diffusion")
    {
      Cerr << "Reading and typing of the diffusion operator : " << finl;

      if (with_nu_==0)
        {
          Cerr<<" On associe le champ de diffusion nul afin de faire comme avant !!!!!! " <<finl;
          EChaine tt("Champ_Uniforme 1 0");
          tt>> Champ_don_nul_;
          terme_diffusif.associer_diffusivite(Champ_don_nul_);
        }
      else
        {
          const Fluide_base& fluide_inc = ref_cast(Fluide_base,le_fluide.valeur());
          if (sub_type(Fluide_Quasi_Compressible,fluide_inc))
            terme_diffusif.associer_diffusivite(fluide_inc.viscosite_dynamique());
          else
            terme_diffusif.associer_diffusivite(fluide_inc.viscosite_cinematique());
        }
      is >> terme_diffusif;
      return 1;
    }
  else if (mot=="convection")
    {
      Cerr << "Reading and typing of the convection operator : " << finl;
      const Champ_Inc_base& vitesse_transportante = probleme().equation(0).inconnue();
      associer_vitesse(vitesse_transportante);
      terme_convectif.associer_vitesse(vitesse_transportante);
      is >> terme_convectif;
      return 1;
    }
  else if (mot=="ecrire_fichier_xyz_valeur")
    {
      Cerr << mot << " is not understood by " << que_suis_je() << finl;
      Cerr << "Use this keyword in the Navier Stokes equation, not in KEps equation, please." << finl;
      exit();
    }
  else
    return Transport_K_ou_Eps_base::lire_motcle_non_standard(mot,is);
  return 1;
}

/*! @brief Associe un modele de turbulence K-epsilon a l'equation.
 *
 * @param (Modele_turbulence_hyd_K_Eps& modele) le modele de turbulence K-epsilon a asoocier a l'equation
 */
void Transport_K_ou_Eps::associer_modele_turbulence(const Modele_turbulence_hyd_RANS_Bicephale_base& modele)
{
  const Equation_base& eqn_hydr = modele.equation();
  associer(eqn_hydr);
  associer_milieu_base(eqn_hydr.milieu());
  associer_vitesse(eqn_hydr.inconnue());
  mon_modele = modele;
  discretiser();
}

/*! @brief Renvoie le nombre d'operateurs de l'equation.
 *
 * Ici 2.
 *
 * @return (int) le nombre d'operateurs de l'equation
 */
int Transport_K_ou_Eps::nombre_d_operateurs() const
{
  return 2;
}

/*! @brief Renvoie l'operateur specifie par son index: renvoie terme_diffusif si i = 0
 *
 *      renvoie terme_convectif si i = 1
 *      exit si i>1
 *     (version const)
 *
 * @param (int i) l'index de l'operateur a renvoyer
 * @return (Operateur&) l'operateur specifie
 * @throws l'equation n'a pas plus de 2 operateurs
 */
const Operateur& Transport_K_ou_Eps::operateur(int i) const
{
  switch(i)
    {
    case 0:
      return terme_diffusif;
    case 1:
      return terme_convectif;
    default :
      Cerr << "Error for Transport_K_ou_Eps::operateur("<<i<<") !! " << finl;
      Cerr << "Transport_K_ou_Eps has " << nombre_d_operateurs() <<" operators "<<finl;
      Cerr << "and you are trying to access the " << i <<" th one."<< finl;
      exit();
    }
  // Pour les compilos!!
  return terme_diffusif;
}

/*! @brief Renvoie l'operateur specifie par son index: renvoie terme_diffusif si i = 0
 *
 *      renvoie terme_convectif si i = 1
 *      exit si i>1
 *
 * @param (int i) l'index de l'operateur a renvoyer
 * @return (Operateur&) l'operateur specifie
 * @throws l'equation n'a pas plus de 2 operateurs
 */
Operateur& Transport_K_ou_Eps::operateur(int i)
{
  switch(i)
    {
    case 0:
      return terme_diffusif;
    case 1:
      return terme_convectif;
    default :
      Cerr << "Error for Transport_K_ou_Eps::operateur("<<i<<") !! " << finl;
      Cerr << "Transport_K_ou_Eps has " << nombre_d_operateurs() <<" operators "<<finl;
      Cerr << "and you are trying to access the " << i <<" th one."<< finl;
      exit();
    }
  // Pour les compilos!!
  return terme_diffusif;
}


/*! @brief Renvoie le nom du domaine d'application de l'equation.
 *
 * Ici "Transport_Keps".
 *
 * @return (Motcle&) le nom du domaine d'application de l'equation
 */
const Motcle& Transport_K_ou_Eps::domaine_application() const
{
//   static Motcle domaine = "Transport_K_ou_Eps";
  static Motcle domaine = "Transport_Keps";
  return domaine;
}

DoubleTab& Transport_K_ou_Eps::corriger_derivee_impl(DoubleTab& d)
{
  const Turbulence_paroi_base& loi_paroi=modele_turbulence().loi_paroi();
  loi_paroi.corriger_derivee_impl(d);
  return Transport_K_ou_Eps_base::corriger_derivee_impl(d);
}

