/****************************************************************************
* Copyright (c) 2023, CEA
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
// File:        Transport_K_Omega.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Incompressible/Equations/RANS
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_hyd_K_Omega.h>
#include <Fluide_Quasi_Compressible.h>
#include <Discretisation_base.h>
#include <Transport_K_Omega.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Les_Pb_Turb.h>
#include <TRUSTTrav.h>
#include <EChaine.h>
#include <Param.h>
#include <Debog.h>
#include <Nom.h>

Implemente_instanciable(Transport_K_Omega, "Transport_K_Omega",Transport_K_Omega_base);

/*! @brief Imprime le type de l'equation sur un flot de sortie.
 *
 * @param (Sortie& s) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Transport_K_Omega::printOn(Sortie& s ) const { return s << que_suis_je() << "\n"; }

void Transport_K_Omega::set_param(Param& param)
{
  Transport_K_Omega_base::set_param(param);
  param.ajouter("with_nu", &with_nu_);
  param.dictionnaire("no", 0);
  param.dictionnaire("yes", 1);
}

/*! @brief Lit les specifications d'une equation de transport K-Omega a partir d'un flot d'entree.
 *
 *     Controle dynamique du type du terme source.
 *
 * @param (Entree& s) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Transport_K_Omega::readOn(Entree& s)
{
  with_nu_= 0 ;
  // Lecture des attributs de l'equation
  Transport_K_Omega_base::readOn(s);

  // Ajout automatique du terme source si pas instancie dans le jeu de donnees
  if (les_sources.est_vide())
    {
      Source t;
      Source& so=les_sources.add(t);
      const Probleme_base& pb = probleme();
      Cerr << "Construction and typing for the source term of the Transport_K_Omega equation." << finl;
      REF(Champ_base) ch;
      Nom pbb = probleme().que_suis_je();
      if (sub_type(Pb_Hydraulique_Turbulent,pb) || milieu().que_suis_je()=="Fluide_Quasi_Compressible")
        {
          Nom typ = "Source_Transport_K_Omega";
          Cerr << "TYPAGE DES SOURCES : this " << *this << finl;
          so.typer(typ,*this);
        }
      so->associer_eqn(*this);
    }
  return s;
}


int Transport_K_Omega::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="diffusion")
    {
      Cerr << "Reading and typing of the diffusion operator : " << finl;

      if (with_nu_==0)
        {
          Cerr << "Paramètre with_nu_ nul. Association d'un champ uniforme nul pour la diffusion."
               << finl;
          EChaine tt("Champ_Uniforme 1 0");
          tt >> Champ_don_nul_;
          terme_diffusif.associer_diffusivite(Champ_don_nul_);
        }
      else
        {
          const Fluide_base& fluide_inc = ref_cast(Fluide_base, le_fluide.valeur());
          if (sub_type(Fluide_Quasi_Compressible, fluide_inc))
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
      const Champ_Inc& vitesse_transportante = probleme().equation(0).inconnue();
      associer_vitesse(vitesse_transportante);
      terme_convectif.associer_vitesse(vitesse_transportante);
      is >> terme_convectif;
      return 1;
    }
  else if ((mot == "ecrire_fichier_xyz_valeur") || (mot == "ecrire_fichier_xyz_valeur_bin"))
    {
      Cerr << mot << " is not understood by " << que_suis_je() << finl;
      Cerr << "Use this keyword in the Navier Stokes equation, not in turbulence equation, please." << finl;
      exit();
    }
  else
    return Transport_K_Omega_base::lire_motcle_non_standard(mot,is);
  return 1;
}

/*! @brief Associe un modele de turbulence K-Omega a l'equation.
 *
 * @param (Modele_turbulence_hyd_K_Omega& modele) le modele de turbulence K-Omega a asoocier a l'equation
 */
void Transport_K_Omega::associer_modele_turbulence(const Mod_turb_hyd_RANS_komega& modele)
{
  const Equation_base& eqn_hydr = modele.equation();
  associer(eqn_hydr);
  associer_milieu_base(eqn_hydr.milieu());
  associer_vitesse(eqn_hydr.inconnue());
  mon_modele = ref_cast(Modele_turbulence_hyd_K_Omega, modele);
  discretiser();
}

/*! @brief Renvoie le nombre d'operateurs de l'equation.
 *
 * Ici 2.
 *
 * @return (int) le nombre d'operateurs de l'equation
 */
int Transport_K_Omega::nombre_d_operateurs() const
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
const Operateur& Transport_K_Omega::operateur(int i) const
{
  switch(i)
    {
    case 0:
      return terme_diffusif;
    case 1:
      return terme_convectif;
    default :
      Cerr << "Error for Transport_K_Omega::operateur("<<i<<") !! " << finl;
      Cerr << "Transport_K_Omega has " << nombre_d_operateurs() <<" operators "<<finl;
      Cerr << "and you are trying to access the " << i <<" th one."<< finl;
      exit();
    }
  return terme_diffusif; // Needed for the compiler
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
Operateur& Transport_K_Omega::operateur(int i)
{
  switch(i)
    {
    case 0:
      return terme_diffusif;
    case 1:
      return terme_convectif;
    default :
      Cerr << "Error for Transport_K_Omega::operateur("<<i<<") !! " << finl;
      Cerr << "Transport_K_Omega has " << nombre_d_operateurs() <<" operators "<<finl;
      Cerr << "and you are trying to access the " << i <<" th one."<< finl;
      exit();
    }
  // Pour les compilos!!
  return terme_diffusif;
}


/*! @brief Renvoie le nom du domaine d'application de l'equation.
 *
 * Ici "Transport_Komega".
 *
 * @return (Motcle&) le nom du domaine d'application de l'equation
 */
const Motcle& Transport_K_Omega::domaine_application() const
{
  static Motcle domaine = "Transport_Komega";
  return domaine;
}

DoubleTab& Transport_K_Omega::corriger_derivee_impl(DoubleTab& d)
{
  Nom pbb = probleme().que_suis_je();
  if (pbb.contient("ALE")) corriger_derivee_impl_ALE(d);

  const Turbulence_paroi_base& loi_paroi = modele_turbulence().loi_paroi().valeur();
  loi_paroi.corriger_derivee_impl(d);
  return Transport_K_Omega_base::corriger_derivee_impl(d);
}
