/****************************************************************************
* Copyright (c) 2021, CEA
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
// File:        Convection_Diffusion_Phase_field.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Multiphase/Phase_field/src
// Version:     /main/22
//
//////////////////////////////////////////////////////////////////////////////

#include <Convection_Diffusion_Phase_field.h>
#include <Probleme_base.h>
#include <Discret_Thyd.h>
#include <Matrice_Morse_Sym.h>
#include <Navier_Stokes_phase_field.h>
#include <TRUSTTrav.h>
#include <Frontiere_dis_base.h>
#include <Param.h>
#include <Constituant.h>
#include <Champ_Don.h>

Implemente_instanciable_sans_constructeur(Convection_Diffusion_Phase_field,"Convection_Diffusion_Phase_field",Convection_Diffusion_Concentration);
// XD convection_diffusion_phase_field convection_diffusion_concentration convection_diffusion_phase_field -1 Cahn-Hilliard equation of the Phase Field problem. The unknown of this equation is the concentration C.
// XD attr mu_1 floattant mu_1 1 Dynamic viscosity of the first phase.
// XD attr mu_2 floattant mu_2 1 Dynamic viscosity of the second phase.
// XD attr rho_1 floattant rho_1 1 Density of the first phase.
// XD attr rho_2 floattant rho_2 1 Density of the second phase.
// XD attr potentiel_chimique_generalise chaine(into=["avec_energie_cinetique","sans_energie_cinetique"]) potentiel_chimique_generalise 0 To define (chaine set to avec_energie_cinetique) or not (chaine set to sans_energie_cinetique) if the Cahn-Hilliard equation contains the cinetic energy term.

Convection_Diffusion_Phase_field::Convection_Diffusion_Phase_field()
{
  mutype_=0;
  /*
    Noms& nom=champs_compris_.liste_noms_compris();
    nom.dimensionner(1);
    nom[0]="potentiel_chimique_generalise";
  */
}
/*! @brief Simple appel a: Convection_Diffusion_std::printOn(Sortie&)
 *
 * @param (Sortie& is) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Convection_Diffusion_Phase_field::printOn(Sortie& is) const
{
  return Convection_Diffusion_Concentration::printOn(is);
}


/*! @brief cf Convection_Diffusion_Concentration::readOn(is)
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree& is) le flot d'entree modifie
 */
Entree& Convection_Diffusion_Phase_field::readOn(Entree& is)
{
  Convection_Diffusion_Concentration::readOn(is);
  return is;
}

void Convection_Diffusion_Phase_field::set_param(Param& param)
{
  Convection_Diffusion_Concentration::set_param(param);
  //
  param.ajouter_non_std("potentiel_chimique_generalise",(this),Param::REQUIRED);
}

int Convection_Diffusion_Phase_field::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="potentiel_chimique_generalise")
    {
      Motcle temp_mutype_;
      is >> temp_mutype_;
      if(temp_mutype_=="avec_energie_cinetique")
        {
          mutype_=1;
        }
      else if(temp_mutype_=="sans_energie_cinetique")
        {
          mutype_=0;
        }
      else
        {
          Cerr<<"Allowed keywords for potentiel_chimique_generalise are :"<<finl;
          Cerr<<"avec_energie_cinetique or sans_energie_cinetique"<<finl;
          exit();
        }
      return 1;
    }
  else
    return Convection_Diffusion_Concentration::lire_motcle_non_standard(mot,is);
  return 1;
}

/*! @brief Discretise l'equation.
 *
 */
void Convection_Diffusion_Phase_field::discretiser()
{
  const Discret_Thyd& dis=ref_cast(Discret_Thyd, discretisation());

  Cerr << "Discretisation de mutilde " << finl;
  Cerr << "mutilde discretization" << finl;
  //dis.mutilde(schema_temps(), zone_dis(), ch_mutilde);

  //dis.discretiser_champ("temperature",zone_dis().valeur(),"potentiel_chimique_generalise",".",1,schema_temps().temps_courant(),ch_mutilde);
  dis.discretiser_champ("temperature",zone_dis().valeur(),"potentiel_chimique_generalise",".",constituant().nb_constituants(),schema_temps().temps_courant(),ch_mutilde);
  champs_compris_.ajoute_champ(ch_mutilde);

  const Navier_Stokes_std& eq_ns=ref_cast(Navier_Stokes_std,probleme().equation(0));
  gradient.associer_eqn(eq_ns);
  gradient.typer();
  gradient.l_op_base().associer_eqn(*this);
  Convection_Diffusion_Concentration::discretiser();

  Cerr << "Convection_Diffusion_Phase_field::discretiser() ok" << finl;
}


void Convection_Diffusion_Phase_field::completer()
{
  Equation_base::completer();
  gradient.completer();
}

Operateur_Grad& Convection_Diffusion_Phase_field::operateur_gradient()
{
  return gradient;
}
const Operateur_Grad& Convection_Diffusion_Phase_field::operateur_gradient() const
{
  return gradient;
}


int Convection_Diffusion_Phase_field::preparer_calcul()
{
  Convection_Diffusion_Concentration::preparer_calcul();

  Cerr<<"Convection_Diffusion_Phase_field::preparer_calcul"<<finl;

  // ATTENTION : avec le nouveau modele, on n'utilise plus div_alpha_rho_gradC mais
  // div_alpha_gradC. Si on veut reutiliser l'ancienne forme, il faut remodifier le code.

  // mutilde, div_alpha_gradC, alpha_gradC_carre et pression_thermo
  // ont la meme structure que la concentration

  mutilde = inconnue().valeurs();
  mutilde = 0.;
  // si on traite une variable avec "dis." (voir discretiser()), l'operation "resize" est inutile car "dis." s'en charge.
  // En sequentiel : resize() autorise
  // En parallele : resize() interdit, car alors on ne prend pas en compte les joints

  //div_alpha_rho_gradC = inconnue().valeurs();
  //div_alpha_rho_gradC = 0.;

  div_alpha_gradC = inconnue().valeurs();
  div_alpha_gradC = 0.;


  //alpha_gradC_carre = div_alpha_gradC;
  //pression_thermo = div_alpha_gradC;

  //Cerr<<" mutilde : preparer_calcul = "<< mutilde <<finl;
  //Cerr<<" div_alpha_gradC : preparer_calcul = "<< div_alpha_gradC <<finl;
  //Cerr<<" alpha_gradC_carre : preparer_calcul = "<< alpha_gradC_carre <<finl;
  //Cerr<<" pression_thermo : preparer_calcul = "<< pression_thermo <<finl;

  sources().mettre_a_jour(schema_temps().temps_courant());
  return 1;
}
