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

#include <Modele_turbulence_hyd_RANS_K_Eps_base.h>
#include <Verif_Cl_Turb.h>
#include <Les_mod_turb.h>
#include <Les_Pb_Turb.h>
#include <Verif_Cl.h>

/*! @brief Teste la compatibilite des equations de convection-diffusion et de l'hydraulique.
 *
 * Les tests se fait sur les conditions
 *     aux limites discretisees de chaque equation, ainsi que sur
 *     les modeles de turbulences des equations qui doivent etre
 *     de la meme famille.
 *     Appel la fonction de librairie hors classe:
 *       tester_compatibilite_hydr_concentration(const Domaine_Cl_dis_base&,const Domaine_Cl_dis_base&)
 *
 * @return (int) renvoie toujours 1
 * @throws modeles de turbulence de famille differente pour
 * l'hydraulique et l'equation de concentration
 */
int Pb_Hydraulique_Concentration_Turbulent::verifier()
{
  const Domaine_Cl_dis_base& domaine_Cl_hydr = eq_hydraulique.domaine_Cl_dis();
  const Domaine_Cl_dis_base& domaine_Cl_co = eq_concentration.domaine_Cl_dis();

  // Verification de la compatibilite des conditions aux limites:
  tester_compatibilite_hydr_concentration(domaine_Cl_hydr,domaine_Cl_co);
  if ( sub_type( Modele_turbulence_hyd_RANS_K_Eps_base, eq_hydraulique.get_modele(TURBULENCE).valeur() ))
    {
      const Modele_turbulence_hyd_RANS_K_Eps_base& le_mod_RANS = ref_cast(Modele_turbulence_hyd_RANS_K_Eps_base, eq_hydraulique.get_modele(TURBULENCE).valeur());
      const Transport_K_Eps_base& eqn = ref_cast(Transport_K_Eps_base, le_mod_RANS.eqn_transp_K_Eps());
      const Domaine_Cl_dis_base& domaine_Cl_turb = eqn.domaine_Cl_dis();
      tester_compatibilite_hydr_turb(domaine_Cl_hydr, domaine_Cl_turb);
    }

  // Verification de la compatibilite des modeles de turbulence:
  const Modele_turbulence_hyd_base& le_mod_turb_hyd = eq_hydraulique.modele_turbulence();
  const Modele_turbulence_scal_base& le_mod_turb_co = ref_cast(Modele_turbulence_scal_base,eq_concentration.get_modele(TURBULENCE).valeur());

  if ((sub_type(Modele_turbulence_hyd_LES_base,le_mod_turb_hyd)) || (sub_type(Modele_turbulence_hyd_K_Eps,le_mod_turb_hyd)))
    {
      if (!sub_type(Modele_turbulence_scal_Schmidt,le_mod_turb_co))
        {
          Cerr << "Les modeles de turbulence ne sont pas de la meme famille" << finl;
          Cerr << "pour l'hydraulique et l'equation de concentration" << finl;
          if (sub_type(Modele_turbulence_scal_Prandtl,le_mod_turb_co))
            {
              Cerr << "Pour le modele de turbulence de l'equation de concentration, la syntaxe a changee:" << finl;
              Cerr << "Utiliser le mot cle Schmidt au lieu du mot cle Prandtl." << finl;
            }
          exit();
        }

    }


  return 1;
}

/*! @brief Teste la compatibilite des equations de la fraction massique et de l'hydraulique.
 *
 * Teste la compatibilite des modeles de turbulence
 *
 */
int Pb_Hydraulique_Melange_Binaire_Turbulent_QC::verifier()
{
  const Domaine_Cl_dis_base& domaine_Cl_hydr = eq_hydraulique.domaine_Cl_dis();
  const Domaine_Cl_dis_base& domaine_Cl_fm = eq_frac_mass.domaine_Cl_dis();
  // Verification de la compatibilite des conditions aux limites:
  tester_compatibilite_hydr_fraction_massique(domaine_Cl_hydr,domaine_Cl_fm);

  if ( sub_type(Modele_turbulence_hyd_RANS_K_Eps_base, eq_hydraulique.get_modele(TURBULENCE).valeur() ) )
    {
      const Modele_turbulence_hyd_RANS_K_Eps_base& le_mod_RANS = ref_cast(Modele_turbulence_hyd_RANS_K_Eps_base, eq_hydraulique.get_modele(TURBULENCE).valeur());
      const Transport_K_Eps_base& eqn = ref_cast(Transport_K_Eps_base, le_mod_RANS.eqn_transp_K_Eps());
      const Domaine_Cl_dis_base& domaine_Cl_turb = eqn.domaine_Cl_dis();
      tester_compatibilite_hydr_turb(domaine_Cl_hydr, domaine_Cl_turb);
    }

  // Verification de la compatibilite des modeles de turbulence:
  const Modele_turbulence_hyd_base& le_mod_turb_hyd = eq_hydraulique.modele_turbulence();
  const Modele_turbulence_scal_base& le_mod_turb_th =
    ref_cast(Modele_turbulence_scal_base,eq_frac_mass.get_modele(TURBULENCE).valeur());

  if  (sub_type(Modele_turbulence_hyd_K_Eps,le_mod_turb_hyd))
    {
      if (!sub_type(Modele_turbulence_scal_Prandtl,le_mod_turb_th))
        {
          Cerr << "Les modeles de turbulence ne sont pas de la meme famille" << finl;
          Cerr << "pour l'hydraulique et la fraction massique" << finl;
          Process::exit();
        }
    }

  return 1;
}

/*! @brief Teste la compatibilite des equations de convection-diffusion et de l'hydraulique.
 *
 * Le test se fait sur les conditions
 *     aux limites discretisees de chaque equation.
 *     Appel la fonction de librairie hors classe:
 *       tester_compatibilite_hydr_turb(const Domaine_Cl_dis_base&)
 *
 * @return (int) code de retour propage
 */
int Pb_Hydraulique_Turbulent::verifier()
{
  const Domaine_Cl_dis_base& domaine_Cl_hydr = eq_hydraulique.domaine_Cl_dis();

  if ( sub_type(Modele_turbulence_hyd_RANS_K_Eps_base, eq_hydraulique.get_modele(TURBULENCE).valeur() ))
    {
      const Modele_turbulence_hyd_RANS_K_Eps_base& le_mod_RANS = ref_cast(Modele_turbulence_hyd_RANS_K_Eps_base, eq_hydraulique.get_modele(TURBULENCE).valeur());
      const Transport_K_Eps_base& eqn = ref_cast(Transport_K_Eps_base, le_mod_RANS.eqn_transp_K_Eps());
      const Domaine_Cl_dis_base& domaine_Cl_turb = eqn.domaine_Cl_dis();
      tester_compatibilite_hydr_turb(domaine_Cl_hydr, domaine_Cl_turb);
    }

  return 1;
}

/*! @brief Teste la compatibilite des equations de convection-diffusion et de l'hydraulique.
 *
 * Le test se fait sur les conditions
 *     aux limites discretisees de chaque equation ainsi que sur
 *     les modeles de turbulence des equations qui doivent etre
 *     de la meme famille.
 *     Appels aux fonctions de librairie hors classe:
 *       tester_compatibilite_hydr_thermique(const Domaine_Cl_dis_base&,const Domaine_Cl_dis_base&)
 *       tester_compatibilite_hydr_concentration(const Domaine_Cl_dis_base&,const Domaine_Cl_dis_base&)
 *
 * @return (int) renvoie toujours 1
 * @throws les modeles de turbulence ne sont pas de la meme
 * famille pour l'hydraulique et la thermique
 * @throws les modeles de turbulence ne sont pas de la meme
 * famille pour l'hydraulique et l'equation de concentration
 */
int Pb_Thermohydraulique_Concentration_Turbulent::verifier()
{
  const Domaine_Cl_dis_base& domaine_Cl_hydr = eq_hydraulique.domaine_Cl_dis();
  const Domaine_Cl_dis_base& domaine_Cl_th = eq_thermique.domaine_Cl_dis();
  const Domaine_Cl_dis_base& domaine_Cl_co = eq_concentration.domaine_Cl_dis();

  // Verification de la compatibilite des conditions aux limites
  tester_compatibilite_hydr_thermique(domaine_Cl_hydr,domaine_Cl_th);
  tester_compatibilite_hydr_concentration(domaine_Cl_hydr,domaine_Cl_co);
  if ( sub_type(Modele_turbulence_hyd_RANS_K_Eps_base, eq_hydraulique.get_modele(TURBULENCE).valeur() ))
    {
      const Modele_turbulence_hyd_RANS_K_Eps_base& le_mod_RANS = ref_cast(Modele_turbulence_hyd_RANS_K_Eps_base, eq_hydraulique.get_modele(TURBULENCE).valeur());
      const Transport_K_Eps_base& eqn = ref_cast(Transport_K_Eps_base, le_mod_RANS.eqn_transp_K_Eps());
      const Domaine_Cl_dis_base& domaine_Cl_turb = eqn.domaine_Cl_dis();
      tester_compatibilite_hydr_turb(domaine_Cl_hydr, domaine_Cl_turb);
    }

  // Verification de la compatibilite des modeles de turbulence
  const Modele_turbulence_hyd_base& le_mod_turb_hyd = eq_hydraulique.modele_turbulence();
  const Modele_turbulence_scal_base& le_mod_turb_th = ref_cast(Modele_turbulence_scal_base,eq_thermique.get_modele(TURBULENCE).valeur());
  const Modele_turbulence_scal_base& le_mod_turb_co = ref_cast(Modele_turbulence_scal_base,eq_concentration.get_modele(TURBULENCE).valeur());
  if ((sub_type(Modele_turbulence_hyd_LES_base,le_mod_turb_hyd))  || (sub_type(Modele_turbulence_hyd_K_Eps,le_mod_turb_hyd)))
    {
      if (!sub_type(Modele_turbulence_scal_Prandtl,le_mod_turb_th))
        {
          Cerr << "Les modeles de turbulence ne sont pas de la meme famille" << finl;
          Cerr << "pour l'hydraulique et la thermique" << finl;
          exit();
        }
      if (!sub_type(Modele_turbulence_scal_Schmidt,le_mod_turb_co))
        {
          Cerr << "Les modeles de turbulence ne sont pas de la meme famille" << finl;
          Cerr << "pour l'hydraulique et l'equation de concentration" << finl;
          if (sub_type(Modele_turbulence_scal_Prandtl,le_mod_turb_co))
            {
              Cerr << "Pour le modele de turbulence de l'equation de concentration, la syntaxe a changee:" << finl;
              Cerr << "Utiliser le mot cle Schmidt au lieu du mot cle Prandtl." << finl;
            }
          exit();
        }

    }
  return 1;
}

/*! @brief Teste la compatibilite des equations de la fraction massique et de l'hydraulique.
 *
 * Teste la compatibilite des modeles de turbulence
 *
 */
int Pb_Thermohydraulique_Turbulent_QC::verifier()
{
  const Domaine_Cl_dis_base& domaine_Cl_hydr = eq_hydraulique.domaine_Cl_dis();
  const Domaine_Cl_dis_base& domaine_Cl_th = eq_thermique.domaine_Cl_dis();

  // Verification de la compatibilite des conditions aux limites:
  tester_compatibilite_hydr_thermique(domaine_Cl_hydr,domaine_Cl_th);
  if ( sub_type(Modele_turbulence_hyd_RANS_K_Eps_base, eq_hydraulique.get_modele(TURBULENCE).valeur() ) )
    {
      const Modele_turbulence_hyd_RANS_K_Eps_base& le_mod_RANS = ref_cast(Modele_turbulence_hyd_RANS_K_Eps_base, eq_hydraulique.get_modele(TURBULENCE).valeur());
      const Transport_K_Eps_base& eqn = ref_cast(Transport_K_Eps_base, le_mod_RANS.eqn_transp_K_Eps());
      const Domaine_Cl_dis_base& domaine_Cl_turb = eqn.domaine_Cl_dis();
      tester_compatibilite_hydr_turb(domaine_Cl_hydr, domaine_Cl_turb);
    }

  // Verification de la compatibilite des modeles de turbulence:
  const Modele_turbulence_hyd_base& le_mod_turb_hyd = eq_hydraulique.modele_turbulence();
  const Modele_turbulence_scal_base& le_mod_turb_th = ref_cast(Modele_turbulence_scal_base,eq_thermique.get_modele(TURBULENCE).valeur());

  if  (sub_type(Modele_turbulence_hyd_K_Eps,le_mod_turb_hyd))
    {
      if (!sub_type(Modele_turbulence_scal_Prandtl,le_mod_turb_th))
        {
          Cerr << "Les modeles de turbulence ne sont pas de la meme famille" << finl;
          Cerr << "pour l'hydraulique et la thermique" << finl;
          exit();
        }
    }
  return 1;
}

/*! @brief Teste la compatibilite des equations de la thermique et de l'hydraulique.
 *
 * Les tests se font sur les conditions
 *     aux limites discretisees de chaque equation et sur les
 *     modeles de turbulences respectifs des equations
 *     de l'hydraulique et de la thermique (qui doivent etre de la meme famille).
 *     Appel la fonction de librairie hors classe:
 *       tester_compatibilite_hydr_thermique(const Domaine_Cl_dis_base&,const Domaine_Cl_dis_base&)
 *
 * @return (int) renvoie toujours 1
 * @throws modeles de turbulence de famille differente pour
 * l'hydraulique et la thermique
 */
int Pb_Thermohydraulique_Turbulent::verifier()
{
  const Domaine_Cl_dis_base& domaine_Cl_hydr = eq_hydraulique.domaine_Cl_dis();
  const Domaine_Cl_dis_base& domaine_Cl_th = eq_thermique.domaine_Cl_dis();

  // Verification de la compatibilite des conditions aux limites:
  tester_compatibilite_hydr_thermique(domaine_Cl_hydr,domaine_Cl_th);
  if ( sub_type(Modele_turbulence_hyd_RANS_K_Eps_base, eq_hydraulique.get_modele(TURBULENCE).valeur() ))
    {
      const Modele_turbulence_hyd_RANS_K_Eps_base& le_mod_RANS = ref_cast(Modele_turbulence_hyd_RANS_K_Eps_base, eq_hydraulique.get_modele(TURBULENCE).valeur());
      const Transport_K_Eps_base& eqn = ref_cast(Transport_K_Eps_base, le_mod_RANS.eqn_transp_K_Eps());
      const Domaine_Cl_dis_base& domaine_Cl_turb = eqn.domaine_Cl_dis();
      tester_compatibilite_hydr_turb(domaine_Cl_hydr, domaine_Cl_turb);
    }
  return 1;
}
