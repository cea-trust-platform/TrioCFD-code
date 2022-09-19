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

#include <Verif_Cl_Turb.h>
#include <Zone_Cl_dis.h>
#include <Periodique.h>
#include <Dirichlet_paroi_fixe.h>
#include <Dirichlet_paroi_defilante.h>
#include <Entree_fluide_vitesse_imposee.h>
#include <Entree_fluide_K_Eps_impose.h>
#include <Neumann_paroi.h>
#include <Neumann_paroi_flux_nul.h>
#include <Symetrie.h>
#include <Neumann_sortie_libre.h>
#include <Motcle.h>
#include <Frontiere_dis_base.h>
#include <Probleme_base.h>
#include <Equation_base.h>
#include <Nom.h>

/*! @brief Teste la compatibilite des conditions aux limites hydrauliques et turbulentes.
 *
 *     La liste des compatibilites est la suivante:
 *     -----------------------------------------------------------------------
 *     Hydraulique                      |       Turbulence
 *     -----------------------------------------------------------------------
 *     Dirichlet_paroi_fixe  ou
 *     Dirichlet_paroi_defilante =======> Neumann_paroi_flux_nul
 *     =================================> Dirichlet_paroi_fixe (pour bas_reynolds)
 *     -----------------------------------------------------------------------
 *     Symetrie ========================> Symetrie
 *     =================================> Neumann_paroi_flux_nul
 *     -----------------------------------------------------------------------
 *     Neumann_sortie_libre ============> Neumann_sortie_libre
 *     =================================> Entree_fluide_K_Eps_impose
 *     -----------------------------------------------------------------------
 *     Entree_fluide_vitesse_imposee ===> Entree_fluide_K_Eps_impose
 *     =================================> Neumann_sortie_libre
 *     =================================> Symetrie
 *     -----------------------------------------------------------------------
 *     Periodique ======================> Periodique
 *     -----------------------------------------------------------------------
 *
 * @param (Zone_Cl_dis& zone_Cl_hydr)
 * @param (Zone_Cl_dis& zone_Cl_turb)
 * @return (int) renvoie toujours 1
 * @throws nombres de conditions aux limites differents
 * @throws conditions aux limite hydraulique et turbulentes incompatible
 */
int tester_compatibilite_hydr_turb(const Zone_Cl_dis& zone_Cl_hydr, const Zone_Cl_dis& zone_Cl_turb)
{

  int nb_Cl = zone_Cl_hydr.nb_cond_lim();

  if (zone_Cl_turb.nb_cond_lim() != nb_Cl)
    {
      Cerr << "The two objects of Zone_Cl_dis type don't have" << finl;
      Cerr << "the same number of boundary conditions." << finl;
      Process::exit();
    }

  for (int num_Cl=0; num_Cl<nb_Cl; num_Cl++)
    {
      const Cond_lim& la_cl_turb = zone_Cl_turb.les_conditions_limites(num_Cl);
      const Cond_lim& la_cl_hydr = zone_Cl_hydr.les_conditions_limites(num_Cl);

      // on teste si une des deux CL est periodique
      // et s'il y en a une, on verifie que les CL des deux eq sont periodiques
      if ( (sub_type(Periodique,la_cl_hydr.valeur()) || sub_type(Periodique,la_cl_turb.valeur())) &&
           ( ! ( sub_type(Periodique,la_cl_hydr.valeur()) && sub_type(Periodique,la_cl_turb.valeur()) ) ) )
        {
          message_erreur_turb( la_cl_hydr, la_cl_turb, num_Cl);
        }

      // hyd symetrie et turb (symetrie ou Neumann_paroi_flux_nul("paroi"))
      if ( sub_type(Symetrie,la_cl_hydr.valeur())  &&
           ! ( sub_type(Symetrie,la_cl_turb.valeur()) ||
               sub_type(Neumann_paroi_flux_nul,la_cl_turb.valeur()) ) )
        {
          message_erreur_turb( la_cl_hydr, la_cl_turb, num_Cl);
        }

      // hyd (paroi_fixe ou defilante) et turb (paroi ou paroi_fixe)
      if ( (sub_type(Dirichlet_paroi_fixe,la_cl_hydr.valeur()) || sub_type(Dirichlet_paroi_defilante,la_cl_hydr.valeur())) &&
           ! ( sub_type(Neumann_paroi_flux_nul,la_cl_turb.valeur()) || sub_type(Dirichlet_paroi_fixe,la_cl_turb.valeur()) ) )
        {
          message_erreur_turb( la_cl_hydr, la_cl_turb, num_Cl);
        }

      // hyd (sortie_libre) et turb (sortie_libre ou k_eps_impose)
      if ( sub_type(Neumann_sortie_libre,la_cl_hydr.valeur()) &&
           ! ( sub_type(Neumann_sortie_libre,la_cl_turb.valeur()) ||
               sub_type(Entree_fluide_K_Eps_impose,la_cl_turb.valeur()) ) )
        {
          message_erreur_turb( la_cl_hydr, la_cl_turb, num_Cl);
        }

      // hyd (Entree_fluide_vitesse_imposee ou Frontiere_ouverte_vitesse_imposee)
      // et turb (Entree_fluide_K_Eps_impose ou Frontiere_ouverte_K_Eps_impose)
      Nom pbb = zone_Cl_hydr->equation().probleme().que_suis_je();
      if ( sub_type(Entree_fluide_vitesse_imposee,la_cl_hydr.valeur()) &&
           ! ( sub_type(Neumann_sortie_libre,la_cl_turb.valeur()) ||
               sub_type(Entree_fluide_K_Eps_impose,la_cl_turb.valeur()) ||
               sub_type(Symetrie,la_cl_turb.valeur()) ) &&
           +	       !pbb.contient("ALE"))
        {
          message_erreur_turb( la_cl_hydr, la_cl_turb, num_Cl);
        }

    }
  return 1;

}


/*! @brief Affiche un message d'erreur pour la fonction precedente
 *
 * @param (Zone_Cl_dis& zone_Cl_hydr)
 * @param (Zone_Cl_dis& zone_Cl_turb)
 * @param (int num_Cl) numero de la CL
 * @return (int) renvoie toujours 1
 */
int message_erreur_turb(const Cond_lim& la_cl_hydr, const Cond_lim& la_cl_turb, int& num_Cl)
{
  Cerr << "The hydraulic and turbulent boundary conditions are not consitent on border:" << finl;
  Cerr << "Boundary conditions number " << num_Cl << " \"" << la_cl_turb.frontiere_dis().le_nom() << "\" have been assigned to : " << finl;
  Cerr << la_cl_hydr.valeur().que_suis_je() << " and " << la_cl_turb.valeur().que_suis_je() << " !! " << finl;
  Process::exit();
  return 1;
}

