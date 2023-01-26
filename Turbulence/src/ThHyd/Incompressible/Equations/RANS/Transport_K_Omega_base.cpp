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
// File:        Transport_K_Omega_base.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Incompressible/Equations/RANS
//
//////////////////////////////////////////////////////////////////////////////

#include <Transport_K_Omega_base.h>
#include <Schema_Temps_base.h>
#include <Champ_Inc_P0_base.h>
#include <communications.h>
#include <Probleme_base.h>
#include <Discret_Thyd.h>
#include <Zone_VF.h>
#include <Param.h>
#include <Debog.h>

Implemente_base_sans_constructeur(Transport_K_Omega_base,
                                  "Transport_K_Omega_base",
                                  Transport_2eq_base);

Transport_K_Omega_base::Transport_K_Omega_base()
{
  champs_compris_.ajoute_nom_compris("K_Omega_residu");
}

/*! @brief
 *
 * @param (Sortie& is) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Transport_K_Omega_base::printOn(Sortie& is) const
{ return is << que_suis_je() << "\n"; }

/*! @brief Simple appel a Equation_base::readOn(Entree&)
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Transport_K_Omega_base::readOn(Entree& is)
{
  Equation_base::readOn(is);
  return is;
}

void Transport_K_Omega_base::discretiser()
{
  if (!sub_type(Discret_Thyd, discretisation()))
    {
      Cerr << " Transport_K_Omega_base::discretiser " << finl;
      Cerr << "Discretization " << discretisation().que_suis_je() << " not recognized." << finl;
      exit();
    }

  Cerr << "K-Omega transport equation (" << que_suis_je() << ") discretization" << finl;
  discretiser_K_Omega(schema_temps(), zone_dis(), le_champ_K_Omega);
  champs_compris_.ajoute_champ(le_champ_K_Omega);
  if (modele_turbulence().equation().calculate_time_derivative())
    set_calculate_time_derivative(1);

  Equation_base::discretiser();
}

void Transport_K_Omega_base::discretiser_K_Omega(const Schema_Temps_base& sch,
                                                 Zone_dis& z, Champ_Inc& ch) const
{
  Cerr << "K_Omega field discretization" << finl;
  Noms noms(2);
  Noms unit(2);
  noms[0] = "K";
  noms[1] = "omega";
  unit[0] = "m2/s2";
  unit[1] = "1/s1";

  // cAlan : possibilité de mutualiser ça dans Transport_RANS_2eq
  const Discretisation_base& dis = discretisation();
  dis.discretiser_champ("temperature", z.valeur(), multi_scalaire,
                        noms, unit, 2, sch.nb_valeurs_temporelles(),
                        sch.temps_courant(), ch);
  ch.valeur().nommer("K_Omega");
}

// cAlan : Mutualiser dans 2eq. Possibilité de faire un templace pour
// avoir zone_vf.xv ou zone_vf.xp ?

//For VEF-like scheme
void Transport_K_Omega_base::get_position_faces(Nom& position, int& n)
{
  const Zone_VF& zone_vf = ref_cast(Zone_VF, zone_dis().valeur());

  position = "x=";
  position += (Nom)zone_vf.xv(n, 0);
  position += " y=";
  position += (Nom)zone_vf.xv(n, 1);
  if (dimension == 3)
    {
      position += " z=";
      position += (Nom)zone_vf.xv(n, 2);
    }
}

// For VDF-like scheme
void Transport_K_Omega_base::get_position_cells(Nom& position, int& n)
{
  const Zone_VF& zone_vf = ref_cast(Zone_VF, zone_dis().valeur());

  position = "x=";
  position += (Nom)zone_vf.xp(n, 0);
  position += " y=";
  position += (Nom)zone_vf.xp(n, 1);
  if (dimension == 3)
    {
      position += " z=";
      position += (Nom)zone_vf.xp(n, 2);
    }
}

/*! @brief Controle le champ inconnue K-Omega en forcant a zero les valeurs du champ
 *
 *     inferieurs a 1.e-10.
 *
 * @return (int) renvoie toujours 1
 */
int Transport_K_Omega_base::controler_K_Omega()
{
  DoubleTab& K_Omega = le_champ_K_Omega.valeurs();
  int size = K_Omega.dimension(0);
  if (size < 0)
    {
      if (!sub_type(Champ_Inc_P0_base, le_champ_K_Omega.valeur()))
        {
          Cerr << "Unsupported K_Omega field in Transport_K_Omega_base::controler_K_Omega()" << finl;
          Process::exit(-1);
        }
      size = le_champ_K_Omega.valeur().equation().zone_dis().zone().nb_elem();
    }

  //int size_tot=mp_sum(size);
  // On estime pour eviter un mp_sum toujours couteux:
  int size_tot = size * Process::nproc();
  ArrOfInt neg(3);
  neg = 0;
  int control = 1;
  int lquiet = modele_turbulence().get_lquiet(); // cAlan remonter ce lquiet dans modele_turbu

  // cAlan, le 20/01/2023 : on force les valeurs au min et max comme pour le K_Eps.
  const Zone_VF& zone_vf = ref_cast(Zone_VF, zone_dis().valeur());
  double OMEGA_MIN = modele_turbulence().get_OMEGA_MIN();
  double OMEGA_MAX = modele_turbulence().get_OMEGA_MAX();
  double K_MIN = modele_turbulence().get_K_MIN();
  const IntTab& face_voisins = zone_vf.face_voisins();
  const IntTab& elem_faces = zone_vf.elem_faces();
  // PL on ne fixe au seuil minimum que si negatifs
  // car la loi de paroi peut fixer a des valeurs tres petites
  // et le rapport K*K/eps est coherent

  // Changement: 13/12/07: en cas de valeurs negatives pour k OU omega
  // on fixe k ET omega a une valeur moyenne des 2 elements voisins
  // cAlan, le 20/01/2023. Traitement identique que pour k-epsilon

  Nom position;
  Debog::verifier("Transport_K_Omega_base::controler_K_Omega K_Omega before", K_Omega);

  for (int n = 0; n < size; n++)
    {
      double& k = K_Omega(n, 0);
      double& omega = K_Omega(n, 1);
      if (k < 0 || omega < 0)
        {
          neg[0] += (  k<0 ? 1 : 0);
          neg[1] += (omega<0 ? 1 : 0);

          get_position_faces(position, n);

          // On impose une valeur plus physique (moyenne des elements voisins)
          k = 0;
          omega = 0;
          int nk = 0;
          int nomega = 0;
          int nb_faces_elem = elem_faces.line_size();
          if (size == face_voisins.dimension(0))
            {
              // cAlan : faire une fonction dans Transport_RANS_2eq qui fait la même chose ?
              // K-Eps on faces (eg:VEF)
              for (int i = 0; i < 2; i++)
                {
                  int elem = face_voisins(n, i);
                  if (elem != -1)
                    for (int j = 0; j < nb_faces_elem; j++)
                      if (j != n)
                        {
                          double& k_face = K_Omega(elem_faces(elem, j), 0);
                          if (k_face > K_MIN)
                            {
                              k += k_face;
                              nk++;
                            }
                          double& o_face = K_Omega(elem_faces(elem, j), 1);
                          if (o_face > OMEGA_MIN)
                            {
                              omega += o_face;
                              nomega++;
                            }
                        }
                }
            }
          else // (size != face_voisins.dimension(0))
            {
              get_position_cells(position, n);
              nk = 0;   // k -> k_min
              nomega = 0; // omega -> omega_min
            } // fin de (size != face_voisins.dimension(0))

          if (nk != 0) k /= nk;
          else k = K_MIN;
          if (nomega != 0) omega /= nomega;
          else omega = OMEGA_MIN;

          if (schema_temps().limpr() && !lquiet)
            {
              // Warnings printed:
              Cerr << (control ? "***Warning***: " : "***Error***: ");
              Cerr << "k forced to " << k << " on node " << n << " : " << position << finl;
              Cerr << (control ? "***Warning***: " : "***Error***: ");
              Cerr << "omega forced to " << omega << " on node " << n << " : " << position << finl;
            }
        } // fin (k < 0 || omega < 0)
      else if (omega > OMEGA_MAX)
        {
          neg[2] += 1;

          if (size == face_voisins.dimension(0))
            get_position_faces(position, n);
          else
            get_position_cells(position, n);

          omega = OMEGA_MAX;
          if (schema_temps().limpr() && !lquiet)
            {
              // Warnings printed:
              Cerr << (control ? "***Warning***: " : "***Error***: ");
              Cerr << "omega forced to " << omega << " on node " << n << " : " << position << finl;
            }
        }
    }

  K_Omega.echange_espace_virtuel();
  if (schema_temps().limpr())
    {
      mp_sum_for_each_item(neg);
      if (neg[0] || neg[1])
        {
          if (Process::je_suis_maitre())
            {
              const double time = le_champ_K_Omega.temps();
              Cerr << "Values forced for k and omega because:" << finl;
              if (neg[0])
                Cerr << "Negative values found for k on " << neg[0] << "/" << size_tot << " nodes at time " << time << finl;
              if (neg[1])
                Cerr << "Negative values found for omega on " << neg[1] << "/" << size_tot << " nodes at time " << time << finl;
              // Warning if more than 0.01% of nodes are values fixed
              // cAlan : mettre une variable "experte" dans le jdd pour ajuster ce seuil ?
              double ratio_k = 100. * neg[0] / size_tot;
              double ratio_omega = 100. * neg[1] / size_tot;
              if ((ratio_k > 0.01 || ratio_omega > 0.01) && !lquiet)
                {
                  // cAlan : adapter le texte pour omega
                  Cerr << "It is possible your initial and/or boundary conditions on k and/or omega are wrong." << finl;
                  Cerr << "Check the initial and boundary values for k and eps by using:" << finl;
                  Cerr << "k~" << (dimension == 2 ? "" : "3/2*") << "(t*U)^2 (t turbulence rate, U mean velocity) ";
                  Cerr << "and eps~Cmu^0.75 k^1.5/l with l turbulent length scale and Cmu a k-eps model parameter whose value is typically given as 0.09." << finl;
                  Cerr << "See explanations here: http://www.esi-cfd.com/esi-users/turb_parameters/" << finl;
                  Cerr << "Remark : by giving the velocity field (u) and the hydraulic diameter (d), by using boundary_field_uniform_keps_from_ud and field_uniform_keps_from_ud,  " << finl;
                  Cerr << "respectively for boudnaries and initial conditions, TRUST will determine automatically values for k and eps." << finl;
                  if (probleme().is_dilatable() == 1)
                    {
                      Cerr << "Please, don't forget (sorry for this TRUST syntax weakness) that when using Quasi-Compressible module" << finl;
                      Cerr << "the unknowns for which you define initial and boundary conditions are rho*k and rho*eps." << finl;
                    }
                }
            }
          if (!control && !lquiet)
            {
              // On quitte en postraitant pour trouver les noeuds
              // qui posent probleme
              Cerr << "The problem is postprocessed in order to find the nodes where K or Omega values go below 0." << finl;
              probleme().postraiter(1);
              exit();
            };
        }
      if (neg[2])
        {
          if (Process::je_suis_maitre())
            {
              const double time = le_champ_K_Omega.temps();
              Cerr << "Values forced for omega because:" << finl;
              Cerr << "Maximum values found for omega on " << neg[2] << "/" << size_tot << " nodes at time " << time << finl;
            }
        }
    }
  Debog::verifier("Transport_K_Omega_base::controler_K_Omega K_Omega after", K_Omega);
  return 1;
}

/*! @brief on remet omega et K positifs
 *
 */
void Transport_K_Omega_base::valider_iteration()
{
  controler_K_Omega();
}
