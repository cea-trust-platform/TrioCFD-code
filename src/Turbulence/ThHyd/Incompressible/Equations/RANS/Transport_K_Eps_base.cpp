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
// File:        Transport_K_Eps_base.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Incompressible/Equations/RANS
//
//////////////////////////////////////////////////////////////////////////////

#include <Transport_K_Eps_base.h>
#include <Schema_Temps_base.h>
#include <Domaine_VF.h>
#include <Champ_Inc_P0_base.h>
#include <communications.h>
#include <Probleme_base.h>
#include <Discret_Thyd.h>
#include <Param.h>
#include <Debog.h>

Implemente_base(Transport_K_Eps_base, "Transport_K_Eps_base", Transport_2eq_base);


/*! @brief
 *
 * @param (Sortie& is) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Transport_K_Eps_base::printOn(Sortie& is) const
{ return is << que_suis_je() << "\n"; }

/*! @brief Simple appel a Equation_base::readOn(Entree&)
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Transport_K_Eps_base::readOn(Entree& is)
{
  Equation_base::readOn(is);
  return is;
}

void Transport_K_Eps_base::discretiser()
{
  if (sub_type(Discret_Thyd,discretisation()))
    {
      Cerr << "K,Eps transport equation ("<< que_suis_je() <<") discretization" << finl;
      discretiser_K_Eps(schema_temps(),domaine_dis(),le_champ_K_Eps);
      champs_compris_.ajoute_champ(le_champ_K_Eps);
      if (modele_turbulence().equation().calculate_time_derivative())
        {
          set_calculate_time_derivative(1);
        }

      Equation_base::discretiser();
    }
  else
    {
      Cerr<<" Transport_K_Eps_base::discretiser "<<finl;
      Cerr<<"Discretization "<<discretisation().que_suis_je()<<" not recognized."<<finl;
      exit();
    }
}

void Transport_K_Eps_base::discretiser_K_Eps(const Schema_Temps_base& sch,
                                             Domaine_dis& z, Champ_Inc& ch) const
{
  Cerr << "K_Eps field discretization" << finl;
  Noms noms(2);
  Noms unit(2);
  noms[0]="K";
  noms[1]="eps";
  unit[0]="m2/s2";
  unit[1]="m2/s3";

  const Discretisation_base& dis = discretisation();
  dis.discretiser_champ("temperature",z.valeur(),multi_scalaire,noms,unit,2,sch.nb_valeurs_temporelles(),sch.temps_courant(),ch);
  ch->nommer("K_Eps");
}

// For VEF-like scheme
void Transport_K_Eps_base::get_position_faces(Nom& position, int& n)
{
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());

  position = "x=";
  position += (Nom)domaine_vf.xv(n, 0);
  position += " y=";
  position += (Nom)domaine_vf.xv(n, 1);
  if (dimension == 3)
    {
      position += " z=";
      position += (Nom)domaine_vf.xv(n, 2);
    }
}

// For VDF-like scheme
void Transport_K_Eps_base::get_position_cells(Nom& position, int& n)
{
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF, domaine_dis().valeur());

  position = "x=";
  position += (Nom)domaine_vf.xp(n, 0);
  position += " y=";
  position += (Nom)domaine_vf.xp(n, 1);
  if (dimension == 3)
    {
      position += " z=";
      position += (Nom)domaine_vf.xp(n, 2);
    }
}

/*! @brief Controle le champ inconnue K-epsilon en forcant a zero les valeurs du champ
 *
 *     inferieurs a 1.e-10.
 *
 * @return (int) renvoie toujours 1
 */
int Transport_K_Eps_base::controler_K_Eps()
{
  DoubleTab& K_Eps = le_champ_K_Eps.valeurs();
  int size = K_Eps.dimension(0);
  if (size < 0)
    {
      if (sub_type(Champ_Inc_P0_base, le_champ_K_Eps.valeur()))
        size = le_champ_K_Eps->equation().domaine_dis().domaine().nb_elem();
      else
        {
          Cerr << "Unsupported K_Eps field in Transport_K_Eps_base::controler_K_Eps()" << finl;
          Process::exit(-1);
        }
    }

  //int size_tot=mp_sum(size);
  // On estime pour eviter un mp_sum toujours couteux:
  int size_tot = size * Process::nproc();
  ArrOfInt neg(3);
  neg = 0;
  int control = 1;
  int lquiet = modele_turbulence().get_lquiet();
  // On interdit K-Eps negatif pour le K-Eps seulement
  // Les autres modeles (2 couches, Launder, ne sont pas assez valides)
  /* 1.6.3 : on renonce en debug a stopper le code quand kEps<0
     #ifndef NDEBUG
     if (this->que_suis_je()=="Transport_K_Eps") control=0;
     #endif
  */
  const Domaine_VF& domaine_vf = ref_cast(Domaine_VF,domaine_dis().valeur());
  double LeEPS_MIN = modele_turbulence().get_EPS_MIN();
  double LeEPS_MAX = modele_turbulence().get_EPS_MAX();
  double LeK_MIN = modele_turbulence().get_K_MIN();
  const IntTab& face_voisins = domaine_vf.face_voisins();
  const IntTab& elem_faces = domaine_vf.elem_faces();
  // PL on ne fixe au seuil minimum que si negatifs
  // car la loi de paroi peut fixer a des valeurs tres petites
  // et le rapport K*K/eps est coherent
  // Changement: 13/12/07: en cas de valeurs negatives pour k OU eps
  // on fixe k ET eps a une valeur moyenne des 2 elements voisins
  Nom position;
  Debog::verifier("Transport_K_Eps_base::controler_K_Eps K_Eps before", K_Eps);

  for (int n = 0; n < size; n++)
    {
      double& k   = K_Eps(n, 0);
      double& eps = K_Eps(n, 1);
      if (k < 0 || eps < 0)
        {
          neg[0] += (  k<0 ? 1 : 0);
          neg[1] += (eps<0 ? 1 : 0);

//          position="x=";
//          position+=(Nom)domaine_vf.xv(n,0);
//          position+=" y=";
//          position+=(Nom)domaine_vf.xv(n,1);
//          if (dimension==3)
//            {
//              position+=" z=";
//              position+=(Nom)domaine_vf.xv(n,2);
//            }

          get_position_faces(position, n);

          // On impose une valeur plus physique (moyenne des elements voisins)
          k = 0;
          eps = 0;
          int nk = 0;
          int neps = 0;
          int nb_faces_elem = elem_faces.line_size();
          if (size == face_voisins.dimension(0))
            {
              // K-Eps on faces (eg:VEF)
              for (int i = 0; i < 2; i++)
                {
                  int elem = face_voisins(n, i);
                  if (elem != -1)
                    for (int j = 0; j < nb_faces_elem; j++)
                      if (j != n)
                        {
                          double& k_face = K_Eps(elem_faces(elem, j), 0);
                          if (k_face > LeK_MIN)
                            {
                              k += k_face;
                              nk++;
                            }
                          double& e_face = K_Eps(elem_faces(elem, j), 1);
                          if (e_face > LeEPS_MIN)
                            {
                              eps += e_face;
                              neps++;
                            }
                        }
                }
            }
          else // (size != face_voisins.dimension(0))
            {
              // K-Eps on cells (eg:VDF)
//              position="x=";
//              position+=(Nom)domaine_vf.xp(n,0);
//              position+=" y=";
//              position+=(Nom)domaine_vf.xp(n,1);
//              if (dimension==3)
//                {
//                  position+=" z=";
//                  position+=(Nom)domaine_vf.xp(n,2);
//                }

              get_position_cells(position, n);

              nk = 0;   // k -> k_min
              neps = 0; // eps -> eps_min
              /* Error in algorithm ?
              for (int j=0;j<nb_faces_elem;j++)
              {
                 int face = elem_faces(n,j);
              for (int i=0; i<2; i++)
              {
               int elem = face_voisins(face,i);
               if (elem!=-1 && elem!=n)
               {
                  double& k_elem = K_Eps(elem,0);
                              if (k_elem > LeK_MIN)
                                {
                                  k += k_elem;
                                  nk++;
                                }
                              double& e_elem = K_Eps(elem,1);
                              if (e_elem > LeEPS_MIN)
                                {
                                  eps += e_elem;
                                  neps++;
                                }
               }
              }
              }*/
            } // fin de (size != face_voisins.dimension(0))

          if (nk != 0) k /= nk;
          else k = LeK_MIN;
          if (neps != 0) eps /= neps;
          else eps = LeEPS_MIN;
          if (schema_temps().limpr() && !lquiet)
            {
              // Warnings printed:
              Cerr << (control ? "***Warning***: " : "***Error***: ");
              Cerr << "k forced to " << k << " on node " << n << " : " << position << finl;
              Cerr << (control ? "***Warning***: " : "***Error***: ");
              Cerr << "eps forced to " << eps << " on node " << n << " : " << position << finl;
            }
        } // fin (k < 0 || eps < 0)
      else if (eps > LeEPS_MAX)
        {
          neg[2] += 1;

          if (size == face_voisins.dimension(0))
            {
              // K-Eps on faces (eg:VEF)

//              position="x=";
//              position+=(Nom)domaine_vf.xv(n,0);
//              position+=" y=";
//              position+=(Nom)domaine_vf.xv(n,1);
//              if (dimension==3)
//                {
//                  position+=" z=";
//                  position+=(Nom)domaine_vf.xv(n,2);
//                }

              get_position_faces(position, n);
            }
          else
            {
              // K-Eps on cells (eg:VDF)
//              position="x=";
//              position+=(Nom)domaine_vf.xp(n,0);
//              position+=" y=";
//              position+=(Nom)domaine_vf.xp(n,1);
//              if (dimension==3)
//                {
//                  position+=" z=";
//                  position+=(Nom)domaine_vf.xp(n,2);
//                }

              get_position_cells(position, n);
            }

          eps = LeEPS_MAX;
          if (schema_temps().limpr() && !lquiet)
            {
              // Warnings printed:
              Cerr << (control ? "***Warning***: " : "***Error***: ");
              Cerr << "eps forced to " << eps << " on node " << n << " : " << position << finl;
            }
        }
    }

  K_Eps.echange_espace_virtuel();
  if (schema_temps().limpr())
    {
      mp_sum_for_each_item(neg);
      if (neg[0] || neg[1])
        {
          if (Process::je_suis_maitre())
            {
              const double time = le_champ_K_Eps.temps();
              Cerr << "Values forced for k and eps because:" << finl;
              if (neg[0])
                Cerr << "Negative values found for k on " << neg[0] << "/" << size_tot << " nodes at time " << time << finl;
              if (neg[1])
                Cerr << "Negative values found for eps on " << neg[1] << "/" << size_tot << " nodes at time " << time << finl;
              // Warning if more than 0.01% of nodes are values fixed
              double ratio_k = 100. * neg[0] / size_tot;
              double ratio_eps = 100. * neg[1] / size_tot;
              if ((ratio_k > 0.01 || ratio_eps > 0.01) && !lquiet)
                {
                  Cerr << "It is possible your initial and/or boundary conditions on k and/or eps are wrong." << finl;
                  Cerr << "Check the initial and boundary values for k and eps by using:" << finl;
                  Cerr << "k~" << (dimension == 2 ? "" : "3/2*") << "(t*U)^2 (t turbulence rate, U mean velocity) ";
                  Cerr << "and eps~Cmu^0.75 k^1.5/l with l turbulent length scale and Cmu a k-eps model parameter whose value is typically given as 0.09." << finl;
                  Cerr << "Remark : by giving the velocity field (u) and the hydraulic diameter (d), by using boundary_field_uniform_keps_from_ud and field_uniform_keps_from_ud,  " << finl;
                  Cerr << "respectively for boudnaries and initial conditions, TrioCFD will determine automatically values for k and eps." << finl;
                  if (probleme().is_dilatable() == 1)
                    {
                      Cerr << "Please, don't forget (sorry for this TrioCFD syntax weakness) that when using Quasi-Compressible module" << finl;
                      Cerr << "the unknowns for which you define initial and boundary conditions are rho*k and rho*eps." << finl;
                    }
                }
            }
          if (!control && !lquiet)
            {
              // On quitte en postraitant pour trouver les noeuds
              // qui posent probleme
              Cerr << "The problem is postprocessed in order to find the nodes where K or Eps values go below 0." << finl;
              probleme().postraiter(1);
              exit();
            };
        }
      if (neg[2])
        {
          if (Process::je_suis_maitre())
            {
              const double time = le_champ_K_Eps.temps();
              Cerr << "Values forced for eps because:" << finl;
              Cerr << "Maximum values found for eps on " << neg[2] << "/" << size_tot << " nodes at time " << time << finl;
            }
        }
    }
  Debog::verifier("Transport_K_Eps_base::controler_K_Eps K_Eps after", K_Eps);
  return 1;
}

/*! @brief on remet eps et K positifs
 *
 */
void Transport_K_Eps_base::valider_iteration()
{
  controler_K_Eps();
}
