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
// File:        Transport_K_ou_Eps_base.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Incompressible/Equations/RANS
//
//////////////////////////////////////////////////////////////////////////////

#include <Transport_K_ou_Eps_base.h>
#include <Discret_Thyd.h>
#include <Probleme_base.h>
#include <Schema_Temps_base.h>
#include <Zone_VF.h>
#include <Param.h>
#include <Debog.h>
#include <communications.h>
#include <Champ_Inc_P0_base.h>

Implemente_base_sans_constructeur(Transport_K_ou_Eps_base,"Transport_K_ou_Eps_base",Equation_base);

Transport_K_ou_Eps_base::Transport_K_ou_Eps_base()
{
  champs_compris_.ajoute_nom_compris("residu");
}
/*! @brief
 *
 * @param (Sortie& is) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Transport_K_ou_Eps_base::printOn(Sortie& is) const
{
  return is << que_suis_je() << "\n";

}

/*! @brief Simple appel a Equation_base::readOn(Entree&)
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Transport_K_ou_Eps_base::readOn(Entree& is)
{
  Equation_base::readOn(is);
  return is;
}

void Transport_K_ou_Eps_base::set_param(Param& param)
{
  Equation_base::set_param(param);
  param.ajouter_non_std("diffusion",(this));
  param.ajouter_non_std("convection",(this));
  param.ajouter_condition("is_read_diffusion","The diffusion operator must be read, select negligeable type if you want to neglect it.");
  param.ajouter_condition("is_read_convection","The convection operator must be read, select negligeable type if you want to neglect it.");
}

void Transport_K_ou_Eps_base::discretiser()
{
  if (sub_type(Discret_Thyd,discretisation()))
    {
      Cerr << "K,Eps transport equation ("<< que_suis_je() <<") discretization" << finl;
      discretiser_K_Eps(schema_temps(),zone_dis(),le_champ_);
      champs_compris_.ajoute_champ(le_champ_);
      if (modele_turbulence().equation().calculate_time_derivative())
        {
          set_calculate_time_derivative(1);
        }

      Equation_base::discretiser();
    }
  else
    {
      Cerr<<" Transport_K_ou_Eps_base::discretiser "<<finl;
      Cerr<<"Discretization "<<discretisation().que_suis_je()<<" not recognized."<<finl;
      exit();
    }
  creer_champ( "residu" );

}

void Transport_K_ou_Eps_base::discretiser_K_Eps(const Schema_Temps_base& sch,
                                                Zone_dis& z, Champ_Inc& ch) const
{
  Cerr << "K_or_Eps field discretization" << finl;
  Noms noms(1);
  Noms unit(1);

  if ( transporte_K_ )
    {
      Cerr << "K field discretization" << finl;
      noms[0]="K";
    }
  else
    {
      Cerr << "Epsilon field discretization" << finl;
      noms[0]="Eps";
    }

  unit[0]="m2/s2";

  const Discretisation_base& dis = discretisation();
  dis.discretiser_champ("temperature",z.valeur(),multi_scalaire,noms,unit,1,sch.nb_valeurs_temporelles(),sch.temps_courant(),ch);

  if ( transporte_K_ )
    {
      ch.valeur().nommer("K");
    }
  else
    {
      ch.valeur().nommer("Eps");
    }
}

/*! @brief Associe un milieu physique a l'equation.
 *
 * @param (Milieu_base& un_milieu) le milieu physique a associer a l'equation
 */
void Transport_K_ou_Eps_base::associer_milieu_base(const Milieu_base& un_milieu)
{
  le_fluide =  un_milieu;
}
/*! @brief Renvoie le milieu (fluide) associe a l'equation.
 *
 * @return (Milieu_base&) le milieu (fluide) associe a l'equation
 */
Milieu_base& Transport_K_ou_Eps_base::milieu()
{
  if(!le_fluide.non_nul())
    {
      Cerr << "No fluid has been associated to the Transport K_Epsilon"
           << que_suis_je()<< " equation." << finl;
      exit();
    }
  return le_fluide.valeur();
}


/*! @brief Renvoie le milieu (fluide) associe a l'equation.
 *
 * (version const)
 *
 * @return (Milieu_base&) le milieu (fluide) associe a l'equation
 */
const Milieu_base& Transport_K_ou_Eps_base::milieu() const
{
  if(!le_fluide.non_nul())
    {
      Cerr << "No fluid has been associated to the Transport K_Epsilon"
           << que_suis_je() <<" equation." << finl;
      exit();
    }
  return le_fluide.valeur();
}

void Transport_K_ou_Eps_base::associer(const Equation_base& eqn_hydr)
{
  Equation_base::associer_pb_base(eqn_hydr.probleme());
  Equation_base::associer_sch_tps_base(eqn_hydr.schema_temps());
  Equation_base::associer_zone_dis(eqn_hydr.zone_dis());
}

/*! @brief Controle le champ inconnue K-epsilon en forcant a zero les valeurs du champ
 *
 *     inferieurs a 1.e-10.
 *
 * @return (int) renvoie toujours 1
 */
int Transport_K_ou_Eps_base::controler_variable()
{
  DoubleTab& K_ou_Eps = le_champ_.valeurs();
  int size=K_ou_Eps.dimension(0);
  if (size<0)
    {
      if (sub_type(Champ_Inc_P0_base, le_champ_.valeur()))
        size = le_champ_.valeur().equation().zone_dis().zone().nb_elem();
      else
        {
          Cerr << "Unsupported K_ou_Eps field in Transport_K_ou_Eps_base::controler_variable()" << finl;
          Process::exit(-1);
        }
    }
  //int size_tot=mp_sum(size);
  // On estime pour eviter un mp_sum toujours couteux:
  int size_tot = size * Process::nproc();
  ArrOfInt neg(2);
  neg=0;
  int control=1;
  int lquiet = modele_turbulence().get_lquiet();
  // On interdit K-Eps negatif pour le K-Eps seulement
  // Les autres modeles (2 couches, Launder, ne sont pas assez valides)
  /* 1.6.3 : on renonce en debug a stopper le code quand kEps<0
     #ifndef NDEBUG
     if (this->que_suis_je()=="Transport_K_Eps") control=0;
     #endif
  */
  const Zone_VF& zone_vf = ref_cast(Zone_VF,zone_dis().valeur());

  double Le_MIN = modele_turbulence().get_LeEPS_MIN();

  double LeEPS_MAX = modele_turbulence().get_LeEPS_MAX();

  if ( transporte_K_ )
    {
      Le_MIN = modele_turbulence().get_LeK_MIN();
    }

  const IntTab& face_voisins = zone_vf.face_voisins();
  const IntTab& elem_faces = zone_vf.elem_faces();
  // PL on ne fixe au seuil minimum que si negatifs
  // car la loi de paroi peut fixer a des valeurs tres petites
  // et le rapport K*K/eps est coherent
  // Changement: 13/12/07: en cas de valeurs negatives pour k OU eps
  // on fixe k ET eps a une valeur moyenne des 2 elements voisins
  Nom position;
  Debog::verifier("Transport_K_ou_Eps_base::controler_variable before",K_ou_Eps);
  for (int n=0; n<size; n++)
    {
      double& var   = K_ou_Eps(n);
      if (var < 0)
        {
          neg[0] += (  var<0 ? 1 : 0);

          position="x=";
          position+=(Nom)zone_vf.xv(n,0);
          position+=" y=";
          position+=(Nom)zone_vf.xv(n,1);
          if (dimension==3)
            {
              position+=" z=";
              position+=(Nom)zone_vf.xv(n,2);
            }
          // On impose une valeur plus physique (moyenne des elements voisins)
          var = 0;
          int nvar = 0;
          int nb_faces_elem = elem_faces.line_size();
          if (size==face_voisins.dimension(0))
            {
              // K or Eps on faces (eg:VEF)
              for (int i=0; i<2; i++)
                {
                  int elem = face_voisins(n,i);
                  if (elem!=-1)
                    for (int j=0; j<nb_faces_elem; j++)
                      if (j != n)
                        {
                          double& var_face = K_ou_Eps(elem_faces(elem,j));
                          if (var_face > Le_MIN)
                            {
                              var += var_face;
                              nvar++;
                            }
                        }
                }
            }
          else
            {
              // K-Eps on cells (eg:VDF)
              position="x=";
              position+=(Nom)zone_vf.xp(n,0);
              position+=" y=";
              position+=(Nom)zone_vf.xp(n,1);
              if (dimension==3)
                {
                  position+=" z=";
                  position+=(Nom)zone_vf.xp(n,2);
                }
              nvar = 0;   // var -> var_min
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
            }
          if (nvar!=0) var /= nvar;
          else var = Le_MIN;
          if (schema_temps().limpr() && !lquiet)
            {
              // Warnings printed:
              Cerr << (control ? "***Warning***: " : "***Error***: ");
              Cerr << "k or eps forced to " << var << " on node " << n << " : " << position << finl;
            }
        }
      else if ( ( var > LeEPS_MAX ) and ( ! transporte_K_ ) )
        {
          neg[1] += 1;

          if (size==face_voisins.dimension(0))
            {
              // K-Eps on faces (eg:VEF)

              position="x=";
              position+=(Nom)zone_vf.xv(n,0);
              position+=" y=";
              position+=(Nom)zone_vf.xv(n,1);
              if (dimension==3)
                {
                  position+=" z=";
                  position+=(Nom)zone_vf.xv(n,2);
                }
            }
          else
            {
              // K-Eps on cells (eg:VDF)
              position="x=";
              position+=(Nom)zone_vf.xp(n,0);
              position+=" y=";
              position+=(Nom)zone_vf.xp(n,1);
              if (dimension==3)
                {
                  position+=" z=";
                  position+=(Nom)zone_vf.xp(n,2);
                }
            }

          var = LeEPS_MAX;
          if (schema_temps().limpr() && !lquiet)
            {
              // Warnings printed:
              Cerr << (control ? "***Warning***: " : "***Error***: ");

              Cerr << "eps forced to " << var << " on node " << n << " : " << position << finl;
            }
        }
    }

  K_ou_Eps.echange_espace_virtuel();
  if (schema_temps().limpr())
    {
      mp_sum_for_each_item(neg);
      if (neg[0] )
        {
          if (Process::je_suis_maitre())
            {
              const double time = le_champ_.temps();
              Cerr << "Values forced for k or eps because:" << finl;
              if (neg[0])
                Cerr << "Negative values found for k or eps on " << neg[0] << "/" << size_tot << " nodes at time " << time
                     << finl;
              // Warning if more than 0.01% of nodes are values fixed
              double ratio_var = 100. * neg[0] / size_tot;
              if ( ( ratio_var > 0.01 ) && !lquiet)
                {
                  Cerr << "It is possible your initial and/or boundary conditions on k or eps are wrong." << finl;
                  Cerr << "Check the initial and boundary values for k and eps by using:" << finl;
                  Cerr << "k~" << (dimension == 2 ? "" : "3/2*") << "(t*U)^2 (t turbulence rate, U mean velocity) ";
                  Cerr
                      << "and eps~Cmu^0.75 k^1.5/l with l turbulent length scale and Cmu a k-eps model parameter whose value is typically given as 0.09."
                      << finl;
                  Cerr << "See explanations here: http://www.esi-cfd.com/esi-users/turb_parameters/" << finl;
                  Cerr
                      << "Remark : by giving the velocity field (u) and the hydraulic diameter (d), by using boundary_field_uniform_keps_from_ud and field_uniform_keps_from_ud,  "
                      << finl;
                  Cerr
                      << "respectively for boudnaries and initial conditions, TRUST will determine automatically values for k and eps."
                      << finl;
                  if (probleme().is_dilatable() == 1)
                    {
                      Cerr
                          << "Please, don't forget (sorry for this TRUST syntax weakness) that when using Quasi-Compressible module"
                          << finl;
                      Cerr
                          << "the unknowns for which you define initial and boundary conditions are rho*k and rho*eps."
                          << finl;
                    }
                }
            }
          if (!control && !lquiet)
            {
              // On quitte en postraitant pour trouver les noeuds
              // qui posent probleme
              Cerr << "The problem is postprocessed in order to find the nodes where K or Eps values go below 0."
                   << finl;
              probleme().postraiter(1);
              exit();
            };
        }
      if ( neg[1] )
        {
          if (Process::je_suis_maitre())
            {
              const double time = le_champ_.temps();
              Cerr << "Values forced for eps because:" << finl;
              Cerr << "Maximum values found for eps on " << neg[1] << "/" << size_tot << " nodes at time " << time
                   << finl;
            }
        }
    }
  Debog::verifier("Transport_K_ou_Eps_base::controler_variable after",K_ou_Eps);
  return 1;
}

/*! @brief on remet eps et K positifs
 *
 */
void Transport_K_ou_Eps_base::valider_iteration()
{
  controler_variable();
}
double Transport_K_ou_Eps_base::calculer_pas_de_temps() const
{
  // on prend le pas de temps de l'eq de NS.
  return probleme().equation(0).calculer_pas_de_temps();
}

const Champ_base& Transport_K_ou_Eps_base::get_champ( const Motcle& nom ) const
{
  REF( Champ_base ) ref_champ;

  double temps_init = schema_temps().temps_init();

  if (nom=="residu")
    {
      Champ_Fonc_base& ch=ref_cast_non_const(Champ_Fonc_base,residu_.valeur());
      if ( (( ch.temps()!=le_champ_->temps() ) || ( ch.temps()==temps_init))  )
        {
          ch.mettre_a_jour( le_champ_->temps() );
        }
      return champs_compris_.get_champ( nom );
    }

  try
    {
      return Equation_base::get_champ(nom);
    }
  catch (Champs_compris_erreur)
    {
    }

  throw Champs_compris_erreur();

  return ref_champ;
}
