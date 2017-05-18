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
// File:        Schema_Euler_Implicite_Stationnaire.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Schema_Euler_Implicite_Stationnaire/src/New
// Version:     /main/57
//
//////////////////////////////////////////////////////////////////////////////

#include <Equation_base.h>
#include <Probleme_base.h>
#include <Probleme_Couple.h>
#include <Milieu_base.h>
#include <Param.h>
#include <Debog.h>
#include <LecFicDiffuse.h>
#include <communications.h>
#include <string>
#include <Schema_Euler_Implicite_Stationnaire.h>

Implemente_instanciable(Schema_Euler_Implicite_Stationnaire,"Schema_Euler_Implicite_Stationnaire|Implicit_Euler_Steady_Scheme",Schema_Euler_Implicite);
// XD schema_euler_implicite_stationnaire schema_implicite_base schema_euler_implicite_stationnaire -1
//     printOn()
/////

Sortie& Schema_Euler_Implicite_Stationnaire::printOn(Sortie& s) const
{
  return  Schema_Implicite_base::printOn(s);
}

//// readOn
//

Entree& Schema_Euler_Implicite_Stationnaire::readOn(Entree& s)
{
  steady_security_facteur_ = 0.5;
  steady_global_dt_ = 1.e+2;
  max_dt_locaux_=0.;
  nb_steps_div0_imposed_=1;
  div0_imposed_=0;
  nb_ite_max=200;
  residu_old_=0;
  // steady_=1;
  Schema_Implicite_base::readOn(s);
  if(!le_solveur.non_nul())
    {
      Cerr << "A solver must be selected." << finl;
      Cerr << "Syntax : " << finl
           << "Solveur solver_name [ solver parameters ] " << finl;
      exit();
    }
  if (diffusion_implicite())
    {
      Cerr << "diffusion_implicite option cannot be used with an implicit time scheme." << finl;
      exit();
    }

  if (le_solveur.le_nom()!="Implicit_steady")
    {
      Cerr << "You can only select the Implicit_steady solver" << finl;
      Cerr << "Syntax : " << finl
           << "Solveur Implicit_steady [ solver parameters ] " << finl;
      exit();
    }


  return s;
}



void Schema_Euler_Implicite_Stationnaire::set_param(Param& param)
{
  param.ajouter("max_iter_implicite",&nb_ite_max); // XD_ADD_P int not_set
  param.ajouter( "steady_security_facteur",&steady_security_facteur_); // XD_ADD_P double not_set
  param.ajouter( "steady_global_dt",&steady_global_dt_);  // XD_ADD_P double not_set

  Schema_Implicite_base::set_param(param);
}


int Schema_Euler_Implicite_Stationnaire::mettre_a_jour()
{

  return Schema_Temps_base::mettre_a_jour();
}

int Schema_Euler_Implicite_Stationnaire::reprendre(Entree&)
{

  return 1;
}

void Schema_Euler_Implicite_Stationnaire::ajouter_inertie(Matrice_Base& mat_morse,DoubleTab& secmem,const Equation_base& eqn) const
{
  // codage pour euler_implicite
  DoubleVect dt_locaux = pas_de_temps_locaux();

  // On ne penalise pas les matrices et les secmems meme dans les cas
  // dirichlet , symetrie
  int pen=0;
  if (dt_locaux.size()!= secmem.size()) //c'est le cas pour la vitesse en VEF car taille = dimension*nb_faces alors que taille(dt_locaux)=nb_faces
    {
      int nb_dim=eqn.inconnue().dimension;
      if(dt_locaux.size()* nb_dim == secmem.size())
        {
          DoubleVect  dt_locaux_taille_vitesse(secmem);
          int j = 0;
          for(int i=0; i<dt_locaux_taille_vitesse.size(); i++)
            {
              if(i != 0 && i%nb_dim == 0) j++;
              dt_locaux_taille_vitesse[i]= dt_locaux[j];
            }
          dt_locaux_taille_vitesse.echange_espace_virtuel();
          eqn.solv_masse().ajouter_masse_dt_local(dt_locaux_taille_vitesse,secmem,eqn.inconnue().passe(),pen);
          eqn.solv_masse().ajouter_masse_dt_local(dt_locaux_taille_vitesse,mat_morse,pen);
        }
      else
        {
          Cerr<<"The size of the 'dt_locaux' does not match in 'Schema_Temps_base::ajouter_inertie_pas_temps_locaux' !"<<finl;
          Cerr << "Please, contact TRUST support." << finl;
          exit();
        }
    }
  else
    {
      eqn.solv_masse().ajouter_masse_dt_local(dt_locaux,secmem,eqn.inconnue().passe(),pen);
      eqn.solv_masse().ajouter_masse_dt_local(dt_locaux,mat_morse,pen);
    }
}

void Schema_Euler_Implicite_Stationnaire::calculer_pas_de_temps_local_pb()
{
  // pb_base().calculer_pas_de_temps_local(dt_locaux_);

  pb_base().equation(0).calculer_pas_de_temps_locaux(dt_locaux_);
  for(int i=1; i< pb_base().nombre_d_equations(); i++)
    {
      DoubleTab dt_temp (dt_locaux_);
      pb_base().equation(i).calculer_pas_de_temps_locaux(dt_temp);
      for (int count=0; count<dt_locaux_.size(); count++)
        {
          if (dt_locaux_[count] > dt_temp[count])
            dt_locaux_[count] =  dt_temp[count];
        }
    }

  dt_locaux_*=steady_security_facteur_;
}
void Schema_Euler_Implicite_Stationnaire::initialize()
{
  Schema_Implicite_base::initialize();
  calculer_pas_de_temps_local_pb();
}

void Schema_Euler_Implicite_Stationnaire::mettre_a_jour_dt_stab()
{
  imprimer(Cout);
  if(nb_pas_dt_ <1)
    dt_stab_=pb_base().calculer_pas_de_temps();
  calculer_pas_de_temps_local_pb();
  max_dt_locaux_= dt_locaux_.mp_max_abs_vect();
}
bool Schema_Euler_Implicite_Stationnaire::corriger_dt_calcule(double& dt_calc) const
{
  // Print the RAM
  if (limpr())
    {
      Cout << finl;
      imprimer_ram_totale();
    }
  // Compute the time step dt as the minimal value between dt_max_ and stability time step (dt_stab) * security factor (facsec)
  double dt = min (dt_max_, dt_stab_ * facsec_);

  //Definition d'un pas de temps global dans la procedure du pas de temps local
  if ( max_dt_locaux_ > 1.e-8 && nb_pas_dt_ >0)
    {
      dt = steady_global_dt_*max_dt_locaux_;
    }
  if ((dt - dt_min_)/(dt+DMINFLOAT) < -1.e-6)
    {
      // Calculation stops if time step dt is less than dt_min
      Cerr << "---------------------------------------------------------" << finl;
      Cerr << "Problem with the time step " << dt << " which is less than dt_min " << dt_min_ << finl;
      Cerr << "Lower dt_min value or check why the time step decreases..." << finl;
      Cerr << "Results are saved to provide help." << finl;
      Cerr << "---------------------------------------------------------" << finl;
      Probleme_base& pb = ref_cast_non_const(Probleme_base, mon_probleme.valeur());
      pb.postraiter(1);
      pb.sauver();
      Process::exit();
    }
  dt_calc=dt;
  return 1;
}
