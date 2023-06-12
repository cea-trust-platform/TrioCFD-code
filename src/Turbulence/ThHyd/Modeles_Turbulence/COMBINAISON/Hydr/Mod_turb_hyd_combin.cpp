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
// File:        Mod_turb_hyd_combin.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/COMBINAISON/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Mod_turb_hyd_combin.h>
#include <Equation_base.h>
#include <Probleme_base.h>
#include <Domaine_VF.h>
#include <TRUSTTabs.h>
#include <Modifier_nut_pour_fluide_dilatable.h>
#include <stat_counters.h>
#include <Param.h>
#include <Champ_Inc_P0_base.h>
#include <Champ_Inc_P1_base.h>
#include <Champ_Inc_Q1_base.h>
#include <Champ_P1NC.h>
#include <Champ_Q1NC.h>
#include <Domaine_EF.h>
#include <Domaine_VEF.h>

Implemente_instanciable_sans_constructeur(Mod_turb_hyd_combin,"Modele_turbulence_hyd_combinaison",Mod_turb_hyd_base);

Mod_turb_hyd_combin::Mod_turb_hyd_combin()
{
  nb_var_=0;
  les_var.dimensionner(0);
}

/*! @brief Simple appel a Mod_turb_hyd_base::printOn(Sortie&)
 *
 * @param (Sortie& is) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Mod_turb_hyd_combin::printOn(Sortie& is) const
{
  return Mod_turb_hyd_base::printOn(is);
}


/*! @brief Simple appel a Mod_turb_hyd_base::readOn(Entree&)
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Mod_turb_hyd_combin::readOn(Entree& is)
{
  Mod_turb_hyd_base::readOn(is);
  return is;
}

void Mod_turb_hyd_combin::set_param(Param& param)
{
  Mod_turb_hyd_base::set_param(param);
  param.ajouter("nb_var",&les_var);
  param.ajouter("fonction",&la_fct,Param::REQUIRED);
}

void Mod_turb_hyd_combin::discretiser()
{
  Mod_turb_hyd_base::discretiser();
  discretiser_K(mon_equation->schema_temps(),mon_equation->domaine_dis(),energie_cinetique_turb_);
}

// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour:
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
void Mod_turb_hyd_combin::completer()
{
  nb_var_ = les_var.size();
  fxyz.dimensionner(1);
  fxyz[0].setNbVar(4+nb_var_);
  fxyz[0].setString(la_fct);
  fxyz[0].addVar("x");
  fxyz[0].addVar("y");
  fxyz[0].addVar("z");
  fxyz[0].addVar("t");
  for (int i=0; i<nb_var_; i++)
    fxyz[0].addVar(les_var[i]);

  fxyz[0].parseString();
}

int Mod_turb_hyd_combin::preparer_calcul()
{
  Mod_turb_hyd_base::preparer_calcul();
  mettre_a_jour(0);
  return 1;
}

void Mod_turb_hyd_combin::mettre_a_jour(double temps)
{
  statistiques().begin_count(nut_counter_);
  calculer_viscosite_turbulente();
  loipar->calculer_hyd(la_viscosite_turbulente,energie_cinetique_turbulente());
  limiter_viscosite_turbulente();
  if (equation().probleme().is_dilatable()) correction_nut_et_cisaillement_paroi_si_qc(*this);
  la_viscosite_turbulente->valeurs().echange_espace_virtuel();
  statistiques().end_count(nut_counter_);
}

Champ_Fonc& Mod_turb_hyd_combin::calculer_viscosite_turbulente()
{
  const Domaine_VF& domaine_VF = ref_cast(Domaine_VF,equation().domaine_dis().valeur());
  const DoubleTab& xp = domaine_VF.xp();
  DoubleTab& viscosite_valeurs =  la_viscosite_turbulente.valeurs();
  const Probleme_base& mon_pb = equation().probleme();
  double temps = equation().inconnue().temps();
  int nb_ddl = domaine_VF.nb_elem();
  double x=0;
  double y=0;
  double z=0;
  DoubleTabs sources_val(nb_var_);

  IntTab conv_to_elem(nb_var_); // conv_to_elem() = 0 <=> champ P0, = 1 <=> champ P1 ou Q1, = 2 <=> champ P1NC ou Q1NC
  int go_conv = 0; // go_conv = 1 <=> domaine EF, go_conv = 2 <=> domaine VEF
  IntTab nb_dim_so(nb_var_);
  IntTab dim_2_so(nb_var_);

  for (int so=0; so<nb_var_; so++)
    {
      REF(Champ_base) ch_ref;
      ch_ref = mon_pb.get_champ(les_var[so]);
      const DoubleTab& source_so_val = ch_ref->valeurs();
      sources_val[so] = source_so_val;

      // type de champ ? doit on convertir en champ par element ?
      if ( sub_type(Champ_Inc_P0_base,mon_pb.get_champ(les_var[so])) )
        {
          conv_to_elem(so) = 0; // on ne fait rien
        }
      else if ( sub_type(Champ_Inc_P1_base,mon_pb.get_champ(les_var[so])) || sub_type(Champ_Inc_Q1_base,mon_pb.get_champ(les_var[so])) )
        {
          conv_to_elem(so) = 1; // champ nodal a convertir en champ element
          go_conv = 1;
        }
      else if ( sub_type(Champ_P1NC,mon_pb.get_champ(les_var[so])) || sub_type(Champ_Q1NC,mon_pb.get_champ(les_var[so])) )
        {
          conv_to_elem(so) = 2; // champ face a convertir en champ element
          go_conv = 2;
        }
      else
        {
          Cerr<<"Mod_turb_hyd_combin::calculer_viscosite_turbulente : error "<< les_var[so] <<" not Champ_P0_XX"<<finl;
          exit();
        }

      if ( conv_to_elem(so) != 0 )
        Cerr<<"Mod_turb_hyd_combin::calculer_viscosite_turbulente : " << les_var[so] <<" Conversion " << " => Champ_P0_XX"<<finl;

      // champ scalaire ? si non, dim seconde dimension ?
      nb_dim_so(so) = source_so_val.nb_dim();
      if (nb_dim_so(so) < 1 || nb_dim_so(so) > 2) // nb_dim in [1; 2] only
        {
          Cerr<<"Mod_turb_hyd_combin::calculer_viscosite_turbulente : "<<les_var[so]<<" nb dimension = "<<source_so_val.nb_dim()<<" != 1 or 2 "<<finl;
          exit();
        }
      dim_2_so(so) = 0;
      if ( nb_dim_so(so) != 1)
        {
          dim_2_so(so) = source_so_val.dimension(1);
          Cerr<<"Mod_turb_hyd_combin::calculer_viscosite_turbulente : "<<les_var[so]<<" nb dimension = "<<nb_dim_so(so)<<" and second dimension = "<<" "<<dim_2_so(so)<<finl;
        }

    }

  int nb_loop_contr = 0;
  IntTab* elem_contr = nullptr;
  if ( go_conv == 1) // Specialization domaine EF
    {
      nb_loop_contr = ref_cast(Domaine_EF,equation().domaine_dis().valeur()).domaine().nb_som_elem();
      elem_contr = &ref_cast(Domaine_EF,equation().domaine_dis().valeur()).domaine().les_elems();
    }
  else if ( go_conv == 2) // Specialization domaine VEF
    {
      nb_loop_contr = ref_cast(Domaine_VEF,equation().domaine_dis().valeur()).domaine().nb_faces_elem();
      elem_contr = &ref_cast(Domaine_VEF,equation().domaine_dis().valeur()).elem_faces();
    }

  for (int i=0; i<nb_ddl; i++) // boucle sur les elements
    {
      x=xp(i,0);
      y=xp(i,1);
      z=0;
      if (dimension >2 )
        z = xp(i,2);

      fxyz[0].setVar("x",x);
      fxyz[0].setVar("y",y);
      fxyz[0].setVar("z",z);
      fxyz[0].setVar("t",temps);

      double vale2;
      for (int so=0; so<nb_var_; so++)
        {
          const Nom nom_source = les_var[so];;
          const DoubleTab& source_so_val = sources_val[so];

          if ( conv_to_elem(so) == 0) // champ elem
            {
              if ( nb_dim_so(so) == 1 || dim_2_so(so) == 1 )
                fxyz[0].setVar(nom_source,source_so_val(i));
              else
                {
                  vale2 = 0.;
                  for (int i2=0; i2<(dim_2_so(so)-1); i2++) vale2 += source_so_val(i,i2)*source_so_val(i,i2);
                  fxyz[0].setVar(nom_source,sqrt(vale2));
                }
            }
          else // not champ elem
            {
              double vale = 0.;
              for (int icontrl=0; icontrl<nb_loop_contr; icontrl++)
                {
                  int contrl = ( elem_contr!=nullptr ? (*elem_contr)(i,icontrl) : 0);
                  if ( nb_dim_so(so) == 1 || dim_2_so(so) == 1 )
                    vale += source_so_val(contrl)/nb_loop_contr;
                  else
                    {
                      vale2 = 0.;
                      for (int i2=0; i2<(dim_2_so(so)-1); i2++) vale2 += source_so_val(contrl,i2)*source_so_val(contrl,i2);
                      vale += sqrt(vale2)/nb_loop_contr;
                    }
                }
              fxyz[0].setVar(nom_source,vale);
            }
        }
      viscosite_valeurs(i) = fxyz[0].eval();
    }
  la_viscosite_turbulente.changer_temps(temps);

  return la_viscosite_turbulente;
}
