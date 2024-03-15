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
// File:        Op_Diff_K_Omega_VDF_base.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Operateurs/Op_Dift
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_hyd_K_Omega.h>
#include <Navier_Stokes_Turbulent.h>
#include <Op_Diff_K_Omega_VDF_base.h>
#include <Fluide_Dilatable_base.h>
#include <Operateur_base.h>
#include <Probleme_base.h>
#include <Champ_P0_VDF.h>
#include <Statistiques.h>

extern Stat_Counter_Id diffusion_counter_;

Implemente_base(Op_Diff_K_Omega_VDF_base,"Op_Diff_K_Omega_VDF_base",Op_Diff_K_Omega_base);

Sortie& Op_Diff_K_Omega_VDF_base::printOn(Sortie& s) const { return s << que_suis_je(); }
Entree& Op_Diff_K_Omega_VDF_base::readOn(Entree& s) { return s; }

void Op_Diff_K_Omega_VDF_base::completer()
{
  Operateur_base::completer();
  iter->completer_();
  iter->associer_champ_convecte_ou_inc(equation().inconnue(), nullptr);
  iter->set_name_champ_inco(equation().inconnue().le_nom().getString());
  iter->set_convective_op_pb_type(false /* diff op */, 0 /* pas pb_multiphase */);

  // diffuse_k_seul   = false;
  // diffuse_eps_seul = false;

  if(sub_type(Transport_K_Omega,mon_equation.valeur()))
    {
      const Transport_K_Omega& eqn_transport = ref_cast(Transport_K_Omega,mon_equation.valeur());
      if(sub_type( Modele_turbulence_hyd_K_Omega,eqn_transport.modele_turbulence()))
        {
          const Modele_turbulence_hyd_K_Omega& mod_turb = ref_cast(Modele_turbulence_hyd_K_Omega,eqn_transport.modele_turbulence());
          Eval_Diff_K_Omega_VDF_Elem& eval_diff = static_cast<Eval_Diff_K_Omega_VDF_Elem&> (iter->evaluateur());
          eval_diff.associer_Pr_K_Omega(mod_turb.get_Prandtl_K(),mod_turb.get_Prandtl_Omega());
        }
    }
  else
    {
      Cerr << "Bad equation type, a Transport_K_Omega type is needed." << finl;
      Process::exit();
    }
}

void Op_Diff_K_Omega_VDF_base::associer_diffusivite(const Champ_base& ch_diff)
{
  Eval_Diff_K_Omega_VDF& eval_diff_turb = dynamic_cast<Eval_Diff_K_Omega_VDF&> (iter->evaluateur());
  eval_diff_turb.associer(ch_diff);
}

const Champ_base& Op_Diff_K_Omega_VDF_base::diffusivite() const
{
  const Eval_Diff_K_Omega_VDF& eval_diff_turb = dynamic_cast<const Eval_Diff_K_Omega_VDF&> (iter->evaluateur());
  return eval_diff_turb.diffusivite();
}

void Op_Diff_K_Omega_VDF_base::associer_diffusivite_turbulente()
{
  assert(mon_equation.non_nul());
  if(sub_type(Transport_K_Omega_base,mon_equation.valeur()))
    {
      const Transport_K_Omega_base& eqn_transport = ref_cast(Transport_K_Omega_base,mon_equation.valeur());
      {
        const Modele_turbulence_hyd_RANS_komega_base& mod_turb = ref_cast(Modele_turbulence_hyd_RANS_komega_base,eqn_transport.modele_turbulence());
        const Champ_Fonc& diff_turb = mod_turb.viscosite_turbulente();
        Eval_Diff_K_Omega_VDF& eval_diff = dynamic_cast<Eval_Diff_K_Omega_VDF&> (iter->evaluateur());
        eval_diff.associer_diff_turb(diff_turb);
      }
    }
  else
    {
      Cerr<<" erreur dans Op_Diff_K_Omega_VDF_base::associer_diffusivite ... cas non prevu "<<finl;
      Process::exit();
    }
  // association de la masse_volumique (pour le QC)
  {
    Eval_Diff_K_Omega_VDF& eval_diff = dynamic_cast<Eval_Diff_K_Omega_VDF&> (iter->evaluateur());
    const Fluide_base& mil = ref_cast(Fluide_base,mon_equation->milieu());
    const Champ_base& mvol = mil.masse_volumique().valeur();
    eval_diff.associer_mvolumique(mvol);
  }
}

const Champ_Fonc& Op_Diff_K_Omega_VDF_base::diffusivite_turbulente() const
{
  const Eval_Diff_K_Omega_VDF& eval_diff = dynamic_cast<const Eval_Diff_K_Omega_VDF&> (iter->evaluateur());
  return eval_diff.diffusivite_turbulente();
}

void Op_Diff_K_Omega_VDF_base::modifier_pour_Cl(Matrice_Morse& matrice, DoubleTab& secmem) const
{
  Op_VDF_Elem::modifier_pour_Cl(iter->domaine(), iter->domaine_Cl(), matrice, secmem);

  const Navier_Stokes_Turbulent& eqn_hydr = ref_cast(Navier_Stokes_Turbulent,equation().probleme().equation(0) ) ;
  const Modele_turbulence_hyd& mod_turb = eqn_hydr.modele_turbulence();
  const Turbulence_paroi& loipar = mod_turb.loi_paroi();

  Nom type_loi = loipar.valeur().que_suis_je();

  if ( !(type_loi.debute_par("negligeable")) )
    {
      const IntVect& tab1=matrice.get_tab1();
      DoubleVect& coeff = matrice.get_set_coeff();

      const IntTab& face_voisins = iter->domaine().face_voisins();

      if(sub_type(Transport_K_Omega,mon_equation.valeur()))
        {
          const Transport_K_Omega& eqn_k_omega=ref_cast(Transport_K_Omega,equation());
          const DoubleTab& val=eqn_k_omega.inconnue().valeurs();

          const Domaine_Cl_dis_base& domaine_Cl_dis_base = ref_cast(Domaine_Cl_dis_base,eqn_hydr.domaine_Cl_dis().valeur());

          const Conds_lim& les_cl = domaine_Cl_dis_base.les_conditions_limites();
          int nb_cl=les_cl.size();
          for (int num_cl=0; num_cl<nb_cl; num_cl++)
            {
              //Boucle sur les bords
              const Cond_lim& la_cl = les_cl[num_cl];
              const Front_VF& la_front_dis = ref_cast(Front_VF,la_cl.frontiere_dis());
              int nbfaces = la_front_dis.nb_faces();
              int ndeb = la_front_dis.num_premiere_face();
              int nfin = ndeb + nbfaces;

              if ((sub_type(Dirichlet_paroi_fixe,la_cl.valeur()))||(sub_type(Dirichlet_paroi_defilante,la_cl.valeur())))
                for (int num_face=ndeb; num_face<nfin; num_face++)
                  {
                    int elem= face_voisins(num_face,0);
                    if (elem==-1) elem= face_voisins(num_face,1);
                    int nb_comp=2;
                    for (int comp=0; comp<nb_comp; comp++)
                      {
                        // on doit remettre la ligne a l'identite et le secmem a l'inconnue
                        int idiag = tab1[elem*nb_comp+comp]-1;
                        coeff[idiag]=1;
                        // pour les voisins
                        int nbvois = tab1[elem*nb_comp+1+comp] - tab1[elem*nb_comp+comp];
                        for (int k=1; k < nbvois; k++) coeff[idiag+k]=0;

                        secmem(elem,comp)=val(elem,comp);
                      }
                  }
            }
        }
      else
        {
          Cerr<<" erreur dans Op_Diff_K_Omega_VDF_base::modifier_pour_Cl ... cas non prevu "<<finl;
          Process::exit();
        }
    }
}

void Op_Diff_K_Omega_VDF_base::dimensionner_blocs(matrices_t matrices, const tabs_t& semi_impl) const
{
  const std::string& nom_inco = equation().inconnue().le_nom().getString();
  Matrice_Morse *mat = matrices.count(nom_inco) ? matrices.at(nom_inco) : nullptr, mat2;
  Op_VDF_Elem::dimensionner(iter->domaine(), iter->domaine_Cl(), mat2);
  mat->nb_colonnes() ? *mat += mat2 : *mat = mat2;
}


void Op_Diff_K_Omega_VDF_base::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  statistiques().begin_count(diffusion_counter_);
  const std::string& nom_inco = equation().inconnue().le_nom().getString();
  Matrice_Morse* mat = matrices.count(nom_inco) ? matrices.at(nom_inco) : nullptr;
  const DoubleTab& inco = semi_impl.count(nom_inco) ? semi_impl.at(nom_inco) : equation().inconnue().valeur().valeurs();

  if(mat) iter->ajouter_contribution(inco, *mat);
  mettre_a_jour_diffusivite();
  iter->ajouter(inco,secmem);
  statistiques().end_count(diffusion_counter_);

}
