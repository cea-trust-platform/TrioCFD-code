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
// File:        Op_Diff_K_Eps_VDF_base.cpp
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_hyd_K_Eps_Bicephale.h>
#include <Modele_turbulence_hyd_K_Eps_Realisable.h>
#include <Modele_turbulence_hyd_K_Eps.h>
#include <Navier_Stokes_Turbulent.h>
#include <Op_Diff_K_Eps_VDF_base.h>
#include <Fluide_Dilatable_base.h>
#include <Operateur_base.h>
#include <Champ_P0_VDF.h>

Implemente_base(Op_Diff_K_Eps_VDF_base,"Op_Diff_K_Eps_VDF_base",Op_Diff_K_Eps_base);

Sortie& Op_Diff_K_Eps_VDF_base::printOn(Sortie& s) const { return s << que_suis_je(); }
Entree& Op_Diff_K_Eps_VDF_base::readOn(Entree& s) { return s; }

void Op_Diff_K_Eps_VDF_base::completer()
{
  Operateur_base::completer();
  iter->completer_();

  diffuse_k_seul   = false;
  diffuse_eps_seul = false;

  if(sub_type(Transport_K_Eps,mon_equation.valeur()))
    {
      const Transport_K_Eps& eqn_transport = ref_cast(Transport_K_Eps,mon_equation.valeur());
      if(sub_type( Modele_turbulence_hyd_K_Eps,eqn_transport.modele_turbulence()))
        {
          const Modele_turbulence_hyd_K_Eps& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps,eqn_transport.modele_turbulence());
          Eval_Diff_K_Eps_VDF_Elem& eval_diff = static_cast<Eval_Diff_K_Eps_VDF_Elem&> (iter->evaluateur());
          eval_diff.associer_Pr_K_Eps(mod_turb.get_Prandtl_K(),mod_turb.get_Prandtl_Eps());
        }
    }
  else if(sub_type(Transport_K_Eps_Realisable,mon_equation.valeur()))
    {
      const Transport_K_Eps_Realisable& eqn_transport = ref_cast(Transport_K_Eps_Realisable,mon_equation.valeur());
      if(sub_type( Modele_turbulence_hyd_K_Eps_Realisable,eqn_transport.modele_turbulence()))
        {
          const Modele_turbulence_hyd_K_Eps_Realisable& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Realisable,eqn_transport.modele_turbulence());
          Eval_Diff_K_Eps_VDF_Elem& eval_diff = static_cast<Eval_Diff_K_Eps_VDF_Elem&> (iter->evaluateur());
          eval_diff.associer_Pr_K_Eps(mod_turb.get_Prandtl_K(),mod_turb.get_Prandtl_Eps());
        }
    }
  else if(sub_type(Transport_K_ou_Eps,mon_equation.valeur())) // Bicephale
    {
      const Transport_K_ou_Eps& eqn_transport = ref_cast(Transport_K_ou_Eps,mon_equation.valeur());

      if ( eqn_transport.transporte_t_il_K( ) ) diffuse_k_seul = true;
      else diffuse_eps_seul = true;

      if(sub_type( Modele_turbulence_hyd_K_Eps_Bicephale,eqn_transport.modele_turbulence()))
        {
          const Modele_turbulence_hyd_K_Eps_Bicephale& mod_turb = ref_cast(Modele_turbulence_hyd_K_Eps_Bicephale,eqn_transport.modele_turbulence());
          if ( diffuse_k_seul )
            {
              Eval_Diff_K_VDF_Elem& eval_diff = static_cast<Eval_Diff_K_VDF_Elem&> (iter->evaluateur());
              eval_diff.associer_Pr_K_Eps(mod_turb.get_Prandtl_K(),mod_turb.get_Prandtl_Eps());
            }
          else
            {
              Eval_Diff_Eps_VDF_Elem& eval_diff = static_cast<Eval_Diff_Eps_VDF_Elem&> (iter->evaluateur());
              eval_diff.associer_Pr_K_Eps(mod_turb.get_Prandtl_K(),mod_turb.get_Prandtl_Eps());
            }
        }
    }
}

void Op_Diff_K_Eps_VDF_base::associer_diffusivite(const Champ_base& ch_diff)
{
  Eval_Diff_K_Eps_VDF& eval_diff_turb = dynamic_cast<Eval_Diff_K_Eps_VDF&> (iter->evaluateur());
  eval_diff_turb.associer(ch_diff);
}

const Champ_base& Op_Diff_K_Eps_VDF_base::diffusivite() const
{
  const Eval_Diff_K_Eps_VDF& eval_diff_turb = dynamic_cast<const Eval_Diff_K_Eps_VDF&> (iter->evaluateur());
  return eval_diff_turb.diffusivite();
}

void Op_Diff_K_Eps_VDF_base::associer_diffusivite_turbulente()
{
  assert(mon_equation.non_nul());
  if(sub_type(Transport_K_Eps_base,mon_equation.valeur()))
    {
      const Transport_K_Eps_base& eqn_transport = ref_cast(Transport_K_Eps_base,mon_equation.valeur());
      {
        const Mod_turb_hyd_RANS& mod_turb = ref_cast(Mod_turb_hyd_RANS,eqn_transport.modele_turbulence());
        const Champ_Fonc& diff_turb = mod_turb.viscosite_turbulente();
        Eval_Diff_K_Eps_VDF& eval_diff = dynamic_cast<Eval_Diff_K_Eps_VDF&> (iter->evaluateur());
        eval_diff.associer_diff_turb(diff_turb);
      }
    }
  else if(sub_type(Transport_K_ou_Eps_base,mon_equation.valeur()))
    {
      const Transport_K_ou_Eps_base& eqn_transport = ref_cast(Transport_K_ou_Eps_base,mon_equation.valeur());
      {
        const Mod_turb_hyd_RANS_Bicephale& mod_turb = ref_cast(Mod_turb_hyd_RANS_Bicephale,eqn_transport.modele_turbulence());
        const Champ_Fonc& diff_turb = mod_turb.viscosite_turbulente();
        if ( eqn_transport.transporte_t_il_K( ) )
          {
            Eval_Diff_K_VDF_Elem& eval_diff = static_cast<Eval_Diff_K_VDF_Elem&> (iter->evaluateur());
            eval_diff.associer_diff_turb(diff_turb);
          }
        else
          {
            Eval_Diff_Eps_VDF_Elem& eval_diff = static_cast<Eval_Diff_Eps_VDF_Elem&> (iter->evaluateur());
            eval_diff.associer_diff_turb(diff_turb);
          }
      }
    }
  else
    {
      Cerr<<" erreur dans Op_Diff_K_Eps_VDF_base::associer_diffusivite ... cas non prevu "<<finl;
      Process::exit();
    }
  // association de la masse_volumique (pour le QC)
  {
    Eval_Diff_K_Eps_VDF& eval_diff = dynamic_cast<Eval_Diff_K_Eps_VDF&> (iter->evaluateur());
    const Fluide_base& mil = ref_cast(Fluide_base,mon_equation->milieu());
    const Champ_base& mvol = mil.masse_volumique();
    eval_diff.associer_mvolumique(mvol);
  }
}

const Champ_Fonc& Op_Diff_K_Eps_VDF_base::diffusivite_turbulente() const
{
  const Eval_Diff_K_Eps_VDF& eval_diff = dynamic_cast<const Eval_Diff_K_Eps_VDF&> (iter->evaluateur());
  return eval_diff.diffusivite_turbulente();
}

void Op_Diff_K_Eps_VDF_base::modifier_pour_Cl(Matrice_Morse& matrice, DoubleTab& secmem) const
{
  Op_VDF_Elem::modifier_pour_Cl(iter.zone(), iter.zone_Cl(), matrice, secmem);

  const Navier_Stokes_Turbulent& eqn_hydr = ref_cast(Navier_Stokes_Turbulent,equation().probleme().equation(0) ) ;
  const Mod_turb_hyd& mod_turb = eqn_hydr.modele_turbulence();
  const Turbulence_paroi& loipar = mod_turb.loi_paroi();

  Nom type_loi = loipar.valeur().que_suis_je();

  if ( !(type_loi.debute_par("negligeable")) )
    {
      const IntVect& tab1=matrice.get_tab1();
      DoubleVect& coeff = matrice.get_set_coeff();

      const IntTab& face_voisins = iter.zone().face_voisins();

      if(sub_type(Transport_K_Eps,mon_equation.valeur()))
        {
          const Transport_K_Eps& eqn_k_eps=ref_cast(Transport_K_Eps,equation());
          const DoubleTab& val=eqn_k_eps.inconnue().valeurs();

          const Zone_Cl_dis_base& zone_Cl_dis_base = ref_cast(Zone_Cl_dis_base,eqn_hydr.zone_Cl_dis().valeur());

          const Conds_lim& les_cl = zone_Cl_dis_base.les_conditions_limites();
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
      else if(sub_type(Transport_K_Eps_Realisable,mon_equation.valeur()))
        {
          const Transport_K_Eps_Realisable& eqn_k_eps=ref_cast(Transport_K_Eps_Realisable,equation());
          const DoubleTab& val=eqn_k_eps.inconnue().valeurs();
          const Zone_Cl_dis_base& zone_Cl_dis_base = ref_cast(Zone_Cl_dis_base,eqn_hydr.zone_Cl_dis().valeur());
          const Conds_lim& les_cl = zone_Cl_dis_base.les_conditions_limites();
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
      else if(sub_type(Transport_K_ou_Eps,mon_equation.valeur()))
        {
          const Transport_K_ou_Eps& eqn=ref_cast(Transport_K_ou_Eps,equation());
          const DoubleTab& val=eqn.inconnue().valeurs();
          const Zone_Cl_dis_base& zone_Cl_dis_base = ref_cast(Zone_Cl_dis_base,eqn_hydr.zone_Cl_dis().valeur());
          const Conds_lim& les_cl = zone_Cl_dis_base.les_conditions_limites();
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

                    // on doit remettre la ligne a l'identite et le secmem a l'inconnue
                    int idiag = tab1[elem]-1;
                    coeff[idiag]=1;
                    // pour les voisins
                    int nbvois = tab1[elem+1] - tab1[elem];
                    for (int k=1; k < nbvois; k++) coeff[idiag+k]=0;

                    secmem(elem)=val(elem);
                  }
            }
        }
      else
        {
          Cerr<<" erreur dans Op_Diff_K_Eps_VDF_base::modifier_pour_Cl ... cas non prevu "<<finl;
          Process::exit();
        }
    }
}
