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
// File:        Modele_turbulence_hyd_K_Eps_Realisable.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_hyd_K_Eps_Realisable.h>
#include <Probleme_base.h>
#include <Schema_Temps_base.h>
#include <Fluide_Incompressible.h>
#include <Champ_Uniforme.h>
#include <Schema_Temps.h>
#include <Debog.h>
#include <stat_counters.h>
#include <Param.h>
#include <Modele_turbulence_hyd_K_Eps_Realisable.h>
#include <Modele_Shih_Zhu_Lumley_VEF.h>
#include <Modifier_pour_fluide_dilatable.h>
#include <Modele_turbulence_scal_base.h>
#include <TRUSTTrav.h>
#include <communications.h>
#include <Champ_Inc_P0_base.h>

Implemente_instanciable(Modele_turbulence_hyd_K_Eps_Realisable, "Modele_turbulence_hyd_K_Epsilon_Realisable", Modele_turbulence_hyd_RANS_K_Eps_base);

// XD K_Eps_Realisable mod_turb_hyd_rans K_Epsilon_Realisable -1 Realizable K-Epsilon Turbulence Model.

Sortie& Modele_turbulence_hyd_K_Eps_Realisable::printOn(Sortie& s) const
{
  return s << que_suis_je() << " " << le_nom();
}

Entree& Modele_turbulence_hyd_K_Eps_Realisable::readOn(Entree& is)
{
  return Modele_turbulence_hyd_RANS_K_Eps_base::readOn(is);
}

void Modele_turbulence_hyd_K_Eps_Realisable::set_param(Param& param)
{
  Modele_turbulence_hyd_RANS_K_Eps_base::set_param(param);
  param.ajouter_non_std("Transport_K_Epsilon_Realisable", (this), Param::REQUIRED); // XD_ADD_P chaine Keyword to define the realisable (k-eps) transportation equation.
  param.ajouter_non_std("Modele_Fonc_Realisable", (this), Param::REQUIRED); // XD_ADD_P Modele_Fonc_Realisable_base This keyword is used to set the model used
  param.ajouter("PRANDTL_K", &Prandtl_K_, Param::REQUIRED); // XD_ADD_P double Keyword to change the Prk value (default 1.0).
  param.ajouter("PRANDTL_EPS", &Prandtl_Eps_, Param::REQUIRED); // XD_ADD_P double Keyword to change the Pre value (default 1.3)
}

int Modele_turbulence_hyd_K_Eps_Realisable::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot == "Transport_K_Epsilon_Realisable")
    {
      eqn_transp_K_Eps().associer_modele_turbulence(*this);
      is >> eqn_transp_K_Eps();
      Cerr << "Realizable K_Epsilon equation type " << eqn_transp_K_Eps().que_suis_je() << finl;
      return 1;
    }
  else if (mot == "Modele_Fonc_Realisable")
    {
      get_modele_fonction().associer_eqn(eqn_transp_K_Eps());
      is >> mon_modele_fonc_;
      get_modele_fonction().discretiser();
      Cerr << "Realizable K_Epsilon model type " << get_modele_fonction().que_suis_je() << finl;
      return 1;
    }
  else
    return Modele_turbulence_hyd_RANS_K_Eps_base::lire_motcle_non_standard(mot, is);
}

Champ_Fonc& Modele_turbulence_hyd_K_Eps_Realisable::calculer_viscosite_turbulente(double temps)
{
  const Champ_base& chK_Eps = eqn_transp_K_Eps().inconnue().valeur();
  Nom type = chK_Eps.que_suis_je();
  const DoubleTab& tab_K_Eps = chK_Eps.valeurs();
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Realisable::calculer_viscosite_turbulente K_Eps", tab_K_Eps);
  DoubleTab& visco_turb = la_viscosite_turbulente_.valeurs();

  int n = tab_K_Eps.dimension(0);

  const DoubleTab& Cmu = get_modele_fonction().get_Cmu(); // attention : il faut qu'il soit deja calcule!

  Debog::verifier("Modele_turbulence_hyd_K_Eps_Realisable::calculer_viscosite_turbulente Cmu", Cmu);

  // dans le cas d'un domaine nul on doit effectuer le dimensionnement
  double non_prepare = 1;
  if (visco_turb.size() == n)
    non_prepare = 0.;
  non_prepare = mp_max(non_prepare);

  if (non_prepare == 1)
    {
      Champ_Inc visco_turb_au_format_K_eps_Rea;
      visco_turb_au_format_K_eps_Rea.typer(type);
      Champ_Inc_base& ch_visco_turb_K_eps_Rea = visco_turb_au_format_K_eps_Rea.valeur();
      ch_visco_turb_K_eps_Rea.associer_domaine_dis_base(eqn_transp_K_Eps().domaine_dis().valeur());
      ch_visco_turb_K_eps_Rea.nommer("diffusivite_turbulente");
      ch_visco_turb_K_eps_Rea.fixer_nb_comp(1);
      ch_visco_turb_K_eps_Rea.fixer_nb_valeurs_nodales(n);
      ch_visco_turb_K_eps_Rea.fixer_unite("inconnue");
      ch_visco_turb_K_eps_Rea.changer_temps(0.);

      DoubleTab& visco_turb_K_eps_Rea = ch_visco_turb_K_eps_Rea.valeurs();
      if (visco_turb_K_eps_Rea.size() != n)
        {
          Cerr << "visco_turb_K_eps_Rea size is " << visco_turb_K_eps_Rea.size() << " instead of " << n << finl;
          exit();
        }

      for (int i = 0; i < n; i++)
        {
          if (tab_K_Eps(i, 1) <= EPS_MIN_)
            visco_turb_K_eps_Rea[i] = 0;
          else
            visco_turb_K_eps_Rea[i] = Cmu(i) * tab_K_Eps(i, 0) * tab_K_Eps(i, 0) / tab_K_Eps(i, 1);

        }

      la_viscosite_turbulente_->affecter(visco_turb_au_format_K_eps_Rea.valeur());

    }
  else
    {
      for (int i = 0; i < n; i++)
        {
          if (tab_K_Eps(i, 1) <= EPS_MIN_)
            visco_turb[i] = 0;
          else
            visco_turb[i] = Cmu(i) * tab_K_Eps(i, 0) * tab_K_Eps(i, 0) / tab_K_Eps(i, 1);
        }
    }
  la_viscosite_turbulente_.changer_temps(temps);
  return la_viscosite_turbulente_;
}

void Modele_turbulence_hyd_K_Eps_Realisable::imprimer_evolution_keps_realisable(int avant) const
{
  const Schema_Temps_base& sch = eqn_transp_K_Eps().schema_temps();
  const Champ_Inc& le_champ_K_Eps = K_Eps();

  if (sch.nb_pas_dt() == 0 || sch.limpr())
    {
      const DoubleTab& tabKEps = K_Eps().valeurs();
      double k_min = DMAXFLOAT;
      double eps_min = DMAXFLOAT;
      double nut_min = DMAXFLOAT;
      double k_max = 0;
      double eps_max = 0;
      double nut_max = 0;
      int loc_k_min = -1;
      int loc_eps_min = -1;
      int loc_nut_min = -1;
      int loc_k_max = -1;
      int loc_eps_max = -1;
      int loc_nut_max = -1;
      int size = tabKEps.dimension(0);
      if (size < 0)
        {
          if (sub_type(Champ_Inc_P0_base, le_champ_K_Eps.valeur()))
            size = le_champ_K_Eps.valeur().equation().domaine_dis().domaine().nb_elem();
          else
            {
              Cerr << "Unsupported K_Eps field in Modele_turbulence_hyd_K_Eps_realisable::imprimer_evolution_keps_realisable()" << finl;
              Process::exit(-1);
            }
        }
      //ConstDoubleTab_parts parts(le_champ_K_Eps.valeurs());
      for (int n = 0; n < size; n++)
        {
          const double k = tabKEps(n, 0);
          const double eps = tabKEps(n, 1);
          double nut = 0;
          if (eps > 0)
            nut = LeCmu_ * k * k / eps;
          if (k < k_min)
            {
              k_min = k;
              loc_k_min = n;
            }
          else if (k > k_max)
            {
              k_max = k;
              loc_k_max = n;
            }
          if (eps < eps_min)
            {
              eps_min = eps;
              loc_eps_min = n;
            }
          else if (eps > eps_max)
            {
              eps_max = eps;
              loc_eps_max = n;
            }
          if (nut < nut_min)
            {
              nut_min = nut;
              loc_nut_min = n;
            }
          else if (nut > nut_max)
            {
              nut_max = nut;
              loc_nut_max = n;
            }
        }
      ArrOfDouble values(3);

      values[0] = k_min;
      values[1] = eps_min;
      values[2] = nut_min;
      mp_min_for_each_item(values);
      k_min = values[0];
      eps_min = values[1];
      nut_min = values[2];

      values[0] = k_max;
      values[1] = eps_max;
      values[2] = nut_max;
      mp_max_for_each_item(values);
      k_max = values[0];
      eps_max = values[1];
      nut_max = values[2];
      if (Process::je_suis_maitre())
        {
          Cout << finl << "K_Eps evolution (" << (avant ? "before" : "after") << " law of the wall applies) at time " << le_champ_K_Eps.temps() << ":" << finl;
          Cout << "std::min(k)=" << k_min;
          if (Process::nproc() == 1)
            Cout << " located at node " << loc_k_min;
          Cout << finl;
          Cout << "std::min(eps)=" << eps_min;
          if (Process::nproc() == 1)
            Cout << " located at node " << loc_eps_min;
          Cout << finl;
          Cout << "std::min(nut)=" << nut_min;
          if (Process::nproc() == 1)
            Cout << " located at node " << loc_nut_min;
          Cout << finl;
          Cout << "std::max(k)=" << k_max;
          if (Process::nproc() == 1)
            Cout << " located at node " << loc_k_max;
          Cout << finl;
          Cout << "std::max(eps)=" << eps_max;
          if (Process::nproc() == 1)
            Cout << " located at node " << loc_eps_max;
          Cout << finl;
          Cout << "std::max(nut)=" << nut_max;
          if (Process::nproc() == 1)
            Cout << " located at node " << loc_nut_max;
          Cout << finl;
        }
    }
}

int Modele_turbulence_hyd_K_Eps_Realisable::preparer_calcul()
{
  eqn_transp_K_Eps().preparer_calcul();
  Modele_turbulence_hyd_base::preparer_calcul();
  // GF pour initialiser la loi de paroi thermique en TBLE
//   if (equation().probleme().nombre_d_equations()>1)
//     {
//       const RefObjU& modele_turbulence = equation().probleme().equation(1).get_modele(TURBULENCE);
//       if (sub_type(Modele_turbulence_scal_base,modele_turbulence.valeur()))
//         {
//           Turbulence_paroi_scal& loi_paroi_T = ref_cast_non_const(Modele_turbulence_scal_base,modele_turbulence.valeur()).loi_paroi();
//           loi_paroi_T.init_lois_paroi();
//         }
//     }
  // GF quand on demarre un calcul il est bon d'utliser la ldp
  // encore plus quand on fait une reprise !!!!!!!!
  Champ_Inc& ch_K_Eps = K_Eps();

  const Milieu_base& mil = equation().probleme().milieu();
  if (equation().probleme().is_dilatable())
    diviser_par_rho_si_dilatable(ch_K_Eps.valeurs(), mil);
  imprimer_evolution_keps_realisable(1);
  loipar_.calculer_hyd(ch_K_Eps);
  eqn_transp_K_Eps().controler_K_Eps();
  calculer_viscosite_turbulente(K_Eps().temps());
  limiter_viscosite_turbulente();
  // on remultiplie K_eps par rho
  if (equation().probleme().is_dilatable())
    {
      multiplier_par_rho_si_dilatable(ch_K_Eps.valeurs(), mil);
      correction_nut_et_cisaillement_paroi_si_qc(*this);
    }
  la_viscosite_turbulente_.valeurs().echange_espace_virtuel();
  imprimer_evolution_keps_realisable(0);
  return 1;

}

void Modele_turbulence_hyd_K_Eps_Realisable::mettre_a_jour(double temps)
{
  Champ_Inc& ch_K_Eps = K_Eps();
  Schema_Temps_base& sch = eqn_transp_K_Eps().schema_temps();
  // Voir Schema_Temps_base::faire_un_pas_de_temps_pb_base
  eqn_transp_K_Eps().domaine_Cl_dis().mettre_a_jour(temps);
  if (!eqn_transp_K_Eps().equation_non_resolue())
    sch.faire_un_pas_de_temps_eqn_base(eqn_transp_K_Eps());
  eqn_transp_K_Eps().mettre_a_jour(temps);

  statistiques().begin_count(nut_counter_);
  const Milieu_base& mil = equation().probleme().milieu();
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Realisable::mettre_a_jour la_viscosite_turbulente before", la_viscosite_turbulente_.valeurs());
  // on divise K_eps par rho en QC pour revenir a K et Eps
  if (equation().probleme().is_dilatable())
    diviser_par_rho_si_dilatable(ch_K_Eps.valeurs(), mil);
  imprimer_evolution_keps_realisable(1);
  loipar_.calculer_hyd(ch_K_Eps);
  eqn_transp_K_Eps().controler_K_Eps();
  calculer_viscosite_turbulente(ch_K_Eps.temps());
  limiter_viscosite_turbulente();
  // on remultiplie K_eps par rho
  if (equation().probleme().is_dilatable())
    {
      multiplier_par_rho_si_dilatable(ch_K_Eps.valeurs(), mil);
      correction_nut_et_cisaillement_paroi_si_qc(*this);
    }
  la_viscosite_turbulente_.valeurs().echange_espace_virtuel();
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Realisable::mettre_a_jour apres calculer_viscosite_turbulente la_viscosite_turbulente", la_viscosite_turbulente_.valeurs());
  statistiques().end_count(nut_counter_);
}

bool Modele_turbulence_hyd_K_Eps_Realisable::initTimeStep(double dt)
{
  return eqn_transport_K_Eps_Rea_.initTimeStep(dt);
}

const Equation_base& Modele_turbulence_hyd_K_Eps_Realisable::equation_k_eps(int i) const
{
  assert((i == 0));
  return eqn_transport_K_Eps_Rea_;
}

const Champ_base& Modele_turbulence_hyd_K_Eps_Realisable::get_champ(const Motcle& nom) const
{
  try
    {
      return Modele_turbulence_hyd_RANS_K_Eps_base::get_champ(nom);
    }
  catch (Champs_compris_erreur&)
    {
    }

  if (mon_modele_fonc_.non_nul())
    {
      try
        {
          return mon_modele_fonc_->get_champ(nom);
        }
      catch (Champs_compris_erreur&)
        {
        }
    }

  throw Champs_compris_erreur();
}

void Modele_turbulence_hyd_K_Eps_Realisable::get_noms_champs_postraitables(Noms& nom, Option opt) const
{
  Modele_turbulence_hyd_RANS_K_Eps_base::get_noms_champs_postraitables(nom, opt);
  if (mon_modele_fonc_.non_nul())
    mon_modele_fonc_->get_noms_champs_postraitables(nom, opt);
}

void Modele_turbulence_hyd_K_Eps_Realisable::verifie_loi_paroi()
{
  Nom lp = loipar_.valeur().que_suis_je();
  if (lp == "negligeable_VEF" || lp == "negligeable_VDF")
    if (!associe_modele_fonction().non_nul())
      {
        Cerr << "The turbulence model of type " << que_suis_je() << finl;
        Cerr << "must not be used with a wall law of type negligeable or with a modele_function." << finl;
        Cerr << "Another wall law must be selected with this kind of turbulence model." << finl;
      }
}

void Modele_turbulence_hyd_K_Eps_Realisable::completer()
{
  eqn_transp_K_Eps().completer();
  verifie_loi_paroi();
}
