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
// File:        Modele_turbulence_hyd_K_Eps_Bicephale.cpp
// Directory:   $TURBULENCE_ROOT/src/ThHyd/Modeles_Turbulence/RANS/Hydr
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_turbulence_hyd_K_Eps_Bicephale.h>
#include <Probleme_base.h>
#include <Debog.h>
#include <Modifier_pour_fluide_dilatable.h>
#include <Schema_Temps_base.h>
#include <Schema_Temps.h>
#include <stat_counters.h>
#include <Modele_turbulence_scal_base.h>
#include <Param.h>
#include <communications.h>
#include <Fluide_base.h>
#include <TRUSTTrav.h>
#include <Champ_Uniforme.h>
#include <TRUSTTab_parts.h>
#include <Champ_Inc_P0_base.h>

Implemente_instanciable(Modele_turbulence_hyd_K_Eps_Bicephale,"Modele_turbulence_hyd_K_Epsilon_Bicephale",Modele_turbulence_hyd_RANS_Bicephale_base);
// XD K_Epsilon_Bicephale mod_turb_hyd_rans_bicephale K_Epsilon_Bicephale -1 Turbulence model (k-eps) en formalisation bicephale.

/*! @brief Ecrit le type de l'objet sur un flot de sortie.
 *
 * @param (Sortie& s) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Modele_turbulence_hyd_K_Eps_Bicephale::printOn(Sortie& s ) const
{
  return s << que_suis_je() << " " << le_nom();
}


/*! @brief Simple appel a Modele_turbulence_hyd_RANS_Bicephale_base::readOn(Entree&)
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Modele_turbulence_hyd_K_Eps_Bicephale::readOn(Entree& s )
{
  return Modele_turbulence_hyd_RANS_Bicephale_base::readOn(s);
}

void Modele_turbulence_hyd_K_Eps_Bicephale::set_param(Param& param)
{
  Modele_turbulence_hyd_RANS_Bicephale_base::set_param(param);
  param.ajouter_non_std("Transport_K",(this),Param::REQUIRED); // XD_ADD_P chaine Keyword to define the realisable (k) transportation equation.
  param.ajouter_non_std("Transport_Epsilon",(this),Param::REQUIRED); // XD_ADD_P chaine Keyword to define the realisable (eps) transportation equation.
  param.ajouter_non_std("Modele_Fonc_Bas_Reynolds",(this)); // XD_ADD_P Modele_Fonc_Realisable_base This keyword is used to set the model used
  param.ajouter("CMU",&LeCmu_); // XD_ADD_P double Keyword to modify the Cmu constant of k-eps model : Nut=Cmu*k*k/eps Default value is 0.09
}

int Modele_turbulence_hyd_K_Eps_Bicephale::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="Transport_K")
    {
      eqn_transp_K().transporte_K( );
      eqn_transp_K().associer_modele_turbulence(*this);
      is >> eqn_transp_K();
      return 1;
    }
  else if (mot=="Transport_Epsilon")
    {
      eqn_transp_Eps().transporte_Eps( );
      eqn_transp_Eps().associer_modele_turbulence(*this);
      is >> eqn_transp_Eps();
      return 1;
    }
  else if (mot=="Modele_Fonc_Bas_Reynolds")
    {
      Cerr << "Lecture du modele bas reynolds associe " << finl;
      mon_modele_fonc_.associer_eqn(eqn_transp_K());
      is >> mon_modele_fonc_;
      mon_modele_fonc_.associer_eqn_2(eqn_transp_Eps());
      Cerr << "mon_modele_fonc.que_suis_je() avant discretisation " << mon_modele_fonc_.que_suis_je() << finl;
      mon_modele_fonc_.valeur().discretiser();
      Cerr << "mon_modele_fonc.que_suis_je() " << mon_modele_fonc_.valeur().que_suis_je() << finl;
      mon_modele_fonc_.valeur().lire_distance_paroi();
      return 1;
    }
  else
    return Modele_turbulence_hyd_RANS_Bicephale_base::lire_motcle_non_standard(mot,is);
}

/*! @brief Calcule la viscosite turbulente au temps demande.
 *
 * @param (double temps) le temps auquel il faut calculer la viscosite
 * @return (Champ_Fonc&) la viscosite turbulente au temps demande
 * @throws erreur de taille de visco_turb_K_eps
 */
Champ_Fonc& Modele_turbulence_hyd_K_Eps_Bicephale::calculer_viscosite_turbulente(double temps)
{
  const Champ_base& chK   = eqn_transp_K().inconnue().valeur();
  const Champ_base& chEps = eqn_transp_Eps().inconnue().valeur();
  const Domaine_dis& le_dom_dis = eqn_transp_K().domaine_dis();
  const Domaine_Cl_dis& le_dom_Cl_dis = eqn_transp_K().domaine_Cl_dis();

  Nom type=chK.que_suis_je();
  const DoubleTab& tab_K = chK.valeurs();
  //Nom type_Eps=chEps.que_suis_je();
  const DoubleTab& tab_Eps = chEps.valeurs();

  DoubleTab& visco_turb =  la_viscosite_turbulente_.valeurs();

  DoubleTab Cmu(tab_K.dimension_tot(0)) ;

  int n = tab_K.dimension(0);
  if (n<0)
    {
      if (sub_type(Champ_Inc_P0_base, chK))
        n = eqn_transp_K().domaine_dis().domaine().nb_elem();
      else
        {
          Cerr << "Unsupported K field in Modele_turbulence_hyd_K_Eps_Bicephale::calculer_viscosite_turbulente" << finl;
          Process::exit(-1);
        }
      if (sub_type(Champ_Inc_P0_base, chEps))
        n = eqn_transp_Eps().domaine_dis().domaine().nb_elem();
      else
        {
          Cerr << "Unsupported epsilon field in Modele_turbulence_hyd_K_Eps_Bicephale::calculer_viscosite_turbulente" << finl;
          Process::exit(-1);
        }
    }

  DoubleTrav Fmu,D(tab_K.dimension_tot(0));
  D=0;
  int is_modele_fonc=(mon_modele_fonc_.non_nul());
  // is_modele_fonc=0;
  if (is_modele_fonc)
    {
      // pour avoir nu en incompressible et mu en QC
      // et non comme on a divise K et eps par rho (si on est en QC)
      // on veut toujours nu
      const Champ_Don ch_visco=ref_cast(Fluide_base,eqn_transp_K().milieu()).viscosite_cinematique();
      const Champ_Don& ch_visco_cin =ref_cast(Fluide_base,eqn_transp_K().milieu()).viscosite_cinematique();

      const DoubleTab& tab_visco = ch_visco_cin->valeurs();

      Fmu.resize(tab_K.dimension_tot(0));

      mon_modele_fonc_.Calcul_Fmu_BiK( Fmu,le_dom_dis,le_dom_Cl_dis,tab_K, tab_Eps,ch_visco);

      int is_Cmu_constant = mon_modele_fonc_.Calcul_is_Cmu_constant();
      if (is_Cmu_constant==0)
        {
          Cerr<< " On utilise un Cmu non constant "<< finl;
          const DoubleTab& vitesse = mon_equation_->inconnue().valeurs();
          mon_modele_fonc_.Calcul_Cmu_BiK(Cmu, le_dom_dis, le_dom_Cl_dis,
                                          vitesse, tab_K, tab_Eps, EPS_MIN_);

          /*Paroi*/
          Nom lp=eqn_transp_K().modele_turbulence().loi_paroi().valeur().que_suis_je();
          if (lp!="negligeable_VEF")
            {
              DoubleTab visco_tab(visco_turb.dimension_tot(0));
              assert(sub_type(Champ_Uniforme,ch_visco_cin.valeur()));
              visco_tab = tab_visco(0,0);
              const int idt =  mon_equation_->schema_temps().nb_pas_dt();
              const DoubleTab& tab_paroi = loi_paroi().valeur().Cisaillement_paroi();
              mon_modele_fonc_.Calcul_Cmu_Paroi_BiK(Cmu, le_dom_dis, le_dom_Cl_dis,visco_tab, visco_turb, tab_paroi, idt,
                                                    vitesse, tab_K, tab_Eps, EPS_MIN_);
            }
        }
      else
        Cerr<< " On utilise un Cmu constant "<< finl;

    }


  // dans le cas d'un domaine nul on doit effectuer le dimensionnement
  double non_prepare=1;
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Bicephale::calculer_viscosite_turbulente la_viscosite_turbulente before",la_viscosite_turbulente_.valeurs());
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Bicephale::calculer_viscosite_turbulente tab_K",tab_K);
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Bicephale::calculer_viscosite_turbulente tab_Eps",tab_Eps);
  if (visco_turb.size() == n)
    non_prepare=0.;
  non_prepare=mp_max(non_prepare);

  if (non_prepare==1)
    {
      Champ_Inc visco_turb_au_format_K_eps;

      visco_turb_au_format_K_eps.typer(type);
      Champ_Inc_base& ch_visco_turb_K_eps=visco_turb_au_format_K_eps.valeur();
      ch_visco_turb_K_eps.associer_domaine_dis_base(eqn_transp_K().domaine_dis().valeur());
      ch_visco_turb_K_eps.nommer("diffusivite_turbulente");
      ch_visco_turb_K_eps.fixer_nb_comp(1);
      ch_visco_turb_K_eps.fixer_nb_valeurs_nodales(n);
      ch_visco_turb_K_eps.fixer_unite("inconnue");
      ch_visco_turb_K_eps.changer_temps(0.);
      DoubleTab& visco_turb_K_eps =  ch_visco_turb_K_eps.valeurs();

      if(visco_turb_K_eps.size() != n)
        {
          Cerr << "visco_turb_K_eps size is " << visco_turb_K_eps.size()
               << " instead of " << n << finl;
          exit();
        }
      // A la fin de cette boucle, le tableau visco_turb_K_eps
      // contient les valeurs de la viscosite turbulente
      // au centre des faces du maillage.
      // Debog::verifier("Modele_turbulence_hyd_K_Eps_Bicephale::calculer_viscosite_turbulente visco_turb_K_eps before",visco_turb_K_eps);
      for (int i=0; i<n; i++)
        {
          if (tab_Eps(i) <= EPS_MIN_)
            visco_turb_K_eps[i] = 0;
          else
            {
              if (is_modele_fonc)
                {
                  int is_Cmu_constant = mon_modele_fonc_.Calcul_is_Cmu_constant();
                  if (is_Cmu_constant)
                    visco_turb_K_eps[i] = Fmu(i)*LeCmu_*tab_K(i)*tab_K(i)/(tab_Eps(i)+D(i));
                  else
                    visco_turb_K_eps[i] = Fmu(i)*Cmu(i)*tab_K(i)*tab_K(i)/(tab_Eps(i)+D(i));
                }
              else
                visco_turb_K_eps[i] = LeCmu_*tab_K(i)*tab_K(i)/(tab_Eps(i)+D(i));
            }
        }
      // Debog::verifier("Modele_turbulence_hyd_K_Eps_Bicephale::calculer_viscosite_turbulente visco_turb_K_eps after",visco_turb_K_eps);

      // On connait donc la viscosite turbulente au centre des faces de chaque element
      // On cherche maintenant a interpoler cette viscosite turbulente au centre des
      // elements.
      la_viscosite_turbulente_->affecter(visco_turb_au_format_K_eps.valeur());
      Debog::verifier("Modele_turbulence_hyd_K_Eps_Bicephale::calculer_viscosite_turbulente visco_turb_au_format_K_eps",visco_turb_au_format_K_eps.valeur());
    }
  else
    {
      for (int i=0; i<n; i++)
        {
          if (tab_Eps(i) <= EPS_MIN_)
            {
              visco_turb[i] = 0;
            }
          else
            {
              if (is_modele_fonc)
                {
                  int is_Cmu_constant = mon_modele_fonc_.Calcul_is_Cmu_constant();
                  if (is_Cmu_constant)
                    visco_turb[i] =Fmu(i)*LeCmu_*tab_K(i)*tab_K(i)/(tab_Eps(i)+D(i));
                  else
                    {
                      visco_turb[i] =Fmu(i)*Cmu(i)*tab_K(i)*tab_K(i)/(tab_Eps(i)+D(i));
                    }
                }
              else
                visco_turb[i] = LeCmu_*tab_K(i)*tab_K(i)/(tab_Eps(i)+D(i));
            }
        }
    }

  la_viscosite_turbulente_.changer_temps(temps);
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Bicephale::calculer_viscosite_turbulente la_viscosite_turbulente after",la_viscosite_turbulente_.valeurs());
  return la_viscosite_turbulente_;
}

void imprimer_evolution_keps(const Champ_Inc& le_champ_K,const Champ_Inc& le_champ_Eps, const Schema_Temps_base& sch, double LeCmu, int avant)
{
  if (sch.nb_pas_dt()==0 || sch.limpr())
    {
      const DoubleTab& K   = le_champ_K.valeurs();
      const DoubleTab& Eps = le_champ_Eps.valeurs();
      double k_min=DMAXFLOAT;
      double eps_min=DMAXFLOAT;
      double nut_min=DMAXFLOAT;
      double k_max=0;
      double eps_max=0;
      double nut_max=0;
      int loc_k_min=-1;
      int loc_eps_min=-1;
      int loc_nut_min=-1;
      int loc_k_max=-1;
      int loc_eps_max=-1;
      int loc_nut_max=-1;
      int size=K.dimension(0);
      if (size<0)
        {
          if (sub_type(Champ_Inc_P0_base, le_champ_K.valeur()))
            size = le_champ_K.valeur().equation().domaine_dis().domaine().nb_elem();
          else
            {
              Cerr << "Unsupported K field in Modele_turbulence_hyd_K_Eps_Bicephale::imprimer_evolution_keps()" << finl;
              Process::exit(-1);
            }
          if (sub_type(Champ_Inc_P0_base, le_champ_Eps.valeur()))
            size = le_champ_Eps.valeur().equation().domaine_dis().domaine().nb_elem();
          else
            {
              Cerr << "Unsupported epsilon field in Modele_turbulence_hyd_K_Eps_Bicephale::imprimer_evolution_keps()" << finl;
              Process::exit(-1);
            }
        }
      //ConstDoubleTab_parts parts(le_champ_K.valeurs());
      //ConstDoubleTab_parts parts2(le_champ_Eps.valeurs());
      for (int n=0; n<size; n++)
        {
          const double k = K(n);
          const double eps = Eps(n);
          double nut = 0;
          if (eps > 0) nut = LeCmu*k*k/eps;
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
      /*
      k_min = Process::mp_min(k_min);
      eps_min = Process::mp_min(eps_min);
      nut_min = Process::mp_min(nut_min);
      k_max = Process::mp_max(k_max);
      eps_max = Process::mp_max(eps_max);
      nut_max = Process::mp_max(nut_max);
       */
      ArrOfDouble values(3);

      values[0]=k_min;
      values[1]=eps_min;
      values[2]=nut_min;
      mp_min_for_each_item(values);
      k_min=values[0];
      eps_min=values[1];
      nut_min=values[2];

      values[0]=k_max;
      values[1]=eps_max;
      values[2]=nut_max;
      mp_max_for_each_item(values);
      k_max=values[0];
      eps_max=values[1];
      nut_max=values[2];
      if (Process::je_suis_maitre())
        {
          Cout << finl << "K_Eps evolution (" << (avant?"before":"after") << " law of the wall applies) at time " << le_champ_K.temps() << ":" << finl;
          Cout << "std::min(k)=" << k_min;
          if (Process::nproc()==1) Cout << " located at node " << loc_k_min;
          Cout << finl;
          Cout << "std::min(eps)=" << eps_min;
          if (Process::nproc()==1) Cout << " located at node " << loc_eps_min;
          Cout << finl;
          Cout << "std::min(nut)=" << nut_min;
          if (Process::nproc()==1) Cout << " located at node " << loc_nut_min;
          Cout << finl;
          Cout << "std::max(k)=" << k_max;
          if (Process::nproc()==1) Cout << " located at node " << loc_k_max;
          Cout << finl;
          Cout << "std::max(eps)=" << eps_max;
          if (Process::nproc()==1) Cout << " located at node " << loc_eps_max;
          Cout << finl;
          Cout << "std::max(nut)=" << nut_max;
          if (Process::nproc()==1) Cout << " located at node " << loc_nut_max;
          Cout << finl;
        }
    }
}

int Modele_turbulence_hyd_K_Eps_Bicephale::preparer_calcul()
{
  eqn_transp_K().preparer_calcul();
  eqn_transp_Eps().preparer_calcul();
  Modele_turbulence_hyd_base::preparer_calcul();
  // GF pour initialiser la loi de paroi thermique en TBLE
  if (equation().probleme().nombre_d_equations()>1)
    {
      const RefObjU& modele_turbulence = equation().probleme().equation(1).get_modele(TURBULENCE);
      if (sub_type(Modele_turbulence_scal_base,modele_turbulence.valeur()))
        {
          Turbulence_paroi_scal& loi_paroi_T = ref_cast_non_const(Modele_turbulence_scal_base,modele_turbulence.valeur()).loi_paroi();
          loi_paroi_T.init_lois_paroi();
        }
    }
  // GF quand on demarre un calcul il est bon d'utliser la ldp
  // encore plus quand on fait une reprise !!!!!!!!
  Champ_Inc& ch_K = K();
  Champ_Inc& ch_Eps = Eps();

  const Milieu_base& mil=equation().probleme().milieu();
  if (equation().probleme().is_dilatable())
    {
      diviser_par_rho_si_dilatable(ch_K.valeurs(),mil);
      diviser_par_rho_si_dilatable(ch_Eps.valeurs(),mil);
    }
  imprimer_evolution_keps(ch_K,ch_Eps,eqn_transp_K().schema_temps(),LeCmu_,1);

  loipar_.calculer_hyd_BiK(ch_K,ch_Eps);

  eqn_transp_K().controler_variable();
  eqn_transp_Eps().controler_variable();
  calculer_viscosite_turbulente(ch_K.temps());
  limiter_viscosite_turbulente();

  // on remultiplie K et eps par rho
  if (equation().probleme().is_dilatable())
    {
      multiplier_par_rho_si_dilatable(ch_K.valeurs(),mil);
      multiplier_par_rho_si_dilatable(ch_Eps.valeurs(),mil);
      correction_nut_et_cisaillement_paroi_si_qc(*this);
    }
  la_viscosite_turbulente_.valeurs().echange_espace_virtuel();
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Bicephale::preparer_calcul la_viscosite_turbulente",la_viscosite_turbulente_.valeurs());
  imprimer_evolution_keps(ch_K,ch_Eps,eqn_transp_K().schema_temps(),LeCmu_,0);
  return 1;
}

bool Modele_turbulence_hyd_K_Eps_Bicephale::initTimeStep(double dt)
{
  return (eqn_transport_K_.initTimeStep(dt) and eqn_transport_Eps_.initTimeStep(dt));
}

/*! @brief Effectue une mise a jour en temps du modele de turbulence.
 *
 * Met a jour les equations de transport de K et epsilon,
 *     calcule la viscosite turbulente
 *     au nouveau temps.
 *
 * @param (double temps) le temps de mise a jour
 */
void Modele_turbulence_hyd_K_Eps_Bicephale::mettre_a_jour(double temps)
{
  Champ_Inc& ch_K = K();
  Champ_Inc& ch_Eps = Eps();
  Schema_Temps_base& sch  = eqn_transp_K().schema_temps();
  Schema_Temps_base& sch2 = eqn_transp_Eps().schema_temps();
  // Voir Schema_Temps_base::faire_un_pas_de_temps_pb_base
  eqn_transp_K().domaine_Cl_dis().mettre_a_jour(temps);
  eqn_transp_Eps().domaine_Cl_dis().mettre_a_jour(temps);
  if (!eqn_transp_K().equation_non_resolue())
    sch.faire_un_pas_de_temps_eqn_base(eqn_transp_K());
  eqn_transp_K().mettre_a_jour(temps);
  if (!eqn_transp_Eps().equation_non_resolue())
    sch2.faire_un_pas_de_temps_eqn_base(eqn_transp_Eps());
  eqn_transp_Eps().mettre_a_jour(temps);

  statistiques().begin_count(nut_counter_);
  const Milieu_base& mil=equation().probleme().milieu();
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Bicephale::mettre_a_jour la_viscosite_turbulente before",la_viscosite_turbulente_.valeurs());
  // on divise K_eps par rho en QC pour revenir a K et Eps
  if (equation().probleme().is_dilatable())
    {
      diviser_par_rho_si_dilatable(ch_K.valeurs(),mil);
      diviser_par_rho_si_dilatable(ch_Eps.valeurs(),mil);
    }
  imprimer_evolution_keps(ch_K,ch_Eps,eqn_transp_K().schema_temps(),LeCmu_,1);

  loipar_.calculer_hyd_BiK(ch_K,ch_Eps);

  eqn_transp_Eps().controler_variable();
  eqn_transp_K().controler_variable();
  calculer_viscosite_turbulente(ch_K.temps());
  limiter_viscosite_turbulente();
  // on remultiplie Ket eps par rho
  if (equation().probleme().is_dilatable())
    {
      multiplier_par_rho_si_dilatable(ch_K.valeurs(),mil);
      multiplier_par_rho_si_dilatable(ch_Eps.valeurs(),mil);
      correction_nut_et_cisaillement_paroi_si_qc(*this);
    }
  la_viscosite_turbulente_.valeurs().echange_espace_virtuel();
  Debog::verifier("Modele_turbulence_hyd_K_Eps_Bicephale::mettre_a_jour la_viscosite_turbulente after",la_viscosite_turbulente_.valeurs());
  imprimer_evolution_keps(ch_K,ch_Eps,eqn_transp_K().schema_temps(),LeCmu_,0);
  statistiques().end_count(nut_counter_);
}

const Equation_base& Modele_turbulence_hyd_K_Eps_Bicephale::equation_k_eps(int i) const
{
  assert ((i==0)||(i==1));
  if (i==0)
    {
      return eqn_transport_K_;
    }
  else
    {
      return eqn_transport_Eps_;
    }
}
const Champ_base&  Modele_turbulence_hyd_K_Eps_Bicephale::get_champ(const Motcle& nom) const
{

  try
    {
      return Modele_turbulence_hyd_RANS_Bicephale_base::get_champ(nom);
    }
  catch (Champs_compris_erreur)
    {
    }
  if (mon_modele_fonc_.non_nul())
    {
      try
        {
          return  mon_modele_fonc_.valeur().get_champ(nom);
        }
      catch (Champs_compris_erreur)
        {
        }
    }
  throw Champs_compris_erreur();

}
void Modele_turbulence_hyd_K_Eps_Bicephale::get_noms_champs_postraitables(Noms& nom,Option opt) const
{
  Modele_turbulence_hyd_RANS_Bicephale_base::get_noms_champs_postraitables(nom,opt);
  if (mon_modele_fonc_.non_nul())
    mon_modele_fonc_.valeur().get_noms_champs_postraitables(nom,opt);

}
void Modele_turbulence_hyd_K_Eps_Bicephale::verifie_loi_paroi()
{
  Nom lp=loipar_.valeur().que_suis_je();
  if ( !(lp=="negligeable_VEF" || lp=="negligeable_VDF") )
    {
      Cerr<<"The turbulence model of type "<<que_suis_je()<<finl;
      Cerr<<"must be used with a wall law of type negligeable for now."<<finl;
    }
  if (!associe_modele_fonction().non_nul())
    {
      Cerr<<"The turbulence model of type "<<que_suis_je()<<finl;
      Cerr<<"must be used with a modele fonction for now."<<finl;
    }

}
