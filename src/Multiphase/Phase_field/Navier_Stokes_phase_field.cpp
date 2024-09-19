/****************************************************************************
* Copyright (c) 2021, CEA
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
// File:        Navier_Stokes_phase_field.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Multiphase/Phase_field/src
// Version:     /main/27
//
//////////////////////////////////////////////////////////////////////////////


#include <Navier_Stokes_phase_field.h>
#include <Probleme_base.h>
#include <Fluide_Incompressible.h>
#include <Discret_Thyd.h>
#include <Convection_Diffusion_Phase_field.h>
#include <Assembleur_base.h>
#include <Statistiques.h>
#include <TRUSTTrav.h>
#include <Domaine_VF.h>
#include <Domaine.h>
#include <Param.h>
#include <Champ_Fonc_Tabule.h>
#include <EChaine.h>
#include <Source_Con_Phase_field_base.h>
#include <Operateur_Diff_base.h>
#include <Constituant.h>
#include <Source_Con_Phase_field.h>

extern Stat_Counter_Id temps_total_execution_counter_;
Implemente_instanciable_sans_constructeur_ni_destructeur(Navier_Stokes_phase_field,"Navier_Stokes_phase_field",Navier_Stokes_std);
// XD navier_stokes_phase_field navier_stokes_standard navier_stokes_phase_field -1 Navier Stokes equation for the Phase Field problem.

Navier_Stokes_phase_field::Navier_Stokes_phase_field()
{
  // On initialise les attributs de la classe
  boussi_=-1;
  diff_boussi_=-1;
  g_.resize(0);
  /*
    Noms& nom=champs_compris_.liste_noms_compris();
    nom.dimensionner(2);
    nom[0]="masse_volumique";
    nom[1]="viscosite_dynamique";
  */
}

Navier_Stokes_phase_field::~Navier_Stokes_phase_field() {}

/*! @brief Simple appel a: Navier_Stokes_std::printOn(Sortie&) Ecrit l'equation sur un flot de sortie.
 *
 * @param (Sortie& os) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Navier_Stokes_phase_field::printOn(Sortie& is) const
{
  return Navier_Stokes_std::printOn(is);
}

/*! @brief Appel Navier_Stokes_phase_field::readOn(Entree& is) En sortie verifie que l'on a bien lu:
 *
 *         - le terme diffusif,
 *         - le terme convectif,
 *         - le solveur en pression
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 * @throws terme diffusif non specifie dans jeu de donnees, specifier
 * un type negligeable pour l'operateur si il est a negliger
 * @throws terme convectif non specifie dans jeu de donnees, specifier
 * un type negligeable pour l'operateur si il est a negliger
 * @throws solveur pression non defini dans jeu de donnees
 */
Entree& Navier_Stokes_phase_field::readOn(Entree& is)
{
  Navier_Stokes_std::readOn(is);

  //Creation du terme source de gravite pour discretisation VDF
  if (g_.size()>0)
    {
      Source t;
      Source& so=les_sources.add(t);
      Cerr << "Creation of the buoyancy source term for the Navier_Stokes_phase_field equation"<< finl;

      Nom disc = discretisation().que_suis_je();
      if (disc!="VDF")
        {
          Cerr<<"Only VDF discretization is allowed."<<finl;
          exit();
        }

      Nom type_so = "Source_Gravite_PF_VDF";
      so.typer_direct(type_so);
      so->associer_eqn(*this);
    }
  return is;
}

void Navier_Stokes_phase_field::set_param(Param& param)
{
  param.ajouter_non_std("approximation_de_boussinesq",(this),Param::REQUIRED); // XD_ADD_P approx_boussinesq To use or not the Boussinesq approximation.
  param.ajouter_non_std("viscosite_dynamique_constante",(this),Param::OPTIONAL); // XD_ADD_P visco_dyn_cons To use or not a viscosity which will depends on concentration C (in fact, C is the unknown of Cahn-Hilliard equation).
  param.ajouter("gravite",&g_); // XD_ADD_P list Keyword to define gravity in the case Boussinesq approximation is not used.
  Navier_Stokes_std::set_param(param);
}

// XD approx_boussinesq objet_lecture nul 0 different mass density formulation are available depending if the Boussinesq approximation is made or not
// XD attr yes_or_no chaine(into=["oui","non"]) yes_or_no 0 To use or not the Boussinesq approximation.
// XD attr bloc_bouss bloc_boussinesq bloc_bouss 0 to choose the rho formulation
// XD bloc_boussinesq objet_lecture nul 1 choice of rho formulation
// XD attr probleme ref_Pb_base probleme 1 Name of problem.
// XD attr rho_1 floattant rho_1 1 value of rho
// XD attr rho_2 floattant rho_2 1 value of rho
// XD attr rho_fonc_c bloc_rho_fonc_c rho_fonc_c_ 1 to use for define a general form for rho
// XD bloc_rho_fonc_c objet_lecture nul 0 if rho has a general form
// XD attr Champ_Fonc_Fonction chaine(into=["Champ_Fonc_Fonction"]) Champ_Fonc_Fonction 1 Champ_Fonc_Fonction
// XD attr problem_name ref_Pb_base problem_name 1 Name of problem.
// XD attr concentration chaine(into=["concentration"]) concentration 1 concentration
// XD attr dim entier dim 1 dimension of the problem
// XD attr val chaine val 1 function of rho
// XD attr Champ_Uniforme chaine(into=["Champ_Uniforme"]) Champ_Uniforme 1 Champ_Uniforme
// XD attr fielddim entier fielddim 1 dimension of the problem
// XD attr val2 chaine val2 1 function of rho

// XD visco_dyn_cons objet_lecture nul 0 different treatment of the kinematic viscosity could be done depending of the use of the Boussinesq approximation or the constant dynamic viscosity approximation
// XD attr yes_or_no chaine(into=["oui","non"]) yes_or_no 0 To use or not the constant dynamic viscosity
// XD attr bloc_visco bloc_visco2 bloc_visco 0 to choose the mu formulation
// XD bloc_visco2 objet_lecture nul 1 choice of mu formulation
// XD attr probleme ref_Pb_base probleme 1 Name of problem.
// XD attr mu_1 floattant mu_1 1 value of mu
// XD attr mu_2 floattant mu_2 1 value of mu
// XD attr mu_fonc_c bloc_mu_fonc_c mu_fonc_c_ 1 to use for define a general form for mu
// XD bloc_mu_fonc_c objet_lecture nul 0 if mu has a general form
// XD attr Champ_Fonc_Fonction chaine(into=["Champ_Fonc_Fonction"]) Champ_Fonc_Fonction 1 Champ_Fonc_Fonction
// XD attr problem_name ref_Pb_base problem_name 1 Name of problem.
// XD attr concentration chaine(into=["concentration"]) concentration 1 concentration
// XD attr dim entier dim 1 dimension of the problem
// XD attr val chaine val 1 function of mu

int Navier_Stokes_phase_field::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  Motcle motlu;
  if (mot=="approximation_de_boussinesq")
    {
      Motcle temp_boussi_;
      is >> temp_boussi_;
      bool fonctionsALire = false;
      if(temp_boussi_=="oui")
        {
          boussi_=1;
          diff_boussi_=0;
          if (!betac_->non_nul())
            {
              fonctionsALire = true;
              Cerr << "Approximation de Boussinesq utilisee mais beta_co non defini comme propriete du fluide"<< finl;
              Cerr << "On s'attend donc a lire les fonctions rho(c) et drho/dc"<< finl;
            }
        }
      else
        {
          boussi_=0;
          fonctionsALire = true;
        }
      is >> motlu ;
      if (motlu != "{")
        {
          Cerr << "Error while reading approximation_de_boussinesq" << finl;
          Cerr << "We expected a { instead of " << motlu << finl;
          exit();
        }
      is >> motlu ;
      if (fonctionsALire)
        {
          if (motlu=="rho_fonc_c")
            {
              Cerr << motlu << finl;
              is >> rho_;
              is >> drhodc_;
              is >> motlu;
            }
          else
            {
              double rho1,rho2;
              Nom prob;
              while(motlu!="}")
                {
                  if (motlu=="probleme") is >> prob;
                  if (motlu=="rho_1") is >> rho1;
                  if (motlu=="rho_2") is >> rho2;
                  is >> motlu;
                }
              Nom chaine("Champ_Fonc_Fonction ");
              chaine+=prob;
              chaine+=" concentration 1 ";
              Nom rhoM(0.5*(rho1+rho2));
              Nom drho(rho2-rho1);
              chaine+=rhoM;
              chaine+="+val*(";
              chaine+=drho;
              chaine+=")";
              EChaine echaine(chaine);
              echaine >> rho_;
              Nom chaine2("Champ_Uniforme 1 ");
              chaine2+=drho;
              EChaine echaine2(chaine2);
              Cerr << "@@@ " << chaine2 << finl;
              echaine2 >> drhodc_;
            }
        }
      if (motlu != "}")
        {
          Cerr << "Error while reading approximation_de_boussinesq" << finl;
          Cerr << "We expected a } instead of " << motlu << finl;
          exit();
        }
      return 1;
    }
  else if (mot=="viscosite_dynamique_constante")
    {
      Motcle temp_diffbou_;
      is >> temp_diffbou_;
      if(temp_diffbou_=="oui")
        {
          diff_boussi_=1;
        }
      else
        {
          diff_boussi_=0;
          is >> motlu ;
          assert(motlu=="{");
          is >> motlu ;
          if (motlu=="mu_fonc_c")
            {
              Cerr << motlu << finl;
              is >> mu_;
              is >> motlu;
            }
          else
            {
              double mu1,mu2;
              Nom prob;
              while(motlu!="}")
                {
                  if (motlu=="probleme") is >> prob;
                  if (motlu=="mu_1") is >> mu1;
                  if (motlu=="mu_2") is >> mu2;
                  is >> motlu;
                }
              Nom chaine("Champ_Fonc_Fonction ");
              chaine+=prob;
              chaine+=" concentration 1 ";
              Nom muM(0.5*(mu1+mu2));
              Nom dmu(mu2-mu1);
              chaine+=muM;
              chaine+="+val*(";
              chaine+=dmu;
              chaine+=")";
              EChaine echaine(chaine);
              echaine >> mu_;
            }
          assert(motlu=="}");
          const Discret_Thyd& dis=ref_cast(Discret_Thyd, discretisation());
          dis.nommer_completer_champ_physique(domaine_dis(),"viscosite_dynamique","Pa.s",mu_,probleme());
          champs_compris_.ajoute_champ(mu_);
        }
      if(diff_boussi_==1)
        {
          // RLT: pas sur de tout a fait comprendre si c'est encore d'actualite mais, en tout cas, en lien avec modification
          // commentee dans Navier_Stokes_phase_field::calculer_pas_de_temps
          Cerr << "LE PAS DE TEMPS DE STABILITE (DE DIFFUSION) DOIT ETRE FORCE DANS VDF/Operateurs/OpDifVDFFaCs.cpp : ligne 115" << finl;
          Cerr << "ATTENTION : LE PAS DE TEMPS DE DIFFUSION SERA ALORS FORCE A LA VALEUR DE 1 SECONDE !!!" << finl;
          // Ceci est utile car normalement dans ce cas le pas de temps de diffusion ne depend pas de rho
          // car la viscosite dynamique est constante (on prend alors mu/rho0)
          // Or dans ce cas TRUST semble calculer ce pas temps avec une influence venant de rho (Debug ?)
        }
      return 1;
    }
  else
    return Navier_Stokes_std::lire_motcle_non_standard(mot,is);
}

const Champ_Don& Navier_Stokes_phase_field::diffusivite_pour_transport() const
{
  if(boussi_==0)
    {
      if (diff_boussi_==0)
        return mu_;
      else if (diff_boussi_==1)
        return le_fluide->viscosite_dynamique();
      else
        {
          Cerr << "If the Boussinesq hypothesis is not done," << finl;
          Cerr << "the keyword viscosite_dynamique_constante must be used to indicate" << finl;
          Cerr << "if the dynamical viscosity is constant or c dependent." << finl;
          Cerr << "This selection must be done before reading the diffusion term" << finl;
          Cerr << "of the Navier_Stokes_phase_field equation." << finl;
          exit();
        }
    }
  else if (boussi_==1)
    return le_fluide->viscosite_cinematique();
  else
    {
      Cerr << "The use of Boussinesq hypothesis must be indicated"<<finl;
      Cerr << "by the keyword approximation_de_boussinesq specification."<<finl;
      Cerr << "This selection must be done before reading the diffusion term" << finl;
      Cerr << "of the Navier_Stokes_phase_field equation." << finl;
      exit();
    }
  //Pour compilation
  return Navier_Stokes_phase_field::diffusivite_pour_transport();
}

/*! @brief Dicretise l'equation.
 *
 */
void Navier_Stokes_phase_field::discretiser()
{
  Navier_Stokes_std::discretiser();

  // Selection de l'assembleur en pression specifique au Phase_Field (P_VDF de BM)
  Nom type = "";
  type = "Assembleur_P_VDF_Phase_Field";
  Cerr << "** Pressure assembling tool : " << type << " **" << finl;
  assembleur_pression_.typer(type);
  assembleur_pression_->associer_domaine_dis_base(domaine_dis());
  assembleur_pression_->set_resoudre_increment_pression(1);
  // la discretisation du champ associe
  // - a la masse volumique est "retardee" a la methode creer_champ/completer
  // - a la viscosite dynamique est "retardee" a la methode lire_motcle_non_standard
  //   (et non pas creer_champ/completer comme pour la masse volumique sinon probleme dans terme_diffusif.associer_diffusivite dans Navier_Stokes_std::lire_motcle_non_standard)
  // afin que les options boussi_ et diff_boussi_ soient connues

  const DoubleTab& rho0Tab=le_fluide->masse_volumique()->valeurs();
  int dim=rho0Tab.nb_dim();
  switch(dim)
    {
    case 2:
      rho0_=rho0Tab(0,0);
      break;
    case 3:
      rho0_=rho0Tab(0,0,0);
      break;
    default:
      Cerr <<"N_S_Phase_field : Problem with rho0 dimension:"<<dim<<finl;
      exit();
      break;
    }
  betac_ = le_fluide->beta_c();
}

void Navier_Stokes_phase_field::completer()
{
  Navier_Stokes_std::completer();
  creer_champ("masse_volumique");
  // RLT: cf. commentaire dans Navier_Stokes_phase_field::calculer_pas_de_temps
  /*if (boussi_==0)
    {
      if (diff_boussi_==0)
        {
          terme_diffusif->associer_champ_masse_volumique(rho_);
        }
      else
        {
          Cout << "Navier_Stokes_phase_field::completer rho=" << le_fluide->masse_volumique() << finl;
          terme_diffusif->associer_champ_masse_volumique(le_fluide->masse_volumique());
        }
    }*/
}


void Navier_Stokes_phase_field::creer_champ(const Motcle& motlu)
{
  Navier_Stokes_std::creer_champ(motlu);
  if (motlu == "masse_volumique")
    {
      Cerr << "Navier_Stokes_phase_field::creer_champ " << motlu << " " << (int)rho_.non_nul() << " " << rho_.le_nom() << finl;
      const Discret_Thyd& dis=ref_cast(Discret_Thyd, discretisation());
      if (!rho_.non_nul())
        {
          if (boussi_ != 1)
            {
              Cerr << "Navier_Stokes_phase_field::creer_champ: should not be here (1)"<< finl;
              exit();
            }
          dis.discretiser_champ("CHAMP_ELEM",domaine_dis(),"masse_volumique","kg/m3",1,0.,rho_);
          champs_compris_.ajoute_champ(rho_);
          dis.discretiser_champ("CHAMP_ELEM",domaine_dis(),"derivee_masse_volumique","kg/m3",1,0.,drhodc_);
          champs_compris_.ajoute_champ(drhodc_);
        }
      if (rho_.le_nom()=="anonyme")
        {
          if (boussi_ == 1 && betac_->non_nul())
            {
              Cerr << "Navier_Stokes_phase_field::creer_champ: should not be here (2)"<< finl;
              exit();
            }
          dis.nommer_completer_champ_physique(domaine_dis(),"masse_volumique","kg/m3",rho_,probleme());
          champs_compris_.ajoute_champ(rho_);
          dis.nommer_completer_champ_physique(domaine_dis(),"derivee_masse_volumique","kg/m3",drhodc_,probleme());
          champs_compris_.ajoute_champ(drhodc_);
        }
    }
}

/*! @brief Calcule rho_n+1 pour utilisation dans la pression
 *
 */
void Navier_Stokes_phase_field::calculer_rho(const bool init)
{
  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field, mon_probleme->equation(1));
  Sources& list_sources = ref_cast_non_const(Sources, eq_c.sources());
  Source_Con_Phase_field& source_pf = ref_cast(Source_Con_Phase_field, list_sources(0).valeur());
  int type_systeme_naire = source_pf.get_type_systeme_naire();

  int i = 1;
  double temps = schema_temps().temps_futur(i);
  if (init)
    {
      i = 0;
      temps = schema_temps().temps_courant();
    }
  if (boussi_==0 || !betac_->non_nul())
    {
      // normalement, quelque chose comme 'rho_->mettre_a_jour(schema_temps().temps_futur(i))' devrait faire l'affaire
      // probleme, cette evaluation se fait avec c(n) alors que l'on veut evaluer avec c(n+1)
      // on introduit une methode specifique (qui depend de la discretisation) que l'on met arbitrairement dans Source_Con_Phase_field_base
      Source_Con_Phase_field_base& source_pf_base = ref_cast(Source_Con_Phase_field_base, list_sources(0).valeur());
      source_pf_base.calculer_champ_fonc_c(temps, rho_, eq_c.inconnue().futur(i));
      source_pf_base.calculer_champ_fonc_c(temps, drhodc_, eq_c.inconnue().futur(i));
      //Cerr << "rho = " << rho_->valeurs()[14*48+24] << finl;
    }
  else
    {
      const DoubleTab& c = eq_c.inconnue().futur(i);
      //Mirantsoa 264902
      DoubleTab& rhoTab = rho_->valeurs(); /**/
      DoubleTab& drhodcTab = drhodc_->valeurs();/**/
      DoubleTab rho0Tab(rhoTab);
      //Cerr << "c = "<<c<<finl;
      // DoubleTab& betacTab = betac_->valeurs();
      //      rho_->affecter_(betac_.valeur());


      if (type_systeme_naire==0)
        {
          drhodcTab=rho0_;
          tab_multiply_any_shape(drhodcTab, betac_.valeur()->valeurs());
          rhoTab=c;
          tab_multiply_any_shape(rhoTab, drhodcTab);
          rhoTab+=rho0_;
          //Cerr << "c = "<< c<<finl;

        }
      else if (type_systeme_naire==1)
        {
          DoubleTab beta_(rhoTab);
          DoubleTab betacTab(c);

          //calcul de drhodcTab comme rho0*somme(betac)
          drhodcTab=rho0_;
          for (i=0; i<beta_.dimension(0); i++)
            {
              for (int j=0; j<betacTab.line_size(); j++)
                {
                  betacTab(i,j)=betac_.valeur()->valeurs()(0,j);
                  beta_(i)+=betacTab(i,j);
                }
            }
          tab_multiply_any_shape(drhodcTab, beta_);

          //calcul de rhoTab comme rho0+rho0*somme(betac*c)
          rhoTab=0;
          tab_multiply_any_shape(betacTab, c);
          for (i=0; i<betacTab.dimension(0); i++)
            {
              for (int j=0; j<betacTab.line_size(); j++)
                {
                  rhoTab(i)+=betacTab(i,j);
                }
            }
          rho0Tab=rho0_;
          tab_multiply_any_shape(rhoTab, rho0Tab);
          rhoTab+=rho0_;
          /*Cerr<<"rhoTab+"<<rhoTab<<finl;
          Cerr << "beta_ = "<<beta_<<finl;
          Cerr << "betacTab = "<< betacTab<<finl;
          Cerr << "c = "<< c<<finl;*/


        }
    }
  rho_->valeurs().echange_espace_virtuel();
  drhodc_->valeurs().echange_espace_virtuel();
}

void Navier_Stokes_phase_field::calculer_mu(const bool init)
{
  if (boussi_==0 && diff_boussi_==0)
    {
      const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field, mon_probleme->equation(1));
      int i = 1;
      double temps = schema_temps().temps_futur(i);
      if (init)
        {
          i = 0;
          temps = schema_temps().temps_courant();
        }
      // normalement, quelque chose comme 'mu_->mettre_a_jour(schema_temps().temps_futur(i))' devrait faire l'affaire
      // probleme, cette evaluation se fait avec c(n) alors que l'on veut evaluer avec c(n+1)
      // on introduit une methode specifique (qui depend de la discretisation) que l'on met arbitrairement dans Source_Con_Phase_field_base
      Sources& list_sources = ref_cast_non_const(Sources, eq_c.sources());
      Source_Con_Phase_field_base& source_pf = ref_cast(Source_Con_Phase_field_base, list_sources(0).valeur());
      source_pf.calculer_champ_fonc_c(temps, mu_, eq_c.inconnue().futur(i));
      mu_->valeurs().echange_espace_virtuel();
    }
}

double Navier_Stokes_phase_field::calculer_pas_de_temps() const
{
  const  Navier_Stokes_phase_field& eq_ns=ref_cast(Navier_Stokes_phase_field,probleme().equation(0));
  // Convection_Diffusion_Phase_field& eq_c=ref_cast_non_const(Convection_Diffusion_Phase_field,probleme().equation(1));
  // RLT: plutot que de multiplier dt(diff) par rho ici comme c'est fait dans TrioCFD 1.8.0 (cf. coeff ci-dessous)
  // (en considerant qu'il a ete evalue dans l'operateur de diffusion (Op_Diff_VDF_Face_base) avec la viscosite dynamique - ce qui n'est vrai que pour boussi_==0)
  // on associe le champ masse volumique a l'operateur de diffusion dans le cas boussi_==0 pour que le calcul de son dt soit correct directement
  double coeff = mp_min_vect(eq_ns.rho()->valeurs()); // eq_c.rho1(); // voir d'apres le FT

  double dt=0;
  double dt_op;

  Cerr << "" << finl;
  double tps_courant=eq_ns.schema_temps().temps_courant();
  double num_dt=eq_ns.schema_temps().nb_pas_dt();
  Cerr << "                           *--*" << finl;
  Cerr << "" << finl;
  Cerr << "Calculation at time : " << tps_courant << " s" << finl;
  Cerr << "Time step number : " << num_dt << finl;
  Cerr << "" << finl;

  dt = 1./terme_convectif.calculer_pas_de_temps();
  if (le_schema_en_temps->limpr())
    {
      Cout << "NS_PF : pas de temps de convection : " << 1./dt << finl;
    }
  dt_op = terme_diffusif.calculer_pas_de_temps()*coeff;
  if (le_schema_en_temps->limpr())
    {
      Cout << "NS_PF : pas de temps de diffusion : " << dt_op << finl;
    }
  dt=dt+1./dt_op;
  dt=std::min(1./dt,le_schema_en_temps->pas_temps_max());
  if (le_schema_en_temps->limpr())
    {
      Cout << "NS_PF : pas de temps de calcul : " << dt << finl;
    }

  // Calcul du temps effectue en simulation, ainsi que de l'estimation du temps restant
  double t_total;
  double t_remaining;
  double tps_init=eq_ns.schema_temps().temps_init();
  double tps_max=eq_ns.schema_temps().temps_max();

  double t_passe=statistiques().last_time(temps_total_execution_counter_);
  if (tps_courant!=tps_init)
    {
      t_total=t_passe*((tps_max-tps_init)/(tps_courant-tps_init));
      t_remaining=t_total-t_passe;
      Cerr << "" << finl;
      Cerr << "**************** ELAPSED TIME FOR THE SIMULATION *******************" << finl;
      Cerr << "Elapsed time of calculation until current time : " << floor(t_passe) << " s" << finl;
      Cerr << "Estimated time needed for the remaining part of the calculation : " << floor(t_remaining) << " s" << finl;
      Cerr << "                            --" << finl;
      Cerr << "Estimated time for the complete calculation : " << floor(t_total) << " s" << finl;
      Cerr << "********************************************************************" << finl;
      Cerr << "" << finl;
    }
  return dt;
}

/*! @brief Interpole la masse volumique aux faces (pour l'assembleur pression en particulier)
 *
 */
void Navier_Stokes_phase_field::rho_aux_faces(const DoubleTab& tab_rho, Champ_Don& rho_face_P)
{
  const Domaine_VF& domaineVF = ref_cast(Domaine_VF,domaine_dis());
  const IntTab& face_voisins = domaineVF.face_voisins();
  const DoubleVect& volumes = domaineVF.volumes();
  int nbfaces=domaineVF.nb_faces();
  const int nbfaces_bord = domaineVF.premiere_face_int();
  int el0,el1;
  double vol0,vol1;

  // Calcul de rho_face_P aux faces pour l'assemblage avec l'assembleur pression PF (cf AssembleurPVDF de BM)

  for (int fac=0; fac<nbfaces_bord; fac++)
    {
      el0=face_voisins(fac,0);
      if(el0!=-1)
        {
          rho_face_P->valeurs()(fac) = tab_rho(el0);
        }
      else
        {
          el1=face_voisins(fac,1);
          rho_face_P->valeurs()(fac) = tab_rho(el1);
        }
    }

  for (int fac=nbfaces_bord; fac<nbfaces; fac++)
    {
      el0=face_voisins(fac,0);
      el1=face_voisins(fac,1);
      vol0=volumes(el0);
      vol1=volumes(el1);
      rho_face_P->valeurs()(fac)=(vol0*tab_rho(el0)+vol1*tab_rho(el1))/(vol0+vol1);
    }
}

/*! @brief Dicretise l'equation.
 *
 */
int Navier_Stokes_phase_field::preparer_calcul()
{
  Cerr << "** Navier_Stokes_phase_field::preparer_calcul **" << finl;
  //   Cerr << " - Copie via le Front_Tracking et modifie pour le Phase_field -" << finl;
  Equation_base::preparer_calcul();

  // Calcul de la masse volumique
  rho_->initialiser(schema_temps().temps_courant());
  //Cerr << " rho_"<< rho_<< finl;
  //Cerr << "valeur de rho_"<< rho_.valeur()<< finl;
  //Cerr << "valeur de rho_ valeurs"<< rho_.valeurs()<< finl;

  drhodc_->initialiser(schema_temps().temps_courant());
  //Cerr << " drhodc_"<< drhodc_<< finl;
  //Cerr << "valeur de rho_"<< drhodc_.valeur()<< finl;
  //Cerr << "valeur de rho_ valeurs"<< drhodc_.valeurs()<< finl;

  calculer_rho(true); // RLT: cette minitialisation etait incorrecte dans TrioCFD 1.8.0 car on utilisait c.futur()
  // Calcul de la viscosite dynamique dans le cas boussi_==0 et diff_boussi_==0
  if (mu_.non_nul())
    {
      mu_->initialiser(schema_temps().temps_courant());
      calculer_mu(true); // RLT: cette minitialisation etait incorrecte dans dans TrioCFD 1.8.0 car on utilisait c.futur()
    }

  if(boussi_==0)
    {
      Champ_Don rho_face_P;
      rho_face_P.typer("Champ_Fonc_Face");
      rho_face_P->valeurs()=inconnue().valeurs();
      const DoubleTab& tab_rho=rho_->valeurs();

      rho_aux_faces(tab_rho, rho_face_P);

      assembleur_pression_->assembler_rho_variable(matrice_pression_, rho_face_P.valeur());
    }
  else
    {
      assembleur_pression_->assembler(matrice_pression_);
    }

  // -------------------------------------------------------------------------------------------------------

  //  la_pression.valeurs()=0.;
  la_pression->changer_temps(schema_temps().temps_courant());
  rho_->changer_temps(schema_temps().temps_courant());
  drhodc_->changer_temps(schema_temps().temps_courant());
  if (mu_.non_nul())
    {
      mu_->changer_temps(schema_temps().temps_courant());
    }

  la_pression->valeurs().echange_espace_virtuel();

  Debog::verifier("Navier_Stokes_phase_field::preparer_calcul, gradP av", gradient_P->valeurs());
  Debog::verifier("Navier_Stokes_phase_field::preparer_calcul, la_pression av", la_pression->valeurs());
  //Cerr << "Navier_Stokes_phase_field::preparer_calcul() avant grad" << finl;
  gradient.calculer(la_pression->valeurs(),
                    gradient_P->valeurs());
  Debog::verifier("Navier_Stokes_phase_field::preparer_calcul, gradP ap", gradient_P->valeurs());
  Debog::verifier("Navier_Stokes_phase_field::preparer_calcul, la_pression ap", la_pression->valeurs());

  gradient_P->changer_temps(schema_temps().temps_courant());
  divergence_U->valeurs()=0.;
  divergence_U->changer_temps(schema_temps().temps_courant());
  Cerr << "Projection of initial and boundaries conditions " << finl;

  Debog::verifier("Navier_Stokes_phase_field::preparer_calcul, gradP av projeter", gradient_P->valeurs());
  Debog::verifier("Navier_Stokes_phase_field::preparer_calcul, la_pression av projeter", la_pression->valeurs());
  if (projection_a_faire()) projeter();
  Debog::verifier("Navier_Stokes_phase_field::preparer_calcul, gradP ap projeter", gradient_P->valeurs());
  Debog::verifier("Navier_Stokes_phase_field::preparer_calcul, la_pression ap projeter", la_pression->valeurs());

  Navier_Stokes_std::calculer_la_pression_en_pa();
  Debog::verifier("Navier_Stokes_phase_field::preparer_calcul, vitesse", inconnue());
  Debog::verifier("Navier_Stokes_phase_field::preparer_calcul, pression", la_pression);

  if (le_traitement_particulier.non_nul())
    le_traitement_particulier->preparer_calcul_particulier();

  return 1;
}

void Navier_Stokes_phase_field::mettre_a_jour(double temps)
{
  Navier_Stokes_std::mettre_a_jour(temps);
  calculer_rho();
  rho_->changer_temps(temps);
  drhodc_->changer_temps(temps);
  if (mu_.non_nul())
    {
      calculer_mu(); // RLT: cette mise a jour etait manquante dans TrioCFD 1.8.0
      mu_->changer_temps(temps);
    }

  if(boussi_==0)
    {
      Champ_Don rho_face_P;
      rho_face_P.typer("Champ_Fonc_Face");
      rho_face_P->valeurs()=inconnue().valeurs();
      const DoubleTab& tab_rho=rho_->valeurs();

      rho_aux_faces(tab_rho, rho_face_P);

      assembleur_pression_->assembler_rho_variable(matrice_pression_, rho_face_P.valeur());
    }
  else
    {
      assembleur_pression_->assembler(matrice_pression_);
    }

  // -------------------------------------------------------------------------------------------------------

  //   Cerr<<"Navier_Stokes_phase_field::mettre_a_jour "<<temps<<finl;
}

void Navier_Stokes_phase_field::_aff_donnee_P0(const DoubleTab& tab, const Motcle& mot) const
{
  int ielem;
  const Domaine_VF& domaineVF = ref_cast(Domaine_VF,domaine_dis());
  int nbelems = domaineVF.nb_elem_tot();
  const DoubleTab xp = domaineVF.xp();
  Cerr<<mot<<finl;
  for (ielem=0 ; ielem<nbelems ; ielem++)
    {
      //if (xv(iface,1)>0.00215 && xv(iface,2)<0.00005)
      //if (ielem==43)
      {
        Cerr<<"elem "<<ielem<<"  coord= "<<xp(ielem,0)<<" "<<xp(ielem,1);
        if (dimension==3)
          Cerr<<" "<<xp(ielem,2);
        if (tab.nb_dim()==2)
          {
            Cerr<<"  "<<mot<<"= "<<tab(ielem,0)<<" "<<tab(ielem,1);
            if (dimension==3)
              Cerr<<" "<<tab(ielem,2);
            Cerr<<" "<<finl;
          }
        else
          {
            Cerr<<"  "<<mot<<"= "<<tab(ielem)<<finl;
          }
      }
    }
}
void Navier_Stokes_phase_field::_aff_donnee_P1(const DoubleTab& tab, const Motcle& mot) const
{
  int iface;
  const Domaine_VF& domaineVF = ref_cast(Domaine_VF,domaine_dis());
  int nbfaces = domaineVF.nb_faces();
  const DoubleTab xv = domaineVF.xv();
  Cerr<<mot<<finl;
  for (iface=0 ; iface<nbfaces ; iface++)
    {
      //if (xv(iface,1)>0.00215 && xv(iface,2)<0.00005)
      //if (iface==1188 || iface==1189 || iface==1190 || iface ==1191)
      {
        Cerr<<"face "<<iface<<"  coord= "<<xv(iface,0)<<" "<<xv(iface,1);
        if (dimension==3)
          Cerr<<" "<<xv(iface,2);
        if (tab.nb_dim()==2)
          {
            Cerr<<"  "<<mot<<"= "<<tab(iface,0)<<" "<<tab(iface,1);
            if (dimension==3)
              Cerr<<" "<<tab(iface,2);
            Cerr<<" "<<finl;
          }
        else
          {
            Cerr<<"  "<<mot<<"= "<<tab(iface)<<finl;
          }
      }
    }
}
void Navier_Stokes_phase_field::_aff_donnee_P1Bulle(const DoubleTab& tab, const Motcle& mot) const
{
  int ielem,isom;
  const Domaine_VF& domaineVF = ref_cast(Domaine_VF,domaine_dis());
  int nbelems = domaineVF.nb_elem_tot();
  int nbsoms = domaineVF.nb_som();
  const DoubleTab xp = domaineVF.xp();
  const DoubleTab xs = domaineVF.domaine().coord_sommets();
  Cerr<<mot<<finl;
  for (ielem=0 ; ielem<nbelems ; ielem++)
    {
      //if (xv(iface,1)>0.00215 && xv(iface,2)<0.00005)
      //if (ielem==43)
      {
        Cerr<<"elem "<<ielem<<"  coord= "<<xp(ielem,0)<<" "<<xp(ielem,1);
        if (dimension==3)
          Cerr<<" "<<xp(ielem,2);
        if (tab.nb_dim()==2)
          {
            Cerr<<"  "<<mot<<"= "<<tab(ielem,0)<<" "<<tab(ielem,1);
            if (dimension==3)
              Cerr<<" "<<tab(ielem,2);
            Cerr<<" "<<finl;
          }
        else
          {
            Cerr<<"  "<<mot<<"= "<<tab(ielem)<<finl;
          }
      }
    }
  for (isom=0 ; isom<nbsoms ; isom++)
    {
      //if (xv(iface,1)>0.00215 && xv(iface,2)<0.00005)
      //if (isom==25 || isom==219 || isom==226 || isom ==240)
      {
        Cerr<<"som  "<<isom<<"  coord= "<<xs(isom,0)<<" "<<xs(isom,1);
        if (dimension==3)
          Cerr<<" "<<xs(isom,2);
        if (tab.nb_dim()==2)
          {
            Cerr<<"  "<<mot<<"= "<<tab(nbelems+isom,0)<<" "<<tab(nbelems+isom,1);
            if (dimension==3)
              Cerr<<" "<<tab(nbelems+isom,2);
            Cerr<<" "<<finl;
          }
        else
          {
            Cerr<<"  "<<mot<<"= "<<tab(nbelems+isom)<<finl;
          }
      }
    }
}

DoubleTab& Navier_Stokes_phase_field::derivee_en_temps_inco(DoubleTab& vpoint)
{

  if(boussi_==0)
    {

      //   Cerr << "** Navier_Stokes_phase_field::derivee_en_temps_inco **" << finl;
      //   Cerr << " - Copie via Navier_Stokes_Front_Tracking et modifie pour le Phase_field - " << finl;
      DoubleTrav secmem(la_pression->valeurs());
      DoubleTrav gradP(la_vitesse->valeurs());
      double dt=le_schema_en_temps->pas_de_temps();

      //   const Probleme_base& pb=ref_cast(Probleme_base, probleme());
      const Domaine_VF& zvf=ref_cast(Domaine_VF, domaine_dis());
      const DoubleTab& tab_rho=rho_->valeurs();
      const IntTab& face_voisins = zvf.face_voisins();
      const int nbfaces_bord = zvf.premiere_face_int();
      const int nbfaces = zvf.nb_faces();
      //   const int nbelems = zvf.nb_elem_tot();

      // Calcul de rho sur les faces
      DoubleTab rho_face(nbfaces);
      int face, elem1,elem2, k;
      double rho_moyen;
      for (face=0; face<nbfaces_bord; face++)
        {
          elem1=face_voisins(face, 0);
          if(elem1!=-1)
            {
              rho_face[face] = tab_rho[elem1];
            }
          else
            {
              elem2=face_voisins(face, 1);
              rho_face[face] = tab_rho[elem2];
            }
        }
      for (face=nbfaces_bord; face<nbfaces; face++)
        {
          // On generalise le calcule de la masse volumique
          elem1=face_voisins(face, 0);
          elem2=face_voisins(face, 1);
          rho_moyen = (zvf.volumes(elem1)*tab_rho[elem1]+zvf.volumes(elem2)*tab_rho[elem2]);
          rho_moyen /= zvf.volumes(elem1) + zvf.volumes(elem2);
          rho_face[face] = rho_moyen;
        }

      //   // Debug de l'affichage
      // #ifdef _AFFDEBUG
      //   {
      //     Cerr<<"RECHERCHER ELEM/FACES/SOMS"<<finl;
      //     for (int ielem=0 ; ielem<zvf.nb_elem_tot() ; ielem++) {
      //       int elem=-1;
      //       {
      //         int nb_som_elem = zvf.domaine().nb_som_elem(), isom,som;
      //         const IntTab& elem_sommets = zvf.domaine().les_elems();
      //         const DoubleTab& coord_sommets = zvf.domaine().domaine().coord_sommets();
      //         int trouve = 0;
      //         for (isom=0 ; isom<nb_som_elem ; isom++) {
      //           som = elem_sommets(ielem,isom);
      //           if (som==240 || som==226)
      //             trouve++;
      //         }
      //         if (trouve==2)
      //           elem = ielem;
      //       }
      //       if (elem!=-1) {
      //         //elem = ielem;
      //         Cerr<<"elem="<<elem<<"  coord="<<zvf.xp(elem,0)<<" "<<zvf.xp(elem,1);
      //         if (dimension==3)
      //           Cerr<<" "<<zvf.xp(elem,2);
      //         Cerr<<finl;
      //         int nb_faces_elem = zvf.domaine().nb_faces_elem(),iface,face;
      //         const IntTab& elem_faces = zvf.elem_faces();
      //         for (iface=0 ; iface<nb_faces_elem ; iface++) {
      //           face = elem_faces(elem,iface);
      //           Cerr<<"  face="<<face<<"  coord="<<zvf.xv(face,0)<<" "<<zvf.xv(face,1);
      //           if (dimension==3)
      //             Cerr<<" "<<zvf.xv(face,2);
      //           Cerr<<finl;
      //         }
      //         int nb_som_elem = zvf.domaine().nb_som_elem(), isom,som;
      //         const IntTab& elem_sommets = zvf.domaine().les_elems();
      //         const DoubleTab& coord_sommets = zvf.domaine().domaine().coord_sommets();
      //         for (isom=0 ; isom<nb_som_elem ; isom++) {
      //           som = elem_sommets(elem,isom);
      //           Cerr<<"  som="<<som<<"  coord="<<coord_sommets(som,0)<<" "<<coord_sommets(som,1);
      //           if (dimension==3)
      //             Cerr<<" "<<coord_sommets(som,2);
      //           Cerr<<finl;
      //         }
      //       }
      //     }
      //   }
      // #endif

      //   // Debug de l'affichage
      // #ifdef _AFFDEBUG
      //   _aff_donnee_P1(rho_face,"RHO_FACE");
      //   _aff_donnee_P1(la_vitesse.valeurs(),"VITESSE_FACE");
      // #endif

      // 1/dt B Un :;
      divergence.calculer(la_vitesse->valeurs(), secmem);
      secmem/=dt;

      //   // Debug de l'affichage
      // #ifdef _AFFDEBUG
      //   if (secmem.size()==zvf.nb_elem_tot())
      //     _aff_donnee_P0(secmem,"SECMEM apres divV");
      //   else
      //     _aff_donnee_P1Bulle(secmem,"SECMEM apres divV");
      // #endif

      vpoint=0.;
      //   if (dimension==3)
      //     Cerr<<"terme_diffusif.ajouter : attention, partie FT de la diffusion non implementee en 3D !!"<<finl;
      // Ajout de Div(tau_d), terme diffusif
      DoubleTrav diff(vpoint);
      terme_diffusif.ajouter(la_vitesse->valeurs(), diff);

      if(diff_boussi_==1)
        {
          // on multiplie par rho :
          if (diff.line_size()==1)
            {
              for (face=0 ; face<nbfaces ; face++)
                {
                  diff(face,0) *= rho_face[face];
                  diff(face,0) /= rho0_;
                }
            }
          else
            {
              double drho;
              for (face=0 ; face<nbfaces ; face++)
                {
                  drho = rho_face[face];
                  for (k=0 ; k<dimension ; k++)
                    {
                      diff(face,k) *= drho;
                      diff(face,k) /= rho0_;
                    }
                }
            }
          // ainsi on a : rho/rho0 div(tau_d), et quand on va diviser par rho ensuite, ne restera que div(tau_d)/rho0
        }

      vpoint+=diff;

      //   // Debug de l'affichage
      // #ifdef _AFFDEBUG
      //   _aff_donnee_P1(vitesse_interpolee_,"VITESSE_INTERP");
      //   _aff_donnee_P1(vpoint,"VPOINT apres diffusion");
      // #endif

      // Modif pour convection :
      DoubleTrav conv(vpoint);
      conv = 0.;
      terme_convectif.ajouter(la_vitesse->valeurs(), conv);

      // on multiplie par rho :
      if (conv.line_size()==1)
        {
          for (face=0 ; face<nbfaces ; face++)
            {
              conv(face,0) *= rho_face[face];
            }
        }
      else
        {
          double drho;
          for (face=0 ; face<nbfaces ; face++)
            {
              drho = rho_face[face];
              for (k=0 ; k<dimension ; k++)
                {
                  conv(face,k) *= drho;
                }
            }
        }
      vpoint+=conv;
      // Ajout de rho u Grad(u), terme convectif

      //   // Debug de l'affichage
      // #ifdef _AFFDEBUG
      //   _aff_donnee_P1(conv,"CONVECTION");
      //   _aff_donnee_P1(vpoint,"VPOINT apres convection");
      // #endif

      // on divise vpoint par rho :
      if (vpoint.line_size()==1)
        {
          for (face=0 ; face<nbfaces ; face++)
            {
              vpoint(face,0) /= rho_face[face];
            }
        }
      else
        {
          double drho;
          for (face=0 ; face<nbfaces ; face++)
            {
              drho = rho_face[face];
              for (k=0 ; k<dimension ; k++)
                {
                  vpoint(face,k) /= drho;
                }
            }
        }

      les_sources.ajouter(vpoint);
      vpoint.echange_espace_virtuel();

      //   // Debug de l'affichage
      // #ifdef _AFFDEBUG
      //   Cerr<<"SOURCES"<<les_sources<<finl;
      //   _aff_donnee_P1(vpoint,"VPOINT apres sources");
      // #endif

      solveur_masse->appliquer(vpoint);
      if (calculate_time_derivative())
        derivee_en_temps().valeurs()=vpoint;
      schema_temps().modifier_second_membre((*this),vpoint);

      // Appliquer le solveur masse => diviser par le volume
      // A la base TRUST gere des variables * volume. Donc ici on retombe sur les VRAIES variables.

      //   // Debug de l'affichage
      // #ifdef _AFFDEBUG
      //   _aff_donnee_P1(vpoint,"VPOINT apres solveur_masse");
      // #endif


      //   // Debug de l'affichage
      // #ifdef _AFFDEBUG
      //   _aff_donnee_P1(vpoint,"VPOINT apres division par Rho");
      // #endif
      // diviser par rho :
      //for (face=0 ; face<nbfaces ; face++) {
      //vpoint[face] /= rho_face[face];
      //}

      //   // Debug de l'affichage
      // #ifdef _AFFDEBUG
      //   if (secmem.size()==zvf.nb_elem_tot())
      //     _aff_donnee_P0(secmem,"SECMEM avant divVp");
      //   else
      //     _aff_donnee_P1Bulle(secmem,"SECMEM avant divVp");
      // #endif

      divergence.ajouter(vpoint, secmem);
      // On calcule Div(vpoint)/dt i-e Div(u_n+1)/dt et on l'ajoute a secmem qui comprend deja Div(u_n)/dt

      //   // Debug de l'affichage
      // #ifdef _AFFDEBUG
      //   if (secmem.size()==zvf.nb_elem_tot())
      //     _aff_donnee_P0(secmem,"SECMEM apres divVp");
      //   else
      //     _aff_donnee_P1Bulle(secmem,"SECMEM apres divVp");
      // #endif

      secmem *= -1;
      secmem.echange_espace_virtuel();

      // Resolution en pression :
      //Assembleur_P_VDF_FT& assembleur_pression=ref_cast(Assembleur_P_VDF_FT, assembleur_pression_.valeur());
      //assembleur_pression.modifier_matrice_pression(matrice_pression_);

      const Domaine_VF& domaine_VF = ref_cast(Domaine_VF,le_dom_dis.valeur());
      const DoubleVect& volumes = domaine_VF.volumes();
      int el0,el1;
      double vol0,vol1;

      // Calcul de rho_face_P aux faces pour l'assemblage avec l'assembleur pression PF (cf AssembleurPVDF de BM)
      //Champ_Fonc_Face rho_face_P;
      Champ_Don rho_face_P;
      rho_face_P.typer("Champ_Fonc_Face");
      rho_face_P->valeurs()=inconnue().valeurs();

      for (int fac=0; fac<nbfaces_bord; fac++)
        {
          el0=face_voisins(fac,0);
          if(el0!=-1)
            {
              rho_face_P->valeurs()(fac) = tab_rho(el0);
            }
          else
            {
              el1=face_voisins(fac,1);
              rho_face_P->valeurs()(fac) = tab_rho(el1);
            }
        }

      for (int fac=nbfaces_bord; fac<nbfaces; fac++)
        {
          el0=face_voisins(fac,0);
          el1=face_voisins(fac,1);
          vol0=volumes(el0);
          vol1=volumes(el1);
          rho_face_P->valeurs()(fac)=(vol0*tab_rho(el0)+vol1*tab_rho(el1))/(vol0+vol1);
        }

      assembleur_pression_->assembler_rho_variable(matrice_pression_, rho_face_P.valeur());
      assembleur_pression_->modifier_secmem(secmem);

      solveur_pression_->reinit();
      secmem.echange_espace_virtuel();
      solveur_pression_.resoudre_systeme(matrice_pression_.valeur(),secmem,la_pression->valeurs());
      assembleur_pression_->modifier_solution(la_pression->valeurs());

      la_pression->valeurs().echange_espace_virtuel();

      //   // Debug de l'affichage
      // #ifdef _AFFDEBUG
      //   if (secmem.size()==zvf.nb_elem_tot()) {
      //     _aff_donnee_P0(secmem,"SECMEM Pression");
      //     _aff_donnee_P0(la_pression.valeurs(),"PRESIONN");
      //   } else {
      //     _aff_donnee_P1Bulle(secmem,"SECMEM Pression");
      //     _aff_donnee_P1Bulle(la_pression.valeurs(),"PRESIONN");
      //   }
      // #endif

      // M-1 Bt P
      gradient.calculer(la_pression->valeurs(), gradP);
      gradP.echange_espace_virtuel();

      //   // Debug de l'affichage
      // #ifdef _AFFDEBUG
      //   _aff_donnee_P1(gradP,"GRADP avant SolvMasse");
      // #endif

      solveur_masse->appliquer(gradP);
      // Idem qu'avant

      //   // Debug de l'affichage
      // #ifdef _AFFDEBUG
      //   _aff_donnee_P1(gradP,"GRADP apres SolvMasse");
      // #endif

      // Pression a diviser par rho :
      if (gradP.line_size()==1)
        {
          for (face=0 ; face<nbfaces ; face++)
            {
              gradP(face,0) /= rho_face[face];
            }
        }
      else
        {
          double drho;
          for (face=0 ; face<nbfaces ; face++)
            {
              drho = rho_face[face];
              for (k=0 ; k<dimension ; k++)
                {
                  gradP(face,k) /= drho;
                }
            }
        }
      //   Cerr << "************ TEST d'evolution de rho ************ : rho(100) = " << rho_face[100] << finl;

      //   // Debug de l'affichage
      // #ifdef _AFFDEBUG
      //   _aff_donnee_P1(gradP,"GRADP apres /rho");
      // #endif

      // Correction en pression
      vpoint -= gradP;
      vpoint.echange_espace_virtuel();

      //   // Debug de l'affichage
      // #ifdef _AFFDEBUG
      //   _aff_donnee_P1(vpoint,"VPOINT apres -gradP");
      // #endif

    }

  else
    {
      Navier_Stokes_std::derivee_en_temps_inco(vpoint);
    }

  return vpoint;

}
