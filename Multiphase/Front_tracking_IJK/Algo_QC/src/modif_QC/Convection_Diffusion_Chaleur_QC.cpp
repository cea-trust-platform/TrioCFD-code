/****************************************************************************
* Copyright (c) 2015, CEA
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
// File:        Convection_Diffusion_Chaleur_QC.cpp
// Directory:   $TRIO_U_ROOT/src/ThHyd/Quasi_Compressible
// Version:     /main/52
//
//////////////////////////////////////////////////////////////////////////////

#include <Convection_Diffusion_Chaleur_QC.h>
#include <Probleme_base.h>
#include <Discret_Thyd.h>
#include <Domaine.h>
#include <Avanc.h>
#include <Debog.h>
#include <Frontiere_dis_base.h>
#include <EcritureLectureSpecial.h>
#include <Champ_Uniforme.h>
#include <Matrice_Morse.h>
#include <Navier_Stokes_std.h>
#include <DoubleTrav.h>
#include <Neumann_sortie_libre.h>
#include <Op_Conv_negligeable.h>
#include <Loi_Etat_GP.h>
#include <Dirichlet.h>
#include <EChaine.h>
#include <Navier_Stokes_QC_impl.h>
#include <Param.h>

#define old_forme

Implemente_instanciable_sans_constructeur(Convection_Diffusion_Chaleur_QC,"Convection_Diffusion_Chaleur_QC",Convection_Diffusion_std);

Convection_Diffusion_Chaleur_QC::Convection_Diffusion_Chaleur_QC():mode_convection_(0)
{
}

// Description:
//    Simple appel a: Convection_Diffusion_std::printOn(Sortie&)
// Precondition:
// Parametre: Sortie& is
//    Signification: un flot de sortie
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Sortie&
//    Signification: le flot de sortie modifie
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
Sortie& Convection_Diffusion_Chaleur_QC::printOn(Sortie& is) const
{
  return Convection_Diffusion_std::printOn(is);
}


// Description:
//    Verifie si l'equation a une inconnue et un fluide associe
//    et appelle Convection_Diffusion_std::readOn(Entree&).
// Precondition: l'objet a une inconnue associee
// Precondition: l'objet a un fluide associe
// Parametre: Entree& is
//    Signification: un flot d'entree
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: Entree& is
//    Signification: le flot d'entree modifie
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
Entree& Convection_Diffusion_Chaleur_QC::readOn(Entree& is)
{
  assert(l_inco_ch.non_nul());
  assert(le_fluide.non_nul());
  Convection_Diffusion_std::readOn(is);
  Nom unite;
  if (dimension+bidim_axi==2) unite="[W/m]";
  else unite="[W]";
  terme_convectif.set_fichier("Convection_chaleur");
  terme_convectif.set_description((Nom)"Convective heat transfer rate=Integral(-rho*cp*T*u*ndS) "+unite);
  terme_diffusif.set_fichier("Diffusion_chaleur");
  terme_diffusif.set_description((Nom)"Conduction heat transfer rate=Integral(lambda*grad(T)*ndS) "+unite);

  //On modifie le nom ici pour que le champ puisse etre reconnu si une sonde d enthalpie est demandee
  if (le_fluide->type_fluide()=="Gaz_Reel")
    l_inco_ch->nommer("enthalpie");

  return is;
}

void Convection_Diffusion_Chaleur_QC::set_param(Param& param)
{
  Convection_Diffusion_std::set_param(param);
  param.ajouter_non_std("mode_calcul_convection",(this));
}

int Convection_Diffusion_Chaleur_QC::lire_motcle_non_standard(const Motcle& motlu, Entree& is)
{

  Motcles les_mots(3);
  {
    les_mots[0] = "diffusion";
    les_mots[1] = "convection";
    les_mots[2] = "mode_calcul_convection";
  }
  {

    int rang=les_mots.search(motlu);
    switch(rang)
      {
      case 0:
        {
          Cerr << "Lecture de l'operateur de diffusion " << finl;
          Cerr << "et typage: ";
          //associe la conductivite dans la diffusivite
          terme_diffusif.associer_diffusivite(le_fluide->conductivite());
          is >> terme_diffusif;
          terme_diffusif.associer_diffusivite_pour_pas_de_temps(le_fluide->diffusivite());
          Cerr<<"diffusivite="<<terme_diffusif.diffusivite().le_nom()<<finl;
          break;
        }
      case 1 :
        {
          const Probleme_base& pb = probleme();
          const Navier_Stokes_std& eqn_hydr = ref_cast(Navier_Stokes_std,pb.equation(0));
          const Champ_Inc& vitessetransportante =( mode_convection_==2?eqn_hydr.rho_la_vitesse():eqn_hydr.inconnue());
          associer_vitesse(vitessetransportante);
          Cerr << "Lecture de l'operateur de convection" << finl;
          Cerr << "Association de la vitesse de convection de Crank Nicholson" << finl;
          terme_convectif.associer_vitesse(vitessetransportante);

          // Il ne faut pas rajouter le terme en dpth/dt ....
          if (0)
            {
              //l'equation de la chaleur en quasi compressible contient un terme source dP/dt
              Cerr << "Creation du terme source de l'equation de la chaleur :"<< finl;
              Source t;
              Source& so=les_sources.add(t);
              Nom type_so = "Source_Quasi_Compressible_Chaleur_";
              Nom disc = discretisation().que_suis_je();
              if (disc=="VEFPreP1B")
                disc = "VEF";
              type_so+=disc;
              so.typer_direct(type_so);
              so->associer_eqn(*this);
              Cerr<<so->que_suis_je()<<finl;
              Cerr << "et typage de l'operateur de convection: ";
            }

          is >> terme_convectif;
          break;
        }
      case 2:
        {
          if (terme_convectif.non_nul())
            {
              Cerr<<" l'option "<<motlu <<" doit etre mise avant la lecture de la convection. modifier votre jeu de donnees"<<finl;
              exit();
            }
          Motcles modes(3);
          modes[0]="ancien";
          modes[1]="divuT_moins_Tdivu";
          modes[2]="divrhouT_moins_Tdivrhou";
          Motcle mot2;
          is>>mot2;
          mode_convection_=modes.search(mot2);
          if (mode_convection_==-1)
            {
              Cerr<<" On ne comprend pas "<<mot2<<" apres mode_convection "<<finl;
              Cerr<<"on comprend: "<<modes<<finl;
              exit();
            }
          break;
        }
      default :
        {
          return Convection_Diffusion_std::lire_motcle_non_standard(motlu,is);

        }
      }
  }
  return 1;
}

// Description:
//    Associe un milieu physique a l'equation,
//    le milieu est en fait caste en Fluide_Incompressible ou en Fluide_Ostwald.
// Precondition:
// Parametre: Milieu_base& un_milieu
//    Signification:
//    Valeurs par defaut:
//    Contraintes: reference constante
//                 doit pourvoir etre force au type "Fluide_Incompressible"
//    Acces: entree
// Retour:
//    Signification:
//    Contraintes:
// Exception: les proprietes physiques du fluide ne sont pas toutes specifiees
// Effets de bord:
// Postcondition:
void Convection_Diffusion_Chaleur_QC::associer_milieu_base(const Milieu_base& un_milieu)
{
  const Fluide_Quasi_Compressible& un_fluideQC = ref_cast(Fluide_Quasi_Compressible,un_milieu);
  associer_fluide(un_fluideQC);
}

const Champ_Don& Convection_Diffusion_Chaleur_QC::diffusivite_pour_transport()
{
  return milieu().conductivite();
}

const Champ_base& Convection_Diffusion_Chaleur_QC::diffusivite_pour_pas_de_temps()
{
  return milieu().diffusivite();
}

const Champ_base& Convection_Diffusion_Chaleur_QC::vitesse_pour_transport()
{
  const Probleme_base& pb = probleme();
  const Navier_Stokes_std& eqn_hydr = ref_cast(Navier_Stokes_std,pb.equation(0));
  const Champ_Inc& vitessetransportante =( mode_convection_==2?eqn_hydr.rho_la_vitesse():eqn_hydr.inconnue());
  return vitessetransportante;
}


// Description:
//    Discretise l'equation.
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
// Postcondition: l'equation est discretisee
void Convection_Diffusion_Chaleur_QC::discretiser()
{
  const Discret_Thyd& dis=ref_cast(Discret_Thyd, discretisation());
  Cerr << "Energy equation discretization " << finl;
  dis.temperature(schema_temps(), zone_dis(), l_inco_ch);
  //Cerr<<finl<<"On va effectuer proprietes_physiques_fluide_Incompressible"<<finl;
  //dis.proprietes_physiques_fluide_incompressible(zone_dis(), fluide(), l_inco_ch);
  champs_compris_.ajoute_champ(l_inco_ch);
  Equation_base::discretiser();
  Cerr << "Convection_Diffusion_Chaleur_QC::discretiser() ok" << finl;
}


// Description:
//    Renvoie le milieu physique de l'equation.
//    (un Fluide_Incompressible upcaste en Milieu_base)
//    (version const)
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Milieu_base&
//    Signification: le Fluide_Incompressible upcaste en Milieu_base
//    Contraintes: reference constante
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
const Milieu_base& Convection_Diffusion_Chaleur_QC::milieu() const
{
  return fluide();
}


// Description:
//    Renvoie le milieu physique de l'equation.
//    (un Fluide_Incompressible upcaste en Milieu_base)
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Milieu_base&
//    Signification: le Fluide_Incompressible upcaste en Milieu_base
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition:
Milieu_base& Convection_Diffusion_Chaleur_QC::milieu()
{
  return fluide();
}


// Description:
//    Renvoie le fluide incompressible associe a l'equation.
//    (version const)
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Fluide_Incompressible&
//    Signification: le fluide incompressible associe a l'equation
//    Contraintes: reference constante
// Exception: pas de fluide associe a l'eqaution
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
const Fluide_Incompressible& Convection_Diffusion_Chaleur_QC::fluide() const
{
  if (!le_fluide.non_nul())
    {
      Cerr << "You forgot to associate the fluid to the problem named " << probleme().le_nom() << finl;
      Process::exit();
    }
  return le_fluide.valeur();
}


// Description:
//    Renvoie le fluide incompressible associe a l'equation.
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Fluide_Incompressible&
//    Signification: le fluide incompressible associe a l'equation
//    Contraintes:
// Exception: pas de fluide associe a l'eqaution
// Effets de bord:
// Postcondition:
Fluide_Incompressible& Convection_Diffusion_Chaleur_QC::fluide()
{
  assert(le_fluide.non_nul());
  return le_fluide.valeur();
}



void Convection_Diffusion_Chaleur_QC::calculer_div_u_ou_div_rhou(DoubleTab& derivee2) const
{
  // autant tester tout de suite
  if (sub_type(Op_Conv_negligeable,operateur(1).l_op_base()))
    {
      derivee2=0;
      return;
    }
  if (mode_convection_!=0)
    {
      const DoubleTab& T=inconnue().valeurs();


      DoubleTrav unite(T);
      unite=1;
      ref_cast_non_const(Operateur_base,operateur(1).l_op_base()).associer_zone_cl_dis(zcl_modif_.valeur());

      operateur(1).ajouter(unite,derivee2);
      ref_cast_non_const(Operateur_base,operateur(1).l_op_base()).associer_zone_cl_dis(zone_Cl_dis().valeur());
    }
  else
    {
      Cerr<<"YB - 'C_D_Chaleur_QC' - le mode de convection n'est pas bon, on ne renvoie pas div(rhou)"<<finl;
      abort();
      /*
        const Navier_Stokes_std& eqn_hydr = ref_cast(Navier_Stokes_std,probleme().equation(0));
        const DoubleTab& pression = eqn_hydr.pression().valeurs();

        const Champ_Inc& vitesse = eqn_hydr.inconnue();
        // if (pression.dimension(0)!=derivee2.dimension(0))
        if (vitesse.valeurs().nb_dim()==2)
        {
          DoubleTab secmem1(pression);
          DoubleTab secmem2(vitesse.valeurs());
          eqn_hydr.operateur_divergence().calculer(vitesse.valeurs(),secmem1);
          secmem1.echange_espace_virtuel();
          ref_cast_non_const(Fluide_Quasi_Compressible,le_fluide.valeur()).divu_discvit(secmem1,secmem2);
          int i, nsom = secmem2.dimension_tot(0);
          for (i=0 ; i<nsom ; i++)
            derivee2(i)=secmem2(i,0);
        }
        else
        {
        eqn_hydr.operateur_divergence().calculer(vitesse.valeurs(),derivee2);
        }
        derivee2*=-1;
      */
    }
}


// Description:
//     Renvoie la derivee en temps de l'inconnue de l'equation.
//     Le calcul est le suivant:
//         d(inconnue)/dt = M^{-1} * (sources - somme(Op_{i}(inconnue))) / rho
// Precondition:
// Parametre: DoubleTab& derivee
//    Signification: le tableau des valeurs de la derivee en temps du champ inconnu
//    Valeurs par defaut:
//    Contraintes: ce parametre est remis a zero des l'entree de la methode
//    Acces: sortie
// Retour: DoubleTab&
//    Signification: le tableau des valeurs de la derivee en temps du champ inconnu
//    Contraintes:
// Exception:
// Effets de bord: des communications (si version parallele) sont generees pas cet appel
// Postcondition:
DoubleTab& Convection_Diffusion_Chaleur_QC::derivee_en_temps_inco(DoubleTab& derivee)
{
  derivee=0;
  derivee_en_temps_inco_sans_solveur_masse(derivee);
  if (!schema_temps().diffusion_implicite())
    {
      solveur_masse.appliquer(derivee);
      derivee.echange_espace_virtuel();
    }
  return derivee;
}
DoubleTab& Convection_Diffusion_Chaleur_QC::derivee_en_temps_inco_sans_solveur_masse(DoubleTab& derivee)
{
  {
    // test
    const Navier_Stokes_std& eqns=ref_cast(Navier_Stokes_std,probleme().equation(0));
    DoubleTab uu;
    if (0)
      {
        uu=(probleme().get_champ("pression").valeurs());

        eqns.operateur_divergence().calculer(terme_convectif.vitesse().valeurs(),uu);

        //  Cerr<<" iii2 "<<mp_max_abs_vect(uu)<<finl;
      }
    if (mode_convection_==2)
      {
        DoubleTab& rhou=ref_cast_non_const(DoubleTab,probleme().get_champ("rho_u").valeurs());
        //rho_vitesse_impl(inconnue().valeurs(),eqns.inconnue().valeurs(),rhou);
        Fluide_Quasi_Compressible& fluide_QC=ref_cast(Fluide_Quasi_Compressible,le_fluide.valeur());

        rhou=fluide_QC.rho_face_np1();
        /*
          DoubleTab test(rhou);
          test-=fluide_QC.rho_face_np1();
          double toto=mp_max_abs_vect(test);
          if(toto!=0)
          {
          Cerr<<" GF je ne sais pas si il faut mettre rhon ou rhonp1 "<<toto<<finl;
          exit();
          }
        */
        rho_vitesse_impl(rhou,eqns.inconnue().valeurs(),rhou);
        //Cerr<<mpf_min_abs_vect(inconnue().valeurs())<<"max_value"<<mp_max_abs_vect(inconnue().valeurs())<<finl;
        if (0)
          {
            eqns.operateur_divergence().calculer(terme_convectif.vitesse().valeurs(),uu);

            Cerr<<" iii3 "<<mp_max_abs_vect(uu)<<finl;
          }


      }

  }
  {
    //        Fluide_Quasi_Compressible& fluide_QC=ref_cast(Fluide_Quasi_Compressible,le_fluide.valeur());  modif YB 28/08/09
    const DoubleTab& tab_rho = le_fluide->masse_volumique().valeurs();
    //Cerr<<"YB - 'C_D_Chaleur_QC' - masse volumique = "<<min(tab_rho)<<" - "<<max(tab_rho)<<finl;
    int i,n=tab_rho.dimension(0);
    for ( i=0 ; i<n ; i++ )
      {
        derivee(i)=tab_rho(i);
        derivee(i)=0;
      }
    if (mode_convection_==0)
      {
        la_zone_Cl_dis.les_conditions_limites().set_modifier_val_imp(1);
        operateur(1).ajouter(derivee);
        la_zone_Cl_dis.les_conditions_limites().set_modifier_val_imp(0);
      }
    else
      calculer_div_u_ou_div_rhou(derivee);

    //Cerr<<"YB - 'C_D_Chaleur_QC' - d(rhou)/dt = "<<min(derivee)<<" - "<<max(derivee)<<finl;
    /* modif YB 28/08/09
       const Loi_Etat_GP& la_loi = ref_cast(Loi_Etat_GP,fluide_QC.loi_etat().valeur());
       const double R = la_loi.R();
       const double Pth = fluide_QC.pression_th();
       for ( i=0 ; i<n ; i++ )
       derivee(i) *= -Pth/(R*tab_rho(i)*tab_rho(i));
       Cerr<<"YB - 'C_D_Chaleur_QC' - d(T)/dt = "<<min(derivee)<<" - "<<max(derivee)<<finl;
       Cerr<<"YB - 'C_D_Chaleur_QC' - R = "<<R<<finl;
       Cerr<<"YB - 'C_D_Chaleur_QC' - Pth = "<<Pth<<finl;
    */
    /*
      const Schema_Temps_base& sch=schema_temps();
      int diffusion_implicite=sch.diffusion_implicite();
      const DoubleTab& tab_rho = le_fluide->masse_volumique().valeurs();
      int n = tab_rho.dimension(0);
      la_zone_Cl_dis.les_conditions_limites().set_modifier_val_imp(0);
      // Calcul de derivee=Div(lambda*gradient(temperature))
      if (le_fluide->type_fluide()=="Gaz_Parfait")
      {
      if (!diffusion_implicite)
      operateur(0).ajouter(derivee);
      }
      else
      {
      operateur(0).ajouter(le_fluide->temperature(),derivee);
      assert(diffusion_implicite==0);
      }

      la_zone_Cl_dis.les_conditions_limites().set_modifier_val_imp(1);
      derivee.echange_espace_virtuel();
      Debog::verifier("Convection_Diffusion_Chaleur_QC::derivee_en_temps_inco_sans_solveur_masse derivee apres convection",derivee);
      // On ajoute le terme source a derivee
      les_sources.ajouter(derivee);


      // on stocke rho_cp
      const Champ_base& champ_rho_cp=get_champ("rho_cp");
      DoubleTab& rho_cp=ref_cast_non_const(DoubleTab,champ_rho_cp.valeurs());
      // On divise derivee par rho*Cp dans le cas d'un gaz parfait, par rho dans le cas d'un gaz reel
      if (le_fluide->type_fluide()=="Gaz_Parfait") {
    {
      double Cp;
      if (sub_type(Champ_Uniforme,le_fluide->capacite_calorifique().valeur()))
      {
      ////double Cp = le_fluide->capacite_calorifique().valeurs()(0,0);
      Cp = le_fluide->capacite_calorifique().valeurs()(0,0);
      for (int som=0 ; som<n ; som++) {
            {
      double rhocp=tab_rho(som)*Cp;

      derivee(som) /= rhocp;
      rho_cp(som)=rhocp;
      }
      }
      else {
        {
          rho_cp= le_fluide->capacite_calorifique().valeurs();
          for (int som=0 ; som<n ; som++)
            {
      rho_cp(som) *= tab_rho(som);
      derivee(som) /= rho_cp(som);
            }
      }

      }
      else {
      for (int som=0 ; som<n ; som++) {
      for (int som=0 ; som<n ; som++)
        {
          derivee(som) /= tab_rho(som);
      }
      }
      derivee.echange_espace_virtuel();
      Debog::verifier("Convection_Diffusion_Chaleur_QC::derivee_en_temps_inco_sans_solveur_masse derivee",derivee);

      //ajout de la convection (avec la vitesse de Crank Nicholson)
      // derivee=(Div(lambda*grad(T))+S)/rho/Cp-ugrad(T))
      //  if (!diffusion_implicite)

      {
      // ajouter test sur pas de temps convectif....

      DoubleTrav derivee2(derivee);


      // ajout de -T div(rho *u * 1)
      const DoubleTab& T=inconnue().valeurs();
      calculer_div_u_ou_div_rhou(derivee2);

      for (int som=0 ; som<n ; som++) {
      {
        derivee2(som) *= (-T(som));
      }
      // ajout de div(rho u T)
      operateur(1).ajouter(derivee2);
      // division par rho
      if (mode_convection_==2)
      for (int som=0 ; som<n ; som++)
        {
      derivee2(som) /= tab_rho(som);
        }

      derivee+=derivee2;
      }
      // else
      if (diffusion_implicite)
      {
      //  abort();
      DoubleTrav secmem(derivee);
      //derivee_en_temps_conv(secmem,T);
      secmem=derivee;
      solv_masse().appliquer(secmem);
      derivee=(Tfutur);
      solveur_masse->set_name_of_coefficient_temporel("rho_cp");
      double seuil = le_schema_en_temps->seuil_diffusion_implicite();
      int niter = Equation_base::Gradient_conjugue_diff_impl( secmem, derivee,seuil ) ;

      solveur_masse->set_name_of_coefficient_temporel("no_coeff");
      Cout << "Viscose Fluxes Implicite: Gradient conjugue convergence in " << niter << " iterationen " << finl ;

      double dt=le_schema_en_temps->pas_de_temps();
      //derivee=delta_T;
      derivee/=dt;
    */
  }
  return derivee;
}

// derivee_en_temps_inco non std donc assembler non plus
void Convection_Diffusion_Chaleur_QC::assembler( Matrice_Morse& matrice,const DoubleTab& inco, DoubleTab& resu)
{

  int test_op=0;
  {
    char* theValue = getenv("TRIO_U_TEST_OPERATEUR_IMPLICITE");
    if (theValue != NULL) test_op=2;
  }
  {
    char* theValue = getenv("TRIO_U_TEST_OPERATEUR_IMPLICITE_BLOQUANT");
    if (theValue != NULL) test_op=1;
  }


  const DoubleTab& tab_rho = le_fluide->masse_volumique().valeurs();
  int n = tab_rho.dimension(0);
  //ajout diffusion (avec la conductivite et la temperature)
  operateur(0).l_op_base().contribuer_a_avec(inco, matrice );

  les_sources.contribuer_a_avec(inco,matrice);


  const IntVect& tab1= matrice.get_tab1();

  DoubleVect& coeff=matrice.get_set_coeff();
  DoubleVect coeff_diffusif(coeff);
  // on calcule les coefficients de l'op de convection on obtient les coeff de div(p*u*T) , il faudrait multiplier par cp puis divisier par rho cp on le fera d'un coup...
  coeff=0;

  operateur(1).l_op_base().contribuer_a_avec(inco, matrice );


  // on calcule div(rhou)
  DoubleTrav derivee2(resu);

  calculer_div_u_ou_div_rhou(derivee2);


  int som;
  if (le_fluide->type_fluide()=="Gaz_Parfait")
    {
      double Cp=-5;
      int is_cp_unif=sub_type(Champ_Uniforme,le_fluide->capacite_calorifique().valeur());
      const DoubleTab& tab_cp =le_fluide->capacite_calorifique().valeurs();
      if (is_cp_unif)
        Cp=tab_cp(0,0);
      for (som=0 ; som<n ; som++)
        {
          if (!is_cp_unif)
            Cp=tab_cp(som);
          double inv_rho=1./tab_rho(som);
          if (mode_convection_!=2)
            inv_rho=1;
          double rapport=1./(tab_rho(som)*Cp);
          // il faut multiplier toute la ligne de la matrice par rapport
          //derivee(som) /= tab_rho(som)*Cp;
          for (int k=tab1(som)-1; k<tab1(som+1)-1; k++)
            coeff(k)= (coeff(k)*inv_rho+coeff_diffusif(k)*rapport);
          // ajout de Tdiv(rhou )/rho

          matrice(som,som)+=derivee2(som)*inv_rho;
        }


    }

  else
    {
      Cerr<<"The implicit algorithm is available only for perfect gas."<<finl;
      exit();
    }

  //  operateur(1).l_op_base().contribuer_a_avec(inco, matrice );
  // on a la matrice approchee on recalcule le residu;
  resu=0;
  derivee_en_temps_inco_sans_solveur_masse(resu);
  //double max_resu=max_abs(resu);Cerr<<__LINE__<< "resu "<<max_resu<<finl;
  matrice.ajouter_multvect(inco,resu);
  //max_resu=max_abs(resu);Cerr<<__LINE__<< "resu "<<max_resu<<finl;
  //if (max_resu>1e5) exit();

  if (test_op)
    {
      DoubleTrav diff(resu);
      DoubleTrav conv(resu);
      operateur(0).l_op_base().contribuer_au_second_membre(diff);
      operateur(1).l_op_base().contribuer_au_second_membre(conv);
      les_sources.ajouter(diff);
      double Cp=-5;
      int is_cp_unif=sub_type(Champ_Uniforme,le_fluide->capacite_calorifique().valeur());
      const DoubleTab& tab_cp =le_fluide->capacite_calorifique().valeurs();
      if (is_cp_unif)
        Cp=tab_cp(0,0);
      for (som=0 ; som<n ; som++)
        {
          if (!is_cp_unif)
            Cp=tab_cp(som);
          double inv_rho=1./tab_rho(som);
          if (mode_convection_!=2)
            inv_rho=1;
          double rapport=1./(tab_rho(som)*Cp);
          diff(som)=resu(som)-conv(som)*inv_rho-diff(som)*rapport;
        }
      solv_masse().appliquer(diff);
      double err=mp_max_abs_vect(diff);
      Cerr<<que_suis_je()<<" erreur assemblage "<<err<<finl;;

      if (err>1e-5)
        {
          {
            DoubleVect& diff_=diff;
            Cerr<<" size "<< diff_.size()<<finl;
            for (int i=0; i<diff_.size(); i++)
              if (std::fabs(diff_(i))>1e-5)
                {
                  Cerr<<i << " "<< diff_(i)<< " "<<finl;
                }
          }
          if (test_op==1)
            {
              Cerr<<" pb max case "<<imin_array(diff)<<" ou " <<imax_array(diff)<<finl;
              exit();
            }
        }

    }
}

int Convection_Diffusion_Chaleur_QC::sauvegarder(Sortie& os) const
{
  int bytes=0;
  //  Equation_base::sauvegarder(os);
// bytes += Equation_base::sauvegarder(os);
  // en mode ecriture special seul le maitre ecrit
  int a_faire,special;
  /*
  //Tentative pour sauvegarder la temperature et non la masse volumique
  //(devenue l'inconnue de l'equation convetion_diffusion_chaleur)
  //AT 14/04/2010
  const Nom temp("temperature_qc");
  l_inco_ch->nommer(temp);
  */
  EcritureLectureSpecial::is_ecriture_special(special,a_faire);
  const Champ_Fonc_base& ch=ref_cast(Champ_Fonc_base,get_champ("temperature_qc"));
  /*
    if (a_faire)
    {
    Nom mon_ident(ch.le_nom());
    mon_ident += ch.que_suis_je();
    mon_ident += ch.equation().probleme().domaine().le_nom();
    mon_ident += Nom(ch.temps(),"%e");
    os << mon_ident << finl;
    os << ch.que_suis_je() << finl;
    os << ch.temps() << finl;
    }
    if (special)
    EcritureLectureSpecial::ecriture_special(ch,os);
    else
    ch.valeurs().ecrit(os);
    if (a_faire)
    {
    // fich << flush ; Ne flushe pas en binaire !
    os.flush();
    }
    //  if (Process::je_suis_maitre())
    //  Cerr << "Sauvegarde du champ " << nom_ << " effectuee au temps : " << Nom(temps_,"%e") << finl;

  */
  bytes+=ch.sauvegarder(os);

  if (a_faire)
    {
      // CORRECTION BOGUE 26/10/2000
      Nom ident_Pth("pression_thermo");
      // FIN CORRECTION BOGUE 26/10/2000
      ident_Pth += probleme().domaine().le_nom();
      double temps = inconnue().temps();
      ident_Pth += Nom(temps,"%e");
      os << ident_Pth<<finl;
      os << "constante"<<finl;
      os << le_fluide->pression_th();
      os << flush ;
      Cerr << "Saving thermodynamic pressure at time : " <<  Nom(temps,"%e") << finl;
    }
  return bytes;
}

// Description:
//     Effectue une reprise a partir d'un flot d'entree.
//     Appelle Equation_base::reprendre()
//     et reprend l'inconnue de la chaleur et la pression thermodynamique
// Precondition:
// Parametre: Entree& is
//    Signification: un flot d'entree
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: int
//    Signification: renvoie toujours 1
//    Contraintes:
// Exception: la reprise a echoue, identificateur de la pression non trouve
// Effets de bord:
// Postcondition:
int Convection_Diffusion_Chaleur_QC::reprendre(Entree& is)
{
  if (le_fluide->type_fluide() != "Gaz_Parfait")
    {
      l_inco_ch->nommer("enthalpie");
    }

  //  Equation_base::reprendre(is);
  double temps = schema_temps().temps_courant();
  Champ_Fonc_base& ch=ref_cast_non_const(Champ_Fonc_base,get_champ("temperature_qc"));
  Nom ident_lu,type_inco;
  is >> ident_lu;
  is >> type_inco;
  Cerr<<"ooooo" <<ident_lu<<" "<<type_inco<<finl;
  ch.reprendre(is);
  inconnue().valeurs()=ch.valeurs();

  Nom ident_Pth("pression_thermo");
  ident_Pth += probleme().domaine().le_nom();
  ident_Pth += Nom(temps,probleme().reprise_format_temps());

  avancer_fichier(is, ident_Pth);
  double pth;
  is>>pth;
  le_fluide->set_pression_th(pth);

  return 1;

}

// Description:
//    Impression des flux sur les bords sur un flot de sortie.
//    Appelle Equation_base::impr(Sortie&)
// Precondition: Sortie&
// Parametre: Sortie& os
//    Signification: un flot de sortie
//    Valeurs par defaut:
//    Contraintes:
//    Acces: entree/sortie
// Retour: int
//    Signification: code de retour propage
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
int Convection_Diffusion_Chaleur_QC::impr(Sortie& os) const
{
  return Equation_base::impr(os);
}


// Description:
//    Renvoie le nom du domaine d'application de l'equation.
//    Ici "Thermique".
// Precondition:
// Parametre:
//    Signification:
//    Valeurs par defaut:
//    Contraintes:
//    Acces:
// Retour: Motcle&
//    Signification: le nom du domaine d'application de l'equation
//    Contraintes: toujours egal a "Thermique"
// Exception:
// Effets de bord:
// Postcondition: la methode ne modifie pas l'objet
const Motcle& Convection_Diffusion_Chaleur_QC::domaine_application() const
{
  static Motcle domaine = "Thermique_H";
  if (le_fluide->type_fluide()=="Gaz_Parfait")
    {
      domaine = "Thermique";
    }
  return domaine;
}

// Description:
//    Associe un fluide incompressible a l'equation.
// Precondition:
// Parametre: Fluide_Incompressible& un_fluide
//    Signification: le milieu fluide incompressible a associer a l'equation
//    Valeurs par defaut:
//    Contraintes: reference constante
//    Acces: entree
// Retour:
//    Signification:
//    Contraintes:
// Exception:
// Effets de bord:
// Postcondition: l'equation a un milieu associe
void Convection_Diffusion_Chaleur_QC::associer_fluide(const Fluide_Incompressible& un_fluide)
{
  assert(sub_type(Fluide_Quasi_Compressible,un_fluide));
  le_fluide = ref_cast(Fluide_Quasi_Compressible,un_fluide);
}

int Convection_Diffusion_Chaleur_QC::preparer_calcul()
{
  Convection_Diffusion_std::preparer_calcul();

  // afin de dimensionner flux_fords
  DoubleTrav tmp(inconnue().valeurs());
  operateur(1).ajouter(tmp,tmp);
  if (mode_convection_==0)
    return 1;
  // remplissage de la zone cl modifiee avec 1 partout au bord...
  zcl_modif_=(zone_Cl_dis());
  Conds_lim& condlims=zcl_modif_.valeur().les_conditions_limites();
  int nb=condlims.size();
  for (int i=0; i<nb; i++)
    {
      // pour chaque condlim on recupere le champ_front et on met 1
      // meme si la cond lim est un flux (dans ce cas la convection restera
      // nullle.)

      if (sub_type(Neumann_sortie_libre,condlims[i].valeur()))
        ref_cast(Neumann_sortie_libre,condlims[i].valeur()).tab_ext()=1;
      if (sub_type(Dirichlet,condlims[i].valeur()))
        {
          const Frontiere_dis_base& frdis=condlims[i].valeur().frontiere_dis();
          EChaine toto(" Champ_front_uniforme 1 1");
          toto>> condlims[i].valeur();
          condlims[i].valeur().associer_fr_dis_base(frdis);
        }
      DoubleTab& T=condlims[i].valeur().champ_front().valeurs();
      T=1.;

    }
  zcl_modif_.les_conditions_limites().set_modifier_val_imp(0);
  return 1;
}
