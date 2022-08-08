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
// File:        Source_Con_Phase_field.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Multiphase/Phase_field/src/VDF
// Version:     /main/25
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Con_Phase_field.h>
#include <Zone_VDF.h>
#include <Zone_Cl_VDF.h>
#include <Probleme_base.h>
#include <Milieu_base.h>
#include <Navier_Stokes_phase_field.h>
#include <Dimension.h>
#include <Source_Qdm_VDF_Phase_field.h>
#include <Check_espace_virtuel.h>
#include <Convection_Diffusion_Phase_field.h>
#include <Champ_Fonc_Tabule.h>
#include <Champ_Uniforme.h>
#include <Convection_Diffusion_Concentration.h> //Mirantsoa 264902
#include <Equation_base.h>
#include <SolveurSys.h>
#include <EChaine.h>
#include <MD_Vector_composite.h>
#include <ConstDoubleTab_parts.h>

Implemente_instanciable(Source_Con_Phase_field,"Source_Con_Phase_field_VDF_P0_VDF",Source_Con_Phase_field_base);
// XD source_con_phase_field source_base source_con_phase_field 1 Keyword to define the source term of the Cahn-Hilliard equation.
// XD attr temps_d_affichage entier temps_d_affichage 0 Time during the caracteristics of the problem are shown before calculation.
// XD attr alpha floattant alpha 0 Internal capillary coefficient alfa.
// XD attr beta floattant beta 0 Parameter beta of the model.
// XD attr kappa floattant kappa 0 Mobility coefficient kappa0.
// XD attr kappa_variable bloc_kappa_variable kappa_variable 0 To define a mobility which depends on concentration C.
// XD attr moyenne_de_kappa chaine moyenne_de_kappa 0 To define how mobility kappa is calculated on faces of the mesh according to cell-centered values (chaine is arithmetique/harmonique/geometrique).
// XD attr multiplicateur_de_kappa floattant multiplicateur_de_kappa 0 To define the parameter of the mobility expression when mobility depends on C.
// XD attr couplage_NS_CH chaine couplage_NS_CH 0 Evaluating time choosen for the term source calculation into the Navier Stokes equation (chaine is mutilde(n+1/2)/mutilde(n), in order to be conservative, the first choice seems better).
// XD attr implicitation_CH chaine(into=["oui","non"]) implicitation_CH 0 To define if the Cahn-Hilliard will be solved using a implicit algorithm or not.
// XD attr gmres_non_lineaire chaine(into=["oui","non"]) gmres_non_lineaire 0 To define the algorithm to solve Cahn-Hilliard equation (oui: Newton-Krylov method, non: fixed point method).
// XD attr seuil_cv_iterations_ptfixe floattant seuil_cv_iterations_ptfixe 0 Convergence threshold (an option of the fixed point method).
// XD attr seuil_residu_ptfixe floattant seuil_residu_ptfixe 0 Threshold for the matrix inversion used in the method (an option of the fixed point method).
// XD attr seuil_residu_gmresnl floattant seuil_residu_gmresnl 0 Convergence threshold (an option of the Newton-Krylov method).
// XD attr dimension_espace_de_krylov entier dimension_espace_de_krylov 0 Vector numbers used in the method (an option of the Newton-Krylov method).
// XD attr nb_iterations_gmresnl entier nb_iterations_gmresnl 0 Maximal iteration (an option of the Newton-Krylov method).
// XD attr residu_min_gmresnl floattant residu_min_gmresnl 0 Minimal convergence threshold (an option of the Newton-Krylov method).
// XD attr residu_max_gmresnl floattant residu_max_gmresnl 0 Maximal convergence threshold (an option of the Newton-Krylov method).
// XD attr potentiel_chimique bloc_potentiel_chim potentiel_chimique 1 chemical potential function
// XD bloc_kappa_variable objet_lecture nul 0 if the parameter of the mobility, kappa, depends on C
// XD attr expr bloc_lecture expr 0 choice for kappa_variable
// XD bloc_potentiel_chim objet_lecture nul 0 if the chemical potential function is an univariate function
// XD attr expr bloc_lecture expr 0 choice for potentiel_chimique

Sortie& Source_Con_Phase_field::printOn(Sortie& s ) const { return s << que_suis_je() ; }

Entree& Source_Con_Phase_field::readOn(Entree& is )
{
  Cerr<<"Source_Con_Phase_field::readOn"<<finl;
  Motcles les_mots(23);
  les_mots[0]="nb_equation_CH";
  les_mots[1]="alpha";
  les_mots[2]="potentiel_chimique";
  les_mots[3]="beta";
  les_mots[4]="kappa";
  les_mots[5]="kappa_variable";
  les_mots[6]="moyenne_de_kappa";
  les_mots[7]="multiplicateur_de_kappa";
  les_mots[8]="couplage_NS_CH";
  les_mots[9]="implicitation_CH";
  les_mots[10]="gmres_non_lineaire";
  les_mots[11]="seuil_cv_iterations_ptfixe";
  les_mots[12]="seuil_residu_ptfixe";
  les_mots[13]="seuil_residu_gmresnl";
  les_mots[14]="dimension_espace_de_krylov";
  les_mots[15]="nb_iterations_gmresnl";
  les_mots[16]="residu_min_gmresnl";
  les_mots[17]="residu_max_gmresnl";
  les_mots[18]="Temps_d_affichage";
  les_mots[19]="kappa_auto_diffusion";
  les_mots[20]="alpha_rotation";
  les_mots[21]="min_x";
  les_mots[22]="max_x";

  Motcle motlu;
  is >> motlu;
  if (motlu!="{")
    {
      Cerr<<"Source_Con_Phase_field::readOn: We are expecting { at the beginning of Source_Con_Phase_field block instead of " << motlu <<finl;
      exit();
    }
  int cpt = 0;
  is >> motlu;
  if (motlu!="systeme_naire")
    {
      Cerr<<"Source_Con_Phase_field::readOn: We are expecting 'systeme_naire' instead of " << motlu <<finl;
      exit ();
    }
  else
    {
      Motcle temp_systeme_naire;
      is >> temp_systeme_naire;
      if (temp_systeme_naire=="non")
        {
          type_systeme_naire_=0;
          is >> motlu;
          if (motlu!="{")
            {
              Cerr<<"Source_Con_Phase_field::readOn: We are expecting { after 'non' instead of "<< motlu <<finl;
              exit();
            }
          is >>motlu;
          int cpt0=0;
          while(motlu!="}")
            {
              int rang=les_mots.search(motlu);
              switch(rang)
                {
                case 1:
                  {
                    cpt0++;
                    is >> alpha;
                    break;
                  }
                case 2:
                  {
                    cpt0++;
                    is >> motlu;
                    if (motlu!="{")
                      {
                        Cerr<<"Source_Con_Phase_field::readOn: We are expecting { after potentiel_chimique instead of "<< motlu <<finl;
                        exit();
                      }
                    Motcle temp_potentiel_chimique_;
                    is >> temp_potentiel_chimique_;
                    if (temp_potentiel_chimique_=="defaut")
                      {
                        dWdc=&Source_Con_Phase_field::dWdc_defaut;
                      }
                    else if (temp_potentiel_chimique_=="fonction")
                      {
                        potentiel_chimique_expr_.lire_f(is,1);
                        dWdc=&Source_Con_Phase_field::dWdc_general;
                      }
                    else
                      {
                        Cerr<<"Source_Con_Phase_field::readOn: "<<temp_potentiel_chimique_<<" is not a valid keyword in potentiel_chimique block"<<finl;
                        exit();
                      }
                    is >> motlu;
                    if (motlu!="}")
                      {
                        Cerr<<"Source_Con_Phase_field::readOn: We are expecting } at the end of potentiel_chimique block instead of "<< motlu <<finl;
                        exit();
                      }
                    break;
                  }
                case 3:
                  {
                    cpt0++;
                    is >> beta;
                    break;
                  }
                case 4:
                  {
                    cpt0++;
                    is >> kappa;
                    break;
                  }
                case 5:
                  {
                    cpt0++;
                    is >> motlu;
                    if (motlu!="{")
                      {
                        Cerr<<"Source_Con_Phase_field::readOn: We are expecting { at the end of kappa_variable block instead of " << motlu <<finl;
                        exit();
                      }
                    Motcle temp_type_kappa_;
                    is >> temp_type_kappa_;
                    if(temp_type_kappa_=="non")
                      {
                        type_kappa_=0;
                      }
                    else
                      {
                        type_kappa_=1;
                        if (temp_type_kappa_=="defaut")
                          {
                            kappa_func_c=&Source_Con_Phase_field::kappa_func_c_defaut;
                          }
                        else if (temp_type_kappa_=="fonction")
                          {
                            kappa_forme_expr_.lire_f(is,1);
                            kappa_func_c=&Source_Con_Phase_field::kappa_func_c_general;
                          }
                        else
                          {
                            Cerr<<"Source_Con_Phase_field::readOn: "<<temp_type_kappa_<<" is not a valid keyword in kappa_variable block"<<finl;
                            Cerr<<"You can specify only non, defaut or fonction expression for kappa_variable"<<finl;
                            exit();
                          }
                      }
                    is >> motlu;
                    if (motlu!="}")
                      {
                        Cerr<<"Source_Con_Phase_field::readOn: We are expecting } at the end of kappa_variable block instead of " << motlu <<finl;
                        exit();
                      }
                    break;
                  }
                default :
                  {
                    Cerr << "Source_Con_Phase_field::readOn: Error while reading systeme_naire" << finl;
                    Cerr << motlu << " is not understood."<< finl;
                    Cerr << "We are expecting a keyword among ";
                    for (int i=1; i<8; i++)
                      Cerr << les_mots(i) << " ";
                    Cerr<< finl;
                    exit();
                  }
                  if(cpt0 != 4)
                    {
                      Cerr << "Source_Con_Phase_field::readOn: Error while reading Source_Con_Phase_field: wrong number of parameters" << finl;
                      Cerr << "You should specify all these parameters: " << les_mots << finl;
                      exit();
                    }
                }
              is >>motlu;

            }
        }
      else if (temp_systeme_naire=="oui")
        {
          type_systeme_naire_=1;
          is >> motlu;
          if (motlu!="{")
            {
              Cerr<<"Source_Con_Phase_field::readOn: We are expecting { after 'oui' instead of "<< motlu <<finl;
              exit();
            }
          is >>motlu;
          int cpt0=0;
          while(motlu!="}")
            {
              int rang=les_mots.search(motlu);
              switch(rang)
                {
                case 0:
                  {
                    cpt0++;
                    is >> nb_equation_CH;
                    break;
                  }
                case 20:
                  {
                    cpt0++;
                    is >> motlu;
                    if (motlu=="oui")
                      {
                        DoubleVect temp_alpha(nb_equation_CH*nb_equation_CH); //define temp_alpha to avoid resize here
                        alphaMatrix = temp_alpha;
                        type_alpha_rotation_=1;
                        is >> motlu;
                        if (motlu!="{")
                          {
                            Cerr<<"Source_Con_Phase_field::readOn: We are expecting { after 'oui' instead of "<< motlu <<finl;
                            exit();
                          }
                        if (nb_equation_CH==2)
                          {
                            is >> motlu;
                            int cpt1=0;
                            Motcles param_alpha(3);
                            param_alpha(0)="alpha_ref";
                            param_alpha(1)="rotation_angle";
                            param_alpha(2)="diagonal_coefficient";
                            while(motlu!="}")
                              {
                                int rang0 = param_alpha.search(motlu);
                                switch(rang0)
                                  {
                                  case 0:
                                    {
                                      cpt1++;
                                      is >> alpha_ref;
                                      break;
                                    }
                                  case 1:
                                    {
                                      cpt1++;
                                      is >> angle_alphaMatrix;
                                      break;
                                    }
                                  case 2:
                                    {
                                      cpt1++;
                                      is >> diagonal_coeff;
                                      break;
                                    }
                                  default :
                                    {
                                      Cerr << "Source_Con_Phase_field::readOn: Error while reading alpha_rotation " << finl;
                                      Cerr << motlu << " is not understood."<< finl;
                                      Cerr << "We are expecting a keyword among " << param_alpha << finl;
                                      exit();
                                    }
                                    if(cpt1 != 3)
                                      {
                                        Cerr << "Source_Con_Phase_field::readOn: Error while reading Source_Con_Phase_field: wrong number of parameters" << finl;
                                        Cerr << "You should specify all these parameters: " << param_alpha << finl;
                                        exit();
                                      }
                                  }
                                is>> motlu;
                              }
                            double coeff_mult;
                            coeff_mult = alpha_ref/(2*pow(cos(angle_alphaMatrix* PI / 180.0),2)+diagonal_coeff*pow(sin(angle_alphaMatrix* PI / 180.0),2));
                            alphaMatrix(0)=2*pow(cos(angle_alphaMatrix* PI / 180.0),2)+diagonal_coeff*pow(sin(angle_alphaMatrix* PI / 180.0),2);
                            alphaMatrix(1)=2*cos(angle_alphaMatrix* PI / 180.0)*sin(angle_alphaMatrix* PI / 180.0)-diagonal_coeff*sin(angle_alphaMatrix* PI / 180.0)*cos(angle_alphaMatrix* PI / 180.0);
                            alphaMatrix(2)=alphaMatrix(1);
                            alphaMatrix(3)=2*pow(sin(angle_alphaMatrix* PI / 180.0),2)+diagonal_coeff*pow(cos(angle_alphaMatrix* PI / 180.0),2);
                            alphaMatrix*=coeff_mult;
                          }
                        else
                          {
                            Cerr<<"alpha rotation not yet implemented for nb_equation_CH>3"<<finl;
                          }
                      }
                    else if (motlu=="non")
                      {
                        type_alpha_rotation_=0;
                        is >> motlu;
                        if (motlu!="{")
                          {
                            Cerr<<"Source_Con_Phase_field::readOn: We are expecting { after 'non' instead of "<< motlu <<finl;
                            exit();
                          }
                        is >> motlu;
                        while(motlu!="}")
                          {
                            int rang1 = les_mots.search(motlu);
                            switch(rang1)
                              {
                              case 1:
                                {
                                  DoubleVect temp_alpha(nb_equation_CH*nb_equation_CH); //define temp_alpha to avoid resize here
                                  alphaMatrix = temp_alpha;
                                  for(int i=0; i< temp_alpha.size(); i++)
                                    is >> alphaMatrix(i);
                                  break;
                                }
                              default :
                                {
                                  Cerr << "Source_Con_Phase_field::readOn: Error while reading alpha_rotation " << finl;
                                  Cerr << motlu << " is not understood."<< finl;
                                  Cerr << "We are expecting a keyword among " << les_mots(1) << finl;
                                  exit();
                                }
                              }
                            is >> motlu;
                          }
                      }
                    break; //test break
                  }
                case 2:
                  {
                    cpt0++;
                    Motcle temp_potentiel_chimique_;
                    is >> temp_potentiel_chimique_;
                    if (temp_potentiel_chimique_=="defaut")
                      {
                        type_potentiel_analytique_=0;
                        is >> motlu;
                        if (motlu!="{")
                          {
                            Cerr<<"Source_Con_Phase_field::readOn: We are expecting { after 'defaut' instead of "<< motlu <<finl;
                            exit();
                          }
                        is >> motlu;
                        int cpt1=0;
                        Motcles param_potentiel_chimique(2);
                        param_potentiel_chimique(0)="xEq_phase1";
                        param_potentiel_chimique(1)="xEq_phase2";
                        while(motlu!="}")
                          {
                            int rang0 = param_potentiel_chimique.search(motlu);
                            switch(rang0)
                              {
                              case 0:
                                {
                                  cpt1++;
                                  DoubleVect temp_equilibre_phase1(nb_equation_CH); //define temp to avoid resize here
                                  eq1 = temp_equilibre_phase1;
                                  for(int i=0; i< temp_equilibre_phase1.size(); i++)
                                    is >> eq1(i);
                                  break;
                                }
                              case 1:
                                {
                                  cpt1++;
                                  DoubleVect temp_equilibre_phase2(nb_equation_CH); //define temp to avoid resize here
                                  eq2 = temp_equilibre_phase2;
                                  for(int i=0; i< temp_equilibre_phase2.size(); i++)
                                    is >> eq2(i);
                                  break;
                                }
                              default :
                                {
                                  Cerr << "Source_Con_Phase_field::readOn: Error while reading potentiel_chimique " << finl;
                                  Cerr << motlu << " is not understood."<< finl;
                                  Cerr << "We are expecting a keyword among " << param_potentiel_chimique << finl;
                                  exit();
                                }
                                if(cpt1 != 2)
                                  {
                                    Cerr << "Source_Con_Phase_field::readOn: Error while reading Source_Con_Phase_field: wrong number of parameters" << finl;
                                    Cerr << "You should specify all these parameters: " << param_potentiel_chimique << finl;
                                    exit();
                                  }
                              }
                            is>> motlu;
                          }
                        dWdc_naire=&Source_Con_Phase_field::dWdc_naire_defaut;
                      }
                    else if (temp_potentiel_chimique_=="fonction_analytique")
                      {
                        type_potentiel_analytique_=1;
                        if (nb_equation_CH==2)
                          {
                            is >> motlu;
                            if (motlu!="{")
                              {
                                Cerr<<"Source_Con_Phase_field::readOn: We are expecting { after 'fonction_analytique' instead of "<< motlu <<finl;
                                exit();
                              }
                            is >> motlu;
                            int cpt1=0;
                            Motcles param_potentiel_analytique_ternaire(5);
                            param_potentiel_analytique_ternaire(0)="psi";
                            param_potentiel_analytique_ternaire(1)="xEqComp0";
                            param_potentiel_analytique_ternaire(2)="xEqComp1";
                            param_potentiel_analytique_ternaire(3)="aEqComp0";
                            param_potentiel_analytique_ternaire(4)="aEqComp1";
                            while(motlu!="}")
                              {
                                int rang0 = param_potentiel_analytique_ternaire.search(motlu);
                                switch(rang0)
                                  {
                                  case 0:
                                    {
                                      cpt1++;
                                      DoubleVect temp_psi(2); //define temp to avoid resize here
                                      angle_psi=temp_psi;
                                      for(int i=0; i< 2; i++)
                                        is >> angle_psi(i);
                                      break;
                                    }
                                  case 1:
                                    {
                                      cpt1++;
                                      DoubleVect temp_x0Eq(2); //define temp to avoid resize here
                                      x0Eq = temp_x0Eq;
                                      for(int i=0; i< 2; i++)
                                        is >> x0Eq(i);
                                      break;
                                    }
                                  case 2:
                                    {
                                      cpt1++;
                                      DoubleVect temp_x1Eq(2); //define temp to avoid resize here
                                      x1Eq = temp_x1Eq;
                                      for(int i=0; i< 2; i++)
                                        is >> x1Eq(i);
                                      break;
                                    }
                                  case 3:
                                    {
                                      cpt1++;
                                      DoubleVect temp_a0Eq(2); //define temp to avoid resize here
                                      a0Eq=temp_a0Eq;
                                      for(int i=0; i< 2; i++)
                                        is >> a0Eq(i);
                                      break;
                                    }
                                  case 4:
                                    {
                                      cpt1++;
                                      DoubleVect temp_a1Eq(2); //define temp to avoid resize here
                                      a1Eq = temp_a1Eq;
                                      for(int i=0; i< 2; i++)
                                        is >> a1Eq(i);
                                      break;
                                    }
                                  default :
                                    {
                                      Cerr << "Source_Con_Phase_field::readOn: Error while reading potentiel_chimique " << finl;
                                      Cerr << motlu << " is not understood."<< finl;
                                      Cerr << "We are expecting a keyword among " << param_potentiel_analytique_ternaire << finl;
                                      exit();
                                    }
                                    if(cpt1 != 5)
                                      {
                                        Cerr << "Source_Con_Phase_field::readOn: Error while reading Source_Con_Phase_field: wrong number of parameters" << finl;
                                        Cerr << "You should specify all these parameters: " << param_potentiel_analytique_ternaire << finl;
                                        exit();
                                      }
                                  }
                                is>> motlu;
                              }
                            dWdc_naire_analytique_ter=&Source_Con_Phase_field::dWdc_naire_analytique_ternaire;
                          }
                        else if (nb_equation_CH==3)
                          {
                            is >> motlu;
                            if (motlu!="{")
                              {
                                Cerr<<"Source_Con_Phase_field::readOn: We are expecting { after 'fonction_analytique' instead of "<< motlu <<finl;
                                exit();
                              }
                            is >> motlu;
                            int cpt1=0;
                            Motcles param_potentiel_analytique_quaternaire(8);
                            param_potentiel_analytique_quaternaire(0)="psi";
                            param_potentiel_analytique_quaternaire(1)="phi";
                            param_potentiel_analytique_quaternaire(2)="xEqComp0";
                            param_potentiel_analytique_quaternaire(3)="xEqComp1";
                            param_potentiel_analytique_quaternaire(4)="xEqComp2";
                            param_potentiel_analytique_quaternaire(5)="aEqComp0";
                            param_potentiel_analytique_quaternaire(6)="aEqComp1";
                            param_potentiel_analytique_quaternaire(7)="aEqComp2";

                            while(motlu!="}")
                              {
                                int rang0 = param_potentiel_analytique_quaternaire.search(motlu);
                                switch(rang0)
                                  {
                                  case 0:
                                    {
                                      cpt1++;
                                      DoubleVect temp_psi(2); //define temp to avoid resize here
                                      angle_psi=temp_psi;
                                      for(int i=0; i< 2; i++)
                                        is >> angle_psi(i);
                                      break;
                                    }
                                  case 1:
                                    {
                                      cpt1++;
                                      DoubleVect temp_phi(2); //define temp to avoid resize here
                                      angle_phi=temp_phi;
                                      for(int i=0; i< 2; i++)
                                        is >> angle_phi(i);
                                      break;
                                    }
                                  case 2:
                                    {
                                      cpt1++;
                                      DoubleVect temp_x0Eq(2); //define temp to avoid resize here
                                      x0Eq = temp_x0Eq;
                                      for(int i=0; i< 2; i++)
                                        is >> x0Eq(i);
                                      break;
                                    }
                                  case 3:
                                    {
                                      cpt1++;
                                      DoubleVect temp_x1Eq(2); //define temp to avoid resize here
                                      x1Eq = temp_x1Eq;
                                      for(int i=0; i< 2; i++)
                                        is >> x1Eq(i);
                                      break;
                                    }
                                  case 4:
                                    {
                                      cpt1++;
                                      DoubleVect temp_x2Eq(2); //define temp to avoid resize here
                                      x2Eq = temp_x2Eq;
                                      for(int i=0; i< 2; i++)
                                        is >> x2Eq(i);
                                      break;
                                    }
                                  case 5:
                                    {
                                      cpt1++;
                                      DoubleVect temp_a0Eq(2); //define temp to avoid resize here
                                      a0Eq=temp_a0Eq;
                                      for(int i=0; i< 2; i++)
                                        is >> a0Eq(i);
                                      break;
                                    }
                                  case 6:
                                    {
                                      cpt1++;
                                      DoubleVect temp_a1Eq(2); //define temp to avoid resize here
                                      a1Eq = temp_a1Eq;
                                      for(int i=0; i< 2; i++)
                                        is >> a1Eq(i);
                                      break;
                                    }
                                  case 7:
                                    {
                                      cpt1++;
                                      DoubleVect temp_a2Eq(2); //define temp to avoid resize here
                                      a2Eq=temp_a2Eq;
                                      for(int i=0; i< 2; i++)
                                        is >> a2Eq(i);
                                      break;
                                    }
                                  default :
                                    {
                                      Cerr << "Source_Con_Phase_field::readOn: Error while reading potentiel_chimique " << finl;
                                      Cerr << motlu << " is not understood."<< finl;
                                      Cerr << "We are expecting a keyword among " << param_potentiel_analytique_quaternaire << finl;
                                      exit();
                                    }
                                    if(cpt1 != 8)
                                      {
                                        Cerr << "Source_Con_Phase_field::readOn: Error while reading Source_Con_Phase_field: wrong number of parameters" << finl;
                                        Cerr << "You should specify all these parameters: " << param_potentiel_analytique_quaternaire << finl;
                                        exit();
                                      }
                                  }
                                is>> motlu;
                              }
                            dWdc_naire_analytique_quater=&Source_Con_Phase_field::dWdc_naire_analytique_quaternaire;
                          }
                      }
                    else if (temp_potentiel_chimique_=="fonction")
                      {
                        //Mirantsoa 264902 generalisation of lire_f (is, 1) tp lire_f(is,nb_comp)
                        const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
                        const int& nb_comp =eq_c.constituant().nb_constituants();
                        Cerr << "nb_comp = "<<nb_comp<<finl;
                        potentiel_chimique_expr_.lire_f(is,nb_comp);
                        dWdc=&Source_Con_Phase_field::dWdc_general;
                      }
                    else
                      {
                        Cerr<<"Source_Con_Phase_field::readOn: "<<temp_potentiel_chimique_<<" is not a valid keyword in potentiel_chimique block"<<finl;
                        exit();
                      }
                    break;
                  }
                case 3:
                  {
                    cpt0++;
                    DoubleVect temp_beta(nb_equation_CH);
                    betaMatrix = temp_beta;
                    for(int i=0; i< nb_equation_CH; i++)
                      is >> betaMatrix(i);
                    break;
                  }
                case 21:
                  {
                    cpt0++;
                    DoubleVect temp_min_x(nb_equation_CH);
                    minX = temp_min_x;
                    for(int i=0; i< nb_equation_CH; i++)
                      is >> minX(i);
                    break;
                  }
                case 22:
                  {
                    cpt0++;
                    DoubleVect temp_max_x(nb_equation_CH);
                    maxX = temp_max_x;
                    for(int i=0; i< nb_equation_CH; i++)
                      is >> maxX(i);
                    break;
                  }
                case 19:
                  {
                    cpt0++;
                    is >> motlu;
                    if (motlu=="oui")
                      {
                        type_kappa_auto_diffusion_=1;
                        is >> motlu;
                        if (motlu!="{")
                          {
                            Cerr<<"Source_Con_Phase_field::readOn: We are expecting { after 'oui' instead of "<< motlu <<finl;
                            exit();
                          }
                        is >> motlu;
                        int cpt1=0;
                        Motcles param_mobilite(3);
                        param_mobilite(0)="coefficient_auto_diffusion";
                        param_mobilite(1)="temperature";
                        param_mobilite(2)="volume_molaire";
                        while(motlu!="}")
                          {
                            int rang0 = param_mobilite.search(motlu);
                            switch(rang0)
                              {
                              case 0:
                                {
                                  cpt1++;
                                  DoubleVect temp_coeff_diffusion(nb_equation_CH+1); //define temp_coeff_diffusion to avoid resize here
                                  coeff_auto_diffusion = temp_coeff_diffusion;
                                  for(int i=0; i< temp_coeff_diffusion.size(); i++)
                                    is >> coeff_auto_diffusion(i);
                                  break;
                                }
                              case 1:
                                {
                                  cpt1++;
                                  is >> temperature;
                                  break;
                                }
                              case 2:
                                {
                                  cpt1++;
                                  is >> molarVolume;
                                  break;
                                }
                              default :
                                {
                                  Cerr << "Source_Con_Phase_field::readOn: Error while reading kappa_auto_diffusion " << finl;
                                  Cerr << motlu << " is not understood."<< finl;
                                  Cerr << "We are expecting a keyword among " << param_mobilite << finl;
                                  exit();
                                }
                                if(cpt1 != 3)
                                  {
                                    Cerr << "Source_Con_Phase_field::readOn: Error while reading Source_Con_Phase_field: wrong number of parameters" << finl;
                                    Cerr << "You should specify all these parameters: " << param_mobilite << finl;
                                    exit();
                                  }
                              }
                            is>> motlu;
                          }
                        //is >> motlu;
                        kappaMatrix_func_c=&Source_Con_Phase_field::kappa_func_auto_diffusion;
                      }
                    else if (motlu=="non")
                      {
                        type_kappa_auto_diffusion_=0;
                        is >> motlu;
                        if (motlu!="{")
                          {
                            Cerr<<"Source_Con_Phase_field::readOn: We are expecting { after 'non' instead of "<< motlu <<finl;
                            exit();
                          }
                        is >> motlu;
                        while(motlu!="}")
                          {
                            int rang1 = les_mots.search(motlu);
                            switch(rang1)
                              {
                              case 4:
                                {
                                  DoubleVect temp_kappa(nb_equation_CH*nb_equation_CH); //define temp_kappa to avoid resize here
                                  kappaMatrix = temp_kappa;
                                  for(int i=0; i< temp_kappa.size(); i++)
                                    is >> kappaMatrix(i);
                                  break;
                                }
                              default :
                                {
                                  Cerr << "Source_Con_Phase_field::readOn: Error while reading kappa_auto_diffusion " << finl;
                                  Cerr << motlu << " is not understood."<< finl;
                                  Cerr << "We are expecting a keyword among " << les_mots(4) << finl;
                                  exit();
                                }
                              }
                            is >> motlu;
                          }
                      }
                    else
                      {
                        Cerr << "Source_Con_Phase_field::readOn: Error while reading systeme_naire" << finl;
                        Cerr << motlu << " is not understood."<< finl;
                        Cerr << "We are expecting a keyword among " << les_mots << finl;
                        exit();
                      }

                    break; //test break
                  }
                  if(cpt0 != 7)
                    {
                      Cerr << "Source_Con_Phase_field::readOn: Error while reading Source_Con_Phase_field: wrong number of parameters" << finl;
                      Cerr << "You should specify all these parameters: " << les_mots << finl;
                      exit();
                    }
                }

              is >> motlu;
            }
        }
    }
  is >> motlu;
  while (motlu!="}")
    {
      int rang=les_mots.search(motlu);
      switch(rang)
        {
        case 18:
          {
            cpt++;
            is >> tpsaff;
            if(tpsaff<0. || tpsaff >100.)
              {
                Cerr << "Source_Con_Phase_field::readOn: Temps_d_affichage should be in the range 0 - 100 seconds." << finl;
                exit();
              }
            break;
          }
        case 6:
          {
            cpt++;

            Motcle temp_kappa_moy_;
            is >> temp_kappa_moy_;
            if(temp_kappa_moy_=="arithmetique")
              {
                kappa_moy_=0;
              }
            else if(temp_kappa_moy_=="harmonique")
              {
                kappa_moy_=1;
              }
            else
              {
                kappa_moy_=2;
              }
            break;
          }
        case 7:
          {
            cpt++;
            is >> mult_kappa;
            break;
          }
        case 8:
          {
            cpt++;
            Motcle temp_couplage_;
            is >> temp_couplage_;
            if(temp_couplage_=="mutilde(n)")
              {
                couplage_=0;
              }
            else
              {
                couplage_=1;
              }
            break;
          }
        case 9:
          {
            cpt++;
            Motcle temp_implicitation_;
            is >> temp_implicitation_;
            if(temp_implicitation_=="oui")
              {
                implicitation_=1;
              }
            else
              {
                implicitation_=0;
              }
            break;
          }
        case 10:
          {
            cpt++;

            Motcle temp_gmres_;
            is >> temp_gmres_;
            if(temp_gmres_=="oui")
              {
                gmres_=1;
              }
            else
              {
                gmres_=0;
              }
            break;
          }
        case 11:
          {
            cpt++;

            is >> epsilon_;
            break;
          }
        case 12:
          {
            cpt++;
            is >> eps_;
            break;
          }
        case 13:
          {
            cpt++;
            is >> epsGMRES;
            break;
          }
        case 14:
          {
            cpt++;
            is >> nkr;
            break;
          }
        case 15:
          {
            cpt++;
            is >> nit;
            break;
          }

        case 16:
          {
            cpt++;
            is >> rec_min;
            break;
          }
        case 17:
          {
            cpt++;
            is >> rec_max;
            break;
          }
        default :
          {
            Cerr << "Source_Con_Phase_field::readOn: Error while reading Source_Con_Phase_field" << finl;
            Cerr << motlu << " is not understood."<< finl;
            Cerr << "We are expecting a keyword among " << les_mots << finl;
            exit();
          }
        }
      is >>motlu;

    }
  if(cpt != 13)
    {
      Cerr << "Source_Con_Phase_field::readOn: Error while reading Source_Con_Phase_field: wrong number of parameters" << finl;
      Cerr << "You should specify all these parameters: " << les_mots << finl;
      exit();
    }

  if (type_systeme_naire_==0)
    if(kappa>0)
      {
        Cerr <<" theoretical time step = ?     dx4*"<< 1./alpha/kappa/2.<<" dx2*"<<1./2./kappa/beta<<finl;
      }
  return is ;

}

void Source_Con_Phase_field::associer_pb(const Probleme_base& pb)
{
  le_probleme2=pb;
  Navier_Stokes_phase_field& eq_ns=ref_cast(Navier_Stokes_phase_field,le_probleme2->equation(0));
  Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  rho0 = eq_ns.rho0();
  drhodc_ = eq_ns.drhodc();
  // Dans le modele cette derivee est constante - On la calcule au debut.

  boussi_=eq_ns.get_boussi_();
  if (boussi_!=1 && boussi_!=0)
    {
      Cerr << "Erreur dans le choix du parametre boussi_" << finl;
      exit();
    }

  diff_boussi_=eq_ns.get_diff_boussi_();
  if (diff_boussi_!=1 && diff_boussi_!=0)
    {
      Cerr << "Erreur dans le choix du parametre diff_boussi_" << finl;
      exit();
    }

  g_=eq_ns.get_g_();

  mutype_=eq_c.get_mutype_();
  if (mutype_!=0 && mutype_!=1)
    {
      Cerr << "Erreur dans le choix du parametre mutype_" << finl;
      exit();
    }

  if(implicitation_==1)
    {
      if(gmres_!=0 && gmres_!=1)
        {
          Cerr << "** Erreur du choix de la methode implicite GMRES - L'execution doit stopper **" << finl;
          Cerr << "Choix possibles : 0, 1" << finl;
          exit();
        }
    }
  else if(implicitation_!=0)
    {
      Cerr<<"================================================="<<finl;
      Cerr<<"Choix de l'implicitation incorrect ! (ni 0, ni 1)"<<finl;
      Cerr<<"================================================="<<finl;
      exit();
    }
  if (type_systeme_naire_==0)
    {
      if(type_kappa_==1)
        {
          kappa_ind=1;
        }
      else if(type_kappa_==0)
        {
          kappa_ind=0;
        }
      else
        {
          Cerr << "Erreur dans le choix de kappa !" << finl;
          exit();
        }
    }
  if (type_systeme_naire_==1)
    {
      if(type_kappa_auto_diffusion_==1)
        {
          kappa_ind=1;
        }
      else if(type_kappa_auto_diffusion_==0)
        {
          kappa_ind=0;
        }
      else
        {
          Cerr << "Erreur dans le choix de kappa !" << finl;
          exit();
        }
    }

  if(implicitation_==1)
    {
      if(kappa_moy_!=0 && kappa_moy_!=1 && kappa_moy_!=2)
        {
          Cerr << "Erreur dans le choix de la moyenne de kappa !" << finl;
          exit();
        }
    }

  if (type_systeme_naire_==0)
    {
      if(mult_kappa<=0)
        {
          Cerr << "Erreur dans le choix du multiplicateur de kappa !" << finl;
        }
    }


  // Recapitulatif des parametres

  Nom choix_boussi;
  if(boussi_==1)
    {
      choix_boussi="oui";
    }
  else
    {
      choix_boussi="non";
    }

  Nom choix_diff_boussi;
  if(diff_boussi_==1)
    {
      choix_diff_boussi="oui";
    }
  else
    {
      choix_diff_boussi="non";
    }
  if(boussi_==1 && diff_boussi_==1)
    {
      Cerr << "ATTENTION : les options 'approximation de Boussinesq' et" << "'approximation de Boussinesq dans le terme de diffusion' ne devrait pas etre utilisees simultanement'." << "(voir Source_Con_Phase_field.cpp)" << finl;
    }

  const int terme_source=eq_ns.getset_terme_source_();

  Nom choix_source;
  if(terme_source==1)
    {
      choix_source="c grad(mutilde)";
    }
  else if(terme_source==2)
    {
      choix_source="c grad(laplacien(c))";
    }
  else if(terme_source==3)
    {
      choix_source="c grad(laplacien(c))-div((grad(c))^2)/2";
    }
  else if(terme_source==4)
    {
      choix_source="-laplacien(c) grad(c)";
    }

  Nom choix_implicite;
  if(implicitation_==1)
    {
      choix_implicite="oui";
    }
  else
    {
      choix_implicite="non";
    }

  Nom type_implicite;
  if(gmres_==1)
    {
      type_implicite="GMRES non lineaire";
    }
  else
    {
      type_implicite="point fixe";
    }

  Nom mutilde_d;
  if(mutype_==1)
    {
      mutilde_d="avec le terme d'Ec";
    }
  else
    {
      mutilde_d="sans le terme d'Ec";
    }

  Nom type_couplage;
  if(couplage_==1)
    {
      type_couplage="potentiel au temps n+1/2";
    }
  else
    {
      type_couplage="potentiel au temps n";
    }
  Nom mobilite_variable;
  Nom moyenne_kappa;

  if(type_kappa_==1 || type_kappa_auto_diffusion_==1)
    {
      mobilite_variable="oui";
    }
  else
    {
      mobilite_variable="non";
    }

  if (type_systeme_naire_==0) //implicitation ??
    {
      if(kappa_moy_==0)
        {
          moyenne_kappa="arithmetique (+)";
        }
      else if(kappa_moy_==1)
        {
          moyenne_kappa="harmonique (/)";
        }
      else if(kappa_moy_==2)
        {
          moyenne_kappa="geometrique (*)";
        }
    }

  Cerr << "" << finl;
  Cerr << "" << finl;
  Cerr << "******************************************************************************" << finl;
  Cerr << "*********************** RECAPITULATIF DU PARAMETRAGE *************************" << finl;
  Cerr << "" << finl;
  Cerr << "1) Equation de Navier-Stokes" << finl;
  Cerr << "  - approximation de Boussinesq                     : " << choix_boussi << finl;
  if(boussi_==0)
    {
      Cerr << "  - approximation de Boussinesq dans la diffusion   : " << choix_diff_boussi << finl;
    }
  else
    {
      Cerr << "  - masse volumique de reference                    : " << rho0 << " kg/m3" << finl;
    }
  Cerr << "  Dans le cas ou la viscosite dynamique est variable : " << finl;
  Cerr << "  - terme source                                    : " << choix_source << finl;
  Cerr << "  - intensite du champ de gravite                   : " << norme_array(g_)  << " m/s^2" << finl;
  Cerr << "" << finl;
  Cerr << "2) Equation de Cahn-Hilliard" << finl;
  Cerr << "  - discretisation implicite                        : " << choix_implicite << finl;
  if(implicitation_==1)
    {
      Cerr << "  - algorithme implicite                            : " << type_implicite << finl;
    }
  Cerr << "  - potentiel chimique generalise                   : " << mutilde_d << finl;
  Cerr << "  - couplage NS / CH via le potentiel chimique      : " << type_couplage << finl;

  Cerr << "  - mobilite variable                               : " << mobilite_variable << finl;
  if(type_kappa_==1 || type_kappa_auto_diffusion_==1)
    {
      Cerr << "  - type de moyenne pour la mobilite                : " << moyenne_kappa << finl;
      Cerr << "  - multiplicateur de la mobilite                   : " << mult_kappa << finl;
    }

  if(implicitation_==1 && gmres_==0)
    {
      Cerr << "  - dans le cas du point fixe ... " << finl;
      Cerr << "        seuil de convergence entre deux iterations  : " << epsilon_ << finl;
      Cerr << "        seuil du residu a convergence               : " << eps_ << finl;
    }
  else if(implicitation_==1 && gmres_==1)
    {
      Cerr << "  - dans le cas du GMRES non lineaire ... " << finl;
      Cerr << "        seuil du residu a convergence               : " << epsGMRES << finl;
      Cerr << "        dimension de l'espace de Krylov             : " << nkr << finl;
      Cerr << "        nombre d'iterations maximal de l'algorithme : " << nit << finl;
      Cerr << "        residu minimum de convergence               : " << rec_min << finl;
      Cerr << "        residu maximum de convergence               : " << rec_max << finl;
    }
  Cerr << "" << finl;
  Cerr << "********************** Ce message persistera " << tpsaff << " secondes *********************" << finl;
  Cerr << "*************** Pour continuer avant les " << tpsaff << " secondes : Ctrl + C **************" << finl;
  Cerr << "******************************************************************************" << finl;
  Cerr << "" << finl;
  Cerr << "" << finl;

  Nom tps_sleep="sleep ";
  tps_sleep+=Nom(tpsaff);
  Cerr << (int) system(tps_sleep) << finl;
}

void Source_Con_Phase_field::associer_zones(const Zone_dis& zone_dis,
                                            const Zone_Cl_dis& zone_Cl_dis)
{
  la_zone_VDF = ref_cast(Zone_VDF, zone_dis.valeur());
  la_zone_Cl_VDF = ref_cast(Zone_Cl_VDF, zone_Cl_dis.valeur());
}


inline double mobilite(const double c)
{
  return (1.);
  //   return (c*c);
  //   return (0.5+2*c*c);
  //  return (std::min(1.,4.*c*c));
  //   return(0.);
  //   const double clim1 = -0.4;
  //   const double clim2 =  0.4;
  //   if (c<=clim1)
  //     return(std::min(1.,-(c-clim1)/(0.5+clim1)));
  //   else if (c>=clim2)
  //     return(std::min(1.,(c-clim2)/(0.5-clim2)));
  //   else return(0.);
}

DoubleTab& Source_Con_Phase_field::laplacien(const DoubleTab& F, DoubleTab& resu) const
{
  const Navier_Stokes_std& eq_ns=ref_cast(Navier_Stokes_std,le_probleme2->equation(0));
  const Operateur_Grad& opgrad=eq_ns.operateur_gradient();
  const Operateur_Div& opdiv=eq_ns.operateur_divergence();

  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const DoubleTab& c=eq_c.inconnue().valeurs();
  const int nb_comp = eq_c.constituant().nb_constituants();


  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const IntTab& face_voisins = zone_VDF.face_voisins();
  const DoubleVect& volumes = zone_VDF.volumes();

  DoubleTab& prov_face=ref_cast_non_const( DoubleTab, prov_face_);
  if (prov_face.size()==0)
    prov_face=eq_ns.inconnue().valeurs();
  prov_face=0.;
  resu=0.;

  if (type_systeme_naire_==0)
    {
      // Grad(F)
      opgrad.calculer(F,prov_face);
      // M*Grad(F)
      int ndeb=zone_VDF.premiere_face_int();
      int nbfaces=zone_VDF.nb_faces();
      int el0,el1;
      double cface,vol0,vol1;
      for (int fac=ndeb; fac<nbfaces; fac++)
        {
          el0=face_voisins(fac,0);
          el1=face_voisins(fac,1);
          vol0=volumes(el0);
          vol1=volumes(el1);
          cface=(vol0*c(el0)+vol1*c(el1))/(vol0+vol1);
          prov_face(fac)=mobilite(cface)*prov_face(fac);
        }
      //   int taille=prov_face.size();
      //   Cerr << "taille : " << taille << finl;
      //   for (int i=0;i<taille;i++)
      //     {
      //       prov_face(i)*=mobilite(c(i));
      //     }
      //   Cerr << "Fin multiplication mobilite" << finl;
      // Application solveur masse
      eq_ns.solv_masse().appliquer(prov_face);

      // Div(M*Grad(F))
      opdiv.calculer(prov_face,resu);
    }
  else if (type_systeme_naire_==1)
    {
      // Grad(F)
      DoubleTab temp_resu(resu.dimension_tot(0),1);
      temp_resu=0;
      DoubleTab temp_F(F.dimension_tot(0),1);
      for (int j=0; j<nb_comp; j++)
        {
          prov_face=0.;
          for (int i=0; i<temp_F.dimension(0); i++)
            {
              temp_F(i,0)=F(i,j);
            }
          opgrad.calculer(temp_F,prov_face);

          // Application solveur masse
          eq_ns.solv_masse().appliquer(prov_face);

          // Laplacien = Div(Grad(F))
          opdiv.calculer(prov_face,temp_resu);
          for (int i=0; i<temp_resu.dimension(0); i++)
            {
              resu(i,j)=temp_resu(i,0);
            }
        }

    }
  return resu;

}

DoubleTab& Source_Con_Phase_field::div_kappa_grad(const DoubleTab& F, const DoubleTab& kappa_var, DoubleTab& resu) const
{
  const Navier_Stokes_std& eq_ns=ref_cast(Navier_Stokes_std,le_probleme2->equation(0));
  const Operateur_Grad& opgrad=eq_ns.operateur_gradient();
  const Operateur_Div& opdiv=eq_ns.operateur_divergence();

  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const DoubleTab& c=eq_c.inconnue().valeurs();
  const int nb_comp = eq_c.constituant().nb_constituants();

  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const IntTab& face_voisins = zone_VDF.face_voisins();
  const DoubleVect& volumes = zone_VDF.volumes();


  DoubleTab& prov_face=ref_cast_non_const( DoubleTab, prov_face_);
  if (prov_face.size()==0)
    prov_face=eq_ns.inconnue().valeurs();
  prov_face=0.;
  resu=0.;

  if (type_systeme_naire_==0)
    {
      // Grad(F)
      opgrad.calculer(F,prov_face);

      // M*Kappa*Grad(F)
      int ndeb=zone_VDF.premiere_face_int();
      int nbfaces=zone_VDF.nb_faces();
      int el0,el1;
      double cface,kappa_face,vol0,vol1;

      for (int fac=ndeb; fac<nbfaces; fac++)
        {
          el0=face_voisins(fac,0);
          el1=face_voisins(fac,1);
          vol0=volumes(el0);
          vol1=volumes(el1);
          kappa_face=(vol0*kappa_var(el0)+vol1*kappa_var(el1))/(vol0+vol1);
          cface=(vol0*c(el0)+vol1*c(el1))/(vol0+vol1);
          prov_face(fac)=kappa_face*mobilite(cface)*prov_face(fac);
        }

      // Application solveur masse (/M)
      eq_ns.solv_masse().appliquer(prov_face);

      // Div(Kappa*Grad(F))
      opdiv.calculer(prov_face,resu);

    }
  else if (type_systeme_naire_==1)
    {
      // Grad(F)
      DoubleTab temp_resu(resu.dimension_tot(0),1);
      temp_resu=0;
      DoubleTab temp_F(F.dimension_tot(0),1);
      DoubleTab temp_prov_face(prov_face.dimension_tot(0),F.line_size());
      for (int j=0; j<nb_comp; j++)
        {
          for (int i=0; i<temp_F.dimension(0); i++)
            {
              temp_F(i,0)=F(i,j);
            }
          opgrad.calculer(temp_F,prov_face);

          for (int i=0; i<prov_face.dimension(0); i++)
            {
              temp_prov_face(i,j)=prov_face(i,0);
            }
        }

      // Determine kappa_face et cface (faces interieures)
      int ndeb=zone_VDF.premiere_face_int();
      int nbfaces=zone_VDF.nb_faces();
      int el0,el1;
      double vol0,vol1;
      DoubleTab kappa_face(prov_face.dimension(0),nb_comp*nb_comp);
      DoubleTab cface(prov_face.dimension(0),nb_comp);
      for (int j=0; j<nb_comp*nb_comp; j++)
        {
          for (int fac=ndeb; fac<nbfaces; fac++)
            {
              el0=face_voisins(fac,0);
              el1=face_voisins(fac,1);
              vol0=volumes(el0);
              vol1=volumes(el1);
              // on calcule kappa_face comme la moyenne des kappa des voisins
              kappa_face(fac,j)=(vol0*kappa_var(el0,j)+vol1*kappa_var(el1,j))/(vol0+vol1);
            }
        }

      // kappa*Grad(F) - sur les faces
      DoubleTab temp_prov_face2 = temp_prov_face;
      temp_prov_face2 = 0;
      for (int j=0; j<nb_comp; j++)
        {
          for (int i=0; i<prov_face.dimension(0); i++)
            {
              prov_face(i,0)=temp_prov_face(i,j);
              for (int k=0; k<nb_comp; k++)
                {
                  temp_prov_face2(i,k)+=prov_face(i,0)*kappa_face(i,j+k*nb_comp);
                }
            }
        }

      // Application solveur masse (/M)
      eq_ns.solv_masse().appliquer(temp_prov_face2);

      // Div(Kappa*Grad(F))
      for (int j=0; j<nb_comp; j++)
        {
          for (int i=0; i<prov_face.dimension(0); i++)
            {
              prov_face(i,0)=temp_prov_face2(i,j);
            }
          opdiv.calculer(prov_face,temp_resu);
          for (int i=0; i<temp_resu.dimension(0); i++)
            {
              resu(i,j)=temp_resu(i,0);
            }
        }
    }

  return resu;
}


DoubleTab& Source_Con_Phase_field::ajouter(DoubleTab& resu) const
{
  //const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));

  resu+=accr;
  return resu;
}

DoubleTab& Source_Con_Phase_field::calculer(DoubleTab& resu) const
{
  resu = 0.;
  return ajouter(resu);
}

void Source_Con_Phase_field::mettre_a_jour(double temps)
{
  Cerr << "Temps : " << temps << " s" << finl;
  Cerr << "" << finl;

  Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  Champ_Fonc& ch=eq_c.set_mutilde_();
  ch.mettre_a_jour(temps);
}


/*! @brief Calcule le premier demi pas de temps dans le cas implicite Calcule le pas de temps dans le cas explicite
 *
 */
void Source_Con_Phase_field::premier_demi_dt()
{
  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const int nb_elem = zone_VDF.nb_elem_tot();

  Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  DoubleTab& c=eq_c.inconnue().valeurs();

  DoubleTab& mutilde=eq_c.set_mutilde_().valeurs();
  DoubleTab& mutilde_demi=eq_c.set_mutilde_demi();
  DoubleTab& c_demi=eq_c.set_c_demi();

  // Utilise dans l'ancien modele de Didier Jamet

  /* commente par mr264902 car pas utilise dans la version présente
  //DoubleTab& alpha_gradC_carre=eq_c.set_alpha_gradC_carre();
  //calculer_alpha_gradC_carre(alpha_gradC_carre);
  //DoubleTab& div_alpha_rho_gradC=eq_c.set_div_alpha_rho_gradC();
  //
  //calculer_div_alpha_rho_gradC(div_alpha_rho_gradC);
  //
  // ---
   */
  DoubleTab& div_alpha_gradC=eq_c.set_div_alpha_gradC();
  calculer_div_alpha_gradC(div_alpha_gradC);
  // commente par mr264902 car pas utilise dans la version presente
  //DoubleTab& pression_thermo=eq_c.set_pression_thermo();
  //calculer_pression_thermo(pression_thermo);

  accr = c;
  accr = 0.;
  c_demi = c;
  c_demi = 0.;

  if(implicitation_==1)
    {
      // Cerr<<"Schema implicite"<<finl;
      Matrice_Morse matrice_diffusion_CH;
      mutilde_demi = c;
      mutilde_demi = 0.;

      // Pour le GMRES NL ----------------------------

      // Commente par DJ
      //----------------
      //       DoubleTab x1(2*nb_elem);
      //       x1=0.;
      //       if(gmres_==1)
      //         {

      //           // Ecriture de x1(0)
      //           for(int n_elem=0;n_elem<nb_elem;n_elem++)
      //             {
      //               x1(n_elem)=c(n_elem);
      //               x1(n_elem+nb_elem)=mutilde(n_elem);
      //             }
      //         }
      //----------------
      // ---------------------------------------------

      calculer_u2_elem(u_carre_);
      assembler_matrice_point_fixe(matrice_diffusion_CH);

      if(gmres_==0)
        {
          calculer_point_fixe(c, mutilde, matrice_diffusion_CH, c_demi, mutilde_demi);
        }
      else if(gmres_==1)
        {
          // Modifie par DJ
          //---------------
          //           non_lin_gmres(c, mutilde, matrice_diffusion_CH, x1);

          if (getenv("TRUST_GMRES"))
            {
              /*
                   // Creation  d'un solveur et typage:
                   SolveurSys solveur_;
                   //Nom str = "Petsc Gmres { precond diag { } seuil 1.e-6 impr }";
                   Nom str = "Petsc Cholesky { impr }";

                   EChaine chl(str);
                   chl >> solveur_;

                   //DoubleVect x,b;
                   // Remplissage guess x et second membre

                   //const int ns = 2*c.size_totale();
                   //int nb_elem_tot = c.size_totale();

                   DoubleTab x1;
                   int nb_inc = 2;
                   ConstDoubleTab_parts tab1(c);
                   ConstDoubleTab_parts tab2(mutilde);

                   DoubleTabs tabs(nb_inc);
                   tabs[0].ref(tab1[0]);
                   tabs[1].ref(tab2[0]);
                   // Autres inconnues

                   MD_Vector_composite mds;
                   int dim = 0;
                   for (int inc=0; inc<nb_inc; inc++)
                     {
                       mds.add_part(tabs[inc].get_md_vector(),0);
                       dim += tabs[inc].dimension_tot(0);
                     }
                   MD_Vector md;
                   md.copy(mds);

                   x1 = DoubleTab(dim ,1);
                   x1.set_md_vector(md);


                   DoubleTab x1(nb_elem_tot,2);
                   x1=0.;
                   zone_VDF.zone().creer_tableau_elements(x1, Array_base::NOCOPY_NOINIT);
                   double* x1_addr = x1.addr();
                   for(int n_elem=0; n_elem<nb_elem; n_elem++)
                     {
                   //                  x1(n_elem,0)=c(n_elem);
                   //                  x1(n_elem,1)=mutilde(n_elem);

                       x1_addr[n_elem]=c(n_elem);
                       x1_addr[n_elem+nb_elem]=mutilde(n_elem);
                     }



                   DoubleVect term_cin(nb_elem);
                   // Ceci correspond a u2/2*drhodc=mutilde_dyn-mutilde

                   DoubleVect second_membre(x1);
                   DoubleTab_parts sm(second_membre); // De taille nombre inconnues nb_inc
                   assert(sm.size()==nb_inc);
                   for(int n_elem=0; n_elem<nb_elem; n_elem++)
                     {
                       sm[0](n_elem)=c(n_elem);
                       term_cin(n_elem)=0.;
                       if (mutype_==1)
                         {
                           term_cin(n_elem)+=(0.5*u_carre_(n_elem))*drhodc(n_elem);
                         }
                       sm[1](n_elem)=term_cin(n_elem)+beta*(this->*dWdc)(c(n_elem));
                     }


                    DoubleTab second_membre(nb_elem, 2);
                   zone_VDF.zone().creer_tableau_elements(second_membre, Array_base::NOCOPY_NOINIT);

                   // Assemblage du second membre
                   double* sm_addr = second_membre.addr();

                   for(int n_elem=0; n_elem<nb_elem; n_elem++)
                     {
                       sm_addr[n_elem]=c(n_elem);

                       term_cin(n_elem)=0.;
                       if (mutype_==1)
                         {
                           term_cin(n_elem)+=(0.5*u_carre_(n_elem))*drhodc(n_elem);
                         }

                       sm_addr[n_elem+nb_elem]=term_cin(n_elem)+beta*(this->*dWdc)(x1_addr[n_elem]);
                     }

                   // Resolution:
                   solveur_.valeur().reinit();
                   solveur_.resoudre_systeme(matrice_diffusion_CH,second_membre,x1);

                   // remplissage des resultats dans mutilde demi et cdemi

                                 for(int n_elem=0; n_elem<nb_elem; n_elem++)
                                   {
                                     c_demi(n_elem)=x1_addr[n_elem];
                                     mutilde_demi(n_elem)=x1_addr[n_elem+nb_elem];
                                   }
                                 c_demi.echange_espace_virtuel();
                                 mutilde_demi.echange_espace_virtuel();
                   DoubleTab_parts x(x1);
                   assert(x.size()==nb_inc);
                   for(int n_elem=0; n_elem<nb_elem; n_elem++)
                     {
                       c_demi(n_elem)=x[0](n_elem);
                       mutilde_demi(n_elem)=x[1](n_elem);
                     }
                   c_demi.echange_espace_virtuel();
                   mutilde_demi.echange_espace_virtuel();*/

              // Creation  d'un solveur et typage:
              SolveurSys solveur_;
              //Nom str = "Petsc Gmres { precond diag { } seuil 1.e-6 impr }";
              Nom str = "Petsc Cholesky { impr }";

              EChaine chl(str);
              chl >> solveur_;

              //DoubleVect x,b;
              // Remplissage guess x et second membre

              //const int ns = 2*c.size_totale();
              int nb_elem_tot = c.size_totale();
              /*
                                DoubleTab x1;
                                int nb_inc = 2;
                                //DoubleTab X_c(nb_equation_CH*nb_elem,1);
                                //DoubleTab X_mutilde(X_c);
                                DoubleTab X_c(c.size_totale(),1);
                                DoubleTab X_mutilde(mutilde.size_totale(),1);

                                for (int j=0; j<nb_equation_CH; j++)
                                  for (int i=0; i<nb_elem; i++)
                                    {
                                      X_c(i+(j*nb_elem),0)=c(i,j);
                                      X_mutilde(i+(j*nb_elem),0)=mutilde(i,j);
                                    }
                                ConstDoubleTab_parts tab_parts_1(X_c);
                                ConstDoubleTab_parts tab_parts_2(X_mutilde);

                                Cerr <<"tab_parts_1(c)"<<tab_parts_1.size()<<finl;
                                Cerr <<"tab_parts_1[0]"<<tab_parts_1[0]<<finl;

                                DoubleTabs tabs(nb_inc);
                                tabs[0].ref(tab_parts_1[0]);
                                tabs[1].ref(tab_parts_2[0]);

                                Cerr <<"tabs[0]"<<tabs[0]<<finl;
                                Cerr <<"tabs[1]"<<tabs[1]<<finl;
                                Cerr <<"tabs"<<tabs<<finl;


                                // Autres inconnues

                                MD_Vector_composite mds;
                                int dim = 0;
                                for (int inc=0; inc<nb_inc; inc++)
                                  {
                                    mds.add_part(tabs[inc].get_md_vector(),0);
                                    dim += tabs[inc].dimension_tot(0);
                                  }
                                Cerr <<"mds"<<mds<<finl;
                                MD_Vector md;
                                md.copy(mds);

                                x1 = DoubleTab(dim ,1);
                                x1.set_md_vector(md);*/
              /*DoubleTab x1;
                  int nb_inc = 2;
                  ConstDoubleTab_parts tab1(c);
                  ConstDoubleTab_parts tab2(mutilde);
                  Cerr <<"tab_parts_1 size"<<tab1.size()<<finl;
                  Cerr <<"tab_parts_1[0]"<<tab1[0]<<finl;

                  DoubleTabs tabs(nb_inc);
                  tabs[0].ref(tab1[0]);
                  tabs[1].ref(tab2[0]);
                  // Autres inconnues
                  Cerr <<"tabs[0]"<<tabs[0]<<finl;
                  Cerr <<"tabs[1]"<<tabs[1]<<finl;
                  Cerr <<"tabs[1]dimension_tot"<<tabs[1].dimension_tot(0)<<finl;

                  Cerr <<"tabs"<<tabs<<finl;
                  MD_Vector_composite mds;
                  int dim = 0;
                  for (int inc=0; inc<nb_inc; inc++)
                    {
                      mds.add_part(tabs[inc].get_md_vector(),0);
                      dim += tabs[inc].dimension_tot(0);
                    }
                  Cerr <<"mds"<<mds<<finl;
                  MD_Vector md;
                  md.copy(mds);

                  x1 = DoubleTab(dim ,1);
                  x1.set_md_vector(md);
               */
              DoubleTab X_c(c.size_totale());
              DoubleTab X_mutilde(mutilde.size_totale());
              DoubleTab x1(nb_elem,2);
              x1=0.;
              zone_VDF.zone().creer_tableau_elements(x1, Array_base::NOCOPY_NOINIT);
              double* x1_addr = x1.addr();
              for (int j=0; j<nb_equation_CH; j++)
                for (int i=0; i<nb_elem; i++)
                  {
                    X_c(i+(j*nb_elem))=c(i,j);
                    X_mutilde(i+(j*nb_elem))=mutilde(i,j);
                  }

              for(int n_elem=0; n_elem<c.size_totale(); n_elem++)
                {
                  //                  x1(n_elem,0)=c(n_elem);
                  //                  x1(n_elem,1)=mutilde(n_elem);

                  x1_addr[n_elem]=X_c(n_elem);
                  x1_addr[n_elem+c.size_totale()]=X_mutilde(n_elem);
                }



              /*
                                DoubleVect term_cin(nb_elem);
                                // Ceci correspond a u2/2*drhodc=mutilde_dyn-mutilde

                                DoubleVect second_membre(x1);
                                DoubleTab_parts sm(second_membre); // De taille nombre inconnues nb_inc
                                assert(sm.size()==nb_inc);
                                for(int n_elem=0; n_elem<nb_elem; n_elem++)
                                  {
                                    sm[0](n_elem)=c(n_elem);
                                    term_cin(n_elem)=0.;
                                    if (mutype_==1)
                                      {
                                        term_cin(n_elem)+=(0.5*u_carre_(n_elem))*drhodc(n_elem);
                                      }
                                    sm[1](n_elem)=term_cin(n_elem)+beta*(this->*dWdc)(c(n_elem));
                                  }*/

              DoubleTab sm_c(nb_elem_tot);
              DoubleTab sm_mutilde(nb_elem_tot);
              DoubleTab second_membre(nb_elem_tot, 2);
              zone_VDF.zone().creer_tableau_elements(second_membre, Array_base::NOCOPY_NOINIT);

              // Assemblage du second membre
              double* sm_addr = second_membre.addr();
              for (int j=0; j<nb_equation_CH; j++)
                for (int i=0; i<nb_elem; i++)
                  {
                    sm_c(i+(j*nb_elem))=c(i,j);

                    sm_mutilde(i+(j*nb_elem))=betaMatrix(j)*(this->*dWdc)(c(i,j));
                  }
              for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
                {
                  sm_addr[n_elem]=sm_c(n_elem);

                  /*term_cin(n_elem)=0.;
                        if (mutype_==1)
                        {
                          term_cin(n_elem)+=(0.5*u_carre_(n_elem))*drhodc(n_elem);
                        }*/

                  sm_addr[n_elem+nb_elem]=sm_mutilde(n_elem);
                }

              // Resolution:
              solveur_.valeur().reinit();
              solveur_.resoudre_systeme(matrice_diffusion_CH,second_membre,x1);

              // remplissage des resultats dans mutilde demi et cdemi

              for(int n_elem=0; n_elem<nb_elem; n_elem++)
                {
                  c_demi(n_elem)=x1_addr[n_elem];
                  mutilde_demi(n_elem)=x1_addr[n_elem+nb_elem];
                }
              c_demi.echange_espace_virtuel();
              mutilde_demi.echange_espace_virtuel();
              /*DoubleTab_parts x(x1);
                  assert(x.size()==nb_inc);
                  for(int n_elem=0; n_elem<nb_elem; n_elem++)
                    {
                      c_demi(n_elem)=x[0](n_elem);
                      mutilde_demi(n_elem)=x[1](n_elem);
                    }
                  c_demi.echange_espace_virtuel();
                  mutilde_demi.echange_espace_virtuel();*/
            }
          else
            non_lin_gmres(c, mutilde, matrice_diffusion_CH, c_demi, mutilde_demi);
          //---------------

          const double theta=0.6;
          // On stocke les nouveaux c et mutilde


          if (type_systeme_naire_==0)
            {
              for(int n_elem=0; n_elem<nb_elem; n_elem++)
                {
                  // Commente par DJ
                  //----------------
                  //               c_demi(n_elem)=x1(n_elem);
                  //----------------
                  c_demi(n_elem)-=(1-theta)*c(n_elem);
                  c_demi(n_elem)/=theta;

                  // Commente par DJ
                  //----------------
                  //               mutilde_demi(n_elem)=x1(n_elem+nb_elem);
                  //----------------
                  mutilde_demi(n_elem)-=(1-theta)*mutilde(n_elem);
                  mutilde_demi(n_elem)/=theta;
                }
            }
          else if (type_systeme_naire_==1)
            {
              for (int j=0; j<nb_equation_CH; j++)
                {
                  for(int n_elem=0; n_elem<nb_elem; n_elem++)
                    {
                      // Commente par DJ
                      //----------------
                      //               c_demi(n_elem)=x1(n_elem);
                      //----------------
                      c_demi(n_elem,j)-=(1-theta)*c(n_elem,j);
                      c_demi(n_elem,j)/=theta;

                      // Commente par DJ
                      //----------------
                      //               mutilde_demi(n_elem)=x1(n_elem+nb_elem);
                      //----------------
                      mutilde_demi(n_elem,j)-=(1-theta)*mutilde(n_elem,j);
                      mutilde_demi(n_elem,j)/=theta;
                    }
                }
            }
          //Cerr <<"c_demi after 1-theta/theta = "<<c_demi<<finl;
          //Cerr <<"mutilde_demi after 1-theta/theta = "<<mutilde_demi<<finl;

        }

      // ATTENTION : A VERIFIER
      //=======================
      //       c_demi.echange_espace_virtuel();
      //       mutilde_demi.echange_espace_virtuel();
      //       c.echange_espace_virtuel();
      //       mutilde.echange_espace_virtuel();

      // Mise a jour

      if (couplage_==0)
        {
          // Si traitement explicite de mutilde dans le cas implicite
          calculer_mutilde(mutilde);
        }
      else if (couplage_==1)
        {
          // Si traitement implicite de mutilde dans le cas implicite
          mutilde=mutilde_demi;
        }
      mutilde.echange_espace_virtuel();

      // Calcul de l'accroissement entre n et n+1/2
      // Utile pour pouvoir utiliser n'importe quel dt
      accr=0.;
    }
  else if(implicitation_==0)
    {
      // Cerr<<"Schema explicite"<<finl;
      calculer_u2_elem(u_carre_);
      calculer_mutilde(mutilde);
      mutilde.echange_espace_virtuel();

      DoubleTab& prov_elem = prov_elem_;
      if (prov_elem.size()==0)
        prov_elem=mutilde;

      if (type_systeme_naire_==0)
        {
          if(kappa_ind==1)
            {
              DoubleTab kappa_var(prov_elem);
              kappa_var=0.;
              for(int ikappa=0; ikappa<nb_elem; ikappa++)
                kappa_var(ikappa)=(this->*kappa_func_c)(c(ikappa));
              // Div(Kappa*Grad(mutilde))
              div_kappa_grad(mutilde, kappa_var, prov_elem);
            }
          else if (kappa_ind==0)
            {
              // Kappa*Laplacien(mutilde)
              laplacien(mutilde,prov_elem);
              prov_elem*=kappa;
            }
        }
      else if (type_systeme_naire_==1)
        {
          if (type_kappa_auto_diffusion_==0)
            {
              DoubleTab temp_mutilde(mutilde.dimension(0),1);
              temp_mutilde=0;
              DoubleTab temp_prov_elem= temp_mutilde;

              // laplacien(mutilde)
              laplacien(mutilde,prov_elem);
              // Cerr << "laplacien mutilde " << prov_elem<<finl;

              //kappa*laplacien(mutilde)
              DoubleTab temp_prov_elem2 = prov_elem;
              temp_prov_elem2 = 0;
              const int nb_comp =eq_c.constituant().nb_constituants();
              for (int j=0; j<prov_elem.line_size(); j++)
                {
                  for (int i=0; i<prov_elem.dimension(0); i++)
                    {
                      temp_prov_elem(i,0)=prov_elem(i,j);
                      for (int k=0; k<nb_comp; k++)
                        {
                          temp_prov_elem2(i,k)+=temp_prov_elem(i,0)*kappaMatrix(j+k*nb_comp);
                        }
                    }
                }
              prov_elem = temp_prov_elem2;
              //Cerr << "kappaMatrix * laplacien(mutilde) " << prov_elem<<finl;
            }
          else if (type_kappa_auto_diffusion_==1)
            {
              DoubleTab kappaMatrix_variable(prov_elem.dimension(0),c.line_size()*c.line_size()); //or (nb_elem,nb_comp*nb_comp)
              //kappaMatrix_variable=0.;
              kappaMatrix_variable = (this->*kappaMatrix_func_c)(c, coeff_auto_diffusion);
              double R=8.314;
              kappaMatrix_variable *= molarVolume/(R*temperature);
              //Cerr <<"kappaMatrix_variable (1/RT)"<<kappaMatrix_variable<<finl;

              //Div_kappa_grad
              div_kappa_grad(mutilde, kappaMatrix_variable, prov_elem);
              //Cerr <<"Div.kappaMatrix_variable *(Vm/RT).grad(mutilde)"<<prov_elem<<finl;

            }
        }

      // Pour equation Allen-Cahn (kappa constant)
      //------------------------------------------
      //   prov_elem=F;
      //   prov_elem*=-kappa;

      accr+=prov_elem;
      //eq_c.solv_masse().appliquer(prov_elem);
      //Cerr<<"masse de la force sur c "<<prov_elem.max_abs()<<finl;

      // Utile pour pouvoir utiliser n'importe quel dt - voir Schema_Phase_Field
      c_demi=c;
    }
  else
    {
      Cerr<<"Type de schema errone !!"<<finl;
    }
}

/*! @brief Calcul de Div(alpha*rho*Grad((C)) au centre des elements
 *
 * @param (DoubleTab& div_alpha_gradC) Div(alpha*Grad((C)) au centre des elements
 */
void Source_Con_Phase_field::calculer_div_alpha_gradC(DoubleTab& div_alpha_gradC) const
{
  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const DoubleTab& c=eq_c.inconnue().valeurs();
  const Navier_Stokes_std& eq_ns=ref_cast(Navier_Stokes_std,le_probleme2->equation(0));
  const Operateur_Grad& opgrad=eq_ns.operateur_gradient();
  const Operateur_Div& opdiv= eq_ns.operateur_divergence();

  if (type_systeme_naire_==0)
    {
      DoubleTab& prov_face=ref_cast_non_const( DoubleTab, prov_face_);
      if (prov_face.size()==0)
        prov_face=eq_ns.inconnue().valeurs();
      prov_face=0.;

      // Grad(C)
      opgrad.calculer(c,prov_face);
      eq_ns.solv_masse().appliquer(prov_face);
      prov_face *= alpha;

      // Div(alpha*Grad(c))
      opdiv.calculer(prov_face,div_alpha_gradC);
      eq_c.solv_masse().appliquer(div_alpha_gradC);
    }
  else if (type_systeme_naire_==1)
    {
      DoubleTab& temp_prov_face= ref_cast_non_const(DoubleTab,prov_face_);
      if (temp_prov_face.size()==0)
        temp_prov_face=eq_ns.inconnue().valeurs();
      temp_prov_face=0.;

      DoubleTab prov_face(temp_prov_face.dimension_tot(0),c.line_size());
      prov_face = 0;
      DoubleTab temp_c(c.dimension_tot(0),1);
      temp_c=0;
      /*
            for (int j=0; j<c.line_size(); j++)
              {
                for (int i=0; i<c.dimension(0); i++)
                  {
                    temp_c(i,0)=c(i,j);
                    // Grad(C)
                    opgrad.calculer(temp_c,temp_prov_face);
                    for (int k=0; k<temp_prov_face.dimension(0); k++)
                      {
                        prov_face(k,j)=temp_prov_face(k,0);
                      }
                  }
              }*/

      for (int j=0; j<c.line_size(); j++)
        {
          for (int i=0; i<temp_c.dimension_tot(0); i++)
            {
              temp_c(i,0)=c(i,j);
            }
          // Grad(C)

          opgrad.calculer(temp_c,temp_prov_face);

          for (int k=0; k<temp_prov_face.dimension_tot(0); k++)
            {
              prov_face(k,j)=temp_prov_face(k,0);
            }
        }
      //eq_ns.solv_masse().appliquer(prov_face); a revoir..

      // alpha*Grad(C)
      DoubleTab temp_prov_face2 = prov_face;
      temp_prov_face2 = 0;
      const int nb_comp =eq_c.constituant().nb_constituants(); //nb_equation_CH
      for (int j=0; j<prov_face.line_size(); j++)
        {
          for (int i=0; i<prov_face.dimension(0); i++)
            {
              temp_prov_face(i,0)=prov_face(i,j);
              for (int k=0; k<nb_comp; k++)
                {
                  temp_prov_face2(i,k)+=temp_prov_face(i,0)*alphaMatrix(j+k*nb_comp);
                }
            }
        }
      prov_face = temp_prov_face2;
      //Cerr <<"alphaMatrix * (gradC apres solv_masse)"<< prov_face<<finl;

      // Div(alpha*Grad(c))
      opdiv.calculer(prov_face,div_alpha_gradC);
      eq_c.solv_masse().appliquer(div_alpha_gradC);
      //Cerr << "div_alpha_gradC apres solv_masse.appliquer "<<div_alpha_gradC<<finl;

    }
}

/*! @brief Calcul de Div(alpha*rho*Grad((C)) au centre des elements
 *
 * @param (DoubleTab& div_alpha_rho_gradC) Div(alpha*rho*Grad((C)) au centre des elements
 */
void Source_Con_Phase_field::calculer_div_alpha_rho_gradC(DoubleTab& div_alpha_rho_gradC) const
{
  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const DoubleTab& c=eq_c.inconnue().valeurs();
  const Navier_Stokes_phase_field& eq_ns=ref_cast(Navier_Stokes_phase_field,le_probleme2->equation(0));
  const Operateur_Grad& opgrad=eq_ns.operateur_gradient();
  const Operateur_Div& opdiv= eq_ns.operateur_divergence();

  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const IntTab& face_voisins = zone_VDF.face_voisins();
  const DoubleVect& volumes = zone_VDF.volumes();
  const int ndeb=zone_VDF.premiere_face_int();
  const int nbfaces=zone_VDF.nb_faces();

  int el0,el1;
  double vol0,vol1;

  DoubleTab& prov_face=ref_cast_non_const( DoubleTab, prov_face_);
  if (prov_face.size()==0)
    prov_face=eq_ns.inconnue().valeurs();
  prov_face=0.;

  // Grad(C)
  opgrad.calculer(c,prov_face);
  eq_ns.solv_masse().appliquer(prov_face);

  const DoubleTab rhoPF=eq_ns.rho().valeurs();
  double rho_face;

  if (boussi_==1)
    {
      prov_face *= alpha*rho0; // Cas approximation de Boussinesq
    }
  else if (boussi_==0)
    {
      for (int fac=ndeb; fac<nbfaces; fac++)
        {
          el0=face_voisins(fac,0);
          el1=face_voisins(fac,1);
          vol0=volumes(el0);
          vol1=volumes(el1);

          rho_face=(vol0*rhoPF(el0)+vol1*rhoPF(el1))/(vol0+vol1);
          prov_face(fac) *= alpha*rho_face;
        }
    }

  // Div(alpha*rho*Grad(c))
  opdiv.calculer(prov_face,div_alpha_rho_gradC);
  eq_c.solv_masse().appliquer(div_alpha_rho_gradC);
}


/*! @brief Assemble la matrice pour le calcul du point fixe
 *
 */
void Source_Con_Phase_field::assembler_matrice_point_fixe(Matrice_Morse& matrice_diffusion_CH)
{
  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const DoubleTab& c=eq_c.inconnue().valeurs();

  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const IntTab& face_voisins = zone_VDF.face_voisins();
  const IntTab& elem_faces = zone_VDF.elem_faces();
  const DoubleTab& positions = zone_VDF.xp();
  const IntVect& ori = zone_VDF.orientation();

  /* Cerr <<"face_ voisins "<< face_voisins<<finl;
   Cerr <<"elem_faces "<< elem_faces<<finl;*/


  int compt=0;
  int compt1=0;
  int compt2=0;

  const int nb_elem = zone_VDF.nb_elem_tot();
  int nb_compo_;
  int f0;
  int min_tri;
  int old_tri;
  int voisin;
  int nb_faces_au_bord;
  int dimensionnement ;
  double dvar2=0.;
  double dvarkeep=0.;
  double dt;
  double valeur_diag=0.;
  const double theta=0.6;
  double kappa_interpolee=0.;
  double kappa_naire_interpolee=0.;

  //   Cerr<<"======================================"<<finl;
  //   Cerr<<"Assemblage de la matrice du point fixe"<<finl;

  // Dimension du pb
  nb_compo_=dimension;

  // Forme de la matrice : ( I A )
  //                       ( B I )

  // Assemblage de la moitie superieure de la matrice : I et A

  // On dimensionne la matrice A et donc en meme temps B (meme nombre de termes non nuls)
  // Pour cela : dimension du vecteur coeff = (2 * dim + 1) * nb_elem - nombre_de_bords

  if (type_systeme_naire_==0)
    {
      DoubleTab kappa_var(c);

      // Cerr << "kappa_ind = " << kappa_ind << finl;
      if(kappa_ind==0)
        {
          // Cerr << "---" << finl;
          // Cerr << "kappa constant pour le point fixe" << finl;
          // Cerr << "---" << finl;
          kappa_var=kappa;
        }
      else
        {
          for(int ikappa=0; ikappa<nb_elem; ikappa++)
            {
              // Cerr << "---" << finl;
              // Cerr << "Calcul de kappa variable pour le point fixe" << finl;
              // Cerr << "---" << finl;
              kappa_var(ikappa)=(this->*kappa_func_c)(c(ikappa));
            }
        }
      //Cerr <<"kappa_var = " << kappa_var<<finl;

      dimensionnement=(2*nb_compo_+1)*nb_elem;

      for(int elem=0; elem<nb_elem; elem++)
        {
          for(int ncomp=0; ncomp<2*nb_compo_; ncomp++)
            {
              f0 = elem_faces(elem,ncomp);
              voisin=face_voisins(f0,0);
              if (face_voisins(f0,0)==elem) voisin=face_voisins(f0,1);

              if (voisin==-1) dimensionnement--;
            }
        }
      // On ajoute le nombre de 1 venant de I, ceci valant pour I et A. Le dimensionnement total est donc le double
      dimensionnement+=nb_elem;
      dimensionnement*=2;

      // Allocation des tableaux specifiques a la matrice morse
      matrice_diffusion_CH.dimensionner(2*nb_elem,dimensionnement);
      DoubleVect& coeff=matrice_diffusion_CH.get_set_coeff();
      IntVect& tab2=matrice_diffusion_CH.get_set_tab2();
      IntVect& tab1=matrice_diffusion_CH.get_set_tab1();
      //Cerr <<"matrice_diffusion_CH = "<<matrice_diffusion_CH<<finl;

      // Boucle sur le nombre d'elements
      for(int elem=0; elem<nb_elem; elem++)
        {
          valeur_diag=0.;
          old_tri=-1;
          nb_faces_au_bord=0;

          // On ajoute le 1 de la diagonale de la sous-matrice I du haut a gauche
          coeff(compt)=1;
          tab2(compt2)=elem;
          compt++;
          compt2++;
          tab1(compt1)=compt2;
          compt1++;

          // Calcul du nombre de bords
          for(int ncomp=0; ncomp<2*nb_compo_; ncomp++)
            {
              f0 = elem_faces(elem,ncomp);
              voisin=face_voisins(f0,0);
              if (face_voisins(f0,0)==elem) voisin=face_voisins(f0,1);

              // On connait le voisin associe a la face en cours. On regarde s'il est au bord
              if (voisin==-1) nb_faces_au_bord++;
            }

          // Boucle sur les faces
          for(int ncomp=0; ncomp<2*nb_compo_-nb_faces_au_bord; ncomp++)
            {

              min_tri=nb_elem+1;
              dt=eq_c.schema_temps().pas_de_temps();    // on fait un calcul sur un demi pas de temps

              for(int ncomp_tri=0; ncomp_tri<2*nb_compo_; ncomp_tri++)
                {
                  f0 = elem_faces(elem,ncomp_tri);

                  // Voisin associe a la face - On s'assure de traiter un voisin et pas l'element resident
                  voisin=face_voisins(f0,0);
                  if (face_voisins(f0,0)==elem) voisin=face_voisins(f0,1);

                  // On va compter le nombre de voisins existants (faces non aux bords de l'element elem)
                  // Ceci sert a calculer la contribution au terme relatif a l'element elem
                  // Cas kappa variable : on utilise la moyenne geometrique des valeurs des mobilites
                  // des elements concernes par la face de calcul ("kappa_interpolee")
                  if (voisin!=-1 && old_tri==-1)
                    {
                      dvar2=pow((positions(elem,ori(f0))-positions(voisin,ori(f0))),2);
                      if(kappa_moy_==2)
                        {
                          kappa_interpolee=pow(dabs(kappa_var(elem)*kappa_var(voisin)),0.5);    // Moyenne geometrique
                        }
                      else if(kappa_moy_==1)
                        {
                          if (kappa_var(elem)==0 || kappa_var(voisin)==0)
                            kappa_interpolee=0;
                          else
                            kappa_interpolee=2./(1./kappa_var(elem)+1./kappa_var(voisin));  // Moyenne harmonique
                        }
                      else if(kappa_moy_==0)
                        {
                          kappa_interpolee=(kappa_var(elem)+kappa_var(voisin))/2.;        // Moyenne arithmetique
                        }
                      valeur_diag+=theta*kappa_interpolee*(dt/dvar2);
                    }

                  // Si le voisin a un numero plus petit que l'ancien minimum, on ne le prend pas en compte
                  // En effet, le nouveau minimum est un numero d'element qui lui sera superieur
                  // Rq : voisin>old_tri => voisin != -1 donc pour une face au bord le if est false
                  if (voisin>old_tri)
                    {
                      min_tri=min(min_tri,voisin);
                      if(min_tri==voisin)
                        {
                          dvarkeep=pow((positions(elem,ori(f0))-positions(voisin,ori(f0))),2);
                        }
                    }
                }

              //Cerr <<" min_tri = min(min_tri,voisin) if voisin>old_tri ====== "<<min_tri<<finl;
              // Si l'ancien minimum est inferieur a elem et le nouveau superieur, il faut traiter elem juste avant le nouveau voisin
              if (old_tri<elem && min_tri>elem)
                {
                  coeff(compt)=valeur_diag;
                  tab2(compt2)=elem+nb_elem;
                  compt++;
                  compt2++;
                }

              // On complete les tableaux avec les valeurs dans les voisins non aux bords du domaine
              if (kappa_moy_==2)
                {
                  kappa_interpolee=pow(dabs(kappa_var(elem)*kappa_var(min_tri)),0.5);    // Moyenne geometrique
                }
              else if(kappa_moy_==1)
                {
                  if (kappa_var(elem)==0 || kappa_var(min_tri)==0)
                    kappa_interpolee=0;
                  else
                    kappa_interpolee=2./(1./kappa_var(elem)+1./kappa_var(min_tri));  // Moyenne harmonique
                }
              else if(kappa_moy_==0)
                {
                  kappa_interpolee=(kappa_var(elem)+kappa_var(min_tri))/2.;        // Moyenne arithmetique
                }
              coeff(compt)=-theta*kappa_interpolee*(dt/dvarkeep);
              tab2(compt2)=min_tri+nb_elem;
              compt++;
              compt2++;

              // Si on est a la fin et que elem n'a pas encore ete mis dans coeff, il faut le faire avant la fin de la boucle sur les faces de l'element
              // Ceci arrive si elem > tous les numeros d'elements de ses voisins
              if (ncomp==2*nb_compo_-nb_faces_au_bord-1 && min_tri<elem)
                {
                  coeff(compt)=valeur_diag;
                  tab2(compt2)=elem+nb_elem;
                  compt++;
                  compt2++;
                }

              // On sauve le numero d'element "minimum"
              old_tri=min_tri;
              //Cerr <<"old_tri = "<<old_tri<<finl;
              //Cerr <<"min_tri fin = "<< min_tri<<finl;
            }

        }

      // Assemblage de la moitie superieure de la matrice : B et I

      // Boucle sur le nombre d'elements
      for(int elem=0; elem<nb_elem; elem++)
        {
          valeur_diag=0;
          old_tri=-1;
          nb_faces_au_bord=0;

          // Calcul du nombre de bords
          for(int ncomp=0; ncomp<2*nb_compo_; ncomp++)
            {
              f0 = elem_faces(elem,ncomp);
              voisin=face_voisins(f0,0);
              if (face_voisins(f0,0)==elem) voisin=face_voisins(f0,1);

              // On connait le voisin associe a la face en cours. On regarde s'il est au bord
              if (voisin==-1) nb_faces_au_bord++;
            }

          // Boucle sur les faces
          for(int ncomp=0; ncomp<2*nb_compo_-nb_faces_au_bord; ncomp++)
            {

              min_tri=nb_elem+1;

              for(int ncomp_tri=0; ncomp_tri<2*nb_compo_; ncomp_tri++)
                {
                  f0 = elem_faces(elem,ncomp_tri);

                  // Voisin associe a la face - On s'assure de traiter un voisin et pas l'element resident
                  voisin=face_voisins(f0,0);
                  if (face_voisins(f0,0)==elem) voisin=face_voisins(f0,1);

                  // On va compter le nombre de voisins existants (faces non aux bords de l'element elem)
                  // Ceci sert a calculer la contribution au terme relatif a l'element elem
                  if (voisin!=-1 && old_tri==-1)
                    {
                      dvar2=pow((positions(elem,ori(f0))-positions(voisin,ori(f0))),2);
                      valeur_diag+=-alpha/dvar2;
                    }

                  // Si le voisin a un numero plus petit que l'ancien minimum, on ne le prend pas en compte
                  // En effet, le nouveau minimum est un numero d'element qui lui sera superieur
                  // Rq : voisin>old_tri => voisin != -1 donc pour une face au bord le if est false
                  if (voisin>old_tri)
                    {
                      min_tri=min(min_tri,voisin);
                      if(min_tri==voisin)
                        {
                          dvarkeep=pow((positions(elem,ori(f0))-positions(voisin,ori(f0))),2);
                        }
                    }

                }

              // Si l'ancien minimum est inferieur a elem et le nouveau superieur, il faut traiter elem juste avant le nouveau voisin
              if (old_tri<elem && min_tri>elem)
                {
                  coeff(compt)=valeur_diag;
                  tab2(compt2)=elem;
                  compt++;
                  compt2++;
                  if (old_tri==-1)
                    {
                      // On ajoute un terme a tab1 dans le cas ou l'element diagonal est le premier terme
                      tab1(compt1)=compt2;
                      compt1++;
                      old_tri=elem;
                    }
                }

              // On complete les tableaux avec les valeurs dans les voisins non aux bords du domaine
              coeff(compt)=alpha/dvarkeep;
              tab2(compt2)=min_tri;
              compt++;
              compt2++;
              if (old_tri==-1)
                {
                  // On ajoute un terme a tab1 dans le cas ou le premier terme est extradiagonal
                  tab1(compt1)=compt2;
                  compt1++;
                }

              // Si on est a la fin et que elem n'a pas encore ete mis dans coeff, il faut le faire avant la fin de la boucle sur les faces de l'element
              // Ceci arrive si elem > tous les numeros d'elements de ses voisins
              // Bien sur dans ce cas la, l'element diagonal est forcement le dernier mais pas le premier - il n'existe pas de maillage avec des mailles seules sans voisins
              if (ncomp==2*nb_compo_-nb_faces_au_bord-1 && min_tri<elem)
                {
                  coeff(compt)=valeur_diag;
                  tab2(compt2)=elem;
                  compt++;
                  compt2++;
                }

              // On sauve le numero d'element "minimum"
              old_tri=min_tri;

            }

          // On ajoute le 1 de la diagonale de la sous-matrice I du bas a droite
          coeff(compt)=1;
          tab2(compt2)=elem+nb_elem;
          compt++;
          compt2++;


        }

      // cf Mat_Morse.h/.cpp, tab1() et tab2() sont a utiliser au sens Fortran
      // EXPLICATION :
      // coeff stocke les valeurs des coefficients
      // tab2 stocke les numeros de colonnes dans la matrice de ces coefficients, au sens FORTRAN (i-e 1 <= tab2[i] <= n)
      // tab1 stocke le rang de l'element de tab2() pour lequel on change de ligne. Ainsi pour tab1[j], on est a la i-eme ligne avec :
      //      ( tab1[i] <= j < tab1[i+1] ) qui implique qu'on considere toujours la i-eme ligne,
      //      et des que cela n'est plus respecte, le passage a la ligne suivante est effectue

      // De par l'algorithme, tab1() est deja au sens FORTRAN. Reste a le faire pour tab2()
      // De plus on doit mettre le dernier terme de tab1() a dimensionnement+1
      tab2+=1;
      tab1(compt1)=dimensionnement+1;
      compt1+=1;

      if (compt!=dimensionnement)
        {
          Cerr << "Erreur lors du calcul de la matrice du point fixe : nombre d'elements non nuls calcules different du nombre d'elements non nuls prevus" << finl;
          exit();
        }
      //   Cerr<<"Nombre d'elements non nuls calcules="<<compt<<finl;
      //   Cerr<<"Nombre d'elements non nuls prevus="<<dimensionnement<<finl;
      //   Cerr<<"Assemblage de la matrice du point fixe : OK"<<finl;
      //   Cerr<<"==========================================="<<finl;

      // Test des tableaux
      Cerr<<"coeff="<<coeff<<finl;
      Cerr<<"tab1="<<tab1<<finl;
      Cerr<<"tab2="<<tab2<<finl;

    }
  else if (type_systeme_naire_==1)
    {
      DoubleTab kappa_Matrix_var(nb_elem,nb_equation_CH*nb_equation_CH);
      // Cerr << "kappa_ind = " << kappa_ind << finl;
      if(kappa_ind==0)
        {
          // Cerr << "---" << finl;
          // Cerr << "kappa constant pour le point fixe" << finl;
          // Cerr << "---" << finl;
          for (int i=0; i<nb_elem; i++)
            {
              for (int j=0; j<nb_equation_CH*nb_equation_CH; j++)
                {
                  kappa_Matrix_var(i,j)=kappaMatrix(j);
                }
            }
        }
      else
        {
          //for(int ikappa=0; ikappa<nb_elem; ikappa++)
          //{
          // Cerr << "---" << finl;
          // Cerr << "Calcul de kappa variable pour le point fixe" << finl;
          // Cerr << "---" << finl;
          kappa_Matrix_var= (this->*kappaMatrix_func_c)(c, coeff_auto_diffusion);
          double R=8.314;
          kappa_Matrix_var *= molarVolume/(R*temperature);
          //kappa_var(ikappa)=(this->*kappa_func_c)(c(ikappa));
          //}
        }
      //Cerr <<"kappaMatrix_Var = "<< kappa_Matrix_var<< finl;

      // Forme de la matrice pour un systeme -naire :
      // ici on montre uniquement pour un systeme quaternaire
      //				  ( 1 0 0 A B C )
      //                  ( 0 1 0 D E F )
      //				  ( 0 0 1 G H J )
      //                  ( K L M 1 0 0 )
      //				  ( N P Q 0 1 0 )
      //                  ( U V W 0 0 1 )

      // Assemblage de la moitie superieure de la matrice : 1 0 0 A B C

      // On dimensionne la matrice A et donc en meme temps toutes les matrices non unitaire B, C, ...(meme nombre de termes non nuls)
      // Pour cela : dimension du vecteur coeff pour chaque matrice A = (2 * dim + 1) * nb_elem - nombre_de_bords

      dimensionnement=(2*nb_compo_+1)*nb_elem;

      for(int elem=0; elem<nb_elem; elem++)
        {
          for(int ncomp=0; ncomp<2*nb_compo_; ncomp++)
            {
              f0 = elem_faces(elem,ncomp);
              voisin=face_voisins(f0,0);
              if (face_voisins(f0,0)==elem) voisin=face_voisins(f0,1);

              if (voisin==-1) dimensionnement--;
            }
        }
      // Sur la meme ligne, on a 1 matrice unitaire et nb_equation_CH* matrice non_unitaire et non nul. Pour cela,
      // On multiplie donc le nombre d'elements non nuls de la matrice A par le nb_equation_CH pour avoir le nombre d'element nnz sur la meme ligne
      // Puis on ajoute le nombre de 1 venant de la matrice unitaire sur la premiere ligne.
      // dimensionnement est le nombre d'element non nuls (nnz) dans la grosse matrice_diffusion_CH
      // Le dimensionnement total est donc le nnz de la premiere ligne multiplié par 2*nb_equation_CH
      dimensionnement*=nb_equation_CH;
      dimensionnement+=nb_elem;
      dimensionnement*=2*nb_equation_CH;

      // Allocation des tableaux specifiques a la matrice morse
      matrice_diffusion_CH.dimensionner(nb_equation_CH*(2*nb_elem),dimensionnement);
      DoubleVect& coeff=matrice_diffusion_CH.get_set_coeff();
      IntVect& tab2=matrice_diffusion_CH.get_set_tab2();
      IntVect& tab1=matrice_diffusion_CH.get_set_tab1();
      //Cerr <<"matrice_diffusion_CH = "<<matrice_diffusion_CH<<finl;


      // Boucle sur le nombre d'equation pour la generalisation
      // for (int ligne=0; ligne<nb_equation_CH; ligne++)
      //	{
      // en fonction de ligne, la matrice unite va commencer en (ligne*nb_elem)
      //	}


      DoubleVect valeur_diag_naire(nb_equation_CH);
      //double valeur_diag_naire;

      // Boucle sur le nombre de matrice dans la grosse matrice (matrice d'une matrice)
      for (int ligne=0; ligne<nb_equation_CH; ligne++)
        {
          // Boucle sur le nombre d'elements
          for(int elem=0; elem<nb_elem; elem++)
            {
              valeur_diag_naire=0.;
              old_tri=-1;
              nb_faces_au_bord=0;

              // On ajoute le 1 de la diagonale de la sous-matrice I du haut a gauche
              coeff(compt)=1;
              tab2(compt2)=(ligne*nb_elem)+elem;
              compt++;
              compt2++;
              tab1(compt1)=compt2;
              compt1++;

              // Boucle sur le nombre de matrice non nulle correspondant a mutilde dans second_membre
              //for (int mat=0; mat<nb_equation_CH; mat++)
              //{
              //int j=(ligne*nb_equation_CH)+mat;

              // Calcul du nombre de bords
              for(int ncomp=0; ncomp<2*nb_compo_; ncomp++)
                {
                  f0 = elem_faces(elem,ncomp);
                  voisin=face_voisins(f0,0);
                  if (face_voisins(f0,0)==elem) voisin=face_voisins(f0,1);

                  // On connait le voisin associe a la face en cours. On regarde s'il est au bord
                  if (voisin==-1) nb_faces_au_bord++;
                }

              //for (int mat=0; mat<nb_equation_CH; mat++)
              //  {

              // Boucle sur les faces
              //for (int mat=0; mat<nb_equation_CH; mat++)
              //{
              //int j=(ligne*nb_equation_CH)+mat;
              //min_tri=nb_elem+1; //initialisation de min_tri pour chaque face (une sorte de reference car le numero des voisins ne doit pas depasser nb_elem+1)

              for(int ncomp=0; ncomp<2*nb_compo_-nb_faces_au_bord; ncomp++)
                {
                  min_tri=nb_elem+1; //initialisation de min_tri pour chaque face (une sorte de reference car le numero des voisins ne doit pas depasser nb_elem+1)
                  dt=eq_c.schema_temps().pas_de_temps();    // on fait un calcul sur un demi pas de temps
                  for (int mat=0; mat<nb_equation_CH; mat++)
                    {
                      int j=(ligne*nb_equation_CH)+mat;
                      for(int ncomp_tri=0; ncomp_tri<2*nb_compo_; ncomp_tri++)
                        {
                          f0 = elem_faces(elem,ncomp_tri);

                          // Voisin associe a la face - On s'assure de traiter un voisin et pas l'element resident
                          voisin=face_voisins(f0,0);
                          if (face_voisins(f0,0)==elem) voisin=face_voisins(f0,1);

                          // On va compter le nombre de voisins existants (faces non aux bords de l'element elem)
                          // Ceci sert a calculer la contribution au terme relatif a l'element elem
                          // Cas kappa variable : on utilise la moyenne (geometrique, harmonique ou arithmetique) des valeurs des mobilites
                          // des elements concernes par la face de calcul ("kappa_naire_interpolee")

                          //Cerr <<"f0 = "<<f0<<finl;

                          if (voisin!=-1 && old_tri==-1)
                            {
                              dvar2=pow((positions(elem,ori(f0))-positions(voisin,ori(f0))),2);
                              if(kappa_moy_==2)
                                {
                                  kappa_naire_interpolee=pow(dabs(kappa_Matrix_var(elem,j)*kappa_Matrix_var(voisin,j)),0.5);    // Moyenne geometrique
                                }
                              else if(kappa_moy_==1)
                                {
                                  if (kappa_Matrix_var(elem,j)==0 || kappa_Matrix_var(voisin,j)==0)
                                    kappa_naire_interpolee=0;
                                  else
                                    kappa_naire_interpolee=2./(1./kappa_Matrix_var(elem,j)+1./kappa_Matrix_var(voisin,j));  // Moyenne harmonique
                                }
                              else if(kappa_moy_==0)
                                {
                                  kappa_naire_interpolee=(kappa_Matrix_var(elem,j)+kappa_Matrix_var(voisin,j))/2.;        // Moyenne arithmetique
                                }
                              valeur_diag_naire(mat)+=theta*kappa_naire_interpolee*(dt/dvar2);
                            }

                          // Si le voisin a un numero plus petit que l'ancien minimum, on ne le prend pas en compte
                          // En effet, le nouveau minimum est un numero d'element qui lui sera superieur
                          // Rq : voisin>old_tri => voisin != -1 donc pour une face au bord le if est false
                          if (voisin>old_tri)
                            {
                              min_tri=min(min_tri,voisin);
                              if(min_tri==voisin)
                                {
                                  dvarkeep=pow((positions(elem,ori(f0))-positions(voisin,ori(f0))),2);
                                }
                            }
                        }

                      // Si l'ancien minimum est inferieur a elem et le nouveau superieur, il faut traiter elem juste avant le nouveau voisin
                      if (old_tri<elem && min_tri>elem)
                        {
                          coeff(compt)=valeur_diag_naire(mat);
                          tab2(compt2)=nb_elem*(mat+nb_equation_CH)+elem; //la ou commence les matrices A, B, C, etc... donc (nb_elem*nb_equation_CH)+1 ??
                          compt++;
                          compt2++;
                        }

                      // On complete les tableaux avec les valeurs dans les voisins non aux bords du domaine
                      if (kappa_moy_==2)
                        {
                          kappa_naire_interpolee=pow(dabs(kappa_Matrix_var(elem,j)*kappa_Matrix_var(min_tri,j)),0.5);    // Moyenne geometrique
                        }
                      else if(kappa_moy_==1)
                        {
                          if (kappa_Matrix_var(elem,j)==0 || kappa_Matrix_var(min_tri,j)==0)
                            kappa_naire_interpolee=0;
                          else
                            kappa_naire_interpolee=2./(1./kappa_Matrix_var(elem,j)+1./kappa_Matrix_var(min_tri,j));  // Moyenne harmonique
                        }
                      else if(kappa_moy_==0)
                        {
                          kappa_naire_interpolee=(kappa_Matrix_var(elem,j)+kappa_Matrix_var(min_tri,j))/2.;        // Moyenne arithmetique
                        }
                      coeff(compt)=-theta*kappa_naire_interpolee*(dt/dvarkeep);
                      tab2(compt2)=min_tri+nb_elem*(nb_equation_CH+mat);
                      compt++;
                      compt2++;

                      // Si on est a la fin et que elem n'a pas encore ete mis dans coeff, il faut le faire avant la fin de la boucle sur les faces de l'element
                      // Ceci arrive si elem > tous les numeros d'elements de ses voisins
                      if (ncomp==2*nb_compo_-nb_faces_au_bord-1 && min_tri<elem)
                        {
                          coeff(compt)=valeur_diag_naire(mat);
                          tab2(compt2)=elem+(nb_equation_CH+mat)*nb_elem;
                          compt++;
                          compt2++;
                        }
                      //Cerr <<"kappa_naire_interpolee = "<<kappa_naire_interpolee<<finl;
                      // On sauve le numero d'element "minimum"
                    }
                  old_tri=min_tri;
                  //}

                }
            }

        }

      // Test des tableaux
      /*Cerr<<"coeff="<<coeff<<finl;
      Cerr<<"tab1="<<tab1<<finl;
      Cerr<<"tab2="<<tab2<<finl;
      //Cerr <<"matrice_diffusion_CH final = "<<matrice_diffusion_CH<<finl;*/

      // Assemblage de la moitie superieure de la matrice : B et I
      // Boucle sur le nombre de matrice dans la grosse matrice (matrice d'une matrice)
      for (int ligne=0; ligne<nb_equation_CH; ligne++)
        {
          // Boucle sur le nombre d'elements
          for(int elem=0; elem<nb_elem; elem++)
            {
              valeur_diag_naire=0.;
              old_tri=-1;
              nb_faces_au_bord=0;

              // On saute de ligne en prenant la valeur de compte2 -1
              tab1(compt1)=compt2+1;
              compt1++;


              // Boucle sur le nombre de matrice non nulle correspondant a mutilde dans second_membre
              //for (int mat=0; mat<nb_equation_CH; mat++)
              //{
              //int j=(ligne*nb_equation_CH)+mat;

              // Calcul du nombre de bords
              for(int ncomp=0; ncomp<2*nb_compo_; ncomp++)
                {
                  f0 = elem_faces(elem,ncomp);
                  voisin=face_voisins(f0,0);
                  if (face_voisins(f0,0)==elem) voisin=face_voisins(f0,1);

                  // On connait le voisin associe a la face en cours. On regarde s'il est au bord
                  if (voisin==-1) nb_faces_au_bord++;
                }

              //for (int mat=0; mat<nb_equation_CH; mat++)
              //  {

              // Boucle sur les faces
              //for (int mat=0; mat<nb_equation_CH; mat++)
              //{
              //int j=(ligne*nb_equation_CH)+mat;
              //min_tri=nb_elem+1; //initialisation de min_tri pour chaque face (une sorte de reference car le numero des voisins ne doit pas depasser nb_elem+1)

              for(int ncomp=0; ncomp<2*nb_compo_-nb_faces_au_bord; ncomp++)
                {
                  min_tri=nb_elem+1; //initialisation de min_tri pour chaque face (une sorte de reference car le numero des voisins ne doit pas depasser nb_elem+1)
                  for (int mat=0; mat<nb_equation_CH; mat++)
                    {
                      int j=(ligne*nb_equation_CH)+mat;
                      for(int ncomp_tri=0; ncomp_tri<2*nb_compo_; ncomp_tri++)
                        {
                          f0 = elem_faces(elem,ncomp_tri);

                          // Voisin associe a la face - On s'assure de traiter un voisin et pas l'element resident
                          voisin=face_voisins(f0,0);
                          if (face_voisins(f0,0)==elem) voisin=face_voisins(f0,1);

                          // On va compter le nombre de voisins existants (faces non aux bords de l'element elem)
                          // Ceci sert a calculer la contribution au terme relatif a l'element elem
                          // Cas kappa variable : on utilise la moyenne (geometrique, harmonique ou arithmetique) des valeurs des mobilites
                          // des elements concernes par la face de calcul ("kappa_naire_interpolee")

                          if (voisin!=-1 && old_tri==-1)
                            {
                              dvar2=pow((positions(elem,ori(f0))-positions(voisin,ori(f0))),2);
                              valeur_diag_naire(mat)+=-alphaMatrix(j)/dvar2;
                            }

                          // Si le voisin a un numero plus petit que l'ancien minimum, on ne le prend pas en compte
                          // En effet, le nouveau minimum est un numero d'element qui lui sera superieur
                          // Rq : voisin>old_tri => voisin != -1 donc pour une face au bord le if est false
                          if (voisin>old_tri)
                            {
                              min_tri=min(min_tri,voisin);
                              if(min_tri==voisin)
                                {
                                  dvarkeep=pow((positions(elem,ori(f0))-positions(voisin,ori(f0))),2);
                                }
                            }
                        }

                      // Si l'ancien minimum est inferieur a elem et le nouveau superieur, il faut traiter elem juste avant le nouveau voisin
                      if (old_tri<elem && min_tri>elem)
                        {
                          coeff(compt)=valeur_diag_naire(mat);
                          tab2(compt2)=nb_elem*(mat)+elem; //la ou commence les matrices A, B, C, etc... donc (nb_elem*nb_equation_CH)+1 ??
                          compt++;
                          compt2++;
                        }

                      // On complete les tableaux avec les valeurs dans les voisins non aux bords du domaine
                      coeff(compt)=alphaMatrix(j)/dvarkeep;
                      tab2(compt2)=min_tri+nb_elem*(mat);
                      compt++;
                      compt2++;

                      // Si on est a la fin et que elem n'a pas encore ete mis dans coeff, il faut le faire avant la fin de la boucle sur les faces de l'element
                      // Ceci arrive si elem > tous les numeros d'elements de ses voisins
                      if (ncomp==2*nb_compo_-nb_faces_au_bord-1 && min_tri<elem)
                        {
                          coeff(compt)=valeur_diag_naire(mat);
                          tab2(compt2)=elem+(mat)*nb_elem;
                          compt++;
                          compt2++;
                        }
                      // On sauve le numero d'element "minimum"
                    }
                  old_tri=min_tri;
                  //}

                }
              coeff(compt)=1;
              tab2(compt2)=((ligne+nb_equation_CH)*nb_elem)+elem;
              compt++;
              compt2++;
            }

        }

      /*for (int ligne=0; ligne<nb_equation_CH; ligne++)
        {
          // Boucle sur le nombre d'elements
          for(int elem=0; elem<nb_elem; elem++)
            {
              valeur_diag_naire=0;
              old_tri=-1;
              nb_faces_au_bord=0;
              tab1(compt1)=compt2+1;//
              compt1++;//

              // Calcul du nombre de bords
              for(int ncomp=0; ncomp<2*nb_compo_; ncomp++)
                {
                  f0 = elem_faces(elem,ncomp);
                  voisin=face_voisins(f0,0);
                  if (face_voisins(f0,0)==elem) voisin=face_voisins(f0,1);

                  // On connait le voisin associe a la face en cours. On regarde s'il est au bord
                  if (voisin==-1) nb_faces_au_bord++;
                }

              // Boucle sur les faces
              //for (int mat=0; mat<nb_equation_CH; mat++)
              //{
              //int j=(ligne*nb_equation_CH)+mat;
              for(int ncomp=0; ncomp<2*nb_compo_-nb_faces_au_bord; ncomp++)
                {
                  min_tri=nb_elem+1;
                  for (int mat=0; mat<nb_equation_CH; mat++)
                    {
                      int j=(ligne*nb_equation_CH)+mat;
                      for(int ncomp_tri=0; ncomp_tri<2*nb_compo_; ncomp_tri++)
                        {
                          f0 = elem_faces(elem,ncomp_tri);

                          // Voisin associe a la face - On s'assure de traiter un voisin et pas l'element resident
                          voisin=face_voisins(f0,0);
                          if (face_voisins(f0,0)==elem) voisin=face_voisins(f0,1);

                          // On va compter le nombre de voisins existants (faces non aux bords de l'element elem)
                          // Ceci sert a calculer la contribution au terme relatif a l'element elem
                          if (voisin!=-1 && old_tri==-1)
                            {
                              dvar2=pow((positions(elem,ori(f0))-positions(voisin,ori(f0))),2);
                              valeur_diag_naire(mat)+=-alphaMatrix(j)/dvar2;
                            }

                          // Si le voisin a un numero plus petit que l'ancien minimum, on ne le prend pas en compte
                          // En effet, le nouveau minimum est un numero d'element qui lui sera superieur
                          // Rq : voisin>old_tri => voisin != -1 donc pour une face au bord le if est false
                          if (voisin>old_tri)
                            {
                              min_tri=min(min_tri,voisin);
                              if(min_tri==voisin)
                                {
                                  dvarkeep=pow((positions(elem,ori(f0))-positions(voisin,ori(f0))),2);
                                }
                            }

                        }

                      // Si l'ancien minimum est inferieur a elem et le nouveau superieur, il faut traiter elem juste avant le nouveau voisin
                      if (old_tri<elem && min_tri>elem)
                        {
                          coeff(compt)=valeur_diag_naire(mat);
                          tab2(compt2)=elem+(mat*nb_elem);
                          compt++;
                          compt2++;
                           if (old_tri==-1)
                             {
                               // On ajoute un terme a tab1 dans le cas ou l'element diagonal est le premier terme
                               tab1(compt1)=compt2;
                               compt1++;
                               //old_tri=elem+(mat*nb_elem);
                             }
                        }

                      // On complete les tableaux avec les valeurs dans les voisins non aux bords du domaine
                      coeff(compt)=alphaMatrix(j)/dvarkeep;
                      tab2(compt2)=min_tri+(nb_elem*mat);
                      compt++;
                      compt2++;
                      if (old_tri==-1)
                        {
                          // On ajoute un terme a tab1 dans le cas ou le premier terme est extradiagonal
                          tab1(compt1)=compt2;
                          compt1++;
                        }

                      // Si on est a la fin et que elem n'a pas encore ete mis dans coeff, il faut le faire avant la fin de la boucle sur les faces de l'element
                      // Ceci arrive si elem > tous les numeros d'elements de ses voisins
                      // Bien sur dans ce cas la, l'element diagonal est forcement le dernier mais pas le premier - il n'existe pas de maillage avec des mailles seules sans voisins
                      if (ncomp==2*nb_compo_-nb_faces_au_bord-1 && min_tri<elem)
                        {
                          coeff(compt)=valeur_diag_naire(mat);
                          tab2(compt2)=elem+(mat*nb_equation_CH);
                          compt++;
                          compt2++;
                        }

                      // On sauve le numero d'element "minimum"
                    }
                  old_tri=min_tri;

                  //}
                }


              // On ajoute le 1 de la diagonale de la sous-matrice I du bas a droite
              coeff(compt)=1;
              tab2(compt2)=elem+(ligne+nb_equation_CH)*nb_elem;
              compt++;
              compt2++;
              //ajouter par mr264902
              tab1(compt1)=compt2;//

            }
        }*/

      // cf Mat_Morse.h/.cpp, tab1() et tab2() sont a utiliser au sens Fortran
      // EXPLICATION :
      // coeff stocke les valeurs des coefficients
      // tab2 stocke les numeros de colonnes dans la matrice de ces coefficients, au sens FORTRAN (i-e 1 <= tab2[i] <= n)
      // tab1 stocke le rang de l'element de tab2() pour lequel on change de ligne. Ainsi pour tab1[j], on est a la i-eme ligne avec :
      //      ( tab1[i] <= j < tab1[i+1] ) qui implique qu'on considere toujours la i-eme ligne,
      //      et des que cela n'est plus respecte, le passage a la ligne suivante est effectue

      // De par l'algorithme, tab1() est deja au sens FORTRAN. Reste a le faire pour tab2()
      // De plus on doit mettre le dernier terme de tab1() a dimensionnement+1
      tab2+=1;
      tab1(compt1)=dimensionnement+1;
      compt1+=1;

      if (compt!=dimensionnement)
        {
          Cerr << "Erreur lors du calcul de la matrice du point fixe : nombre d'elements non nuls calcules different du nombre d'elements non nuls prevus" << finl;
          exit();
        }
      //   Cerr<<"Nombre d'elements non nuls calcules="<<compt<<finl;
      //   Cerr<<"Nombre d'elements non nuls prevus="<<dimensionnement<<finl;
      //   Cerr<<"Assemblage de la matrice du point fixe : OK"<<finl;
      //   Cerr<<"==========================================="<<finl;


      // Test des tableaux
      /*Cerr<<"coeff="<<coeff<<finl;
      Cerr<<"tab1="<<tab1<<finl;
      Cerr<<"tab2="<<tab2<<finl;*/
      //Cerr <<"matrice_diffusion_CH final = "<<matrice_diffusion_CH<<finl;

    }

}


/*! @brief Calcul du point fixe
 *
 */
void Source_Con_Phase_field::calculer_point_fixe(const DoubleTab& c, const DoubleTab& mutilde, const Matrice_Morse& matrice_diffusion_CH, DoubleTab& c_demi, DoubleTab& mutilde_demi)
{
  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const int nb_elem = zone_VDF.nb_elem_tot();
  const double theta=0.6;
  DoubleVect residu1(nb_elem);
  DoubleVect residu2(nb_elem);
  DoubleVect term_cin(nb_elem);
  // Ceci correspond a u2/2*drhodc=mutilde_dyn-mutilde

  DoubleVect second_membre(2*nb_elem);
  DoubleVect M(nb_elem);
  DoubleVect old_Phi(nb_elem);
  DoubleVect old_M(nb_elem);
  DoubleVect Solu_temp(2*nb_elem);
  DoubleTab Phi(c);

  // Initialisation des variables
  M=mutilde;
  old_Phi=0.;
  old_M=0.;

  // Test de convergence
  for(int n_elem=0; n_elem<nb_elem; n_elem++)
    {
      residu1(n_elem)=Phi(n_elem)-old_Phi(n_elem);
      residu2(n_elem)=M(n_elem)-old_M(n_elem);
    }

  Cerr<<"==========================="<<finl;
  Cerr<<"Assemblage du second membre"<<finl;
  Cerr<<"Inversion du systeme du point fixe"<<finl;
  Cerr<<"(Attention : fonctionne seulement en sequentiel)"<<finl;
  if (Process::nproc()>1)
    exit();

  while(norme_array(residu1) > epsilon_ || norme_array(residu2)>epsilon_)
    {
      // Assemblage du second membre
      for(int n_elem=0; n_elem<nb_elem; n_elem++)
        {
          second_membre(n_elem)=c(n_elem);
          term_cin(n_elem)=0.;
          if (mutype_==1)
            {
              term_cin(n_elem)+=(0.5*u_carre_(n_elem))*drhodc(n_elem);
            }

          second_membre(n_elem+nb_elem)=term_cin(n_elem)+beta*(this->*dWdc)(Phi(n_elem));
        }

      // Inversion de la matrice morse
      matrice_diffusion_CH.inverse(second_membre, Solu_temp, eps_);

      old_Phi=Phi;
      old_M=M;
      for(int n_elem=0; n_elem<nb_elem; n_elem++)
        {
          Phi(n_elem)=Solu_temp(n_elem);
          M(n_elem)=Solu_temp(n_elem+nb_elem);
        }
      Phi.echange_espace_virtuel();
      M.echange_espace_virtuel();

      for(int n_elem=0; n_elem<nb_elem; n_elem++)
        {
          residu1(n_elem)=Phi(n_elem)-old_Phi(n_elem);
          residu2(n_elem)=M(n_elem)-old_M(n_elem);
        }

    }


  // Ecriture de c et mutilde au temps n+1/2, cf le theta-schema:
  // Phi=theta*c_demi+(1.-theta)*c
  // M=theta*mutilde_demi+(1.-theta)*mutilde

  if(theta!=0)
    {
      c_demi=Phi;
      c_demi.ajoute(-(1.-theta),c);
      c_demi/=theta;
      mutilde_demi=M;
      mutilde_demi.ajoute(-(1.-theta),mutilde);
      mutilde_demi/=theta;
    }
  else
    {
      Cerr<<"Attention : le traitement est explicite (theta=0) - De plus pour theta < 0.5 ou theta > 1 le systeme est instable !!"<<finl;
      c_demi=c;
      mutilde_demi=mutilde;
    }
}


/*! @brief Construire le residu du GMRES NL
 *
 */
void Source_Con_Phase_field::construire_systeme(const DoubleTab& c, const Matrice_Morse& matrice_diffusion_CH, DoubleTab& v0, const DoubleTab& x1)
{
  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const int nb_elem = zone_VDF.nb_elem_tot();
  int nb_elem_tot = c.size_totale();

  if (type_systeme_naire_==0)
    {
      DoubleVect term_cin(nb_elem);
      // Ceci correspond a u2/2*drhodc=mutilde_dyn-mutilde

      DoubleVect second_membre(2*nb_elem);
      DoubleTab Ax1(2*nb_elem);

      // Assemblage du second membre

      for(int n_elem=0; n_elem<nb_elem; n_elem++)
        {
          second_membre(n_elem)=c(n_elem);

          term_cin(n_elem)=0.;
          if (mutype_==1)
            {
              term_cin(n_elem)+=(0.5*u_carre_(n_elem))*drhodc(n_elem);
            }

          second_membre(n_elem+nb_elem)=term_cin(n_elem)+beta*(this->*dWdc)(x1(n_elem));
        }
      /*Cerr <<"second_membre = "<<second_membre<<finl;
      Cerr <<"x1 = "<<x1<<finl;*/

      //   // Ajoute par DJ pour Debog
      //   //-------------------------
      //   {
      //     DoubleTab secmem_c(c);
      //     DoubleTab secmem_mutilde(c);
      //     for(int n_elem=0; n_elem<nb_elem; n_elem++)
      //       {
      //         secmem_c(n_elem) = second_membre(n_elem);
      //         secmem_mutilde(n_elem) = second_membre(n_elem+nb_elem);
      //       }
      //     Debog::verifier("Construire Systeme secmem_c : ",secmem_c);
      //     Debog::verifier("Construire Systeme secmem_mutilde : ",secmem_mutilde);
      //   }
      //   //-------------------------

      // Calcul du produit matrice / vecteur utilise
      matrice_diffusion_CH.multvect_(x1,Ax1);
      // Modifie par DJ
      //---------------
      {
        DoubleTab Ax1_c(c);
        DoubleTab Ax1_mutilde(c);
        for(int n_elem=0; n_elem<nb_elem; n_elem++)
          {
            Ax1_c(n_elem) = Ax1(n_elem);
            Ax1_mutilde(n_elem) = Ax1(n_elem+nb_elem);
          }

        Ax1_c.echange_espace_virtuel();
        Ax1_mutilde.echange_espace_virtuel();

        for(int n_elem=0; n_elem<nb_elem; n_elem++)
          {
            Ax1(n_elem) = Ax1_c(n_elem);
            Ax1(n_elem+nb_elem) = Ax1_mutilde(n_elem);
          }
      }
      Cerr <<"Ax1 = "<<Ax1<<finl;

      //---------------
      //   // Ajouter par DJ pour Debog
      //   //--------------------------
      //   {
      //     DoubleTab Ax1_c(c);
      //     DoubleTab Ax1_mutilde(c);
      //     for(int n_elem=0; n_elem<nb_elem; n_elem++)
      //       {
      //         Ax1_c(n_elem) = Ax1(n_elem);
      //         Ax1_mutilde(n_elem) = Ax1(n_elem+nb_elem);
      //       }
      //     Debog::verifier("Construire Systeme Ax1_c : ",Ax1_c);
      //     Debog::verifier("Construire Systemeyes Ax1_mutilde : ",Ax1_mutilde);
      //   }
      //   //-------------------------

      // Calcul de v0 = [A(xn) xn - bn]
      for(int n_elem=0; n_elem<2*nb_elem; n_elem++)
        {
          v0(n_elem)=(Ax1(n_elem)-second_membre(n_elem));
        }
    }
  else if (type_systeme_naire_==1)
    {
      //DoubleVect term_cin(nb_elem_tot);
      // Ceci correspond a u2/2*drhodc=mutilde_dyn-mutilde

      DoubleVect second_membre(2*nb_elem_tot);
      DoubleTab Ax1(2*nb_elem_tot);
      DoubleTab terme_non_lin(c);
      DoubleTab x1_c(c);
      for (int j=0; j<nb_equation_CH; j++)
        {
          for(int n_elem=0; n_elem<nb_elem; n_elem++)
            {
              x1_c(n_elem,j) = x1(n_elem+(j*nb_elem));
            }
        }


      // Assemblage du second membre
      if (type_potentiel_analytique_==0)
        {
          terme_non_lin=(this->*dWdc_naire)(x1_c, eq1, eq2);
        }
      else if (type_potentiel_analytique_==1 && nb_equation_CH==2)
        {
          terme_non_lin=(this->*dWdc_naire_analytique_ter)(x1_c,angle_psi,x0Eq,x1Eq,a0Eq,a1Eq);

        }
      else if (type_potentiel_analytique_==1 && nb_equation_CH==3)
        {
          terme_non_lin=(this->*dWdc_naire_analytique_quater)(x1_c,angle_psi,angle_phi,x0Eq,x1Eq,x2Eq,a0Eq,a1Eq,a2Eq);
        }


      for (int j=0; j<nb_equation_CH; j++)
        {
          for(int n_elem=0; n_elem<nb_elem; n_elem++)
            {
              second_membre(n_elem+(j*nb_elem))=c(n_elem,j);
              if (type_potentiel_analytique_==0)
                {
                  terme_non_lin(n_elem,j)*=betaMatrix(j);
                }
              second_membre(n_elem+nb_elem_tot+(j*nb_elem))=terme_non_lin(n_elem,j);
            }
        }
      /* Cerr <<"c second_membre = "<<c<<finl;
       Cerr <<"x1_c = "<<x1_c<<finl;
       Cerr <<"terme_non_lineaire = "<<terme_non_lin<<finl;
       Cerr <<"second_membre = "<<second_membre<<finl;*/


      //   // Ajoute par DJ pour Debog
      //   //-------------------------
      //   {
      //     DoubleTab secmem_c(c);
      //     DoubleTab secmem_mutilde(c);
      //     for(int n_elem=0; n_elem<nb_elem; n_elem++)
      //       {
      //         secmem_c(n_elem) = second_membre(n_elem);
      //         secmem_mutilde(n_elem) = second_membre(n_elem+nb_elem);
      //       }
      //     Debog::verifier("Construire Systeme secmem_c : ",secmem_c);
      //     Debog::verifier("Construire Systeme secmem_mutilde : ",secmem_mutilde);
      //   }
      //   //-------------------------

      // Calcul du produit matrice / vecteur utilise
      matrice_diffusion_CH.multvect_(x1,Ax1);
      // Modifie par DJ
      //---------------
      {
        DoubleTab Ax1_c(c.size_totale());
        DoubleTab Ax1_mutilde(c.size_totale());
        for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
          {
            Ax1_c(n_elem) = Ax1(n_elem);
            Ax1_mutilde(n_elem) = Ax1(n_elem+nb_elem_tot);
          }

        Ax1_c.echange_espace_virtuel();
        Ax1_mutilde.echange_espace_virtuel();

        for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
          {
            Ax1(n_elem) = Ax1_c(n_elem);
            Ax1(n_elem+nb_elem_tot) = Ax1_mutilde(n_elem);
          }
      }
      //Cerr <<"Ax1 = "<<Ax1<<finl;

      //---------------
      //   // Ajouter par DJ pour Debog
      //   //--------------------------
      //   {
      //     DoubleTab Ax1_c(c);
      //     DoubleTab Ax1_mutilde(c);
      //     for(int n_elem=0; n_elem<nb_elem; n_elem++)
      //       {
      //         Ax1_c(n_elem) = Ax1(n_elem);
      //         Ax1_mutilde(n_elem) = Ax1(n_elem+nb_elem);
      //       }
      //     Debog::verifier("Construire Systeme Ax1_c : ",Ax1_c);
      //     Debog::verifier("Construire Systemeyes Ax1_mutilde : ",Ax1_mutilde);
      //   }
      //   //-------------------------

      // Calcul de v0 = [A(xn) xn - bn]
      for(int n_elem=0; n_elem<2*nb_elem_tot; n_elem++)
        {
          v0(n_elem)=(Ax1(n_elem)-second_membre(n_elem));
        }
      //Cerr <<"v0 = "<<v0<<finl;


    }


  return ;

}


/*! @brief Construire le residu du GMRES NL
 *
 */
void Source_Con_Phase_field::matvect(const DoubleTab& c, const Matrice_Morse& matrice_diffusion_CH, const DoubleTab& v0, const DoubleTab& x1, DoubleTab& v1)
{
  const double delta = 1.e-5;

  DoubleTab v2(v1);
  DoubleTab x1_(v0);

  v1 = 0.;

  x1_ = v0;
  x1_ *= delta;
  x1_ += x1;

  // Construction de la differentielle (dans la direction v0)
  construire_systeme(c, matrice_diffusion_CH, v1, x1_);
  construire_systeme(c, matrice_diffusion_CH, v2, x1);

  v1 -= v2;
  v1/= delta;

  return ;

}


/*! @brief Algorithme GMRES Non Lineaire
 *
 * @param (mutilde) inverse A x = b
 * @return (int)
 */
int Source_Con_Phase_field::non_lin_gmres(const DoubleTab& c, const DoubleTab& mutilde, const Matrice_Morse& matrice_diffusion_CH, DoubleTab& c_demi, DoubleTab& mutilde_demi)
// int Source_Con_Phase_field::non_lin_gmres(DoubleTab& c, DoubleTab& mutilde, Matrice_Morse& matrice_diffusion_CH, DoubleTab& x1)
//---------------
/*--------------------------------------------------------------
  Solves Ax = b by GMRES with nit iterations (reinitializations)
  and nkr Krilov vectors.
*/

// Copied from Sch_Crank_Nicholson.cpp (scalar)//
{

  if (type_systeme_naire_==0)
    {
      int i,j,nk,i0,im,it,ii;
      double tem=1.,res,ccos,ssin ;
      const int ns = 2*c.size_totale();

      // Ajoute par DJ
      //--------------
      int nb_elem_tot = c.size_totale();
      DoubleTab x1(ns);
      x1=0.;
      for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
        {
          x1(n_elem)=c(n_elem);
          x1(n_elem+nb_elem_tot)=mutilde(n_elem);
        }
      //--------------
      /* Cerr <<"c = "<<c<<finl;
       Cerr <<"mutilde = "<<mutilde<<finl;
       Cerr <<"x1 = "<<x1<<finl;*/

      // A present dans le jdd
      //   double epsGMRES=1.e-10;
      //   int nkr=2;                         // dimension de l'espace de Krylov
      //   int nit=10;                        // nombre d'iterations
      //   double rec_min = 1.e-8;
      //   double rec_max = 0.1  ;

      DoubleTab v(ns,nkr);                         // Krilov vectors
      DoubleTab h(nkr+1,nkr);                // Heisenberg maatrix of coefficients
      DoubleVect r(nkr+1);
      DoubleTab v0(x1);
      DoubleTab v1(x1) ;

      // Initialisation
      v = 0. ;
      v0 = 0. ;
      v1 = 0. ;

      // v0 = -1.*construire_systeme(eqn, x1); // DJ

      // Cerr << " gmres : avant construire systeme " << finl ;
      construire_systeme(c, matrice_diffusion_CH, v0, x1) ;
      //   // Ajoute par DJ pour Debog
      //   //-------------------------
      //   {
      //     DoubleTab v0_c(c);
      //     DoubleTab v0_mutilde(c);
      //     for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
      //       {
      //         v0_c(n_elem) = v0(n_elem);
      //         v0_mutilde(n_elem) = v0(n_elem+nb_elem_tot);
      //       }
      //     Debog::verifier("GMRES Non Lineaire v0_c apres construire_systeme initial : ",v0_c);
      //     Debog::verifier("GMRES Non Lineaire v0_mutilde apres construire_systeme initial : ",v0_mutilde);
      //   }
      //   //-------------------------
      v0 *= -1. ;
      // Cerr << " gmres : apres construire systeme " << finl ;
      //   // Ajoute par DJ pour Debog
      //   //-------------------------
      //   {
      //     DoubleTab x1_c(c);
      //     DoubleTab x1_mutilde(c);
      //     DoubleTab systeme_c(c);
      //     DoubleTab systeme_mutilde(c);
      //     for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
      //       {
      //         x1_c(n_elem) = x1(n_elem);
      //         x1_mutilde(n_elem) = x1(n_elem+nb_elem_tot);
      //         systeme_c(n_elem) = v0(n_elem);
      //         systeme_mutilde(n_elem) = v0(n_elem+nb_elem_tot);
      //       }
      //     Debog::verifier("GMRES Non Lineaire construire_systeme x1_c : ", x1_c);
      //     Debog::verifier("GMRES Non Lineaire construire_systeme x1_mutilde : ", x1_mutilde);
      //     Debog::verifier("GMRES Non Lineaire construire_systeme systeme_c : ", systeme_c);
      //     Debog::verifier("GMRES Non Lineaire construire_systeme systeme_mutilde : ", systeme_mutilde);
      //   }
      //   //-------------------------


      res = 0. ;
      // Modifie par DJ
      //---------------
      //   for (ii=0; ii<ns;ii++) res += v0(ii) * v0(ii) ;
      for (ii=0; ii<c.size(); ii++) res += v0(ii) * v0(ii) ;
      for (ii=c.size_totale(); ii<c.size_totale()+c.size(); ii++) res += v0(ii) * v0(ii) ;
      //---------------
      res=mp_sum(res);
      //   Debog::verifier("GMRES Non Lineaire res debut : ",res);
      res = sqrt(res)  ;

      //Cerr<<"initial residual = "<<res<<finl;

      if(res<rec_min)
        return 0; // nothing to do

      rec_min = (rec_min<res*epsGMRES) ? res*epsGMRES : rec_min;
      rec_min = (rec_min<rec_max) ? rec_min : rec_max ;

      Cerr << "Source Concentration Phase Field - GMRES NL" << finl;
      Cerr << "Stopping rule scalar : " << rec_min << finl;

      // iterations
      for(it=0; it<nit; it++)
        {
          nk = nkr;

          //...Orthogonalisation of Arnoldi
          v0 /= res;
          r = 0. ;
          r[0] = res;
          h=0.;

          for(j=0; j<nkr; j++)
            {
              for(ii=0; ii<ns; ii++) v(ii,j) = v0(ii);
              //           // Ajoute par DJ pour Debog
              //           //-------------------------
              //           {
              //             DoubleTab v_c(c);
              //             DoubleTab v_mutilde(c);
              //             for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
              //               {
              //                 v_c(n_elem) = v(n_elem,j);
              //                 v_mutilde(n_elem) = v(n_elem+nb_elem_tot,j);
              //               }
              //             Debog::verifier("GMRES Non Lineaire v_c : ",v_c);
              //             Debog::verifier("GMRES Non Lineaire v_mutilde : ",v_mutilde);
              //           }
              //           //-------------------------
              //          v0 = a * v0; // commente par DJ
              //      Cerr << " x1 avant matvec " << x1 << finl ;
              //           MatVect(...); // Modif pour le PhF
              matvect(c, matrice_diffusion_CH, v0, x1, v1) ;
              //           // Ajoute par DJ pour Debog
              //           //-------------------------
              //           {
              //             DoubleTab v1_c(c);
              //             DoubleTab v1_mutilde(c);
              //             for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
              //               {
              //                 v1_c(n_elem) = v1(n_elem);
              //                 v1_mutilde(n_elem) = v1(n_elem+nb_elem_tot);
              //               }
              //             Debog::verifier("GMRES Non Lineaire v1_c apres matvect : ",v1_c);
              //             Debog::verifier("GMRES Non Lineaire v1_mutilde apres matvect : ",v1_mutilde);
              //           }
              //           //-------------------------
              v0 = v1 ;

              // Modifie par DJ
              //---------------
              for(i=0; i<=j; i++)
                {
                  //               for (ii=0; ii<ns;ii++) h(i,j) += v0(ii) * v(ii,i);
                  for (ii=0; ii<c.size(); ii++) h(i,j) += v0(ii) * v(ii,i);
                  for (ii=c.size_totale(); ii<c.size_totale()+c.size(); ii++) h(i,j) += v0(ii) * v(ii,i);
                  h(i,j)=mp_sum(h(i,j));
                  //               Debog::verifier("GMRES Non Lineaire h(i,j) : ",h(i,j));
                  for (ii=0; ii<ns; ii++) v0(ii) -= h(i,j) * v(ii,i);
                  //               {
                  //                 DoubleTab v0_c(c);
                  //                 DoubleTab v0_mutilde(c);
                  //                 for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
                  //                   {
                  //                     v0_c(n_elem) = v0(n_elem);
                  //                     v0_mutilde(n_elem) = v0(n_elem+nb_elem_tot);
                  //                   }
                  //                 Debog::verifier("GMRES Non Lineaire v0_c apres modif : ",v0_c);
                  //                 Debog::verifier("GMRES Non Lineaire v0_mutilde apres modif : ",v0_mutilde);
                  //               }
                }
              //---------------
              // tem = sqrt(v0 * v0);
              tem = 0. ;
              // Modifie par DJ
              //---------------
              //           for (ii=0; ii<ns;ii++) tem += v0(ii) * v0(ii) ;
              for (ii=0; ii<c.size(); ii++) tem += v0(ii) * v0(ii) ;
              for (ii=c.size_totale(); ii<c.size_totale()+c.size(); ii++) tem += v0(ii) * v0(ii) ;
              //---------------
              tem=mp_sum(tem);
              //           Debog::verifier("GMRES Non Lineaire tem apres mp_sum : ",tem);
              tem = sqrt(tem)  ;
              h(j+1,j) = tem;
              if(tem<rec_min)
                {
                  nk = j+1;
                  // Cerr<<"tem="<<tem<<" nk="<<nk<<finl;
                  goto l5;
                }
              v0 /= tem;
            }

          //...Triangularisation
l5:
          for(i=0; i<nk; i++)
            {
              im = i+1;
              tem = 1./sqrt(h(i,i)*h(i,i) + h(im,i)*h(im,i));
              ccos = h(i,i) * tem;
              ssin = - h(im,i) * tem;
              for(j=i; j<nk; j++)
                {
                  tem = h(i,j);
                  h(i,j) = ccos * tem - ssin * h(im,j);
                  h(im,j) =  ssin * tem + ccos * h(im,j);
                }
              r[im] = ssin * r[i];
              r[i] *= ccos;
            }

          //...Solution of linear system
          for(i=nk-1; i>=0; i--)
            {
              r[i] /= h(i,i);
              for(i0=i-1; i0>=0; i0--)
                r[i0] -= h(i0,i)* r[i];
            }
          for(i=0; i<nk; i++)
            for(ii=0; ii<ns; ii++)  x1(ii) += r[i]*v(ii,i);
          //       // Ajoute par DJ pour Debog
          //       //-------------------------
          //       {
          //         DoubleTab x1_c(c);
          //         DoubleTab x1_mutilde(c);
          //         for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
          //           {
          //             x1_c(n_elem) = x1(n_elem);
          //             x1_mutilde(n_elem) = x1(n_elem+nb_elem_tot);
          //           }
          //         Debog::verifier("GMRES Non Lineaire x1_c apres modif : ", x1_c);
          //         Debog::verifier("GMRES Non Lineaire x1_mutilde apres modif : ", x1_mutilde);
          //       }
          //       //-------------------------
          //       x1.echange_espace_virtuel();

          //Cerr <<" futur in gemres : " << x1 << finl ;
          //New residual and stopping tests
          //v0 = -1.*construire_systeme(eqn, x1); //DJ
          //res = sqrt(v0 * v0);

          construire_systeme(c, matrice_diffusion_CH, v0, x1) ;
          v0 *= -1. ;

          // Cerr << "  v0 in gmres " << v0 << finl ;

          res = 0. ;
          // Modifie par DJ
          //---------------
          //   for (ii=0; ii<ns;ii++) res += v0(ii) * v0(ii) ;
          for (ii=0; ii<c.size(); ii++) res += v0(ii) * v0(ii) ;
          for (ii=c.size_totale(); ii<c.size_totale()+c.size(); ii++) res += v0(ii) * v0(ii) ;
          //---------------
          res=mp_sum(res);
          //       Debog::verifier("GMRES Non Lineaire res fin : ",res);
          res = sqrt(res)  ;

          Cerr<<" - At it = "<< it <<", residu scalar = "<< res << finl;

          if(res<rec_min)
            {
              // Ajoute par DJ
              //--------------
              for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
                {
                  c_demi(n_elem)=x1(n_elem);
                  mutilde_demi(n_elem)=x1(n_elem+nb_elem_tot);
                }
              // L'echange espace virtuel est fait dans premier_demi_dt()
              //--------------
              //         return 1;
              Cerr << "Number of iterations to reach convergence : " << it+1 << finl;
              Cerr << "" << finl;
              Cerr <<"c_demi = "<<c_demi<<finl;
              Cerr <<"mutilde_demi = "<<mutilde_demi<<finl;

              return it;
            }
          else if (it==nit-1)
            {
              // Ajoute par DJ
              //--------------
              // On fait le choix de mettre a jour c_demi et mutilde_demi meme s'il n'y a pas convergence...
              for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
                {
                  c_demi(n_elem)=x1(n_elem);
                  mutilde_demi(n_elem)=x1(n_elem+nb_elem_tot);
                }
              //--------------
              Cerr << "Stopped before convergence" << finl;
              Cerr << "" << finl;
              Cerr <<"c_demi = "<<c_demi<<finl;
              Cerr <<"mutilde_demi = "<<mutilde_demi<<finl;

            }
        }
    }

  else if (type_systeme_naire_==1)
    {

      int i,j,nk,i0,im,it,ii;
      double tem=1.,res,ccos,ssin ;
      const int ns = 2*c.size_totale();
      const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
      const int nb_elem = zone_VDF.nb_elem_tot();

      // Ajoute par DJ
      //--------------
      int nb_elem_tot = c.size_totale();
      DoubleTab x1(ns);
      x1=0.;
      for (int ncomponent=0; ncomponent<nb_equation_CH; ncomponent++)
        {
          for (int nelem=0; nelem<nb_elem; nelem++)
            {
              x1(nelem+(ncomponent*nb_elem))=c(nelem,ncomponent);
              x1(nelem+(ncomponent*nb_elem)+nb_elem_tot)=mutilde(nelem,ncomponent);
            }
        }
      /*Cerr <<"c = "<<c<<finl;
      Cerr <<"mutilde = "<<mutilde<<finl;
      Cerr <<"x1 = "<<x1<<finl;
      */

      //--------------

      // A present dans le jdd
      //   double epsGMRES=1.e-10;
      //   int nkr=2;                         // dimension de l'espace de Krylov
      //   int nit=10;                        // nombre d'iterations
      //   double rec_min = 1.e-8;
      //   double rec_max = 0.1  ;

      DoubleTab v(ns,nkr);                         // Krilov vectors
      DoubleTab h(nkr+1,nkr);                // Heisenberg maatrix of coefficients
      DoubleVect r(nkr+1);
      DoubleTab v0(x1);
      DoubleTab v1(x1) ;

      // Initialisation
      v = 0. ;
      v0 = 0. ;
      v1 = 0. ;

      // v0 = -1.*construire_systeme(eqn, x1); // DJ

      // Cerr << " gmres : avant construire systeme " << finl ;
      construire_systeme(c, matrice_diffusion_CH, v0, x1) ;
      //   // Ajoute par DJ pour Debog
      //   //-------------------------
      //   {
      //     DoubleTab v0_c(c);
      //     DoubleTab v0_mutilde(c);
      //     for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
      //       {
      //         v0_c(n_elem) = v0(n_elem);
      //         v0_mutilde(n_elem) = v0(n_elem+nb_elem_tot);
      //       }
      //     Debog::verifier("GMRES Non Lineaire v0_c apres construire_systeme initial : ",v0_c);
      //     Debog::verifier("GMRES Non Lineaire v0_mutilde apres construire_systeme initial : ",v0_mutilde);
      //   }
      //   //-------------------------
      v0 *= -1. ;
      // Cerr << " gmres : apres construire systeme " << finl ;
      //   // Ajoute par DJ pour Debog
      //   //-------------------------
      //   {
      //     DoubleTab x1_c(c);
      //     DoubleTab x1_mutilde(c);
      //     DoubleTab systeme_c(c);
      //     DoubleTab systeme_mutilde(c);
      //     for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
      //       {
      //         x1_c(n_elem) = x1(n_elem);
      //         x1_mutilde(n_elem) = x1(n_elem+nb_elem_tot);
      //         systeme_c(n_elem) = v0(n_elem);
      //         systeme_mutilde(n_elem) = v0(n_elem+nb_elem_tot);
      //       }
      //     Debog::verifier("GMRES Non Lineaire construire_systeme x1_c : ", x1_c);
      //     Debog::verifier("GMRES Non Lineaire construire_systeme x1_mutilde : ", x1_mutilde);
      //     Debog::verifier("GMRES Non Lineaire construire_systeme systeme_c : ", systeme_c);
      //     Debog::verifier("GMRES Non Lineaire construire_systeme systeme_mutilde : ", systeme_mutilde);
      //   }
      //   //-------------------------


      res = 0. ;
      // Modifie par DJ
      //---------------
      //   for (ii=0; ii<ns;ii++) res += v0(ii) * v0(ii) ;
      /*
      		for (ii=0; ii<c.size(); ii++) res += v0(ii) * v0(ii) ;
      		for (ii=c.size_totale(); ii<c.size_totale()+c.size(); ii++) res += v0(ii) * v0(ii) ;
      */
      for (ii=0; ii<ns; ii++) res += v0(ii) * v0(ii) ;

      //---------------
      res=mp_sum(res);
      //   Debog::verifier("GMRES Non Lineaire res debut : ",res);
      res = sqrt(res)  ;

      //Cerr<<"initial residual = "<<res<<finl;

      if(res<rec_min)
        return 0; // nothing to do

      rec_min = (rec_min<res*epsGMRES) ? res*epsGMRES : rec_min;
      rec_min = (rec_min<rec_max) ? rec_min : rec_max ;

      Cerr << "Source Concentration Phase Field - GMRES NL" << finl;
      Cerr << "Stopping rule scalar : " << rec_min << finl;

      // iterations
      for(it=0; it<nit; it++)
        {
          nk = nkr;

          //...Orthogonalisation of Arnoldi
          v0 /= res;
          r = 0. ;
          r[0] = res;
          h=0.;

          for(j=0; j<nkr; j++)
            {
              for(ii=0; ii<ns; ii++) v(ii,j) = v0(ii);
              //           // Ajoute par DJ pour Debog
              //           //-------------------------
              //           {
              //             DoubleTab v_c(c);
              //             DoubleTab v_mutilde(c);
              //             for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
              //               {
              //                 v_c(n_elem) = v(n_elem,j);
              //                 v_mutilde(n_elem) = v(n_elem+nb_elem_tot,j);
              //               }
              //             Debog::verifier("GMRES Non Lineaire v_c : ",v_c);
              //             Debog::verifier("GMRES Non Lineaire v_mutilde : ",v_mutilde);
              //           }
              //           //-------------------------
              //          v0 = a * v0; // commente par DJ
              //      Cerr << " x1 avant matvec " << x1 << finl ;
              //           MatVect(...); // Modif pour le PhF
              matvect(c, matrice_diffusion_CH, v0, x1, v1) ;
              //           // Ajoute par DJ pour Debog
              //           //-------------------------
              //           {
              //             DoubleTab v1_c(c);
              //             DoubleTab v1_mutilde(c);
              //             for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
              //               {
              //                 v1_c(n_elem) = v1(n_elem);
              //                 v1_mutilde(n_elem) = v1(n_elem+nb_elem_tot);
              //               }
              //             Debog::verifier("GMRES Non Lineaire v1_c apres matvect : ",v1_c);
              //             Debog::verifier("GMRES Non Lineaire v1_mutilde apres matvect : ",v1_mutilde);
              //           }
              //           //-------------------------
              v0 = v1 ;

              // Modifie par DJ
              //---------------
              for(i=0; i<=j; i++)
                {
                  //               for (ii=0; ii<ns;ii++) h(i,j) += v0(ii) * v(ii,i);
                  /*for (ii=0; ii<c.size(); ii++) h(i,j) += v0(ii) * v(ii,i);
                  for (ii=c.size_totale(); ii<c.size_totale()+c.size(); ii++) h(i,j) += v0(ii) * v(ii,i);
                  */
                  for (ii=0; ii<ns; ii++) h(i,j) += v0(ii) * v(ii,i);

                  h(i,j)=mp_sum(h(i,j));
                  //               Debog::verifier("GMRES Non Lineaire h(i,j) : ",h(i,j));
                  for (ii=0; ii<ns; ii++) v0(ii) -= h(i,j) * v(ii,i);
                  //               {
                  //                 DoubleTab v0_c(c);
                  //                 DoubleTab v0_mutilde(c);
                  //                 for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
                  //                   {
                  //                     v0_c(n_elem) = v0(n_elem);
                  //                     v0_mutilde(n_elem) = v0(n_elem+nb_elem_tot);
                  //                   }
                  //                 Debog::verifier("GMRES Non Lineaire v0_c apres modif : ",v0_c);
                  //                 Debog::verifier("GMRES Non Lineaire v0_mutilde apres modif : ",v0_mutilde);
                  //               }
                }
              //---------------
              // tem = sqrt(v0 * v0);
              tem = 0. ;
              // Modifie par DJ
              //---------------
              //           for (ii=0; ii<ns;ii++) tem += v0(ii) * v0(ii) ;
              /*

              				for (ii=0; ii<c.size(); ii++) tem += v0(ii) * v0(ii) ;
              				for (ii=c.size_totale(); ii<c.size_totale()+c.size(); ii++) tem += v0(ii) * v0(ii) ;

              */
              for (ii=0; ii<ns; ii++) tem += v0(ii) * v0(ii) ;

              //---------------
              tem=mp_sum(tem);
              //           Debog::verifier("GMRES Non Lineaire tem apres mp_sum : ",tem);
              tem = sqrt(tem)  ;
              h(j+1,j) = tem;
              if(tem<rec_min)
                {
                  nk = j+1;
                  // Cerr<<"tem="<<tem<<" nk="<<nk<<finl;
                  goto l5naire;
                }
              v0 /= tem;
            }

          //...Triangularisation
l5naire:
          for(i=0; i<nk; i++)
            {
              im = i+1;
              tem = 1./sqrt(h(i,i)*h(i,i) + h(im,i)*h(im,i));
              ccos = h(i,i) * tem;
              ssin = - h(im,i) * tem;
              for(j=i; j<nk; j++)
                {
                  tem = h(i,j);
                  h(i,j) = ccos * tem - ssin * h(im,j);
                  h(im,j) =  ssin * tem + ccos * h(im,j);
                }
              r[im] = ssin * r[i];
              r[i] *= ccos;
            }

          //...Solution of linear system
          for(i=nk-1; i>=0; i--)
            {
              r[i] /= h(i,i);
              for(i0=i-1; i0>=0; i0--)
                r[i0] -= h(i0,i)* r[i];
            }
          for(i=0; i<nk; i++)
            for(ii=0; ii<ns; ii++)  x1(ii) += r[i]*v(ii,i);
          //       // Ajoute par DJ pour Debog
          //       //-------------------------
          //       {
          //         DoubleTab x1_c(c);
          //         DoubleTab x1_mutilde(c);
          //         for(int n_elem=0; n_elem<nb_elem_tot; n_elem++)
          //           {
          //             x1_c(n_elem) = x1(n_elem);
          //             x1_mutilde(n_elem) = x1(n_elem+nb_elem_tot);
          //           }
          //         Debog::verifier("GMRES Non Lineaire x1_c apres modif : ", x1_c);
          //         Debog::verifier("GMRES Non Lineaire x1_mutilde apres modif : ", x1_mutilde);
          //       }
          //       //-------------------------
          //       x1.echange_espace_virtuel();

          //Cerr <<" futur in gemres : " << x1 << finl ;
          //New residual and stopping tests
          //v0 = -1.*construire_systeme(eqn, x1); //DJ
          //res = sqrt(v0 * v0);

          construire_systeme(c, matrice_diffusion_CH, v0, x1) ;
          v0 *= -1. ;

          // Cerr << "  v0 in gmres " << v0 << finl ;

          res = 0. ;
          // Modifie par DJ
          //---------------
          //   for (ii=0; ii<ns;ii++) res += v0(ii) * v0(ii) ;
          /*
          			for (ii=0; ii<c.size(); ii++) res += v0(ii) * v0(ii) ;
          			for (ii=c.size_totale(); ii<c.size_totale()+c.size(); ii++) res += v0(ii) * v0(ii) ;
          */
          for (ii=0; ii<ns; ii++) res += v0(ii) * v0(ii) ;

          //---------------
          res=mp_sum(res);
          //       Debog::verifier("GMRES Non Lineaire res fin : ",res);
          res = sqrt(res)  ;

          Cerr<<" - At it = "<< it <<", residu scalar = "<< res << finl;

          if(res<rec_min)
            {
              // Ajoute par DJ
              //--------------
              for (int ncomponent=0; ncomponent<nb_equation_CH; ncomponent++)
                {
                  for(int n_elem=0; n_elem<nb_elem; n_elem++)
                    {
                      c_demi(n_elem,ncomponent)=x1(n_elem+(ncomponent*nb_elem));
                      mutilde_demi(n_elem,ncomponent)=x1(n_elem+nb_elem_tot+(ncomponent*nb_elem));
                    }
                }
              // L'echange espace virtuel est fait dans premier_demi_dt()
              //--------------
              //         return 1;
              Cerr << "Number of iterations to reach convergence : " << it+1 << finl;
              Cerr << "" << finl;
              //Cerr <<"c_demi = "<<c_demi<<finl;
              //Cerr <<"mutilde_demi = "<<mutilde_demi<<finl;

              return it;
            }
          else if (it==nit-1)
            {
              // Ajoute par DJ
              //--------------
              // On fait le choix de mettre a jour c_demi et mutilde_demi meme s'il n'y a pas convergence...
              for (int ncomponent=0; ncomponent<nb_equation_CH; ncomponent++)
                {
                  for(int n_elem=0; n_elem<nb_elem; n_elem++)
                    {
                      c_demi(n_elem,ncomponent)=x1(n_elem+(ncomponent*nb_elem));
                      mutilde_demi(n_elem,ncomponent)=x1(n_elem+nb_elem_tot+(ncomponent*nb_elem));
                    }
                }

              //--------------
              //Cerr <<"c_demi = "<<c_demi<<finl;
              //Cerr <<"mutilde_demi = "<<mutilde_demi<<finl;
              Cerr << "Stopped before convergence" << finl;
              Cerr << "" << finl;
            }
        }
    }


  return -1;
}


/*! @brief Calcul de mutilde au centre des elements
 *
 * @param (mutilde) mutilde au centre des elements
 */
void Source_Con_Phase_field::calculer_mutilde(DoubleTab& mutilde) const
{
  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const DoubleTab& c=eq_c.inconnue().valeurs();
  const DoubleTab& div_alpha_gradC=eq_c.get_div_alpha_gradC();

  // Calcul de mutilde

  mutilde = div_alpha_gradC;
  mutilde *= -1.;

  //Cerr << "mutilde  = -div_alpha_gradC " << mutilde << finl;

  if (type_systeme_naire_==0)
    {
      const int taille=mutilde.size();
      for (int i=0; i<taille; i++)
        {
          mutilde(i)+=beta*(this->*dWdc)(c(i));
          if(mutype_==1) //avec Ec
            {
              mutilde(i)+=(0.5*u_carre_(i))*drhodc(i);
            }
        }
    }
  else if (type_systeme_naire_==1)
    {
      DoubleTab potent_chimique(mutilde);

      if (type_potentiel_analytique_==1)
        {
          if (nb_equation_CH==2)
            {
              potent_chimique= (this->*dWdc_naire_analytique_ter) (c,angle_psi,x0Eq,x1Eq,a0Eq,a1Eq);
            }
          else if (nb_equation_CH==3)
            {
              potent_chimique= (this->*dWdc_naire_analytique_quater) (c,angle_psi,angle_phi,x0Eq,x1Eq,x2Eq,a0Eq,a1Eq,a2Eq);
            }
          else
            {
              Cerr <<"potentiel analytique for nb_equation>3 not yet implemented. Use defaut potentiel chimique"<<finl;
            }
        }
      else if (type_potentiel_analytique_==0)
        {
          potent_chimique=(this->*dWdc_naire) (c,eq1,eq2);

          for (int j=0; j<c.line_size(); j++)
            {
              for (int i=0; i<c.dimension(0); i++)
                {
                  potent_chimique(i,j)*=betaMatrix(j);
                }
            }

        }
      //Cerr <<"potentiel_chimique"<<potent_chimique<<finl;

      for (int j=0; j<mutilde.line_size(); j++)
        {
          for (int i=0; i<mutilde.dimension(0); i++)
            {
              mutilde(i,j)+= potent_chimique(i,j);
              if (mutype_==1) //avec Ec, n'est pas implemente pour le cas multicomposant
                {
                  Cerr << "potentiel_chimique_generalisee avec_energie_cinetique not implemented for multicomponent system"<<finl;
                  //mutilde(i,j)+=(0.5*u_carre_(i)*drhodc(i));
                }
            }
        }
    }
  //Cerr << "c = " << c << finl;
  //Cerr << "mutilde = " << mutilde << finl;

// L'espace virtuel n'est pas a jour
  assert_invalide_items_non_calcules(mutilde);
}


/*! @brief Calcul de u2 au centre des elements
 *
 */
void Source_Con_Phase_field::calculer_u2_elem(DoubleVect& u_carre)
{
  const Navier_Stokes_std& eq_ns=ref_cast(Navier_Stokes_std,le_probleme2->equation(0));
  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const DoubleTab& u=eq_ns.inconnue().valeurs();
  const DoubleTab& c=eq_c.inconnue().valeurs();

  DoubleVect u2;

  u_carre = c;
  u_carre = 0.;

  u2 = u;
  u2.carre(VECT_ALL_ITEMS);

  // u2 est calcule aux faces. mutildes, pression_thermo et c au centre des elements.
  // On doit donc interpoler u2 pour tout calculer au elements.

  // Interpolation de u2, sortie dans u_carre :
  // ------------------------------------------

  // Interpolation aux elements de u^2 (a partir de Discretisation_SG_VDF::calculer_tenS_capil)
  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const IntTab& elem_faces = zone_VDF.elem_faces();
  const IntTab& face_sommets = zone_VDF.face_sommets();
  const DoubleTab& positions = zone_VDF.xp();
  const Zone& zone_geom = zone_VDF.zone();
  const Domaine& dom=zone_geom.domaine();
  const int nb_elem = zone_VDF.nb_elem_tot();
  int f0,f1, som0,som1;
  double psi,val0,val1;
  DoubleTab u2_elem(nb_elem,dimension);
  const int nb_compo_ = u2_elem.line_size();

  // Boucle sur le nombre d'elements
  for(int elem=0; elem<nb_elem; elem++)
    {
      // Boucle sur le nombre de composantes du vecteur
      for(int ncomp=0; ncomp<nb_compo_; ncomp++)
        {
          f0 = elem_faces(elem,ncomp);
          f1 = elem_faces(elem,dimension+ncomp);

          val0 = u2(f0);
          val1 = u2(f1);

          som0 = face_sommets(f0,0);
          som1 = face_sommets(f1,0);

          psi = ( positions(elem,ncomp) - dom.coord(som0,ncomp) )
                / ( dom.coord(som1,ncomp) - dom.coord(som0,ncomp) ) ;

          if (std::fabs(psi) < 1.e-12)
            u2_elem(elem,ncomp) = val0 ;
          else if (std::fabs(1.-psi) < 1.e-12)
            u2_elem(elem,ncomp) = val1 ;
          else
            u2_elem(elem,ncomp) = val0 + psi * (val1-val0) ;
        }
    }

  // Pour chaque element, on additionne les composantes de u2_elem
  double tmp;
  for (int elem=0; elem<nb_elem; elem++)
    {
      for (int dim=0; dim<dimension; dim++)
        {
          tmp=u2_elem(elem,dim);
          u_carre(elem)+=tmp;
        }
    }

  u_carre.echange_espace_virtuel();

}


/*! @brief Calcul de alpha*(Grad(C))^2 au centre des elements
 *
 * @param (DoubleTab& alpha_gradC_carre) alpha*(Grad(C))^2 au centre des elements
 */
void Source_Con_Phase_field::calculer_alpha_gradC_carre(DoubleTab& alpha_gradC_carre) const
{
  //Cette methode n'est pas utilisee dans cette version systeme multicomposant.(mr264902)

  const Navier_Stokes_std& eq_ns=ref_cast(Navier_Stokes_std,le_probleme2->equation(0));
  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const DoubleTab& c=eq_c.inconnue().valeurs();
  const Operateur_Grad& opgrad=eq_ns.operateur_gradient();

  DoubleTab& prov_face=ref_cast_non_const( DoubleTab, prov_face_);
  if (prov_face.size()==0)
    prov_face=eq_ns.inconnue().valeurs();
  prov_face=0.;

  // On calcule Grad(c) que l'on met dans prov_face
  opgrad.calculer(c,prov_face);
  eq_ns.solv_masse().appliquer(prov_face);

  // On calcule (Grad(c))^2 sur les faces que l'on met dans gradc2
  DoubleVect gradc2;
  gradc2=prov_face;
  gradc2.carre();

  // Interpolation aux elements de (Grad(c))^2 (a partir de Discretisation_SG_VDF::calculer_tenS_capil)
  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const IntTab& elem_faces = zone_VDF.elem_faces();
  const IntTab& face_sommets = zone_VDF.face_sommets();
  const DoubleTab& positions = zone_VDF.xp();
  const Zone& zone_geom = zone_VDF.zone();
  const Domaine& dom=zone_geom.domaine();
  const int nb_elem = zone_VDF.nb_elem();
  int nb_compo_;
  int elem;
  int f0,f1;
  int som0,som1;
  double psi,val0,val1;

  /*
    DoubleTab& temp_prov_face= ref_cast_non_const(DoubleTab,prov_face_);
    if (temp_prov_face.size()==0)
      temp_prov_face=eq_ns.inconnue().valeurs();
    temp_prov_face=0.;

    DoubleTab prov_face(temp_prov_face.dimension(0),c.line_size());
    prov_face = 0;

    DoubleTab temp_c(c.dimension(0),1);
    temp_c=0;

    for (int j=0; j<c.line_size(); j++)
      {
        for (int i=0; i<c.dimension(0); i++)
          {
            temp_c(i,0)=c(i,j);
            opgrad.calculer(temp_c,temp_prov_face);
            for (int k=0; k<temp_prov_face.dimension(0); k++)
              {
                prov_face(k,j)=temp_prov_face(k,0);
              }
          }
      }

    eq_ns.solv_masse().appliquer(prov_face);
    Cerr<<"prov_face apres eq_ns.solv_masse().appliquer"<<prov_face<<finl;


    // On calcule (Grad(c))^2 sur les faces que l'on met dans gradc2
    DoubleTab gradc2;//modif mr264902
    gradc2=prov_face;
    gradc2.carre();
    Cerr<<"prov_face apres gradc2.carre"<<gradc2<<finl;

    // Interpolation aux elements de (Grad(c))^2 (a partir de Discretisation_SG_VDF::calculer_tenS_capil)
    const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
    const IntTab& elem_faces = zone_VDF.elem_faces();
    const IntTab& face_sommets = zone_VDF.face_sommets();
    const DoubleTab& positions = zone_VDF.xp();
    const Zone& zone_geom = zone_VDF.zone();
    const Domaine& dom=zone_geom.domaine();
    const int nb_elem = zone_VDF.nb_elem();
    int elem;
    int f0,f1;
    int som0,som1;
    double psi,val0,val1;
    Cerr<<"elem_faces "<<elem_faces<<finl;
    Cerr<<"face_sommets "<<face_sommets<<finl;
    Cerr<<"positions "<<positions<<finl;
    Cerr<<"nb_elem "<<nb_elem<<finl;*/

  DoubleTab gradc2_elem(nb_elem,dimension);//***a dimensionner avec le nombre de composants des élémentsconcentration

  // Nombre de composantes du vecteur gradc2
  if(gradc2_elem.nb_dim()==1)
    nb_compo_ = 1;
  else
    nb_compo_ = gradc2_elem.dimension(1);

  // Boucle sur le nombre d'elements
  for(elem=0; elem<nb_elem; elem++)
    {
      // Boucle sur le nombre de composantes du vecteur
      //*** composantes issues de dimension (si 2 donc x et y), faudrait-il mettre prob_dimension à la place de nb_compo et ncomp
      //*** pour eviter toute confusion avec le nombre de composants des elements c1 c2 c3 etc//***
      for(int ncomp=0; ncomp<nb_compo_; ncomp++)
        {
          f0 = elem_faces(elem,ncomp);
          f1 = elem_faces(elem,dimension+ncomp);

          val0 = gradc2(f0);
          val1 = gradc2(f1);

          som0 = face_sommets(f0,0);
          som1 = face_sommets(f1,0);

          psi = ( positions(elem,ncomp) - dom.coord(som0,ncomp) )
                / ( dom.coord(som1,ncomp) - dom.coord(som0,ncomp) ) ;

          if (std::fabs(psi) < 1.e-12)
            gradc2_elem(elem,ncomp) = val0 ;
          else if (std::fabs(1.-psi) < 1.e-12)
            gradc2_elem(elem,ncomp) = val1 ;
          else
            gradc2_elem(elem,ncomp) = val0 + psi * (val1-val0) ;
        }
    }

  // Pour chaque element, on additionne les composantes de gradc2_elem
  int dim;
  double tmp;
  for (elem=0; elem<nb_elem; elem++)
    {
      for (dim=0; dim<dimension; dim++)
        {
          tmp=gradc2_elem(elem,dim);
          alpha_gradC_carre(elem)+=tmp;
        }
    }

  alpha_gradC_carre *= alpha;
  // Ajout B.M suite a plantage dans assert_espace_virtuel_vect dans op_conv...
  alpha_gradC_carre.echange_espace_virtuel();
}



/*! @brief Calcul de la pression thermodynamique aux elements
 *
 * @param (DoubleTab& pression_thermo) pression thermodynamique aux elements
 */
void Source_Con_Phase_field::calculer_pression_thermo(DoubleTab& pression_thermo) const
{
  //Cette methode n'est pas utilisee dans cette version systeme multicomposant.(mr264902)

  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const Navier_Stokes_std& eq_ns=ref_cast(Navier_Stokes_std,le_probleme2->equation(0));
  const DoubleTab& c=eq_c.inconnue().valeurs();
  const DoubleTab& div_alpha_gradC=eq_c.get_div_alpha_gradC();

  // Recuperation de la pression en Pascal de l'etape de projection
  //---------------------------------------------------------------
  const DoubleTab& P_Pa=eq_ns.pression_pa().valeurs();

  const int taille = pression_thermo.size();

  for (int i=0; i<taille; i++)
    {
      pression_thermo(i) = P_Pa(i) - c(i)*div_alpha_gradC(i);
    }

//  for (int j=0; j<pression_thermo.line_size(); j++)
//    {
//      for (int i=0; i<pression_thermo.dimension(0); i++)
//        {
//          pression_thermo(i,j)= P_Pa(i)-c(i,j)*div_alpha_gradC(i,j);
//        }
//    }

}


void Source_Con_Phase_field::calculer_champ_fonc_c(const double t, Champ_Don& champ_fonc_c, const DoubleTab& val_c) const
{
  if (sub_type(Champ_Fonc_Tabule,champ_fonc_c.valeur()))
    {
      const Champ_Fonc_Tabule& ch_champ_fonc_c=ref_cast(Champ_Fonc_Tabule, champ_fonc_c.valeur());
      const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
      const Table& table = ch_champ_fonc_c.table();
      const int isfct = table.isfonction();
      DoubleTab& mes_valeurs = champ_fonc_c.valeur().valeurs();
      // code ci-dessous adapte de Champ_Fonc_Tabule_P0_VDF.mettre_a_jour
      if (!(val_c.nb_dim() == mes_valeurs.nb_dim()))
        {
          Cerr << "Erreur a l'initialisation d'un Champ_Fonc_Tabule" << finl;
          Cerr << "Le champ parametre et le champ a initialiser ne sont pas compatibles" << finl;
          Process::exit();
        }
      if (isfct!=2)
        {
          int nb_elems = zone_VDF.nb_elem();
          if (val_c.line_size() == 1)
            for (int num_elem=0; num_elem<nb_elems; num_elem++)
              mes_valeurs(num_elem,0) = table.val(val_c(num_elem));
          else
            {
              int nb_comps = val_c.nb_dim();
              for (int num_elem=0; num_elem<nb_elems; num_elem++)
                for (int ncomp=0; ncomp<nb_comps; ncomp++)
                  mes_valeurs(num_elem,ncomp) = table.val(val_c(num_elem,ncomp));
            }
        }
      else
        {
          const DoubleTab& centres_de_gravites = zone_VDF.xp();
          table.valeurs(val_c,centres_de_gravites,t,mes_valeurs);
        }
    }
}

/*
const DoubleTab& Source_Con_Phase_field::get_terme_non_lineaire(DoubleTab& non_lineaire)
{
  const Convection_Diffusion_Phase_field& eq_c=ref_cast(Convection_Diffusion_Phase_field,le_probleme2->equation(1));
  const DoubleTab& c=eq_c.inconnue().valeurs();

  non_lineaire(c);

  if (type_potentiel_analytique_==1)
    {
      if (nb_equation_CH==2)
        {
          non_lineaire= (this->*dWdc_naire_analytique_ter) (c,angle_psi,x0Eq,x1Eq,a0Eq,a1Eq);
        }
      else if (nb_equation_CH==3)
        {
          non_lineaire= (this->*dWdc_naire_analytique_quater) (c,angle_psi,angle_phi,x0Eq,x1Eq,x2Eq,a0Eq,a1Eq,a2Eq);
        }
      else
        {
          Cerr <<"potentiel analytique for nb_equation>3 not yet implemented. Use defaut potentiel chimique"<<finl;
        }
    }
  else if (type_potentiel_analytique_==0)
    {
      non_lineaire=(this->*dWdc_naire) (c,eq1,eq2);

      for (int j=0; j<c.line_size(); j++)
        {
          for (int i=0; i<c.dimension(0); i++)
            {
              non_lineaire(i,j)*=betaMatrix(j);
            }
        }
    }
  return non_lineaire;
}*/
