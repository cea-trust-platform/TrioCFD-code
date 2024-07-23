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
// File:        Source_Con_Phase_field.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Multiphase/Phase_field/src/VDF
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Source_Con_Phase_field_included
#define Source_Con_Phase_field_included

#include <Source_Con_Phase_field_base.h>
#include <Matrice_Morse.h>
#include <Champ_Don.h>
#include <Table.h>
#include <Convection_Diffusion_Concentration.h>
#include <Equation_base.h>
#include <math.h>
#include <SolveurSys.h>
#include <EChaine.h>
#include <TRUST_Ref.h>

class Domaine_Cl_dis;
class Domaine_Cl_VDF;
class Domaine_VDF;


#define PI 3.14159265

class Source_Con_Phase_field : public Source_Con_Phase_field_base
{

  Declare_instanciable(Source_Con_Phase_field);

public:
  void associer_pb(const Probleme_base& ) override;
  DoubleTab& ajouter(DoubleTab& ) const override;
  DoubleTab& calculer(DoubleTab& ) const override;
  void mettre_a_jour(double ) override;
  void premier_demi_dt() override;

  inline const DoubleVect& get_u_carre() override;
  inline const DoubleVect& get_maxX();
  inline const DoubleVect& get_minX();

  inline const int& get_type_systeme_naire();

  inline double drhodc(const int n_elem) const;

  void calculer_champ_fonc_c(const double t, Champ_Don& champ_fonc_c, const DoubleTab& val_c) const override;

protected:
  int tpsaff;
  double rho0;
  REF(Champ_Don) drhodc_;
  DoubleTab accr;
  DoubleVect u_carre_;
  DoubleTab prov_face_,prov_elem_;
  double alpha, beta;
  DoubleVect alphaMatrix, betaMatrix, minX, maxX;
  DoubleVect angle_psi,angle_phi,x0Eq,x1Eq,x2Eq,a0Eq,a1Eq,a2Eq;
  int nb_equation_CH;
  Table potentiel_chimique_expr_;
  double (Source_Con_Phase_field::*dWdc)(const double) const;
  DoubleTab (Source_Con_Phase_field::*dWdc_naire)(const DoubleTab&, const DoubleVect&, const DoubleVect&) const;
  DoubleTab (Source_Con_Phase_field::*dWdc_naire_analytique_ter)(const DoubleTab& , const DoubleVect&, const DoubleVect& , const DoubleVect&, const DoubleVect& , const DoubleVect& ) const;
  DoubleTab (Source_Con_Phase_field::*dWdc_naire_analytique_quater)(const DoubleTab& , const DoubleVect& , const DoubleVect& ,const DoubleVect& , const DoubleVect& , const DoubleVect& , const DoubleVect& , const DoubleVect&, const DoubleVect&) const ;
  DoubleVect eq1, eq2;
  double kappa;
  DoubleVect kappaMatrix;
  DoubleTab (Source_Con_Phase_field::*kappaMatrix_func_c)(const DoubleTab&, const DoubleVect&) const;
  DoubleVect coeff_auto_diffusion;
  double temperature;
  double molarVolume;
  double prefactor;
  double alpha_ref;
  double angle_alphaMatrix;
  double diagonal_coeff;
  int kappa_ind;
  int type_kappa_ = -1;
  int type_kappa_auto_diffusion_ = -1;
  int type_alpha_rotation_;
  int type_potentiel_analytique_;
  int type_systeme_naire_;
  int kappa_moy_;
  double mult_kappa;
  Table kappa_forme_expr_;
  double (Source_Con_Phase_field::*kappa_func_c)(const double) const;
  double mu1, mu2;
  int implicitation_;
  int gmres_;
  double epsilon_;
  double eps_;
  int boussi_;
  int diff_boussi_;
  DoubleTab g_;
  int mutype_;
  int couplage_;
  int nkr, nit;
  double rec_min, rec_max, epsGMRES;

  void associer_domaines(const Domaine_dis& ,const Domaine_Cl_dis& ) override;
  DoubleTab& laplacien(const DoubleTab&, DoubleTab&) const;
  DoubleTab& div_kappa_grad(const DoubleTab&, const DoubleTab&, DoubleTab&) const;
  void calculer_alpha_gradC_carre(DoubleTab&) const;
  void calculer_div_alpha_rho_gradC(DoubleTab&) const;
  void calculer_div_alpha_gradC(DoubleTab&) const;
  void calculer_pression_thermo(DoubleTab&) const;
  void calculer_mutilde_demi(DoubleTab&, DoubleTab&) const;
  //const DoubleTab& get_terme_non_lineaire(DoubleTab&);
  virtual void assembler_matrice_point_fixe(Matrice_Morse&);
  virtual void calculer_point_fixe(const DoubleTab&, const DoubleTab&, const Matrice_Morse&, DoubleTab&, DoubleTab&);
  virtual void calculer_u2_elem(DoubleVect&);
  virtual void calculer_mutilde(DoubleTab&) const;
  double dWdc_defaut(const double) const;
  double kappa_func_c_defaut(const double) const;
  double dWdc_general(const double) const;
  double kappa_func_c_general(const double) const;
  DoubleTab kappa_func_auto_diffusion(const DoubleTab&, const DoubleVect& ) const;
  DoubleTab dWdc_naire_defaut(const DoubleTab&, const DoubleVect&, const DoubleVect&) const;
  DoubleTab dWdc_naire_analytique_ternaire(const DoubleTab& , const DoubleVect&, const DoubleVect& , const DoubleVect&, const DoubleVect& , const DoubleVect& ) const;
  DoubleTab dWdc_naire_analytique_quaternaire(const DoubleTab& , const DoubleVect& , const DoubleVect& ,const DoubleVect& , const DoubleVect& , const DoubleVect& , const DoubleVect& , const DoubleVect&, const DoubleVect&) const ;

  virtual void construire_systeme(const DoubleTab&, const Matrice_Morse&, DoubleTab&, const DoubleTab&);
  virtual void matvect(const DoubleTab&, const Matrice_Morse&, const DoubleTab&, const DoubleTab&, DoubleTab&);
  // Modifie par DJ
  //---------------
  virtual int non_lin_gmres(const DoubleTab&, const DoubleTab&, const Matrice_Morse&, DoubleTab&, DoubleTab&);
  /*   virtual int non_lin_gmres(const DoubleTab&, const DoubleTab&, const Matrice_Morse&, DoubleTab&); */
  //---------------

  REF(Probleme_base) le_probleme2;
  REF(Domaine_VDF) le_dom_VDF;
  REF(Domaine_Cl_VDF) le_dom_Cl_VDF;
};

inline const DoubleVect& Source_Con_Phase_field::get_u_carre()
{
  return u_carre_;
}
inline const int& Source_Con_Phase_field::get_type_systeme_naire()
{
  return type_systeme_naire_;
}
inline double Source_Con_Phase_field::drhodc(const int n_elem) const
{
  const DoubleTab& tab = drhodc_->valeurs();
  if (tab.dimension(0)==1)
    {
      int dim = tab.nb_dim();
      double value = 0.0;
      switch(dim)
        {
        case 2:
          value = tab(0,0);
          break;
        case 3:
          value = tab(0,0,0);
          break;
        default:
          Cerr <<"Source_Con_Phase_field::drhodc : Problem with drhodc_ dimension:"<<dim<<finl;
          exit();
          break;
        }
      return value;
    }
  else
    {
      return tab(n_elem);
    }
}

inline double Source_Con_Phase_field::dWdc_defaut(const double c) const
{
  return 4*c*(c+0.5)*(c-0.5);
  //return 2*c*(1-c)*(1-2*c);

  //return 2*(c-0)*(c-0.68)*(2*c-0-0.68);

}

inline DoubleTab Source_Con_Phase_field::dWdc_naire_defaut(const DoubleTab& c, const DoubleVect& a, const DoubleVect& b) const
{
  DoubleTab dWdc_(c);
  dWdc_=0;
  for (int j=0; j<c.line_size(); j++)
    {
      for (int i=0; i<c.dimension(0); i++)
        {
          dWdc_(i,j) = 2*(c(i,j)-a(j))*(c(i,j)-b(j))*(2*c(i,j)-a(j)-b(j));//4*c(i,j)*(c(i,j)+0.5)*(c(i,j)-0.5);
        }
    }
  return dWdc_;
}

inline DoubleTab Source_Con_Phase_field::dWdc_naire_analytique_ternaire(const DoubleTab& c, const DoubleVect& psi_, const DoubleVect& cUEq, const DoubleVect& cZrEq, const DoubleVect& aUEq, const DoubleVect& aZrEq) const
{
  DoubleTab dWdc_(c);
  dWdc_=0;
  DoubleTab P(c.dimension(0),2);
  P = 0;
  DoubleTab dPdxU(P);
  DoubleTab dPdxZr(P);


  //here j defines oxide and metal
  for (int j=0; j<2; j++)
    {
      for (int i=0; i<c.dimension(0); i++)
        {
          P(i,j) = pow(((cos(psi_(j)* PI / 180.0)*(c(i,0)-cUEq(j))+sin(psi_(j)* PI / 180.0)*(c(i,1)-cZrEq(j)))/(aUEq(j))),2)+pow(((-sin(psi_(j)* PI / 180.0)*(c(i,0)-cUEq(j))+cos(psi_(j)* PI / 180.0)*(c(i,1)-cZrEq(j)))/(aZrEq(j))),2) ;
        }
    }

  //here j defines oxide and metal
  for (int j=0; j<2; j++)
    {
      for (int i=0; i<c.dimension(0); i++)
        {
          dPdxU(i,j) = 2*(c(i,0)-cUEq(j))*((pow(cos(psi_(j)* PI / 180.0)/(aUEq(j)),2)+pow(sin(psi_(j)* PI / 180.0)/(aZrEq(j)),2)))+ 2*(c(i,1)-cZrEq(j))*((cos(psi_(j)* PI / 180.0)*sin(psi_(j)* PI / 180.0))/(pow(aUEq(j),2))-(sin(psi_(j)* PI / 180.0)*cos(psi_(j)* PI / 180.0))/(pow(aZrEq(j),2)));
          //dPdxZr(i,j) = 2*(c(i,0)-cUEq(j))*((cos(psi_(j)* PI / 180.0)*sin(psi_(j)* PI / 180.0))/(pow(aUEq(j),2))-(sin(psi_(j)* PI / 180.0)*cos(psi_(j)* PI / 180.0))/(pow(aZrEq(j),2)))+ 2*(c(i,1)-cZrEq(j))*((pow(cos(psi_(j)* PI / 180.0)/(aUEq(j)),2)+pow(sin(psi_(j)* PI / 180.0)/(aZrEq(j)),2)));
          dPdxZr(i,j) = 2*(-(c(i,0) - cUEq(j))*sin(psi_(j)* PI / 180.0) + (c(i,1) - cZrEq(j))*cos(psi_(j)* PI / 180.0))*cos(psi_(j)* PI / 180.0)/pow(aZrEq(j),2) + 2*((c(i,0) - cUEq(j))*cos(psi_(j)* PI / 180.0) + (c(i,1) - cZrEq(j))*sin(psi_(j)* PI / 180.0))*sin(psi_(j)* PI / 180.0)/pow(aUEq(j),2);

        }
    }
  for (int i=0; i<c.dimension(0); i++)
    {
      dWdc_(i,0)=P(i,0)*dPdxU(i,1)+P(i,1)*dPdxU(i,0); //mutildeU *** ///// +muUeq a rajouter ??
      dWdc_(i,1)=P(i,0)*dPdxZr(i,1)+P(i,1)*dPdxZr(i,0);//mutildeZr **** //// +muZrEq a rajouter ??
    }

  return dWdc_;
}

inline DoubleTab Source_Con_Phase_field::dWdc_naire_analytique_quaternaire(const DoubleTab& c, const DoubleVect& psi_, const DoubleVect& phi_,const DoubleVect& cUEq, const DoubleVect& cZrEq, const DoubleVect& cFeEq, const DoubleVect& aUEq, const DoubleVect& aZrEq, const DoubleVect& aFeEq) const
{
  DoubleTab dWdc_(c);
  dWdc_=0;
  DoubleTab P(c.dimension(0),2);
  P = 0;
  DoubleTab dPdx0(P);
  DoubleTab dPdx1(P);
  DoubleTab dPdx2(P);

  for (int j=0; j<2; j++)
    {
      for (int i=0; i<c.dimension(0); i++)
        {
          P(i,j) = pow(((cos(psi_(j)* PI / 180.0)*(c(i,0)-cUEq(j))+sin(psi_(j)* PI / 180.0)*cos(phi_(j)* PI /180.0)*(c(i,1)-cZrEq(j))+sin(psi_(j)* PI / 180.0)*sin(phi_(j)* PI / 180.0)*(c(i,2)-cFeEq(j)))/(aUEq(j))),2)+pow(((-sin(psi_(j)* PI / 180.0)*(c(i,0)-cUEq(j))+cos(psi_(j)* PI / 180.0)*cos(phi_(j)* PI / 180.0)*(c(i,1)-cZrEq(j))+cos(psi_(j)* PI / 180.0)*sin(phi_(j)* PI / 180.0)*(c(i,2)-cFeEq(j)))/(aZrEq(j))),2)+pow(((-sin(phi_(j)* PI / 180.0)*(c(i,1)-cZrEq(j))+cos(phi_(j)* PI / 180.0)*(c(i,2)-cFeEq(j)))/(aFeEq(j))),2) ;
        }
    }
  for (int j=0; j<2; j++)
    {
      for (int i=0; i<c.dimension(0); i++)
        {
          dPdx0(i,j) = 2*(c(i,0)-cUEq(j))*((pow(cos(psi_(j)* PI / 180.0)/(aUEq(j)),2)+pow(sin(psi_(j)* PI / 180.0)/(aZrEq(j)),2)))+ 2*(c(i,1)-cZrEq(j))*((cos(psi_(j)* PI / 180.0)*sin(psi_(j)* PI / 180.0)*cos(phi_(j)* PI / 180.0))/(pow(aUEq(j),2))-(cos(psi_(j)* PI / 180.0)*sin(psi_(j)* PI / 180.0)*cos(phi_(j)* PI / 180.0))/(pow(aZrEq(j),2)))+2*(c(i,2)-cFeEq(j))*((cos(psi_(j)* PI / 180.0)*sin(psi_(j)* PI / 180.0)*sin(phi_(j)* PI / 180.0))/(pow(aUEq(j),2))-(cos(psi_(j)* PI / 180.0)*sin(psi_(j)* PI / 180.0)*sin(phi_(j)* PI / 180.0))/(pow(aZrEq(j),2)));
          dPdx1(i,j) = 2*(c(i,0)-cUEq(j))*((sin(psi_(j)* PI / 180.0)*cos(phi_(j)* PI / 180.0)*cos(phi_(j)* PI / 180.0))/(pow(aUEq(j),2))-(sin(psi_(j)* PI / 180.0)*cos(phi_(j)* PI / 180.0)*cos(phi_(j)* PI / 180.0))/(pow(aZrEq(j),2)))+ 2*(c(i,1)-cZrEq(j))*((pow(sin(psi_(j)* PI / 180.0)*cos(phi_(j)* PI / 180.0)/(aUEq(j)),2)+pow(cos(psi_(j)* PI / 180.0)*cos(phi_(j)* PI / 180.0)/(aZrEq(j)),2))+pow(sin(phi_(j)* PI / 180.0)/(aFeEq(j)),2))+2*(c(i,2)-cFeEq(j))*((pow(sin(psi_(j)* PI / 180.0),2)*cos(phi_(j)* PI / 180.0)*sin(phi_(j)* PI / 180.0)/(pow(aUEq(j),2)))+pow(cos(psi_(j)* PI / 180.0),2)*cos(phi_(j)* PI / 180.0)*sin(phi_(j)* PI / 180.0)/(pow(aZrEq(j),2))+(sin(phi_(j)* PI / 180.0)*cos(phi_(j)* PI / 180.0)/(pow(aFeEq(j),2))));
          dPdx2(i,j) = 2*(c(i,0)-cUEq(j))*((sin(psi_(j)* PI / 180.0)*sin(phi_(j)* PI / 180.0)*cos(psi_(j)* PI / 180.0))/(pow(aUEq(j),2))-(sin(psi_(j)* PI / 180.0)*sin(phi_(j)* PI / 180.0)*cos(psi_(j)* PI / 180.0))/(pow(aZrEq(j),2)))+2*(c(i,1)-cZrEq(j))*((pow(sin(psi_(j)* PI / 180.0),2)*sin(phi_(j)* PI / 180.0)*cos(phi_(j)* PI / 180.0)/(pow(aUEq(j),2)))+pow(cos(psi_(j)* PI / 180.0),2)*sin(phi_(j)* PI / 180.0)*cos(phi_(j)* PI / 180.0)/(pow(aZrEq(j),2))+(-sin(phi_(j)* PI / 180.0)*cos(phi_(j)* PI / 180.0)/(pow(aFeEq(j),2))))+ 2*(c(i,2)-cFeEq(j))*(pow(sin(psi_(j)* PI / 180.0)*sin(phi_(j)* PI / 180.0)/(aUEq(j)),2)+pow(cos(psi_(j)* PI / 180.0)*sin(phi_(j)* PI / 180.0)/(aZrEq(j)),2)+pow(cos(phi_(j)* PI / 180.0)/(aFeEq(j)),2));
        }
    }
  for (int i=0; i<c.dimension(0); i++)
    {
      dWdc_(i,0)=P(i,0)*dPdx0(i,1)+P(i,1)*dPdx0(i,0); //+muUeq a rajouter
      dWdc_(i,1)=P(i,0)*dPdx1(i,1)+P(i,1)*dPdx1(i,0);//+muZrEq a rajouter
      dWdc_(i,2)=P(i,0)*dPdx2(i,1)+P(i,1)*dPdx2(i,0);//+muFeEq a rajouter
    }

  return dWdc_;
}
inline double Source_Con_Phase_field::kappa_func_c_defaut(const double c) const
{
  return std::max(mult_kappa*kappa*(c+0.5)*(0.5-c),0.);
}
inline DoubleTab Source_Con_Phase_field::kappa_func_auto_diffusion(const DoubleTab& c, const DoubleVect& coeff_auto_diff) const
{
  DoubleTab temp_c(c.dimension(0),c.line_size()+1);
  DoubleTab temp_kappa(c.dimension(0), nb_equation_CH, nb_equation_CH); //define temp to avoid resize here
  DoubleTab kappaMat=temp_kappa;
  kappaMat = 0;
  DoubleTab kappa_ij (temp_kappa);
  DoubleTab somme_c(c.dimension(0));
  somme_c=0;
  DoubleTab kappaMatrix_(c.dimension(0), nb_equation_CH*nb_equation_CH);

  // creation d'un tableau temp_c qui contiendra toutes les concentrations des n composants non pas les n-1 seulements (n-1 eqs).
  //La derniere colonne correspondra a 1-somme sur (n-1) des concentrations

  for (int k=0; k<c.line_size(); k++)
    {
      for (int l=0; l<c.dimension(0); l++)
        {
          temp_c(l,k)= c(l,k);
          somme_c(l) += c(l,k);
          temp_c(l,c.line_size())=1.0-somme_c(l);
        }
    }

  // calcul de la mobilite kappa dans un doubleTab (n1,n2,n3)
  for(int k=0; k< nb_equation_CH+1; k++)
    {
      for (int j=0; j<nb_equation_CH; j++)
        {
          for (int i=0; i<nb_equation_CH; i++)
            {
              for (int l=0; l<c.dimension(0); l++)
                {
                  //            kappa_ij(l,i,j) = (j==k) ? 1-temp_c(l,j) : -temp_c(l,j);
                  //            kappa_ij(l,i,j) *= (i==k) ? 1-temp_c(l,i) : -temp_c(l,i);
                  //            kappa_ij(l,i,j) *= temp_c(l,k)*coeff_auto_diff(k);
                  //
                  if (j==k)
                    {
                      kappa_ij(l,i,j)=1-temp_c(l,j);
                    }
                  else
                    {
                      kappa_ij(l,i,j)=-temp_c(l,j);
                    }
                  if (i==k)
                    {
                      kappa_ij(l,i,j)*=1-temp_c(l,i);
                    }
                  else
                    {
                      kappa_ij(l,i,j)*=-temp_c(l,i);
                    }
                  kappa_ij(l,i,j) *= temp_c(l,k)*coeff_auto_diff(k);
                }
            }
        }
      kappaMat+=kappa_ij;
    }

  //retour vers un DoubleTab (n1,n2)
  for (int j=0; j<nb_equation_CH; j++)
    for (int i=0; i<nb_equation_CH; i++)
      for (int l=0; l<c.dimension(0); l++)
        {
          kappaMatrix_(l,i+nb_equation_CH*j)=kappaMat(l,i,j);
        }
  return kappaMatrix_;

}

inline double Source_Con_Phase_field::dWdc_general(const double c) const
{
  return potentiel_chimique_expr_.val(c,0);
}

inline double Source_Con_Phase_field::kappa_func_c_general(const double c) const
{
  return std::max(mult_kappa*kappa*kappa_forme_expr_.val(c,0),0.);
}
inline const DoubleVect& Source_Con_Phase_field::get_minX()
{
  return minX;
}
inline const DoubleVect& Source_Con_Phase_field::get_maxX()
{
  return maxX;
}

#endif
