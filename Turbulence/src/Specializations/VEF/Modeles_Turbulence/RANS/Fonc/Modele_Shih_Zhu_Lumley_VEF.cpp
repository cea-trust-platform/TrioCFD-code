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
// File:        Modele_Shih_Zhu_Lumley_VEF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Modeles_Turbulence/RANS/Fonc
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_Shih_Zhu_Lumley_VEF.h>
#include <Zone_VEF.h>
#include <Champ_Uniforme.h>
#include <Champ_P1NC.h>
#include <Periodique.h>
#include <Scatter.h>
#include <Champ_P0_VEF.h>
#include <Discretisation_base.h>
#include <Check_espace_virtuel.h>
#include <LecFicDiffuse.h>
#include <EcritureLectureSpecial.h>

Implemente_instanciable(Modele_Shih_Zhu_Lumley_VEF,"Modele_Shih_Zhu_Lumley_VEF",Modele_Fonc_Realisable_base);


// XD Shih_Zhu_Lumley Modele_Fonc_Realisable_base Shih_Zhu_Lumley -1 Functions necessary to Realizable K-Epsilon Turbulence Model in VEF


///////////////////////////////////////////////////////////////
//   Implementation des fonctions de la classe
///////////////////////////////////////////////////////////////
// printOn et readOn

Sortie& Modele_Shih_Zhu_Lumley_VEF::printOn(Sortie& s ) const
{
  return s;
}

Entree& Modele_Shih_Zhu_Lumley_VEF::readOn(Entree& is )
{
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);
  is_Cmu_constant_ = 0;
  return is;
}

void Modele_Shih_Zhu_Lumley_VEF::set_param(Param& param)
{
  param.ajouter("A0",&A0_);    // XD_ADD_P double value of parameter A0 in U* formula
}


void Modele_Shih_Zhu_Lumley_VEF::Initialisation(const Zone_dis& zone_dis)
{
  const Zone_VEF&  zone_VEF = ref_cast(Zone_VEF,zone_dis.valeur());
  init_tenseur_elem(S_elem_,zone_VEF,2);
  init_tenseur_elem(R_elem_,zone_VEF,2);

  nfaces_ = zone_VEF.nb_faces();
  S_.resize_tab( nfaces_ );
  Cmu_.resize_tab( nfaces_ );
  C1_.resize_tab( nfaces_ );
}

void Modele_Shih_Zhu_Lumley_VEF::Calcul_Tenseurs_S_et_R_elem(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse)
{
  const Zone_VEF&       zone_VEF = ref_cast(Zone_VEF,zone_dis.valeur());
  const Zone_Cl_VEF& zone_Cl_VEF = ref_cast(Zone_Cl_VEF,zone_Cl_dis.valeur());

  DoubleTab gradient_elem;
  init_tenseur_elem(gradient_elem,zone_VEF,2);

  Champ_P1NC::calcul_gradient(vitesse,gradient_elem,zone_Cl_VEF);

  int nelem = S_elem_.dimension(0);

  for (int elem=0; elem<nelem; elem++)
    {
      for (int i=0; i<dimension; i++)
        {
          for (int j=0; j<dimension; j++)
            {
              S_elem_(elem,i,j) =  0.5 * ( gradient_elem(elem,i,j) + gradient_elem(elem,j,i) );
              R_elem_(elem,i,j) =  0.5 * ( gradient_elem(elem,i,j) - gradient_elem(elem,j,i) );
//               if (i==j)
//                 for (int k=0; k<Objet_U::dimension; k++)
//                   {
//                     S_elem_(elem,i,j) -= 1/3 * gradient_elem(elem,k,k);
//                   }
            }
        }
    }

  S_elem_.echange_espace_virtuel();
  R_elem_.echange_espace_virtuel();
}


void Modele_Shih_Zhu_Lumley_VEF::Calcul_Tenseurs_S_et_R_elem_Paroi(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,
                                                                   const DoubleTab& visco_tab, const DoubleTab& visco_turb,
                                                                   const DoubleTab& tab_paroi,const int idt)
{
  const Zone_VEF&       zone_VEF = ref_cast(Zone_VEF,zone_dis.valeur());
  const Zone_Cl_VEF& zone_Cl_VEF = ref_cast(Zone_Cl_VEF,zone_Cl_dis.valeur());

  DoubleTab gradient_elem;
  init_tenseur_elem(gradient_elem,zone_VEF,2);

  Champ_P1NC::calcul_gradient(vitesse,gradient_elem,zone_Cl_VEF);

  if (idt>0)
    Champ_P1NC::calcul_duidxj_paroi(gradient_elem,visco_tab,visco_turb,tab_paroi,zone_Cl_VEF);

  int nelem = S_elem_.dimension(0);

  for (int elem=0; elem<nelem; elem++)
    {
      for (int i=0; i<dimension; i++)
        {
          for (int j=0; j<dimension; j++)
            {
              S_elem_(elem,i,j) =  0.5 * ( gradient_elem(elem,i,j) + gradient_elem(elem,j,i) );
              R_elem_(elem,i,j) =  0.5 * ( gradient_elem(elem,i,j) - gradient_elem(elem,j,i) );
//               if (i==j)
//                 for (int k=0; k<Objet_U::dimension; k++)
//                   {
//                     S_elem_(elem,i,j) -= 1/3 * gradient_elem(elem,k,k);
//                   }
            }
        }
    }

  S_elem_.echange_espace_virtuel();
  R_elem_.echange_espace_virtuel();
}


// Calcul de la norme S SUR LES FACES
void Modele_Shih_Zhu_Lumley_VEF::Calcul_S(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse)
{
  const Zone_VEF&       zone_VEF = ref_cast(Zone_VEF,zone_dis.valeur());
  const Zone_Cl_VEF& zone_Cl_VEF = ref_cast(Zone_Cl_VEF,zone_Cl_dis.valeur());

  DoubleTab S_face;
  init_tenseur_face(S_face,zone_VEF,2);
  calcul_tenseur_face(S_face,S_elem_,zone_VEF,zone_Cl_VEF);


  for (int face=0; face<nfaces_; face++)
    {
      double somme = 0.;
      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          somme += S_face(face,i,j)*S_face(face,i,j);
      S_[face] = sqrt(2.*somme);
    }
}

// Calcul d'un tenseur aux faces a partir d'un tenseur aux elements
DoubleTab& Modele_Shih_Zhu_Lumley_VEF::calcul_tenseur_face(DoubleTab& Tenseur_face, const DoubleTab& Tenseur_elem,
                                                           const Zone_VEF& zone_VEF, const Zone_Cl_VEF& zone_Cl_VEF) const
{
  assert_espace_virtuel_vect(Tenseur_elem);
  const IntTab& face_voisins = zone_VEF.face_voisins();
  int nb_faces = zone_VEF.nb_faces();

  const Conds_lim& les_cl = zone_Cl_VEF.les_conditions_limites();
  int nb_cl=les_cl.size();
  const DoubleVect& volumes = zone_VEF.volumes();

  for (int n_bord=0; n_bord<nb_cl; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
      int ndeb = le_bord.num_premiere_face();
      int nfin = ndeb + le_bord.nb_faces();

      if (sub_type(Periodique,la_cl.valeur()))
        {
          for (int fac=ndeb; fac<nfin; fac++)
            {
              int poly1 = face_voisins(fac,0);
              int poly2 = face_voisins(fac,1);
              double a=volumes(poly1)/(volumes(poly1)+volumes(poly2));
              double b=volumes(poly2)/(volumes(poly1)+volumes(poly2));
              for (int i=0; i<dimension; i++)
                for (int j=0; j<dimension; j++)
                  Tenseur_face(fac,i,j) = a*Tenseur_elem(poly1,i,j) + b*Tenseur_elem(poly2,i,j);
            }
        }
      else
        {
          for (int fac=ndeb; fac<nfin; fac++)
            {
              int poly1 = face_voisins(fac,0);
              for (int i=0; i<dimension; i++)
                for (int j=0; j<dimension; j++)
                  Tenseur_face(fac,i,j) = Tenseur_elem(poly1,i,j);
            }
        }
    }
  int n0 = zone_VEF.premiere_face_int();
  for (int fac = n0; fac<nb_faces; fac++)
    {
      int poly1 = face_voisins(fac,0);
      int poly2 = face_voisins(fac,1);
      double a=volumes(poly1)/(volumes(poly1)+volumes(poly2));
      double b=volumes(poly2)/(volumes(poly1)+volumes(poly2));
      for (int i=0; i<dimension; i++)
        for (int j=0; j<dimension; j++)
          Tenseur_face(fac,i,j) = a*Tenseur_elem(poly1,i,j) + b*Tenseur_elem(poly2,i,j);
    }

  return Tenseur_face;
}

// Initialisation d'une matrice aux elements
void Modele_Shih_Zhu_Lumley_VEF::init_tenseur_elem(DoubleTab& Tenseur, const Zone_VEF& zone_VEF, const int ndim)  const
{
  if(!Tenseur.get_md_vector().non_nul())
    {
      if (ndim==1)
        Tenseur.resize(0, Objet_U::dimension);
      else if (ndim==2)
        Tenseur.resize(0, Objet_U::dimension, Objet_U::dimension);
      zone_VEF.zone().creer_tableau_elements(Tenseur);
    }
  Tenseur = 0.;
}

// Initialisation d'une matrice aux elements
void Modele_Shih_Zhu_Lumley_VEF::init_tenseur_elem(DoubleTab& Tenseur, const Zone_VEF& zone_VEF, const int ndim)
{
  if(!Tenseur.get_md_vector().non_nul())
    {
      if (ndim==1)
        Tenseur.resize(0, Objet_U::dimension);
      else if (ndim==2)
        Tenseur.resize(0, Objet_U::dimension, Objet_U::dimension);
      zone_VEF.zone().creer_tableau_elements(Tenseur);
    }
  Tenseur = 0.;
}


// Initialisation d'une matrice aux faces
void  Modele_Shih_Zhu_Lumley_VEF::init_tenseur_face(DoubleTab& Tenseur,const Zone_VEF& zone_VEF, const int ndim) const
{
  if(!Tenseur.get_md_vector().non_nul())
    {
      if (ndim==1)
        Tenseur.resize(0, Objet_U::dimension);
      else if (ndim==2)
        Tenseur.resize(0, Objet_U::dimension, Objet_U::dimension);
      zone_VEF.creer_tableau_faces(Tenseur);
    }
  Tenseur = 0.;
}

// Initialisation d'une matrice aux faces
void  Modele_Shih_Zhu_Lumley_VEF::init_tenseur_face(DoubleTab& Tenseur,const Zone_VEF& zone_VEF, const int ndim)
{
  if(!Tenseur.get_md_vector().non_nul())
    {
      if (ndim==1)
        Tenseur.resize(0, Objet_U::dimension);
      else if (ndim==2)
        Tenseur.resize(0, Objet_U::dimension, Objet_U::dimension);
      zone_VEF.creer_tableau_faces(Tenseur);
    }
  Tenseur = 0.;
}

void Modele_Shih_Zhu_Lumley_VEF::Calcul_C1 (const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN)
{
  for (int face=0; face<nfaces_; face++)
    {
      double eta;
      // Definition d'un C1 extremum base sur EPS_MIN = 1e-10
//       if (K_Eps(face,1) <= EPS_MIN)
//         eta = S_(face) * K_Eps(face,0)/BR_EPS;
//       else
//         eta = S_(face) * K_Eps(face,0)/K_Eps(face,1);

      eta = S_(face) * K_Eps(face,0)/( K_Eps(face,1) + BR_EPS );

      C1_[face] = std::max( 0.43 , eta / ( 5. + eta ) );
    }

}

void Modele_Shih_Zhu_Lumley_VEF::Calcul_C1_BiK (const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K, const DoubleTab& Eps, const double EPS_MIN)
{
  for (int face=0; face<nfaces_; face++)
    {
      double eta;

      eta = S_(face) * K(face)/( Eps(face) + BR_EPS );

      C1_[face] = std::max( 0.43 , eta / ( 5. + eta ) );
    }

}

void  Modele_Shih_Zhu_Lumley_VEF::Calcul_Cmu_et_S(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse, const DoubleTab& K_Eps, const double EPS_MIN)
{
  const Zone_VEF&       zone_VEF = ref_cast(Zone_VEF,zone_dis.valeur());
  const Zone_Cl_VEF& zone_Cl_VEF = ref_cast(Zone_Cl_VEF,zone_Cl_dis.valeur());

  DoubleTab S_face;
  init_tenseur_face(S_face,zone_VEF,2);
  DoubleTab R_face;
  init_tenseur_face(R_face,zone_VEF,2);

  calcul_tenseur_face(S_face,S_elem_,zone_VEF,zone_Cl_VEF);
  calcul_tenseur_face(R_face,R_elem_,zone_VEF,zone_Cl_VEF);

  DoubleTab U_etoile_face;
  zone_VEF.creer_tableau_faces(U_etoile_face);
  DoubleTab As_face;
  zone_VEF.creer_tableau_faces(As_face);

  for (int face=0; face<nfaces_; face++)
    {
      double somme  = 0.;
      double somme2 = 0.;
      double somme3 = 0.;

      for (int i=0; i<Objet_U::dimension; i++)
        {
          for (int j=0; j<Objet_U::dimension; j++)
            {
              somme  += S_face(face,i,j)*S_face(face,i,j)+R_face(face,i,j)*R_face(face,i,j);
              somme2 += S_face(face,i,j)*S_face(face,i,j);
              for (int k=0; k<Objet_U::dimension; k++)
                {
                  somme3 += S_face(face,i,j)*S_face(face,j,k)*S_face(face,k,i);
                }
            }
        }

      U_etoile_face(face)  = sqrt(somme);
      double S_tilde       = sqrt(somme2);
      S_[face]             = sqrt(2.*somme2);
      double val_cosinus   = sqrt(6.) * somme3 / ( S_tilde * S_tilde * S_tilde +1.e-20 );

      if ( val_cosinus > 1. )
        {
          val_cosinus = 1. ;
        }
      else if ( val_cosinus < -1. )
        {
          val_cosinus = -1. ;
        }

      As_face(face)       = sqrt(6.) * cos ( (1./3.) * acos( val_cosinus ) );

//       // Definition d'un Cmu extremum base sur EPS_MIN = 1e-10
//       if (K_Eps(face,1) <= EPS_MIN)
//         Cmu_[face] = 1./(A0_+As_face(face)*U_etoile_face(face)*K_Eps(face,0)/BR_EPS);
//       else
//         Cmu_[face] = 1./(A0_+As_face(face)*U_etoile_face(face)*K_Eps(face,0)/K_Eps(face,1));

      Cmu_[face] = 1./(A0_+As_face(face)*U_etoile_face(face)*K_Eps(face,0)/( K_Eps(face,1) + BR_EPS));

    }
}

void  Modele_Shih_Zhu_Lumley_VEF::Calcul_Cmu_et_S_BiK(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse, const DoubleTab& K, const DoubleTab& Eps, const double EPS_MIN)
{
  const Zone_VEF&       zone_VEF = ref_cast(Zone_VEF,zone_dis.valeur());
  const Zone_Cl_VEF& zone_Cl_VEF = ref_cast(Zone_Cl_VEF,zone_Cl_dis.valeur());

  DoubleTab S_face;
  init_tenseur_face(S_face,zone_VEF,2);
  DoubleTab R_face;
  init_tenseur_face(R_face,zone_VEF,2);

  calcul_tenseur_face(S_face,S_elem_,zone_VEF,zone_Cl_VEF);
  calcul_tenseur_face(R_face,R_elem_,zone_VEF,zone_Cl_VEF);

  DoubleTab U_etoile_face;
  zone_VEF.creer_tableau_faces(U_etoile_face);
  DoubleTab As_face;
  zone_VEF.creer_tableau_faces(As_face);

  for (int face=0; face<nfaces_; face++)
    {
      double somme  = 0.;
      double somme2 = 0.;
      double somme3 = 0.;

      for (int i=0; i<Objet_U::dimension; i++)
        {
          for (int j=0; j<Objet_U::dimension; j++)
            {
              somme  += S_face(face,i,j)*S_face(face,i,j)+R_face(face,i,j)*R_face(face,i,j);
              somme2 += S_face(face,i,j)*S_face(face,i,j);
              for (int k=0; k<Objet_U::dimension; k++)
                {
                  somme3 += S_face(face,i,j)*S_face(face,j,k)*S_face(face,k,i);
                }
            }
        }

      U_etoile_face(face)  = sqrt(somme);
      double S_tilde       = sqrt(somme2);
      S_[face]             = sqrt(2.*somme2);
      double val_cosinus   = sqrt(6.) * somme3 / ( S_tilde * S_tilde * S_tilde +1.e-20 );

      if ( val_cosinus > 1. )
        {
          val_cosinus = 1. ;
        }
      else if ( val_cosinus < -1. )
        {
          val_cosinus = -1. ;
        }

      As_face(face)       = sqrt(6.) * cos ( (1./3.) * acos( val_cosinus ) );

//       // Definition d'un Cmu extremum base sur EPS_MIN = 1e-10
//       if (K_Eps(face,1) <= EPS_MIN)
//         Cmu_[face] = 1./(A0_+As_face(face)*U_etoile_face(face)*K_Eps(face,0)/BR_EPS);
//       else
//         Cmu_[face] = 1./(A0_+As_face(face)*U_etoile_face(face)*K_Eps(face,0)/K_Eps(face,1));

      Cmu_[face] = 1./(A0_+As_face(face)*U_etoile_face(face)*K(face)/( Eps(face) + BR_EPS));

    }
}


void  Modele_Shih_Zhu_Lumley_VEF::associer(const Zone_dis& zone_dis,
                                           const Zone_Cl_dis& zone_Cl_dis)
{
  la_zone_VEF = ref_cast(Zone_VEF,zone_dis.valeur());
  la_zone_Cl_VEF = ref_cast(Zone_Cl_VEF,zone_Cl_dis.valeur());

  Initialisation( zone_dis );
}

void  Modele_Shih_Zhu_Lumley_VEF::mettre_a_jour(double temps)
{
  ;
}

void Modele_Shih_Zhu_Lumley_VEF::Contributions_Sources(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN)
{
  Calcul_Tenseurs_S_et_R_elem(zone_dis,zone_Cl_dis,vitesse);
  Calcul_Cmu_et_S(zone_dis,zone_Cl_dis,vitesse,K_Eps,EPS_MIN);
  Calcul_C1(zone_dis,zone_Cl_dis,vitesse,K_Eps,EPS_MIN);
}

void Modele_Shih_Zhu_Lumley_VEF::Contributions_Sources_Paroi(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN,
                                                             const DoubleTab& visco, const DoubleTab& visco_turb,const DoubleTab& loi_paroi,const int idt)
{
  Calcul_Tenseurs_S_et_R_elem_Paroi(zone_dis,zone_Cl_dis,vitesse,visco,visco_turb,loi_paroi,idt);
  Calcul_Cmu_et_S(zone_dis,zone_Cl_dis,vitesse,K_Eps,EPS_MIN);
  Calcul_C1(zone_dis,zone_Cl_dis,vitesse,K_Eps,EPS_MIN);
}

void Modele_Shih_Zhu_Lumley_VEF::Contributions_Sources_BiK(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K, const DoubleTab& Eps, const double EPS_MIN)
{
  Calcul_Tenseurs_S_et_R_elem(zone_dis,zone_Cl_dis,vitesse);
  Calcul_Cmu_et_S_BiK(zone_dis,zone_Cl_dis,vitesse,K,Eps,EPS_MIN);
  Calcul_C1_BiK(zone_dis,zone_Cl_dis,vitesse,K,Eps,EPS_MIN);
}

void Modele_Shih_Zhu_Lumley_VEF::Contributions_Sources_Paroi_BiK(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K, const DoubleTab& Eps, const double EPS_MIN,
                                                                 const DoubleTab& visco, const DoubleTab& visco_turb,const DoubleTab& loi_paroi,const int idt)
{
  Calcul_Tenseurs_S_et_R_elem_Paroi(zone_dis,zone_Cl_dis,vitesse,visco,visco_turb,loi_paroi,idt);
  Calcul_Cmu_et_S_BiK(zone_dis,zone_Cl_dis,vitesse,K,Eps,EPS_MIN);
  Calcul_C1_BiK(zone_dis,zone_Cl_dis,vitesse,K,Eps,EPS_MIN);
}



