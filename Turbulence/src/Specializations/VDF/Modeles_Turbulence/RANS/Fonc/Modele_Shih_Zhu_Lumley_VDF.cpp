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
// File:        Modele_Shih_Zhu_Lumley_VDF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Modeles_Turbulence/RANS/Fonc
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_Shih_Zhu_Lumley_VDF.h>
#include <Zone_VDF.h>
#include <Zone_Cl_VDF.h>
#include <TRUSTTrav.h>
#include <Dirichlet.h>
#include <Dirichlet_homogene.h>
#include <Dirichlet_paroi_defilante.h>
#include <Echange_externe_impose.h>
#include <Periodique.h>
#include <Symetrie.h>
#include <Neumann.h>
#include <Neumann_homogene.h>
#include <Champ_Face_VDF.h>
#include <Champ_Uniforme.h>


Implemente_instanciable(Modele_Shih_Zhu_Lumley_VDF,"Modele_Shih_Zhu_Lumley_VDF",Modele_Fonc_Realisable_base);

// XD Modele_Shih_Zhu_Lumley_VDF Modele_Fonc_Realisable_base Modele_Shih_Zhu_Lumley_VDF -1 Functions necessary to Realizable K-Epsilon Turbulence Model in VDF


///////////////////////////////////////////////////////////////
//   Implementation des fonctions de la classe
///////////////////////////////////////////////////////////////
// printOn et readOn

Sortie& Modele_Shih_Zhu_Lumley_VDF::printOn(Sortie& s ) const
{
  return s;
}

Entree& Modele_Shih_Zhu_Lumley_VDF::readOn(Entree& is )
{
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);
  is_Cmu_constant_ = 0;
  return is;
}

void Modele_Shih_Zhu_Lumley_VDF::set_param(Param& param)
{
  param.ajouter("A0",&A0_);    // XD_ADD_P double value of parameter A0 in U* formula
}


void Modele_Shih_Zhu_Lumley_VDF::Initialisation(const Zone_dis& zone_dis)
{
  const Zone_VDF&  zone_VDF = ref_cast(Zone_VDF,zone_dis.valeur());

  nelem_ = zone_VDF.nb_elem();

  S_.resize_tab( nelem_ );
  Cmu_.resize_tab( nelem_ );
  C1_.resize_tab( nelem_ );
}

// Calcul de la norme S SUR LES ELEMENTS
void Modele_Shih_Zhu_Lumley_VDF::Calcul_S(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vit)
{
  const Zone_VDF&       zone_VDF = ref_cast(Zone_VDF,zone_dis.valeur());
  const Zone_Cl_VDF& zone_Cl_VDF = ref_cast(Zone_Cl_VDF,zone_Cl_dis.valeur());

  int nb_elem_tot=zone_VDF.nb_elem_tot();
  DoubleTab gij(nb_elem_tot,dimension,dimension);
  const Champ_Face_VDF& vitesse = ref_cast(Champ_Face_VDF,eq_hydraulique->inconnue().valeur() );

  ref_cast_non_const(Champ_Face_VDF,vitesse).calcul_duidxj( vitesse.valeurs(),gij,zone_Cl_VDF );

  for (int elem=0; elem<nelem_; elem++)
    {
      double somme2 = 0.;

      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          {
            double Sij = 0.5*( gij(elem,i,j) + gij(elem,j,i) ) ;
//             if (i==j)
//               {
//                 for (int l=0; l<Objet_U::dimension; l++)
//                   {
//                     Sij -= 1/3 * gij(elem,l,l);
//                   }
//               }
            somme2 += Sij*Sij;
          }

      S_( elem ) = sqrt(2.*somme2);
    }
}



void Modele_Shih_Zhu_Lumley_VDF::Calcul_C1 (const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN)
{
#ifdef __INTEL_COMPILER
#pragma novector // Desactive vectorisation sur Intel car crash sinon
#endif
  for (int elem=0; elem<nelem_; elem++)
    {
      double eta;
      // Definition d'un C1 extremum base sur EPS_MIN = 1e-10
      if (K_Eps(elem,1) <= EPS_MIN)
        eta = S_(elem) * K_Eps(elem,0)/BR_EPS;
      else
        eta = S_(elem) * K_Eps(elem,0)/K_Eps(elem,1);

      C1_[elem] = std::max( 0.43 , eta / ( 5. + eta ) );
    }

}



void Modele_Shih_Zhu_Lumley_VDF::Calcul_C1_BiK (const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K,const DoubleTab& Eps, const double EPS_MIN)
{
#ifdef __INTEL_COMPILER
#pragma novector // Desactive vectorisation sur Intel car crash sinon
#endif
  for (int elem=0; elem<nelem_; elem++)
    {
      double eta;
      // Definition d'un C1 extremum base sur EPS_MIN = 1e-10
      if (Eps(elem) <= EPS_MIN)
        eta = S_(elem) * K(elem)/BR_EPS;
      else
        eta = S_(elem) * K(elem)/Eps(elem);

      C1_[elem] = std::max( 0.43 , eta / ( 5. + eta ) );
    }

}

void  Modele_Shih_Zhu_Lumley_VDF::Calcul_Cmu_et_S(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vit, const DoubleTab& K_Eps, const double EPS_MIN)
{
  const Zone_VDF&       zone_VDF = ref_cast(Zone_VDF,zone_dis.valeur());
  const Zone_Cl_VDF& zone_Cl_VDF = ref_cast(Zone_Cl_VDF,zone_Cl_dis.valeur());

  int nb_elem_tot=zone_VDF.nb_elem_tot();
  DoubleTab gij(nb_elem_tot,dimension,dimension);
  const Champ_Face_VDF& vitesse = ref_cast(Champ_Face_VDF,eq_hydraulique->inconnue().valeur() );

  ref_cast_non_const(Champ_Face_VDF,vitesse).calcul_duidxj( vitesse.valeurs(),gij,zone_Cl_VDF );

  DoubleTab U_etoile( nelem_);
  DoubleTab As( nelem_);

  for (int elem=0; elem<nelem_; elem++)
    {
      double somme  = 0.;
      double somme2 = 0.;
      double somme3 = 0.;

      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          {
            double Sij = 0.5*( gij(elem,i,j) + gij(elem,j,i) ) ;
            double Rij = 0.5*( gij(elem,i,j) - gij(elem,j,i) ) ;
//             if (i==j)
//               {
//                 for (int l=0; l<Objet_U::dimension; l++)
//                   {
//                     Sij -= 1/3 * gij(elem,l,l);
//                   }
//               }

            somme  += Sij*Sij+Rij*Rij;
            somme2 += Sij*Sij;

            for (int k=0; k<Objet_U::dimension; k++)
              {
                double Sjk = 0.5*( gij(elem,j,k) + gij(elem,k,j) ) ;
                double Ski = 0.5*( gij(elem,k,i) + gij(elem,i,k) ) ;

                somme3 += Sij*Sjk*Ski;
              }
          }

      U_etoile( elem )     = sqrt(somme);
      double S_tilde       = sqrt(somme2);
      S_( elem )           = sqrt(2.*somme2);
      double val_cosinus   = sqrt(6.) * somme3 / ( S_tilde * S_tilde * S_tilde +1.e-20 );

      if ( val_cosinus > 1. )
        {
          val_cosinus = 1. ;
        }
      else if ( val_cosinus < -1. )
        {
          val_cosinus = -1. ;
        }

      As( elem )        = sqrt(6.) * cos( (1./3.) * acos( val_cosinus ) );

      // Definition d'un Cmu extremum base sur EPS_MIN = 1e-10
//       if (K_Eps(elem,1) <= EPS_MIN)
//         Cmu_[elem] = 1./(A0_+As(elem)*U_etoile(elem)*K_Eps(elem,0)/BR_EPS);
//       else
//         Cmu_[elem] = 1./(A0_+As(elem)*U_etoile(elem)*K_Eps(elem,0)/K_Eps(elem,1));

      Cmu_[elem] = 1./(A0_+As(elem)*U_etoile(elem)*K_Eps(elem,0)/( K_Eps(elem,1) + BR_EPS ));

    }
}

void  Modele_Shih_Zhu_Lumley_VDF::Calcul_Cmu_et_S_BiK(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vit, const DoubleTab& K, const DoubleTab&  Eps, const double EPS_MIN)
{
  const Zone_VDF&       zone_VDF = ref_cast(Zone_VDF,zone_dis.valeur());
  const Zone_Cl_VDF& zone_Cl_VDF = ref_cast(Zone_Cl_VDF,zone_Cl_dis.valeur());

  int nb_elem_tot=zone_VDF.nb_elem_tot();
  DoubleTab gij(nb_elem_tot,dimension,dimension);
  const Champ_Face_VDF& vitesse = ref_cast(Champ_Face_VDF,eq_hydraulique->inconnue().valeur() );

  ref_cast_non_const(Champ_Face_VDF,vitesse).calcul_duidxj( vitesse.valeurs(),gij,zone_Cl_VDF );

  DoubleTab U_etoile( nelem_);
  DoubleTab As( nelem_);

  for (int elem=0; elem<nelem_; elem++)
    {
      double somme  = 0.;
      double somme2 = 0.;
      double somme3 = 0.;

      for (int i=0; i<Objet_U::dimension; i++)
        for (int j=0; j<Objet_U::dimension; j++)
          {
            double Sij = 0.5*( gij(elem,i,j) + gij(elem,j,i) ) ;
            double Rij = 0.5*( gij(elem,i,j) - gij(elem,j,i) ) ;
//             if (i==j)
//               {
//                 for (int l=0; l<Objet_U::dimension; l++)
//                   {
//                     Sij -= 1/3 * gij(elem,l,l);
//                   }
//               }

            somme  += Sij*Sij+Rij*Rij;
            somme2 += Sij*Sij;

            for (int k=0; k<Objet_U::dimension; k++)
              {
                double Sjk = 0.5*( gij(elem,j,k) + gij(elem,k,j) ) ;
                double Ski = 0.5*( gij(elem,k,i) + gij(elem,i,k) ) ;

                somme3 += Sij*Sjk*Ski;
              }
          }

      U_etoile( elem )     = sqrt(somme);
      double S_tilde       = sqrt(somme2);
      S_( elem )           = sqrt(2.*somme2);
      double val_cosinus   = sqrt(6.) * somme3 / ( S_tilde * S_tilde * S_tilde +1.e-20 );

      if ( val_cosinus > 1. )
        {
          val_cosinus = 1. ;
        }
      else if ( val_cosinus < -1. )
        {
          val_cosinus = -1. ;
        }

      As( elem )        = sqrt(6.) * cos( (1./3.) * acos( val_cosinus ) );

      // Definition d'un Cmu extremum base sur EPS_MIN = 1e-10
//       if (K_Eps(elem,1) <= EPS_MIN)
//         Cmu_[elem] = 1./(A0_+As(elem)*U_etoile(elem)*K_Eps(elem,0)/BR_EPS);
//       else
//         Cmu_[elem] = 1./(A0_+As(elem)*U_etoile(elem)*K_Eps(elem,0)/K_Eps(elem,1));

      Cmu_[elem] = 1./(A0_+As(elem)*U_etoile(elem)*K(elem)/( Eps(elem) + BR_EPS ));

    }
}


void  Modele_Shih_Zhu_Lumley_VDF::associer(const Zone_dis& zone_dis,
                                           const Zone_Cl_dis& zone_Cl_dis)
{
  la_zone_VDF = ref_cast(Zone_VDF,zone_dis.valeur());
  la_zone_Cl_VDF = ref_cast(Zone_Cl_VDF,zone_Cl_dis.valeur());

  Initialisation( zone_dis );
}

void  Modele_Shih_Zhu_Lumley_VDF::mettre_a_jour(double temps)
{
  ;
}


void Modele_Shih_Zhu_Lumley_VDF::Contributions_Sources(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN)
{
  Calcul_Cmu_et_S(zone_dis,zone_Cl_dis,vitesse,K_Eps,EPS_MIN);
  Calcul_C1(zone_dis,zone_Cl_dis,vitesse,K_Eps,EPS_MIN);
}

void Modele_Shih_Zhu_Lumley_VDF::Contributions_Sources_Paroi(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K_Eps, const double EPS_MIN,
                                                             const DoubleTab& visco_tab, const DoubleTab& visco_turb,const DoubleTab& tab_paroi,const int idt)
{
  // identique a Contributions_Sources pour l'instant
  Calcul_Cmu_et_S(zone_dis,zone_Cl_dis,vitesse,K_Eps,EPS_MIN);
  Calcul_C1(zone_dis,zone_Cl_dis,vitesse,K_Eps,EPS_MIN);
}


void Modele_Shih_Zhu_Lumley_VDF::Contributions_Sources_BiK(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K,const DoubleTab& Eps, const double EPS_MIN)
{
  Calcul_Cmu_et_S_BiK(zone_dis,zone_Cl_dis,vitesse,K,Eps,EPS_MIN);
  Calcul_C1_BiK(zone_dis,zone_Cl_dis,vitesse,K,Eps,EPS_MIN);
}

void Modele_Shih_Zhu_Lumley_VDF::Contributions_Sources_Paroi_BiK(const Zone_dis& zone_dis, const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& vitesse,const DoubleTab& K,const DoubleTab& Eps, const double EPS_MIN,
                                                                 const DoubleTab& visco_tab, const DoubleTab& visco_turb,const DoubleTab& tab_paroi,const int idt)
{
  // identique a Contributions_Sources pour l'instant
  Calcul_Cmu_et_S_BiK(zone_dis,zone_Cl_dis,vitesse,K,Eps,EPS_MIN);
  Calcul_C1_BiK(zone_dis,zone_Cl_dis,vitesse,K,Eps,EPS_MIN);
}


