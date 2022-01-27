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
/////////////////////////////////////////////////////////////////////////////
//
// File      : estimateur_aposteriori_P0_VEF.cpp
// Directory : $APOSTERIORI_FOR_TRIOCFD_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <estimateur_aposteriori_P0_VEF.h>
#include <Zone_VEF.h>
#include <Champ_P1NC.h>
#include <Champ_P0_VEF.h>
#include <Champ_Uniforme.h>
#include <Champ_P1_isoP1Bulle.h>
#include <Zone_Cl_VEF.h>

Implemente_instanciable( estimateur_aposteriori_P0_VEF, "estimateur_aposteriori_P0_VEF", Champ_Fonc_P0_VEF ) ;

Sortie& estimateur_aposteriori_P0_VEF::printOn( Sortie& s ) const
{
  return s << que_suis_je() << " " << le_nom();
}

Entree& estimateur_aposteriori_P0_VEF::readOn( Entree& s )
{
  return s;
}

void estimateur_aposteriori_P0_VEF::associer_champ(const Champ_P1NC& la_vitesse, const Champ_P0_VEF& la_pression, const Champ_Uniforme& la_viscosite_cinematique, const Zone_Cl_dis_base& la_zone_Cl_dis_base)
{
  la_zone_Cl_VEF  = ref_cast(Zone_Cl_VEF, la_zone_Cl_dis_base);
  vitesse_= la_vitesse;
  pression_p0_ = la_pression;
  viscosite_cinematique_ = la_viscosite_cinematique;
}

void estimateur_aposteriori_P0_VEF::associer_champ(const Champ_P1NC& la_vitesse, const Champ_P1_isoP1Bulle& la_pression, const Champ_Uniforme& la_viscosite_cinematique, const Zone_Cl_dis_base& la_zone_Cl_dis_base)
{
  la_zone_Cl_VEF  = ref_cast(Zone_Cl_VEF, la_zone_Cl_dis_base);
  vitesse_= la_vitesse;
  pression_p1isop1b_ = la_pression;
  viscosite_cinematique_ = la_viscosite_cinematique;
}

void estimateur_aposteriori_P0_VEF::mettre_a_jour(double tps)
{
  const Zone_Cl_VEF& zone_cl_VEF = la_zone_Cl_VEF.valeur();
  const Zone_VEF& zone_VEF = zone_cl_VEF.zone_VEF();

  const DoubleTab& face_normales = zone_VEF.face_normales();
  const int nb_faces = zone_VEF.nb_faces_tot();
  const IntTab& face_voisins = zone_VEF.face_voisins();
  int premiere_face_int = zone_VEF.premiere_face_int();
  double sautdirectionnel;
  double estimtot;

  int dimension=Objet_U::dimension;
  int fac,elem1,elem2;
  int i,j,ielem;
  int nb_elem = zone_VEF.nb_elem();
  cout << "estimateur_aposteriori_P0_VEF::mettre_a_jour nb_elem = " << nb_elem << endl;
  int nb_elem_tot = zone_VEF.nb_elem_tot();
  cout << "estimateur_aposteriori_P0_VEF::mettre_a_jour  nb_elem_tot = " << nb_elem_tot << endl;

  DoubleTab gradient_elem(nb_elem_tot,dimension,dimension);
  gradient_elem=0.;

// DoubleVect visc(nb_elem_tot);
  double visc;
// visc = viscosite_cinematique_.valeur().valeurs();
  visc = viscosite_cinematique_.valeur()(0,0);
  Cout << "visc = " << visc;

// const DoubleTab& les_coords = zone_VEF.xv();
  // DoubleTab xv(nb_faces,dimension);
//  zone_VEF.calculer_centres_gravite_aretes(xa);

//   DoubleTab vite_essai(nb_faces,dimension);
//   vite_essai=0.;
//   for (fac=0; fac<nb_faces; fac++)
//     {
//       for (i=0; i<dimension; i++)
//         {
//           for(j=0; j<dimension; j++)
//             {
//               vite_essai(fac,i) += (3.*i+j+1.)*les_coords(fac,j) ;
// //               cout << "face = " << fac << endl;
// //               cout << "j    = " << j << endl;
// //               cout << "coo  = " << les_coords(fac,j) << endl;
// //               cout << "vit  = " << vite_essai(fac,i) << endl;
//             }
//         }
//     }
//
//   Champ_P1NC::calcul_gradient(vite_essai, gradient_elem, la_zone_Cl_VEF.valeur());
//   for (ielem=0; ielem<nb_elem; ielem++)
//     {
//       for (i=0; i<dimension; i++)
//         {
//           for(j=0; j<dimension; j++)
//             {
//               if ( fabs( gradient_elem(ielem,i,j) - (3.*i+j+1.) ) > 1.e-6 )
//                 {
//                   cout << "ielem    =  " << ielem << endl;
//                   cout << "(i,j)    =  " << i << " , " << j << endl;
//                   cout << "gradient =  " << gradient_elem(ielem,i,j) << endl;
//                   cout << "attendu =   " << (3.*i+j+1.) << endl;
//                 }
//             }
//         }
//     }
//
//   for (fac=premiere_face_int; fac<nb_faces; fac++)
//     {
//       elem1=face_voisins(fac,0);
//       elem2=face_voisins(fac,1);
//       if ( (elem1>=0) && (elem2>=0) )
//         {
//           for (i=0; i<dimension; i++)
//             {
//               sautdirectionnel = 0.;
//               for (j=0; j<dimension; j++)
//                 {
//                   sautdirectionnel += ( gradient_elem(elem1,i,j) - gradient_elem(elem2,i,j) ) * face_normales(fac,j);
//                   if( fabs(sautdirectionnel) > 1.e-13 )
//                     {
//                       cout << "saut grad     = " << ( gradient_elem(elem1,i,j) - gradient_elem(elem2,i,j) ) << endl;
//                       cout << "face_normales = " <<  face_normales(fac,j) << endl;
//                     }
//                 }
//               if( fabs(sautdirectionnel) > 1.e-13 )
//                 {
//                   cout << "elem1, elem2 " << elem1 << " , " << elem2 << endl;
//                   cout << "sautdirectionnel avant = " << sautdirectionnel << endl;
//                 }
// //              sautdirectionnel -= ( la_pression(elem1) - la_pression(elem2) ) * face_normales(fac,i);
//               sautdirectionnel *= sautdirectionnel;
// //              cout << "sautdirectionnel apres = " << sautdirectionnel << endl;
// //               cout << "elem1 = " << elem1 << "estimelem1 avant = " << estim(elem1) << endl;
// //               cout << "elem2 = " << elem2 << "estimelem2 avant = " << estim(elem2) << endl;
// //               estim(elem1) += sautdirectionnel;
// //               estim(elem2) += sautdirectionnel;
// //               cout << "elem1 = " << elem1 << "estimelem1 apres = " << estim(elem1) << endl;
// //               cout << "elem2 = " << elem2 << "estimelem2 apres = " << estim(elem2) << endl;
//             }
//
//         }
//     }

  gradient_elem=0.;
  Champ_P1NC::calcul_gradient(vitesse_.valeur().valeurs(), gradient_elem, la_zone_Cl_VEF.valeur());


  DoubleTab& estim = valeurs(); // Estimateur
  estim = 0.;

//  if ( tps > 0.4 )
  {
    for (fac=premiere_face_int; fac<nb_faces; fac++)
      {
        elem1=face_voisins(fac,0);
        elem2=face_voisins(fac,1);
//        cout << "viscosite elem1 " << visc(elem1) << endl;
//        cout << "viscosite elem2 " << visc(elem2) << endl;
        if ( (elem1>=0) && (elem2>=0) )
          {
            for (i=0; i<dimension; i++)
              {
                sautdirectionnel = 0.;
                for (j=0; j<dimension; j++)
                  {
                    sautdirectionnel += ( visc*gradient_elem(elem1,i,j) - visc*gradient_elem(elem2,i,j) ) * face_normales(fac,j);
//                     cout << "saut grad     = " << ( gradient_elem(elem1,i,j) - gradient_elem(elem2,i,j) ) << endl;
//                     cout << "face_normales = " <<  face_normales(fac,j) << endl;
                  }
//                 cout << "sautdirectionnel avant = " << sautdirectionnel << endl;
//              sautdirectionnel -= ( la_pression(elem1) - la_pression(elem2) ) * face_normales(fac,i);
                sautdirectionnel *= sautdirectionnel;
//                   cout << "sautdirectionnel apres = " << sautdirectionnel << endl;
//                   cout << "elem1 = " << elem1 << "estimelem1 avant = " << estim(elem1) << endl;
//                   cout << "elem2 = " << elem2 << "estimelem2 avant = " << estim(elem2) << endl;
                estim(elem1) += sautdirectionnel;
                estim(elem2) += sautdirectionnel;
//                   cout << "elem1 = " << elem1 << "estimelem1 apres = " << estim(elem1) << endl;
//                   cout << "elem2 = " << elem2 << "estimelem2 apres = " << estim(elem2) << endl;
              }

          }
      }
  }
  estimtot = 0.;

  for (ielem=0; ielem<nb_elem; ielem++)
    {
      estimtot += estim(ielem);
      estim(ielem) = sqrt(estim(ielem));
    }
  estimtot = sqrt(estimtot);
  cout << "estimateur_aposteriori_P0_VEF::mettre_a_jour - estimateur = " << estimtot << endl;
// Champ_P1NC::calculer_gradientP1NC(vitesse_.valeur().valeurs(),la_zone_VEF.valeur(),la_zone_Cl_VEF.valeur(),gradient_elem);

//   DoubleVect tmp(nb_elem);
//   vitesse_->calcul_S_barre(vitesse_.valeur().valeurs(),tmp,la_zone_Cl_VEF.valeur());
//   DoubleTab& S = valeurs(); // Shear rate
//   for (int i=0; i<nb_elem; i++)
//   {
//
//   }
  changer_temps(tps);
  Champ_Fonc_base::mettre_a_jour(tps);
}