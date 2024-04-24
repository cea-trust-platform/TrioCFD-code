/****************************************************************************
* Copyright (c) 2015 - 2016, CEA
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
#if 0
#include <Turbulence_Wale_IJK.h>


void Turbulence_Wale_IJK::compute(const IJK_Field_ST& vx, const IJK_Field_ST& vy, const IJK_Field_ST& vz,
                                  IJK_Field_ST& nu_t, IJK_Field_ST& k)
{
  // Reimplemented from Modele_turbulence_hyd_LES_Wale_VDF::calculer_OP1_OP2
  const int dimension = 3;
  // compute duidxj
  // here

  for(i = 0; i < dimension; i++)
    for(j = 0; j < dimension; j++)
      gij2(i,j) = duidxj(i,0) * duidxj(0,j) + duidxj(i,1) * duidxj(1,j) + duidxj(i,2) * duidxj(2,j);

  gkk = gij2(0,0) + gij2(1,1) + gij2(2,2);

  for(i = 0; i < dimension; i++)
    for(j = 0; j < dimension; j++)
      {
        sd(i,j) = 0.5 * (gij2(i,j)+gij2(j,i));
        if(i==j)
          sd(i,j) -= gkk2 * one_third;
      }

  // Calcul de sd2 et Sij2
  sd2=0.;
  Sij2=0.;

  for ( i=0 ; i<dimension ; i++)
    for ( j=0 ; j<dimension ; j++)
      {
        sd2 += sd(i,j)*sd(i,j);
        // Deplacement du calcul de sij
        Sij = 0.5 * (duidxj(i,j) + duidxj(j,i));


        // PQ : 24/01/07 : le stencil de Sij est par contruction de :
        // -  1 maille pour les termes diagonaux Sii
        // - ~2 mailles pour les termes croises Sij
        //
        // Wale s'appuyant a la fois sur sd2 (porte par Sij) et sur Sij2 (porte principalement par Sii)
        // est sensible a cette difference de stencil.
        // En portant le stencil a 3 maille spour le calcul de Sii, on retrouve en THI
        // le bon taux de dissipation ainsi que des spectres possedant la bonne allure en k^-5/3.
        //
        // A traiter : Quid sur canal plan ???

        if(i==j)  //!< augmentation du stencil de Sii
          {
            Sij = duidxj_ijk(mlkqjsdmklcjq);



            face1=elem_faces(elem,i);
            face2=elem_faces(elem,i+dimension);

            elem1=face_voisins(face1,0);
            elem2=face_voisins(face2,1);

            // if(elem1==elem) elem1=face_voisins(face1,1);  // par construction il n'y a pas besoin
            // if(elem2==elem) elem2=face_voisins(face2,0);  // par construction il n'y a pas besoin

            // si pas de bord a proximit\ufffd on passe au stencil de 3 mailles
            // sinon on reste au stencil \ufffd 1 maille

            if( elem1>=0 && elem2>=0 ) Sij = ((duidxj(elem1,i,i)+duidxj(elem,i,i)+duidxj(elem2,i,i)))/3.;
          }


        Sij2+=Sij*Sij;
      }

  // Calcul de OP1 et OP2
  OP1(elem)=pow(sd2,1.5);
  OP2(elem)=pow(Sij2,2.5)+pow(sd2,1.25);

}
#endif
