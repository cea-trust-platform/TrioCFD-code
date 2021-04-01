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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Champ_front_ALE_Beam.cpp
// Directory : $BEAM_MODEL_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Champ_front_ALE_Beam.h>
#include <Beam_model.h>
#include <Domaine.h>
#include <Frontiere_dis_base.h>
#include <Domaine_ALE.h>
#include <Domaine.h>
#include <math.h>
#define pi 3.14159265

using namespace std;


Implemente_instanciable( Champ_front_ALE_Beam, "Champ_front_ALE_Beam", Ch_front_var_instationnaire_dep ) ;

Sortie& Champ_front_ALE_Beam::printOn( Sortie& os ) const
{
  return Champ_front_ALE::printOn(os);
}

Entree& Champ_front_ALE_Beam::readOn(Entree& is)
{
  return Champ_front_ALE::readOn(is);
}
void Champ_front_ALE_Beam::remplir_vit_som_bord_ALE(double tps)
{
  //Cerr<<"Champ_front_ALE::remplir_vit_som_bord_ALE"<<finl;
  const Frontiere& front=la_frontiere_dis->frontiere();
  int nb_faces=front.nb_faces();
  const Zone& zone=front.zone();
  const Faces& faces=front.faces();
  const Domaine& domaine=zone.domaine();
  const Domaine_ALE& dom_ale=ref_cast(Domaine_ALE, zone.domaine());
  double x,y,z;
  int nbsf=faces.nb_som_faces();
  int i,j,k;
  int nb_som_tot=domaine.nb_som_tot();
  vit_som_bord_ALE.resize(nb_som_tot,nb_comp());
  vit_som_bord_ALE=0.;
  for( i=0; i<nb_faces; i++)
    {
      x=y=z=0;
      for( k=0; k<nbsf; k++)
        {
          x=domaine.coord(faces.sommet(i,k),0);
          if(dimension>1)
            y=domaine.coord(faces.sommet(i,k),1);
          if(dimension>2)
            z=domaine.coord(faces.sommet(i,k),2);
          DoubleVect phi(3);
          phi=dom_ale.interpolationOnThe3DSurface(x,y,z);
          /*if(abs(phi[0] - (-1.*0.22*(pi/0.7)*cos(pi*x/0.7)*y))>1.e-2)
            Cerr<<"difference phi x= "<<phi[0] - (-1.*0.22*(pi/0.7)*cos(pi*x/0.7)*y)<<finl;
          if(abs(phi[1] - (0.22*sin(pi*x/0.7)))>1.e-2)
            Cerr<<"difference phi y= "<<phi[1] - (0.22*sin(pi*x/0.7))<<finl;
          if (abs(phi[2]) > 1.e-8)
            Cerr<<"difference phi z= "<<phi[2]<<finl;*/
          for( j=0; j<nb_comp(); j++)
            {
              fxyzt[j].setVar("x",x);
              fxyzt[j].setVar("y",y);
              fxyzt[j].setVar("z",z);
              fxyzt[j].setVar("t",tps);
              vit_som_bord_ALE(faces.sommet(i,k),j)=fxyzt[j].eval()*phi[j];
            }

        }
    }
}
