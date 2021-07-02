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
#include <Zone_Cl_dis_base.h>
#include <Front_VF.h>
#include <Domaine_ALE.h>
#include <Domaine.h>
#include <Zone_VEF.h>
#include <Cond_lim.h>
#include <Navier_Stokes_std.h>
#include <Equation_base.h>
#include <Operateur_Diff.h>
#include <Operateur_Grad.h>
#include <DoubleVect.h>
#include <fstream>
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
void Champ_front_ALE_Beam::initializationBeam_bis(double tps)
{
  Cerr<<"Champ_front_ALE_Beam::initializationBeam"<<finl;
  const Frontiere& front=la_frontiere_dis->frontiere();
  const Zone& zone=front.zone();
  Domaine_ALE& dom_ale=ref_cast_non_const(Domaine_ALE, zone.domaine());
  double x=0,y=0,z=0;
  double velocity=0.;
  int dir=dom_ale.getBeamDirection();
  fxyzt[dir].setVar("x",x);
  fxyzt[dir].setVar("y",y);
  fxyzt[dir].setVar("z",z);
  fxyzt[dir].setVar("t",tps);
  velocity=fxyzt[dir].eval();
  dom_ale.initializationBeam(velocity);
}



void Champ_front_ALE_Beam::remplir_vit_som_bord_ALE(double tps)
{
  Cerr<<"Champ_front_ALE_Beam::remplir_vit_som_bord_ALE "<<finl;
  const Frontiere& front=la_frontiere_dis->frontiere();
  Cerr<<" front nom := "<<front.le_nom()<<finl;
  int nb_faces=front.nb_faces();
  const Zone& zone=front.zone();
  const Faces& faces=front.faces();

  const Domaine& domaine=zone.domaine();
  Domaine_ALE& dom_ale=ref_cast_non_const(Domaine_ALE, zone.domaine());

  double dt = dom_ale.get_dt();
  const int nbModes=dom_ale.getBeamNbModes();

  double x,y,z;
  int nbsf=faces.nb_som_faces();
  int i,j,k;
  int nb_som_tot=domaine.nb_som_tot();
  vit_som_bord_ALE.resize(nb_som_tot,nb_comp());
  vit_som_bord_ALE=0.;
  const DoubleVect& beamVelocity=dom_ale.getBeamVelocity(tps, dt);
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


          DoubleVect value(3);
          value=0.;
          for(int count=0; count<nbModes; count++ )
            {
              const DoubleTab& u=dom_ale.getBeamDisplacement(count);
              const DoubleTab& R=dom_ale.getBeamRotation(count);
              DoubleVect phi(3);
              phi=dom_ale.interpolationOnThe3DSurface(x,y,z, u, R);
              for(int comp=0; comp<nb_comp(); comp++)
                {
                  value[comp] +=beamVelocity[count]*phi[comp];
                }
            }

          for( j=0; j<nb_comp(); j++)
            {
              vit_som_bord_ALE(faces.sommet(i,k),j)=value[j];

            }



          if(front.le_nom()=="Poutre")
            {
              if(i ==2713 && k==1)
                {
                  std::ofstream ofs_P;
                  ofs_P.open ("beam_vel_P_0.5_h.txt", std::ofstream::out | std::ofstream::app);
                  ofs_P<<tps<<" "<<value[0]<<" "<<value[1]<<" "<<value[2]<<endl;
                  ofs_P.close();
                  //Cout<<"time = "<<tps <<" velocity beam Poutre= "<<beamVelocity<<finl;
                  //getchar();

                }
            }
        }

    }
}
void  Champ_front_ALE_Beam::computeFluidForce(const double tps, DoubleVect& fluidForce)
{

  const Frontiere& front=la_frontiere_dis->frontiere();
  const Zone& zone=front.zone();
  Domaine_ALE& dom_ale=ref_cast_non_const(Domaine_ALE, zone.domaine());
  const Equation_base& eqn =  dom_ale.getEquation();
  const Navier_Stokes_std& eqn_hydr = ref_cast(Navier_Stokes_std,eqn);
  int ndeb = front.num_premiere_face();
  int nfin = ndeb + front.nb_faces();
  const Operateur_base& op_grad= eqn_hydr.operateur_gradient().l_op_base();
  const Operateur_base& op_diff= eqn_hydr.operateur_diff().l_op_base();
  const Zone_VEF& la_zone_vef=ref_cast(Zone_VEF,op_grad.equation().zone_dis().valeur());
  const DoubleTab& xv=la_zone_vef.xv();

  DoubleTab& flux_bords_grad=op_grad.flux_bords();
  DoubleTab& flux_bords_diff=op_diff.flux_bords();
  DoubleVect force(3);
  DoubleVect force_p(3);
  DoubleVect force_v(3);
  std::ofstream ofs_k0;
  ofs_k0.open ("force_p.txt", std::ofstream::out | std::ofstream::app);

  std::ofstream ofs_k1;
  ofs_k1.open ("force_v.txt", std::ofstream::out | std::ofstream::app);
  Cout<<" BEAM! ndeb= "<<ndeb<<" nfin = "<<nfin<<finl;
  getchar();
  if(flux_bords_grad.size() == flux_bords_diff.size())
    {

      for (int face=ndeb; face<nfin; face++)
        {
          for(int comp=0; comp<3; comp++)
            {
              force_p[comp] += (flux_bords_grad(face, comp));
              force_v[comp] += (flux_bords_diff(face, comp));
            }
        }
      ofs_k0<<tps<<" "<<force_p[0]<<" "<<force_p[1]<<" "<<force_p[2]<<endl;
      ofs_k0.close();
      ofs_k1<<tps<<" "<<force_v[0]<<" "<<force_v[1]<<" "<<force_v[2]<<endl;
      ofs_k1.close();

      Cout<<"ALE BEAM: temps:= "<<tps<<" force v =  "<<force_v[1]<<finl;
      getchar();
    }

  if(flux_bords_grad.size() == flux_bords_diff.size())
    {
      for(int nbmodes=0; nbmodes<fluidForce.size(); nbmodes++)
        {
          for (int face=ndeb; face<nfin; face++)
            {
              DoubleVect phi(3);
              const DoubleTab& u=dom_ale.getBeamDisplacement(nbmodes);
              const DoubleTab& R=dom_ale.getBeamRotation(nbmodes);
              phi=dom_ale.interpolationOnThe3DSurface(xv(face,0),xv(face,1),xv(face,2), u, R);
              for(int comp=0; comp<3; comp++)
                {
                  force[comp] += (flux_bords_grad(face, comp)+ flux_bords_diff(face, comp))*phi[comp];
                }
            }
          for(int comp=0; comp<3; comp++)
            {
              fluidForce[nbmodes] += force[comp];
            }
        }

    }

  fluidForce=0.;
}
