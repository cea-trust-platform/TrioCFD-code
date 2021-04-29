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
// File      : Beam_model.cpp
// Directory : $BEAM_MODEL_ROOT/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Beam_model.h>
#include <Domaine_ALE.h>
#include <DoubleVect.h>
#include <IntTab.h>
#include <Probleme_base.h>
#include <Zone_VEF.h>
#include <Zone_Cl_VEF.h>
#include <Frontiere.h>
#include <MD_Vector_tools.h>
#include <MD_Vector_std.h>
#include <Ref_Zone.h>
#include <Faces.h>
#include <fstream>
#include <math.h>       /* sin */
#define PI 3.14159265

using namespace std;

Implemente_liste(DoubleTab);

Implemente_instanciable( Beam_model, "Beam_model", Interprete_geometrique_base ) ;

Sortie& Beam_model::printOn( Sortie& os ) const
{
  return Interprete::printOn(os);
}

Entree& Beam_model::readOn( Entree& is )
{
  return Interprete::readOn(is);
}

Entree&  Beam_model::interpreter_(Entree& is)
{
  associer_domaine(is);
  Domaine_ALE& dom=ref_cast(Domaine_ALE, domaine());
  dom.reading_beam_model(is);
  return is;
}

void Beam_model::readInputMassStiffnessFiles(Nom& masse_and_stiffness_file_name)
{
  Cerr << "Read beam model coefficients from "<< masse_and_stiffness_file_name <<" files " << finl;
  mass_.resize(nbModes_);
  stiffness_.resize(nbModes_);
  damping_.resize(nbModes_);
  string const nomFichier(masse_and_stiffness_file_name);

  ifstream monFlux(nomFichier.c_str());

  if(monFlux)
    {
      string line;  // read the first line
      getline(monFlux, line);
      monFlux.ignore();
      double mass, stiffness, damping;
      for(int i=0; i<nbModes_; i++)
        {
          monFlux >> mass >> stiffness >> damping;
          if(abs(mass) > 1.e-16 && abs(stiffness)> 1.e-16)
            {
              mass_[i]=mass;
              stiffness_[i]=stiffness;
              damping_[i]=damping;
            }
          else
            {
              Cerr<<" Beam_model::read_input_files values ​​of mass and stiffness are wrong in the file "<<masse_and_stiffness_file_name<<finl;
              Cerr<<" Value ​​of mass = "<< mass<< " and value of stiffness = "<<stiffness<<finl;
              exit();
            }
        }

      monFlux.close();
    }
  else
    {
      Cerr<< "ERROR: Unable to open the file." <<masse_and_stiffness_file_name<<finl;
    }
}
void Beam_model::readInputAbscFiles(Nom& absc_file_name)
{
  Cerr << "Read beam coordinates from "<< absc_file_name <<" files " << finl;

  string const nomFichier(absc_file_name);
  ifstream monFlux(nomFichier.c_str());
  if(monFlux)
    {
      string line;  // read the first line
      getline(monFlux, line);
      monFlux.ignore();
      double abscissa;
      int i=0;
      while(monFlux)
        {
          monFlux >> abscissa;
          abscissa_.resize(i+1);
          abscissa_[i]=abscissa;
          i++;
        }

      monFlux.close();
    }
  else
    {
      Cerr<< "ERROR: Unable to open the file." <<absc_file_name<<finl;
    }
}



void Beam_model::readInputModalDeformation(Nom& modal_deformation_file_name)
{

  Cerr << "Read beam modal deformation coefficients from "<< modal_deformation_file_name <<" files " << finl;

  string const nomFichier(modal_deformation_file_name);

  ifstream monFlux(nomFichier.c_str());
  int size = abscissa_.size();
  u_.resize(size, 3);
  R_.resize(size, 3);

  if(monFlux)
    {
      string line;  // read the first line
      getline(monFlux, line);
      monFlux.ignore();
      double  ux, uy, uz, rx, ry, rz;
      for(int i=0; i<size; i++)
        {
          monFlux >> ux>> uy>> uz>> rx>> ry>> rz;
          u_(i, 0)=ux;
          u_(i, 1)=uy;
          u_(i, 2)=uz;
          R_(i, 0)=rx;
          R_(i, 1)=ry;
          R_(i, 2)=rz;

        }
      monFlux.close();
    }
  else
    {
      Cerr<< "ERROR: Unable to open the file." <<modal_deformation_file_name<<finl;
    }
}

void Beam_model::initialization(double velocity)
{

  qSpeed_.resize(nbModes_);
  qAcceleration_.resize(nbModes_);
  qDisplacement_.resize(nbModes_);
  qHalfSpeed_.resize(nbModes_);
  qSpeed_=velocity;
  qHalfSpeed_=velocity;
  qAcceleration_=0.;
  qDisplacement_=0.;
}

DoubleVect Beam_model::NewmarkSchemeFD (double dt, double fluidForce)
{

  double halfDt=dt/2;
  for(int j=0; j < nbModes_; j++)
    {
      qHalfSpeed_[j] = qSpeed_[j] + halfDt*qAcceleration_[j];
      qDisplacement_[j] += dt*qHalfSpeed_[j];
      qAcceleration_[j]= (fluidForce - damping_[j]*qHalfSpeed_[j] - stiffness_[j]*qDisplacement_[j])/(mass_[j] + halfDt*damping_[j]);
      qSpeed_[j] = qHalfSpeed_[j] + halfDt*qAcceleration_[j];
      //Cout<<" qHalfSpeed_[j] "<<qHalfSpeed_[j]<<"qDisplacement_[j] "<<qDisplacement_[j]<<"qAcceleration_[j] "<<qAcceleration_[j]<<"qSpeed_[j] "<<qSpeed_[j]<<finl;
      //getchar();
    }

  return qHalfSpeed_;
}

DoubleVect Beam_model::NewmarkSchemeMA (double dt, double fluidForce)
{
  double halfDt=dt/2;
  double squareHalfDt= halfDt*halfDt;
  for(int j=0; j < nbModes_; j++)
    {
      double PreviousqAcceleration= qAcceleration_[j];
      double coeff1 = mass_[j] + halfDt*damping_[j] + squareHalfDt*stiffness_[j];
      double coeff2 = damping_[j]*(qSpeed_[j] + halfDt*qAcceleration_[j]) + stiffness_[j]*(qDisplacement_[j] + dt*qSpeed_[j] + squareHalfDt*qAcceleration_[j]);
      qAcceleration_[j]=(fluidForce - coeff2)/coeff1;
      qDisplacement_[j] += dt*qSpeed_[j] + squareHalfDt*(PreviousqAcceleration + qAcceleration_[j]);
      qSpeed_[j] += halfDt*(PreviousqAcceleration + qAcceleration_[j]);
      //Cout<<"qDisplacement_[j] "<<qDisplacement_[j]<<"qAcceleration_[j] "<<qAcceleration_[j]<<"qSpeed_[j] "<<qSpeed_[j]<<finl;
      //getchar();
    }
  return qSpeed_;
}

DoubleVect Beam_model::interpolationOnThe3DSurface(const double& x, const double& y, const double& z) const
{
  DoubleVect phi(3);
  phi=0.;
  double h = abscissa_[1] -abscissa_[0]; //1d mesh pitch
  int abscissa_size = abscissa_.size();
  double s=0.;
  double xs=x;
  double ys=y;
  double zs=z;
  if (direction_== 0)
    {
      s = x;
      xs=0.;
    }
  else if (direction_== 1)
    {
      s = y;
      ys=0.;
    }

  else
    {
      s = z;
      zs=0.;
    }

  int i, j ;
  i = s/h;
  if((i+1) < abscissa_size)
    {
      j= i+1;
    }
  else
    {
      j=i;
    }
  double ux, uy, uz, Rx, Ry, Rz;

  //linear interpolation between points i and j
  double alpha, betha ;
  if (i==j)
    {
      alpha=1.;
      betha=0.;
    }
  else if(abs(abscissa_[i] - s)< 1.e-5 )
    {
      alpha=1.;
      betha=0.;
    }
  else if (abs(abscissa_[j] - s)< 1.e-5 )
    {
      alpha=0.;
      betha=1.;
    }
  else
    {
      alpha = (abscissa_[j] - s)/h;
      betha = (s - abscissa_[i])/h;
    }

  ux=alpha*u_(i, 0) + betha*u_(j, 0);
  uy=alpha*u_(i, 1) + betha*u_(j, 1);
  uz=alpha*u_(i, 2) + betha*u_(j, 2);
  Rx=alpha*R_(i, 0) + betha*R_(j, 0);
  Ry=alpha*R_(i, 1) + betha*R_(j, 1);
  Rz=alpha*R_(i, 2) + betha*R_(j, 2);

  phi[0] =ux + Ry*zs -Rz*ys;
  phi[1] =uy + Rz*xs -Rx*zs;
  phi[2] =uz + Rx*ys -Ry*xs;

  return phi;
}
/*void Beam_model::interpolationOnThe3DSurface(const Bords& les_bords_ALE)
{
  Cerr<<"Interpolation of the 1d modal deformation to the 3d surface"<<finl;

  if(les_bords_ALE.size()!= 1)
    {
      Cerr<<"Too many beams. For the moment one can treat only a beam with the Beam module!"<<finl;
      exit();
    }
  if(dimension != 3)
    {
      Cerr<<"The dimension of the problem is different from 3"<<finl;
      exit();
    }

  Cerr<<" Interpolation on the "<<les_bords_ALE(0).le_nom()<<" boundary"<<finl;

  const Faces& faces = les_bords_ALE(0).faces();
  int size=faces.nb_faces_tot();
  phi3D_.resize(size, dimension);
  phi3D_ = 0.;
  DoubleTab coord_cg(size, dimension);
  faces.calculer_centres_gravite(coord_cg);
  double h = abscissa_[1] -abscissa_[0]; //1d mesh pitch
  int abscissa_size = abscissa_.size();
  for(int face= 0; face<size; face++)
    {
      if(abs(phi3D_(face, direction_)) <1.e-16) //this node has not yet been treated
        {
          double s = coord_cg(face,direction_);
          double xs= coord_cg(face,0);
          double ys= coord_cg(face,1);
          double zs= coord_cg(face,2);
          if (direction_==0)
            xs=0.;
          else if (direction_==1)
            ys=0.;
          else
            zs=0.;
          //double test = 0.22*sin(PI*s/0.7);
          int i, j ;
          i = s/h;
          if((i+1) < abscissa_size)
            {
              j= i+1;
            }
          else
            {
              j=i;
            }
          double ux, uy, uz, Rx, Ry, Rz;

          //linear interpolation between points i and j
          double	alpha = (abscissa_[j] - s)/h;
          double	betha = (s - abscissa_[i])/h;

          ux=alpha*u_(i, 0) + betha*u_(j, 0);
          uy=alpha*u_(i, 1) + betha*u_(j, 1);
          uz=alpha*u_(i, 2) + betha*u_(j, 2);
          Rx=alpha*R_(i, 0) + betha*R_(j, 0);
          Ry=alpha*R_(i, 1) + betha*R_(j, 1);
          Rz=alpha*R_(i, 2) + betha*R_(j, 2);

          //Cerr<<" test apres interpolation = "<<uy - test<<finl;
          phi3D_(face, 0) =ux + Ry*zs -Rz*ys;
          phi3D_(face, 1) =uy + Rz*xs -Rx*zs;
          phi3D_(face, 2) =uz + Rx*ys -Ry*xs;

        }

    }
}*/
