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
#include <TRUSTVects.h>
#include <TRUSTTabs.h>
#include <Probleme_base.h>
#include <Domaine_VEF.h>
#include <Domaine_Cl_VEF.h>
#include <Frontiere.h>
#include <MD_Vector_tools.h>
#include <MD_Vector_std.h>
#include <Faces.h>
#include <fstream>
#include <iostream>
#include <math.h>       /* sin */
#include <TRUST_Ref.h>
#include <Domaine.h>

#define PI 3.14159265

using namespace std;


Implemente_instanciable_sans_constructeur_ni_destructeur(Beam_model, "Beam_model", Interprete_geometrique_base ) ;

// Syntaxe:
//  Beam_model NOMDOMAINE {
//    nb_beam
//	  Name
//     nb_modes number of modes
//     direction 0|1|2
//     NewmarkTimeScheme MA|FD|HHT
//     Mass_and_stiffness_file_name
//     Absc_file_name
//     Modal_deformation_file_name nb_modes files
//     [ Young_Module ]
//     { Rho_beam ]
//     [ BaseCenterCoordinates position  ]
//     [ CI_file_name ]
//     [ Output_position_1D nb_points  position ]
//     [ Output_position_3D nb_points  position ]
//     [ Restart_file_name file ]
//
//  }
//XD  Beam_model interprete Beam_model 1 Reduced mechanical model: a beam model. Resolution based on a modal analysis. Temporal discretization: Newmark or Hilber-Hughes-Taylor (HHT)
//XD attr nb_beam entier nb 0 Number of beams
//XD attr Name chaine BeamName 0 Name of the beam (the name must match with the name of the edge in the fluid domain)
//XD attr nb_modes entier n 0 Number of modes
//XD attr direction entier dir 0 x=0, y=1, z=2
//XD attr NewmarkTimeScheme chaine NewmarkTimeScheme 0 Solve the beam dynamics. Time integration scheme: choice between MA (Newmark mean acceleration),  FD (Newmark finite differences), and HHT alpha (Hilber-Hughes-Taylor, alpha usually -0.1 )
//XD attr Mass_and_stiffness_file_name chaine  Mass_and_stiffness_file_name 0 Name of the file containing the diagonal modal mass, stiffness, and damping matrices.
//XD attr Absc_file_name chaine Absc_file_name 0 Name of the file containing the coordinates of the Beam
//XD attr  Modal_deformation_file_name chaine  Modal_deformation_file_name 0 Name of the file containing the modal deformation of the Beam (mandatory if different from 0. 0. 0.)
//XD attr Young_Module floattant young 1 Young Module
//XD attr Rho_beam floattant rho 1 Beam density
//XD attr BaseCenterCoordinates floattant pos_center 1 position of the base center coordinates on the Beam
//XD attr CI_file_name chaine CI_file_name 1 Name of the file containing the initial condition of the Beam
//XD attr Restart_file_name chaine Restart_file_name 1 SaveBeamForRestart.txt file to restart the calculation
//XD attr Output_position_1D floattant pt1d 1 nb_points  position Post-traitement of specific points on the Beam
//XD attr Output_position_3D floattant pt3d 1 nb_points  position Post-traitement of specific points on the 3d FSI boundary
Beam_model::Beam_model()
{

  nbModes_=0;;
  direction_=0;
  young_=200.e+9;
  rho_ = 8100.;
  timeScheme_="MA";
  temps_ =0.;
  output_position_1D_.resize(0);
  output_position_1D_=0.;
  output_position_3D_.resize(0);
  output_position_3D_=0.;
  alpha_=0.;
  fluidForceOnBeam_.resize(0);
  fluidForceOnBeam_=0.;
  tempsComputeForceOnBeam_=0.;
  x0_=0.;
  y0_=0.;
  z0_=0.;
}
Beam_model::~Beam_model()
{

}

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
              Cerr<<" Beam_model::read_input_files values of mass and stiffness are wrong in the file "<<masse_and_stiffness_file_name<<finl;
              Cerr<<" Value of mass = "<< mass<< " and value of stiffness = "<<stiffness<<finl;
              exit();
            }
        }

      monFlux.close();
    }
  else
    {
      Cerr<< "ERROR: Unable to open the file." <<masse_and_stiffness_file_name<<finl;
      exit();
    }
  //Cerr<<"mass = "<<mass_<<" stiffnes "<<stiffness_<<" damp = "<<damping_<<finl;
}
void Beam_model::readInputAbscFiles(Nom& absc_file_name)
{
  Cerr << "Read beam coordinates from "<< absc_file_name <<" files " << finl;

  string const nomFichier(absc_file_name);
  ifstream monFlux(nomFichier.c_str());
  if(monFlux)
    {
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
      exit();
    }
}

void Beam_model::readInputCIFile(Nom& CI_file_name)
{

  Cerr << "Read initial condition beam from "<< CI_file_name <<" file " << finl;

  qSpeed_.resize(nbModes_);
  qAcceleration_.resize(nbModes_);
  qDisplacement_.resize(nbModes_);
  qHalfSpeed_.resize(nbModes_);
  fluidForceOnBeam_.resize(nbModes_);
  qSpeed_=0.;
  qHalfSpeed_=0.;
  qAcceleration_=0.;
  qDisplacement_=0.;
  fluidForceOnBeam_=0.;

  string const nomFichier(CI_file_name);

  ifstream monFlux(nomFichier.c_str());

  if(monFlux)
    {
      double displacement;
      for(int i=0; i<nbModes_; i++)
        {
          monFlux >> displacement;
          qDisplacement_[i]=displacement;
        }

      monFlux.close();
    }
  else
    {
      Cerr<< "ERROR: Unable to open the file." <<CI_file_name<<finl;
      exit();
    }

}

void Beam_model::readRestartFile(Nom& Restart_file_name)
{

  Cerr << "Read restart "<< Restart_file_name <<" file " << finl;
  string const nomFichier(Restart_file_name);

  ifstream monFlux(nomFichier.c_str());

  if(monFlux)
    {
      double temps, displacement, speed, acceleration;
      for(int i=0; i<nbModes_; i++)
        {
          monFlux >> temps >> displacement >> speed >>acceleration;
          temps_= temps;
          qDisplacement_[i]=displacement;
          qSpeed_[i] = speed;
          qAcceleration_[i] = acceleration;
        }

      monFlux.close();
      tempsComputeForceOnBeam_=temps_;
    }
  else
    {
      Cerr<< "ERROR: Unable to open the file." <<Restart_file_name<<finl;
    }
}


void Beam_model::readInputModalDeformation(Noms& modal_deformation_file_name)
{

  Cerr << "Read beam modal deformation coefficients from "<< modal_deformation_file_name <<" files " << finl;

  for(int count=0; count<nbModes_; count++)
    {
      string const nomFichier(modal_deformation_file_name[count]);

      ifstream monFlux(nomFichier.c_str());
      int size = abscissa_.size();
      DoubleTab u(size, 3);
      DoubleTab R(size, 3);

      if(monFlux)
        {
          double  ux, uy, uz, rx, ry, rz;
          for(int i=0; i<size; i++)
            {
              monFlux >> ux>> uy>> uz>> rx>> ry>> rz;
              u(i, 0)=ux;
              u(i, 1)=uy;
              u(i, 2)=uz;
              R(i, 0)=rx;
              R(i, 1)=ry;
              R(i, 2)=rz;

            }
          monFlux.close();
          u_.add(u);
          R_.add(R);
        }

      else
        {
          Cerr<< "ERROR: Unable to open the file." <<modal_deformation_file_name[count]<<finl;
          exit();
        }
    }
}

void Beam_model::initialization(double displacement)
{
  qSpeed_.resize(nbModes_);
  qAcceleration_.resize(nbModes_);
  qDisplacement_.resize(nbModes_);
  qHalfSpeed_.resize(nbModes_);
  fluidForceOnBeam_.resize(nbModes_);
  qSpeed_=0.;
  qHalfSpeed_=0.;
  qAcceleration_=0.;
  qDisplacement_=displacement;
  fluidForceOnBeam_=0.;

}
void Beam_model::initialization()
{
  qSpeed_.resize(nbModes_);
  qAcceleration_.resize(nbModes_);
  qDisplacement_.resize(nbModes_);
  qHalfSpeed_.resize(nbModes_);
  fluidForceOnBeam_.resize(nbModes_);
  qSpeed_=0.;
  qHalfSpeed_=0.;
  qAcceleration_=0.;
  qDisplacement_=0.;
  fluidForceOnBeam_=0.;

}
//Solve the beam dynamics. Time integration scheme: Newmark finite differences
DoubleVect& Beam_model::NewmarkSchemeFD (const double& dt)
{
  double halfDt=dt/2.;
  for(int j=0; j < nbModes_; j++)
    {
      qHalfSpeed_[j] = qSpeed_[j] + halfDt*qAcceleration_[j];
      qDisplacement_[j] += dt*qHalfSpeed_[j];
      double coeff1 = mass_[j] + halfDt*damping_[j];
      double coeff2 =	damping_[j]*qHalfSpeed_[j] + stiffness_[j]*qDisplacement_[j];
      qAcceleration_[j]= (fluidForceOnBeam_[j] - coeff2)/coeff1;
      qSpeed_[j] = qHalfSpeed_[j] + halfDt*qAcceleration_[j];
      //qHalfSpeed_[j] = qSpeed_[j] + halfDt*qAcceleration_[j];
    }


  saveBeamForRestart();
  if(output_position_1D_.size()>0) printOutputBeam1D();
  if(output_position_3D_.size()>0) printOutputBeam3D();

  return qSpeed_;
}
//Solve the beam dynamics. Time integration scheme: Newmark mean acceleration
DoubleVect& Beam_model::NewmarkSchemeMA (const double& dt)
{
  double halfDt=dt/2;
  double squareHalfDt= halfDt*halfDt;
  for(int j=0; j < nbModes_; j++)
    {
      double PreviousqAcceleration= qAcceleration_[j];
      double coeff1 = mass_[j] + halfDt*damping_[j] + squareHalfDt*stiffness_[j];
      double coeff2 = damping_[j]*(qSpeed_[j] + halfDt*qAcceleration_[j]) + stiffness_[j]*(qDisplacement_[j] + dt*qSpeed_[j] + squareHalfDt*qAcceleration_[j]);
      qAcceleration_[j]=(fluidForceOnBeam_[j] - coeff2)/coeff1;
      qDisplacement_[j] += dt*qSpeed_[j] + squareHalfDt*(PreviousqAcceleration + qAcceleration_[j]);
      qSpeed_[j] += halfDt*(PreviousqAcceleration + qAcceleration_[j]);
    }

  saveBeamForRestart();
  if(output_position_1D_.size()>0) printOutputBeam1D();
  if(output_position_3D_.size()>0) printOutputBeam3D();

  return qSpeed_;
}

//Solve the beam dynamics. Time integration scheme: HHT
DoubleVect& Beam_model::TimeSchemeHHT (const double& dt)
{
  //Cerr<<" dt = "<<dt<<" temps_ = "<<temps_<<" tempsComputeForceOnBeam_= "<<tempsComputeForceOnBeam_<<finl;
  //Cerr<<" force = "<<fluidForceOnBeam_<<finl;
  double beta = (1-alpha_)*(1-alpha_)/4;
  double gamma= (1- 2*alpha_)/2;
  double squareDt=dt*dt;

  for(int j=0; j < nbModes_; j++)
    {
      double PreviousqAcceleration= qAcceleration_[j];
      double coeff1 = mass_[j] + dt*(1.+ alpha_)*gamma*damping_[j] + squareDt*(1.+ alpha_)*beta*stiffness_[j];
      double coeff2 = (1.+ alpha_)*damping_[j]*(qSpeed_[j] + dt*(1.-gamma)*PreviousqAcceleration);
      double coeff3 = (1.+ alpha_)*stiffness_[j]*(qDisplacement_[j] + dt*qSpeed_[j] + squareDt*(0.5-beta)*PreviousqAcceleration);
      double coeff4 = alpha_*stiffness_[j]*qDisplacement_[j];
      double coeff5 = alpha_*damping_[j]*qSpeed_[j];
      qAcceleration_[j]=(fluidForceOnBeam_[j] - coeff2 - coeff3 + coeff4 + coeff5)/coeff1;
      qDisplacement_[j] += dt*qSpeed_[j] + squareDt*((0.5-beta)*PreviousqAcceleration + beta*qAcceleration_[j]);
      qSpeed_[j] += dt*((1.-gamma)*PreviousqAcceleration + gamma*qAcceleration_[j]);
    }

  saveBeamForRestart();
  if(output_position_1D_.size()>0) printOutputBeam1D();
  if(output_position_3D_.size()>0) printOutputBeam3D();

  return qSpeed_;
}


DoubleVect& Beam_model::getVelocity(const double& tps, const double& dt)
{
  Cerr<<" Beam name : "<<beamName_<<", Time discretization scheme : "<<timeScheme_<<", Nb Modes : "<<nbModes_<<finl;

  if(dt == 0.)
    {
      return qSpeed_;
    }
  else if(temps_!=tps) // update qSpeed_ only once per time step!
    {
      temps_=tps;
      if(timeScheme_=="HHT")
        return TimeSchemeHHT(dt);
      else if(timeScheme_=="FD")
        return NewmarkSchemeFD(dt);
      else
        return NewmarkSchemeMA(dt);
    }
  else
    {
      return qSpeed_;
    }
}

DoubleVect Beam_model::interpolationOnThe3DSurface(const double& x, const double& y, const double& z, const DoubleTab& u, const DoubleTab& R) const
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
      s = xs;
      xs=0.;
    }
  else if (direction_== 1)
    {
      s = ys;
      ys=0.;
    }

  else
    {
      s = zs;
      zs=0.;
    }

  int i, j ;
  i = int(s/h);
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
  else if(abs(abscissa_[i] - s)< 1.e-4)
    {
      alpha=1.;
      betha=0.;
    }
  else if (abs(abscissa_[j] - s)< 1.e-4)
    {
      alpha=0.;
      betha=1.;
    }
  else
    {
      alpha = (abscissa_[j] - s)/h;
      betha = (s - abscissa_[i])/h;
      if(alpha <0.)
        {
          alpha=0.;
          betha=1.;
        }
      else if (betha < 0.)
        {
          alpha=1.;
          betha=0.;
        }

    }

  ux=alpha*u(i, 0) + betha*u(j, 0);
  uy=alpha*u(i, 1) + betha*u(j, 1);
  uz=alpha*u(i, 2) + betha*u(j, 2);
  Rx=alpha*R(i, 0) + betha*R(j, 0);
  Ry=alpha*R(i, 1) + betha*R(j, 1);
  Rz=alpha*R(i, 2) + betha*R(j, 2);

  phi[0] =ux + Ry*(zs-z0_) -Rz*(ys-y0_);
  phi[1] =uy + Rz*(xs-x0_) -Rx*(zs-z0_);
  phi[2] =uz + Rx*(ys-y0_) -Ry*(xs-x0_);



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


void Beam_model::saveBeamForRestart() const
{

  if (je_suis_maitre())
    {
      std::ofstream ofs_sauve;
      ofs_sauve.open (beamName_+"SaveBeamForRestart.txt", std::ofstream::out | std::ofstream::trunc);
      ofs_sauve.precision(32);
      for(int j=0; j < nbModes_; j++)
        {
          ofs_sauve<<temps_<<"  "<<qDisplacement_[j]<<" "<<qSpeed_[j]<<" "<<qAcceleration_[j]<<" "<<endl;
        }
      ofs_sauve.close();
    }

}


void Beam_model::printOutputBeam1D(bool first_writing) const
{

  if (je_suis_maitre())
    {
      int nb_output_points= output_position_1D_.size();
      DoubleTab displacement(nb_output_points,3);
      DoubleTab velocity(nb_output_points,3);
      DoubleTab acceleration(nb_output_points,3);
      displacement=0.;
      velocity=0.;
      acceleration=0.;
      for(int j=0; j < nbModes_; j++)
        {
          const DoubleTab& u=u_(j);
          for(int k=0; k<nb_output_points; k++)
            {
              for(int i=0; i<3; i++)
                {
                  displacement(k, i) += qDisplacement_[j]*u(int(output_position_1D_[k]),i);
                  velocity(k, i) += qSpeed_[j]*u(int(output_position_1D_[k]),i);
                  acceleration(k, i) += qAcceleration_[j]*u(int(output_position_1D_[k]),i);
                }
            }
        }
      Nom filename_disp(beamName_);
      filename_disp+="_Displacement1D.out";
      Nom filename_speed(beamName_);
      filename_speed+="_Velocity1D.out";
      Nom filename_acc(beamName_);
      filename_acc+="_Acceleration1D.out";
      if (!displacement_out_1d_.is_open())
        {
          displacement_out_1d_.ouvrir(filename_disp, (first_writing?ios::out:ios::app));
          displacement_out_1d_.setf(ios::scientific);
        }
      if (!speed_out_1d_.is_open())
        {
          speed_out_1d_.ouvrir(filename_speed, (first_writing?ios::out:ios::app));
          speed_out_1d_.setf(ios::scientific);
        }
      if (!acceleration_out_1d_.is_open())
        {
          acceleration_out_1d_.ouvrir(filename_acc, (first_writing?ios::out:ios::app));
          acceleration_out_1d_.setf(ios::scientific);
        }
      // comments are added to the file header
      if (first_writing)
        {
          displacement_out_1d_<< "# Printing Beam 1D displacement: time  values of x y z -component at points ";
          speed_out_1d_<< "# Printing Beam 1D velocity: time  values of x y z -component at points ";
          acceleration_out_1d_<< "# Printing Beam 1D acceleration: time  values of x y z -component at points ";
          for(int k=0; k<nb_output_points; k++)
            {
              displacement_out_1d_<<output_position_1D_[k]<<" ";
              speed_out_1d_<<output_position_1D_[k]<<" ";
              acceleration_out_1d_<<output_position_1D_[k]<<" ";
            }
          displacement_out_1d_<<finl;
          speed_out_1d_<<finl;
          acceleration_out_1d_<<finl;
        }
      displacement_out_1d_<< temps_<< " ";
      speed_out_1d_<< temps_<< " ";
      acceleration_out_1d_<< temps_<< " ";
      for(int k=0; k<nb_output_points; k++)
        {
          for(int i=0; i<3; i++)
            {
              displacement_out_1d_<<displacement(k, i)<<" ";
              speed_out_1d_<<velocity(k, i)<<" ";
              acceleration_out_1d_<<acceleration(k, i)<<" ";
            }
        }
      displacement_out_1d_<<finl;
      speed_out_1d_<<finl;
      acceleration_out_1d_<<finl;
    }
}
void Beam_model::printOutputBeam3D(bool first_writing) const
{

  if (je_suis_maitre())
    {
      int nb_output_points= output_position_3D_.dimension(0);
      DoubleTab displacement(nb_output_points,3);
      DoubleTab velocity(nb_output_points,3);
      DoubleTab acceleration(nb_output_points,3);
      displacement=0.;
      velocity=0.;
      acceleration=0.;
      DoubleVect phi3D(3);
      for(int j=0; j < nbModes_; j++)
        {
          const DoubleTab& u=u_(j);
          const DoubleTab& R=R_(j);
          for(int k=0; k<nb_output_points; k++)
            {
              phi3D=interpolationOnThe3DSurface(output_position_3D_(k,0),output_position_3D_(k,1),output_position_3D_(k,2), u, R);
              for(int i=0; i<3; i++)
                {
                  displacement(k, i) += qDisplacement_[j]*phi3D[i];
                  velocity(k, i) += qSpeed_[j]*phi3D[i];
                  acceleration(k, i) += qAcceleration_[j]*phi3D[i];
                }
            }
        }

      Nom filename_disp(beamName_);
      filename_disp+="_Displacement3D.out";
      Nom filename_speed(beamName_);
      filename_speed+="_Velocity3D.out";
      Nom filename_acc(beamName_);
      filename_acc+="_Acceleration3D.out";
      if (!displacement_out_3d_.is_open())
        {
          displacement_out_3d_.ouvrir(filename_disp, (first_writing?ios::out:ios::app));
          displacement_out_3d_.setf(ios::scientific);
        }
      if (!speed_out_3d_.is_open())
        {
          speed_out_3d_.ouvrir(filename_speed, (first_writing?ios::out:ios::app));
          speed_out_3d_.setf(ios::scientific);
        }
      if (!acceleration_out_3d_.is_open())
        {
          acceleration_out_3d_.ouvrir(filename_acc, (first_writing?ios::out:ios::app));
          acceleration_out_3d_.setf(ios::scientific);
        }
      // comments are added to the file header
      if (first_writing)
        {
          displacement_out_3d_<< "# Printing Beam 3D displacement: time  values of x y z -component at points ";
          speed_out_3d_<< "# Printing Beam 3D velocity: time  values of x y z -component at points ";
          acceleration_out_3d_<< "# Printing Beam 3D acceleration: time  values of x y z -component at points ";
          for(int k=0; k<nb_output_points; k++)
            {
              displacement_out_3d_<<"( ";
              speed_out_3d_<<"( ";
              acceleration_out_3d_<<"( ";
              for(int i=0; i<3; i++)

                {
                  displacement_out_3d_<<output_position_3D_(k, i)<<" ";
                  speed_out_3d_<<output_position_3D_(k, i)<<" ";
                  acceleration_out_3d_<<output_position_3D_(k, i)<<" ";
                }
              displacement_out_3d_<<")";
              speed_out_3d_<<")";
              acceleration_out_3d_<<")";
            }
          displacement_out_3d_<<finl;
          speed_out_3d_<<finl;
          acceleration_out_3d_<<finl;
        }
      displacement_out_3d_<< temps_<< " ";
      speed_out_3d_<< temps_<< " ";
      acceleration_out_3d_<< temps_<< " ";
      for(int k=0; k<nb_output_points; k++)
        {
          for(int i=0; i<3; i++)
            {
              displacement_out_3d_<<displacement(k, i)<<" ";
              speed_out_3d_<<velocity(k, i)<<" ";
              acceleration_out_3d_<<acceleration(k, i)<<" ";
            }
        }
      displacement_out_3d_<<finl;
      speed_out_3d_<<finl;
      acceleration_out_3d_<<finl;
    }
}

void Beam_model::printOutputFluidForceOnBeam(bool first_writing) const
{

}

void Beam_model::setCenterCoordinates(const double& x0,const double& y0, const double& z0)
{

  x0_=x0;
  y0_=y0;
  z0_=z0;
}

/*void Beam_model::printOutputPosition1D() const
{

  if (je_suis_maitre())
    {
      int nb_output_points= output_position_1D_.size();
      DoubleTab displacement(nb_output_points,3);
      DoubleTab velocity(nb_output_points,3);
      DoubleTab acceleration(nb_output_points,3);
      displacement=0.;
      velocity=0.;
      acceleration=0.;
      for(int j=0; j < nbModes_; j++)
        {
          const DoubleTab& u=u_(j);
          for(int k=0; k<nb_output_points; k++)
            {
              for(int i=0; i<3; i++)
                {
                  displacement(k, i) += qDisplacement_[j]*u(int(output_position_1D_[k]),i);
                  velocity(k, i) += qSpeed_[j]*u(int(output_position_1D_[k]),i);
                  acceleration(k, i) += qAcceleration_[j]*u(int(output_position_1D_[k]),i);
                }
            }
        }
      std::ofstream ofs_1;
      ofs_1.open (beamName_+"Displacement1D.txt", std::ofstream::out | std::ofstream::app);
      std::ofstream ofs_2;
      ofs_2.open (beamName_+"Velocity1D.txt", std::ofstream::out | std::ofstream::app);
      std::ofstream ofs_3;
      ofs_3.open (beamName_+"Acceleration1D.txt", std::ofstream::out | std::ofstream::app);

      ofs_1<<temps_<<" ";
      for(int k=0; k<nb_output_points; k++)
        {
          for(int i=0; i<3; i++)
            {
              ofs_1<< displacement(k, i)<<" ";
            }
        }
      ofs_1<<endl;
      ofs_1.close();

      ofs_2<<temps_<<" ";
      for(int k=0; k<nb_output_points; k++)
        {
          for(int i=0; i<3; i++)
            {
              ofs_2<< velocity(k, i)<<" ";
            }
        }
      ofs_2<<endl;
      ofs_2.close();

      ofs_3<<temps_<<" ";
      for(int k=0; k<nb_output_points; k++)
        {
          for(int i=0; i<3; i++)
            {
              ofs_3<< acceleration(k, i)<<" ";
            }
        }
      ofs_3<<endl;
      ofs_3.close();
    }

}*/
