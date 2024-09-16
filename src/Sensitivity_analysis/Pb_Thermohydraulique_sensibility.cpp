/****************************************************************************
* Copyright (c) 2020, CEA
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
// File      : Pb_Thermohydraulique_sensibility.cpp
// Directory : $Sensitivity_analysis/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Pb_Thermohydraulique_sensibility.h>
#include <Fluide_Ostwald.h>
#include <Verif_Cl.h>
#include <Probleme_base.h>

Implemente_instanciable( Pb_Thermohydraulique_sensibility, "Pb_Thermohydraulique_sensibility", Pb_Fluide_base ) ;
// XD Pb_Thermohydraulique_sensibility pb_thermohydraulique Pb_Thermohydraulique_sensibility -1 Resolution of Resolution of thermohydraulic sensitivity problem
// XD  attr fluide_incompressible fluide_incompressible fluide_incompressible 0 The fluid medium associated with the problem.
// XD  attr Convection_Diffusion_Temperature_Sensibility Convection_Diffusion_Temperature_sensibility convection_diffusion_temperature 0  Convection diffusion temperature sensitivity equation
// XD  attr Navier_Stokes_standard_sensibility Navier_Stokes_standard_sensibility Navier_Stokes_standard_sensibility 0   Navier Stokes sensitivity equation


/*! @brief Simple appel a: Pb_Fluide_base::printOn(Sortie&) Ecrit le probleme sur un flot de sortie.
 *
 * @param (Sortie& os) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Pb_Thermohydraulique_sensibility::printOn( Sortie& os ) const
{
  Pb_Fluide_base::printOn( os );
  return os;
}


/*! @brief Simple appel a: Pb_Fluide_base::readOn(Entree&) Lit le probleme a partir d'un flot d'entree.
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Pb_Thermohydraulique_sensibility::readOn( Entree& is )
{
  Pb_Fluide_base::readOn( is );
  return is;
}

/*! @brief Renvoie le nombre d'equation, Renvoie 2 car il y a 2 equations a un probleme de
 *
 *     thermo-hydraulique standard:
 *         l'equation de Navier Stokes
 *         l' equation de la thermique de type Convection_Diffusion_Temperature
 *
 * @return (int) le nombre d'equation
 */
int Pb_Thermohydraulique_sensibility::nombre_d_equations() const
{
  return 2;
}

/*! @brief Renvoie l'equation d'hydraulique de type Navier_Stokes_std si i=0 Renvoie l'equation de la thermique de type
 *
 *     Convection_Diffusion_Temperature si i=1
 *     (version const)
 *
 * @param (int i) l'index de l'equation a renvoyer
 * @return (Equation_base&) l'equation correspondante a l'index
 */
const Equation_base& Pb_Thermohydraulique_sensibility::equation(int i) const
{
  if ( !( i==0 || i==1 ) )
    {
      Cerr << "\nError in Pb_Thermohydraulique::equation() : Wrong number of equation !" << finl;
      Process::exit();
    }
  if (i == 0)
    return eq_hydraulique;
  else
    return eq_thermique;
}

/*! @brief Renvoie l'equation d'hydraulique de type Navier_Stokes_std si i=0 Renvoie l'equation de la thermique de type
 *
 *     Convection_Diffusion_Temperature si i=1
 *
 * @param (int i) l'index de l'equation a renvoyer
 * @return (Equation_base&) l'equation correspondante a l'index
 */
Equation_base& Pb_Thermohydraulique_sensibility::equation(int i)
{
  if ( !( i==0 || i==1 ) )
    {
      Cerr << "\nError in Pb_Thermohydraulique::equation() : Wrong number of equation !" << finl;
      Process::exit();
    }
  if (i == 0)
    return eq_hydraulique;
  else
    return eq_thermique;
}


/*! @brief Associe le milieu au probleme Le milieu doit etre de type fluide incompressible
 *
 * @param (Milieu_base& mil) le milieu physique a associer au probleme
 * @throws mauvais type de milieu physique
 */
void Pb_Thermohydraulique_sensibility::associer_milieu_base(const Milieu_base& mil)
{
  if (sub_type(Fluide_Incompressible,mil))
    {
      eq_hydraulique.associer_milieu_base(mil);
      eq_thermique.associer_milieu_base(mil);
    }
  else
    {
      Cerr << "Un milieu de type " << mil.que_suis_je() << " ne peut etre associe a "<< finl;
      Cerr << "un probleme de type Pb_Thermohydraulique_sensibility " << finl;
      exit();
    }
}



/*! @brief Teste la compatibilite des equations de la thermique et de l'hydraulique.
 *
 * Le test se fait sur les conditions
 *     aux limites discretisees de chaque equation.
 *     Appel la fonction de librairie hors classe:
 *       tester_compatibilite_hydr_thermique(const Domaine_Cl_dis_base&,const Domaine_Cl_dis_base&)
 *
 * @return (int) code de retour propage
 */
int Pb_Thermohydraulique_sensibility::verifier()
{
  const Domaine_Cl_dis_base& domaine_Cl_hydr = eq_hydraulique.domaine_Cl_dis();
  const Domaine_Cl_dis_base& domaine_Cl_th = eq_thermique.domaine_Cl_dis();
  return tester_compatibilite_hydr_thermique(domaine_Cl_hydr,domaine_Cl_th);
}

/*void Pb_Hydraulique_sensibility::finir()
{

  Probleme_base::finir();

  Navier_Stokes_std& eqnNS = ref_cast(Navier_Stokes_std,equation(0));
  const  DoubleTab& vitesse=  eqnNS.inconnue().valeur();
  if( equation(0).que_suis_je() == "Navier_Stokes_standard")
    {
      Nom   nom_fichier_vitesse("Velocity_state.txt");
      SFichier file_vitesse(nom_fichier_vitesse);
      Nom   nom_fichier_centre_faces("Center_gravity_faces.txt");
      SFichier file_c_faces(nom_fichier_centre_faces);
      const Zone_VEF& la_zone = ref_cast(Zone_VEF,eqnNS.zone_dis().valeur());
      const DoubleTab& xv= la_zone.xv();
      assert(xv.dimension(0)==vitesse.dimension(0));
      for(int face=0; face<vitesse.dimension(0); face++)
        {
          file_vitesse<<vitesse(face, 0)<<"  "<<vitesse(face, 1)<<finl;
          file_c_faces<<xv(face, 0)<<"  "<<xv(face, 1)<<finl;
        }
    }
  else  if( equation(0).que_suis_je() == "Navier_Stokes_standard_sensibility")
    {
      Nom   nom_fichier_vitesse("Velocity_sensibility.txt");
      SFichier file_vitesse(nom_fichier_vitesse);
      for(int face=0; face<vitesse.dimension(0); face++)
        {
          file_vitesse<<vitesse(face, 0)<<"  "<<vitesse(face, 1)<<finl;
        }
    }
  if(nombre_d_equations() == 2 )
    {
      Convection_Diffusion_Temperature& eqnCDT = ref_cast(Convection_Diffusion_Temperature,equation(1));
      const  DoubleTab& temperature=  eqnCDT.inconnue().valeur();
      if (equation(1).que_suis_je() == "Convection_Diffusion_Temperature")
        {
          Nom   nom_fichier_temp("Temperature_state.txt");
          SFichier file_temp(nom_fichier_temp);
          for(int face=0; face<temperature.size(); face++)
            {
              file_temp<<temperature(face)<<finl;
            }

        }
      else if (equation(1).que_suis_je() == "Convection_Diffusion_Temperature_sensibility")
        {
          Nom   nom_fichier_temp("Temperature_sensibility.txt");
          SFichier file_temp(nom_fichier_temp);
          for(int face=0; face<temperature.size(); face++)
            {
              file_temp<<temperature(face)<<finl;
            }
        }
    }

} */
