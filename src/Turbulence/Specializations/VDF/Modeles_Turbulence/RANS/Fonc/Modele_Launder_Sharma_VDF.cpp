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
//////////////////////////////////////////////////////////////////////////////
//
// File:        Modele_Launder_Sharma_VDF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Modeles_Turbulence/RANS/Fonc
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_Launder_Sharma_VDF.h>
#include <Domaine_VDF.h>
#include <Champ_Uniforme.h>

Implemente_instanciable(Modele_Launder_Sharma_VDF,"Modele_Launder_Sharma_VDF",Modele_Jones_Launder_VDF);
// XD Launder_Sharmma modele_fonction_bas_reynolds_base Launder_Sharma -1 Model described in ' Launder, B. E. and Sharma, B. I. (1974), Application of the Energy-Dissipation Model of Turbulence to the Calculation of Flow Near a Spinning Disc, Letters in Heat and Mass Transfer, Vol. 1, No. 2, pp. 131-138.'



///////////////////////////////////////////////////////////////
//   Implementation des fonctions de la classe
///////////////////////////////////////////////////////////////
// printOn et readOn

Sortie& Modele_Launder_Sharma_VDF::printOn(Sortie& s ) const
{
  Modele_Jones_Launder_VDF::printOn(s);
  return s;
}

Entree& Modele_Launder_Sharma_VDF::readOn(Entree& is )
{
  Modele_Jones_Launder_VDF::readOn(is);
  return is;
}

Entree& Modele_Launder_Sharma_VDF::lire(const Motcle& , Entree& is)
{
  return is;
}

void  Modele_Launder_Sharma_VDF::associer(const Domaine_dis_base& domaine_dis,
                                          const Domaine_Cl_dis_base& domaine_Cl_dis)
{
  //  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis);
  //  const Domaine_Cl_VDF& le_dom_Cl = ref_cast(Domaine_Cl_VDF,domaine_Cl_dis);
}

DoubleTab&  Modele_Launder_Sharma_VDF::Calcul_Fmu( DoubleTab& Fmu,const Domaine_dis_base& domaine_dis,const Domaine_Cl_dis_base& domaine_Cl_dis,const DoubleTab& K_eps_Bas_Re,const Champ_Don& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco->valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco.valeur());
  if (is_visco_const)
    visco=tab_visco(0,0);
  // Cerr << " dans Jones Sharma Calc Fmu " << finl;
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis);
  Fmu = 0;
  int nb_elem = le_dom.nb_elem();
  double Re;
  int elem;
  for (elem=0; elem< nb_elem ; elem++)
    {
      if (!is_visco_const)
        visco=tab_visco[elem];
      if (visco>0)
        {
          Re = K_eps_Bas_Re(elem,0)*K_eps_Bas_Re(elem,0)/(K_eps_Bas_Re(elem,1)+DMINFLOAT)/visco;
          //Fmu[elem] = exp(-3.4/((1.+Re/50.)*(1.+Re/50.)*(1+Re/50.)));
          Fmu[elem] = exp(-3.4/((1.+Re/50.)*(1.+Re/50.)));
        }
      else
        Fmu[elem]=1.;
    }
  //Fmu=1;
//  Cerr<<Fmu.mp_min_vect()<<" Fmu "<<Fmu.mp_max_vect()<<finl;


  return Fmu;
}
/*
  DoubleTab&  Modele_Launder_Sharma_VDF::Calcul_Fmu( DoubleTab& Fmu,const Domaine_dis_base& domaine_dis,const DoubleTab& K_eps_Bas_Re,const DoubleTab& tab_visco ) const
  {
  // Cerr << " dans Jones Sharma Calc Fmu " << finl;
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis);
  Fmu = 0;
  int nb_elem = le_dom.nb_elem();
  int elem;
  DoubleTab Re(nb_elem);
  for (elem=0; elem< nb_elem ; elem++) {
  Re(elem) = (K_eps_Bas_Re(elem,0)*K_eps_Bas_Re(elem,0))/(tab_visco(elem)*K_eps_Bas_Re(elem,1));
  Fmu[elem] = exp(-3.4/((1.+Re(elem)/50.)*(1.+Re(elem)/50.)));

  }

  return Fmu;
  }
*/



DoubleTab&  Modele_Launder_Sharma_VDF::Calcul_Fmu_BiK( DoubleTab& Fmu,const Domaine_dis_base& domaine_dis,const Domaine_Cl_dis_base& domaine_Cl_dis,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re,const Champ_Don& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco->valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco.valeur());
  if (is_visco_const)
    visco=tab_visco(0,0);
  // Cerr << " dans Jones Sharma Calc Fmu " << finl;
  const Domaine_VDF& le_dom = ref_cast(Domaine_VDF,domaine_dis);
  Fmu = 0;
  int nb_elem = le_dom.nb_elem();
  double Re;
  int elem;
  for (elem=0; elem< nb_elem ; elem++)
    {
      if (!is_visco_const)
        visco=tab_visco[elem];
      if (visco>0)
        {
          Re = K_Bas_Re(elem)*K_Bas_Re(elem)/(eps_Bas_Re(elem)+DMINFLOAT)/visco;
          //Fmu[elem] = exp(-3.4/((1.+Re/50.)*(1.+Re/50.)*(1+Re/50.)));
          Fmu[elem] = exp(-3.4/((1.+Re/50.)*(1.+Re/50.)));
        }
      else
        Fmu[elem]=1.;
    }
  //Fmu=1;
  Cerr<<Fmu.mp_min_vect()<<" Fmu "<<Fmu.mp_max_vect()<<finl;


  return Fmu;
}

