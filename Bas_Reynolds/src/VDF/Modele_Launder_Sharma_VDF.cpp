/****************************************************************************
* Copyright (c) 2015, CEA
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
// Directory:   $TRUST_ROOT/src/VDF/Turbulence
// Version:     /main/10
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_Launder_Sharma_VDF.h>
#include <Zone_VDF.h>
#include <Champ_Uniforme.h>

Implemente_instanciable(Modele_Launder_Sharma_VDF,"Modele_Launder_Sharma_VDF",Modele_Jones_Launder_VDF);
// XD Launder_Sharmma modele_fonction_bas_reynolds_base Launder_Sharma -1 not_set



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

void  Modele_Launder_Sharma_VDF::associer(const Zone_dis& zone_dis,
                                          const Zone_Cl_dis& zone_Cl_dis)
{
  //  const Zone_VDF& la_zone = ref_cast(Zone_VDF,zone_dis.valeur());
  //  const Zone_Cl_VDF& la_zone_Cl = ref_cast(Zone_Cl_VDF,zone_Cl_dis.valeur());
}

DoubleTab&  Modele_Launder_Sharma_VDF::Calcul_Fmu( DoubleTab& Fmu,const Zone_dis& zone_dis,const Zone_Cl_dis& zone_Cl_dis,const DoubleTab& K_eps_Bas_Re,const Champ_Don& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco.valeur());
  if (is_visco_const)
    visco=tab_visco(0,0);
  // Cerr << " dans Jones Sharma Calc Fmu " << finl;
  const Zone_VDF& la_zone = ref_cast(Zone_VDF,zone_dis.valeur());
  Fmu = 0;
  int nb_elem = la_zone.nb_elem();
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
  Cerr<<Fmu.mp_min_vect()<<" Fmu "<<Fmu.mp_max_vect()<<finl;


  return Fmu;
}
/*
  DoubleTab&  Modele_Launder_Sharma_VDF::Calcul_Fmu( DoubleTab& Fmu,const Zone_dis& zone_dis,const DoubleTab& K_eps_Bas_Re,const DoubleTab& tab_visco ) const
  {
  // Cerr << " dans Jones Sharma Calc Fmu " << finl;
  const Zone_VDF& la_zone = ref_cast(Zone_VDF,zone_dis.valeur());
  Fmu = 0;
  int nb_elem = la_zone.nb_elem();
  int elem;
  DoubleTab Re(nb_elem);
  for (elem=0; elem< nb_elem ; elem++) {
  Re(elem) = (K_eps_Bas_Re(elem,0)*K_eps_Bas_Re(elem,0))/(tab_visco(elem)*K_eps_Bas_Re(elem,1));
  Fmu[elem] = exp(-3.4/((1.+Re(elem)/50.)*(1.+Re(elem)/50.)));

  }

  return Fmu;
  }
*/
