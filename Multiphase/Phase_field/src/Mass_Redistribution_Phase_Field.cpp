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

/*
 * Mass_Redistribution_Phase_Field.cpp
 *
 *  Created on: May 6, 2022
 *      Author: Shambhavi Nandan
 */

#include <Mass_Redistribution_Phase_Field.h>

double Mass_Redistribution_Phase_Field::epsilon_mass_redistribute = 1.0e-7;
DoubleTab Mass_Redistribution_Phase_Field::c_ini;

/* double Mass_Redistribution_Phase_Field::minn = 0.01;
double Mass_Redistribution_Phase_Field::maxx = 0.99;
DoubleTab Mass_Redistribution_Phase_Field::c_ini;
double topvalue = 1;
double minvalue = 0;
 */

/*
 * minvalue: the minimal value that concentration can take
 * topvalue: the top value that concentration can take
 * minnMassredistribution: the minn reference for the massredistribution class
 * maxxMassredistribution: the maxx reference for the massredistribution class
 * maxxMassredistribution-minnMassredistribution should be equal to topvalue-minvalue
 * */


//void Mass_Redistribution_Phase_Field::impose_mass_redistribution(const Zone_VDF& zvdf, DoubleTab& c_present, const DoubleVect minvalue, const DoubleVect topvalue, const DoubleVect minnMassreditribution, const DoubleVect maxxMassredistribution)
void Mass_Redistribution_Phase_Field::impose_mass_redistribution(const Zone_VDF& zvdf, DoubleTab& c_present, DoubleVect minX, DoubleVect maxX)
{
  DoubleTab createMDvect(c_ini.dimension(0));
  const MD_Vector& md = zvdf.zone().md_vector_elements();
  MD_Vector_tools::creer_tableau_distribue(md,createMDvect);
  createMDvect = 0.;

  DoubleTab c_ini_comp = createMDvect;
  DoubleTab c_present_comp = createMDvect;

  for(int nbcomp=0; nbcomp<c_present.line_size(); nbcomp++)
    {
      for(int iElem=0; iElem<c_present.dimension(0); iElem++)
        {
          if(c_present(iElem, nbcomp)>maxX(nbcomp)+epsilon_mass_redistribute)
            c_present(iElem, nbcomp) = maxX(nbcomp);
          else if(c_present(iElem, nbcomp)<minX(nbcomp)+epsilon_mass_redistribute)
            c_present(iElem, nbcomp) = minX(nbcomp);
        }

      int transitionValsElems = 0;
      for(int iElem=0; iElem<c_present.dimension(0); iElem++)
        {
          if(c_present(iElem, nbcomp)<maxX(nbcomp)-0.01 && c_present(iElem, nbcomp)>minX(nbcomp)+0.01)
            transitionValsElems++;
        }
      const int transitionValsElems_parallel = Process::mp_sum(transitionValsElems);

      for(int iElem=0; iElem<c_present.dimension(0); iElem++)
        {
          c_ini_comp(iElem) = c_ini(iElem, nbcomp);
          c_present_comp(iElem) = c_present(iElem, nbcomp);
        }

      const double initial_mass = mp_somme_vect(c_ini_comp);
      double present_mass = mp_somme_vect(c_present_comp);
      const double delta_mass = initial_mass-present_mass;
      const double redistributed_mass = (transitionValsElems_parallel>0)?(delta_mass/transitionValsElems_parallel):0.0;

      for(int iElem=0; iElem<c_present.dimension(0); iElem++)
        {
          if(c_present(iElem, nbcomp)<maxX(nbcomp)-0.01 && c_present(iElem, nbcomp)>minX(nbcomp)+0.01)
            {
              c_present(iElem, nbcomp)+=redistributed_mass;
              if(c_present(iElem, nbcomp)<minX(nbcomp) || c_present(iElem, nbcomp)>maxX(nbcomp))
                {
                  Cerr<<"********************** Unbounded values of concentration in Phase Field even after applying Mass Redistribution strategy **********************"<<finl;
                  Cerr<<"***************** This is due to the complete un-physical simulation resulting in the order of Redistributed Mass >"<<"1.0e-3 *****************"<<finl;
                  Cerr<<"***** Try adjusting the minn and maxx values to the order >"<<"1.0e-3"<<"in Mass_Redistribution_Phase_Field.cpp to continue simulation ********"<<finl;
                  Cerr<<"***************************** However doing that will later cause **** MASS NOT CONSERVED IN PHASE FIELD **** error ***************************"<<finl;
                  abort();
                }
            }
        }

      for(int iElem=0; iElem<c_present.dimension(0); iElem++)
        c_present_comp(iElem) = c_present(iElem, nbcomp);

      present_mass = mp_somme_vect(c_present_comp);

      if(fabs(present_mass-initial_mass)>epsilon_mass_redistribute)
        {
          Cerr<<"**********************"<<"MASS NOT CONSERVED IN PHASE FIELD"<<"***************************"<<finl;
          abort();
        }
    }
  c_present.echange_espace_virtuel();
}
