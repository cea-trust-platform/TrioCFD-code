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
// File:        Restriction_elem_elem.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Operateurs
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Restriction_elem_elem_included
#define Restriction_elem_elem_included


#include <Domaine_VF.h>
#include <Probleme_base.h>
#include <Restriction_base.h>


/*! @brief class Restriction_elem_elem
 *
 */

//////////////////////////////////////////////////////////////////////////////////
//
// CLASS: Restriction_elem_elem
//
/////////////////////////////////////////////////////////////////////////////////

class Restriction_elem_elem : public Restriction_base
{
  Declare_instanciable(Restriction_elem_elem);

public:

  inline IntVect& nb_elemF();
  void restreindre(const Domaine_VF& domaine_VFG,const  Domaine_VF& domaine_VFF,const  IntVect& connect,
                   DoubleTab& incoG,
                   const DoubleTab& incoF,int nbcomp) override;
  inline void restreindre(const Domaine_VF& domaine_VFG, const Domaine_VF& domaine_VFF,
                          const  IntVect& connect,
                          DoubleTab& incoG,
                          const DoubleTab& incoF, int nb_comp,
                          int num_prem_face_frontG) override;
  void calculer(const Domaine_VF& domaine_VFG, const Domaine_VF& domaine_VFF,const  IntVect& connect) override;


private:
  IntVect nb_elemF_;
  /*   Pb_fin probleme_fin; */
  /*   Connectivites connectivites_ff_ee; */
  //Prolongement operateur_P;
  //Restriction operateur_R;

};



inline IntVect& Restriction_elem_elem::nb_elemF()
{
  return nb_elemF_;
}



inline void Restriction_elem_elem::restreindre(const Domaine_VF& domaine_VFG,
                                               const Domaine_VF& domaine_VFF,
                                               const  IntVect& connect,
                                               DoubleTab& incoG,
                                               const DoubleTab& incoF,
                                               int nb_comp,
                                               int num_prem_face_frontG)
{
  Cerr<<"n'est pas codee dans cette classe"<<finl;
  exit();
}


#endif
