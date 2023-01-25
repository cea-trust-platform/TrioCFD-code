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
// File:        Restriction.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Kernel
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////





#ifndef Restriction_included
#define Restriction_included

#include <TRUST_Deriv.h>
#include <Restriction_base.h>


/*! @brief classe Restriction Classe generique de la hierarchie des algorithmes, un objet
 *
 *      Restriction peut referencer n'importe quel objet derivant de
 *      Restriction_base.
 *      La plupart des methodes appellent les methodes de l'objet
 *      Restriction_base sous-jacent via la methode valeur() declaree grace
 * alamacroDeclare_deriv().;
 *
 * @sa Restriction_base
 */
class Restriction : public DERIV(Restriction_base)
{
  Declare_instanciable(Restriction);

public :

  inline void restreindre(const Zone_VF& zone_VFG, const Zone_VF& zone_VFF, const IntVect& connect,
                          DoubleTab& incoG,
                          const DoubleTab& incoF,int nbcomp) ;
  inline void calculer(const Zone_VF& zone_VFG, const Zone_VF& zone_VFF, const IntVect& connect) ;
  inline void restreindre(const Zone_VF& zone_VFG, const Zone_VF& zone_VFF,
                          const  IntVect& connect,
                          DoubleTab& incoG,
                          const DoubleTab& incoF, int nb_comp,
                          int num_prem_face_frontG);
};



void Restriction::restreindre(const Zone_VF& zone_VFG, const Zone_VF& zone_VFF,
                              const IntVect& connect,
                              DoubleTab& incoG,
                              const DoubleTab& incoF,int nbcomp)
{
  valeur().restreindre(zone_VFG, zone_VFF, connect, incoG, incoF,nbcomp);
}



void Restriction::restreindre(const Zone_VF& zone_VFG, const Zone_VF& zone_VFF,
                              const IntVect& connect,
                              DoubleTab& incoG,
                              const DoubleTab& incoF,int nbcomp,
                              int num_prem_face_frontG)
{
  valeur().restreindre(zone_VFG, zone_VFF, connect, incoG, incoF,nbcomp, num_prem_face_frontG);
}


void Restriction::calculer(const Zone_VF& zone_VFG,
                           const Zone_VF& zone_VFF, const IntVect& connect)
{
  valeur().calculer(zone_VFG, zone_VFF, connect);
}

#endif
