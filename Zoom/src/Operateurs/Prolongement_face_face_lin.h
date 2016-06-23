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
// File:        Prolongement_face_face_lin.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Operateurs
// Version:     /main/8
//
//////////////////////////////////////////////////////////////////////////////


#ifndef Prolongement_face_face_lin_included
#define Prolongement_face_face_lin_included


/* #include <Objet_U.h> */
/* #include <IntTab.h> */
/* #include <Zone_VF.h> */
/* #include <Pb_1G.h> */
/* #include <Pb_fin.h> */
/* #include <Probleme_base.h> */
/* #include <DoubleTab.h> */
/* #include <Zone_dis.h> */
/* #include <Connectivites.h> */
/* #include <Champ_Inc_base.h> */
#include <Prolongement_base.h>
#include <Zone_VDF.h>
#include <Zone_VEF.h>


class Pb_1G;

//
// .DESCRIPTION class Prolongement_face_face
//
// .SECTION voir aussi


//////////////////////////////////////////////////////////////////////////////////
//
// CLASS: Prolongement_face_face_lin
//
/////////////////////////////////////////////////////////////////////////////////

class Prolongement_face_face_lin : public Prolongement_base
{
  Declare_instanciable(Prolongement_face_face_lin);

public:

  void prolonger(Zone_VF& zone_VFG, Zone_VF& zone_VFF,
                 const Frontiere& frontF,IntVect& connect,
                 const DoubleTab& incoG,
                 DoubleTab& tab, int nb_comp);
  void calculer(Zone_VF& zonef,
                Zone_VF& zoneg,
                IntVect& connect_ff);

private:


};




#endif
