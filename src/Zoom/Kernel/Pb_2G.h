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
// File:        Pb_2G.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Zoom/src/Kernel
// Version:     /main/12
//
//////////////////////////////////////////////////////////////////////////////

// Tout ca n'est pas tres beau, juste pour que ca continue a marcher...


#ifndef Pb_2G_included
#define Pb_2G_included


#include <TRUSTTabs_forward.h>
#include <Connectivites_base.h>
#include <TRUST_List.h>
#include <Prolongement.h>
#include <Restriction.h>
#include <TRUST_Ref.h>
#include <Domaine_forward.h>

class Pb_MG;
class Prolongement;
class Restriction;
class Pb_grossier;
class Champ_front_zoom;
class Pb_1G;
class Probleme_base;
/*! @brief class Pb_2G
 *
 */

//////////////////////////////////////////////////////////////////////////////////
//
// CLASS: Pb_2G
//
/////////////////////////////////////////////////////////////////////////////////

class Pb_2G : public Objet_U
{
  Declare_instanciable(Pb_2G);

public:
  //inline Pb_2G(Pb_fin&);
  inline void set_nb_prol(int i);
  inline void set_nb_rest(int i);
  inline OWN_PTR(Connectivites_base)& connectivites();
  Probleme_base& pb_Fin();
  const Probleme_base& pb_Fin() const ;
  inline Prolongement& mon_prolongement(int i);
  inline LIST(Prolongement)& mon_prolongement();
  inline Restriction& ma_restriction(int i);
  inline LIST(Restriction)& ma_restriction();
  Probleme_base& pbG();
  inline int nb_prolongement();
  inline int nb_restriction();
  inline void typer_Connectivites(const char* type)
  {
    connectivites_ff_ee.typer(type);
  }
  void preparer_calcul_pbFin();
  void postraiter_pbFin();
  void calculer_connectivites_2G();
  void prolonger_et_resoudre_pbFin();
  void prolonger_2G(IntVect& connect, DoubleTab& valG, int nb_compo,
                    const Frontiere& front, DoubleTab& tab, int num_prolongement);
  void restreindre_2G(IntVect& connect, DoubleTab& val_incoG, const DoubleTab& val_incoF, int nb_comp_incoG, int num_rest);
  void restreindre_2G(IntVect& connect, DoubleTab& val_incoG, const DoubleTab& val_incoF, int nb_comp_incoG, int num_prem_face_frontG, int num_rest);
  void associer_probleme_fin(Pb_MG& mg, int index);

private:
  REF(Pb_MG) pb_MG;
  int index_pb_fin;
  OWN_PTR(Connectivites_base) connectivites_ff_ee;
  LIST(Prolongement) mon_prolongement_;
  LIST(Restriction) ma_restriction_;
  int nb_prolongement_;
  int nb_restriction_;
};



inline OWN_PTR(Connectivites_base)& Pb_2G::connectivites()
{
  return connectivites_ff_ee;
}


inline Prolongement& Pb_2G::mon_prolongement(int i)
{
  return mon_prolongement_(i);
}


inline Restriction& Pb_2G::ma_restriction(int i)
{
  return ma_restriction_(i);
}


inline int Pb_2G::nb_prolongement()
{
  return nb_prolongement_;
}


inline int Pb_2G::nb_restriction()
{
  return nb_restriction_;
}


inline LIST(Prolongement)& Pb_2G::mon_prolongement()
{
  return mon_prolongement_;
}


inline LIST(Restriction)& Pb_2G::ma_restriction()
{
  return ma_restriction_;
}


inline void Pb_2G::set_nb_prol(int i)
{
  nb_prolongement_ = i;
}

inline void Pb_2G::set_nb_rest(int i)
{
  nb_restriction_ = i;
}

#endif
