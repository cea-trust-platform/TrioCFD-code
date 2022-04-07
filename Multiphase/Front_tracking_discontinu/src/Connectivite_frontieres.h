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
// File:        Connectivite_frontieres.h
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Front_tracking_discontinu/src
// Version:     /main/7
//
//////////////////////////////////////////////////////////////////////////////
#ifndef Connectivite_frontieres_included
#define Connectivite_frontieres_included

#include <TRUSTTab.h>
#include <Ref_Zone_VF.h>

class Zone_VF;

class Connectivite_frontieres : public Objet_U
{
  Declare_instanciable(Connectivite_frontieres);
public:
  virtual void     associer_zone_vf(const Zone_VF& zone_vf);
  const IntTab&    faces_voisins() const;
  const ArrOfInt& face_cond_lim() const;
  const IntTab&    def_face_aretes() const;

protected:
  void remplir_def_face_aretes(const Zone_VF& zone_vf);
  void remplir_faces_voisins(const Zone_VF& zone_vf);

  REF(Zone_VF) refzone_vf_;

  // Tableau faces_voisins_ :
  //  dimension(0) = zone.nb_faces_frontiere()
  //  dimension(1) = nombres d'aretes d'une face de bord
  //  Contenu : faces_voisins_(i,j) est l'indice dans la zone
  //            de la face voisine de la face frontiere i par l'arete j.
  //            L'indice i est le numero d'une face frontiere reelle.
  //            L'indice j est un numero d'arete tel qu'il est defini dans
  //            def_face_aretes_
  //            Il existe toujours une facette voisine, eventuellement
  //            c'est une facette virtuelle (indice superieur a nb_faces()).
  IntTab faces_voisins_;

  // Tableau face_cond_lim_ :
  //  size_array() = nb_faces_frontiere()
  //  Contenu : numero de la condition aux limites de la face
  //            ( entre 0 et zone.nb_front_Cl() ) pour les faces reelles.
  ArrOfInt face_cond_lim_;

  // Tableau def_face_aretes_ :
  //  dimension(0) = nombre d'aretes par face
  //  dimension(1) = 2
  //  Contenu :
  //   Ce tableau donne la definition des aretes d'une face :
  //   L'arete i va du sommet def_face_arete_(i,0) au sommet def_face_arete_(i,1) en 3D.
  //   En 2D, c'est les aretes sont reduites a des points, def_face_arete_(i,1) = -1
  //   Le numero du sommet est un numero local sur la face (compris entre 0
  //   et le nombre de sommets par face). Exemple:
  //    i_sommet   = def_face_arete_(i_arete, 0); // numero du sommet sur la face
  //    num_sommet = zone_vf.face_sommets(num_face, i_sommet); // numero du sommet dans la zone
  IntTab def_face_aretes_;

private:
  Connectivite_frontieres(const Connectivite_frontieres&);  // Interdit
  const Connectivite_frontieres& operator=(const Connectivite_frontieres&);   // Interdit
};

inline const IntTab& Connectivite_frontieres::faces_voisins() const
{
  assert(refzone_vf_.non_nul());
  return faces_voisins_;
}

inline const IntTab& Connectivite_frontieres::def_face_aretes() const
{
  assert(refzone_vf_.non_nul());
  return def_face_aretes_;
}

#endif
