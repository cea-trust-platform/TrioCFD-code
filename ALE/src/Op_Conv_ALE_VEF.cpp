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
// File:        Op_Conv_ALE_VEF.cpp
// Directory:   $TRUST_ROOT/src/ALE
// Version:     /main/14
//
//////////////////////////////////////////////////////////////////////////////

#include <Op_Conv_ALE_VEF.h>
#include <Op_Conv_VEF_base.h>
#include <Domaine_ALE.h>

Implemente_instanciable(Op_Conv_ALE_VEF,"Op_Conv_ALE_VEF",Op_Conv_ALE);


Sortie& Op_Conv_ALE_VEF::printOn(Sortie& os) const
{
  return os;
}


Entree& Op_Conv_ALE_VEF::readOn(Entree& is)
{
  Cerr << "OpConvALE_VEF ::readOn " << finl;
  op_conv.associer_eqn(equation());
  op_conv.associer_vitesse(la_vitesse.valeur());
  op_conv.lire(is);
  Cerr << "OpConvALE_VEF : " << op_conv.que_suis_je() << finl;
  return is;
}

void Op_Conv_ALE_VEF::associer (const Zone_dis& zone_dis ,
                                const Zone_Cl_dis& zone_cl_dis,
                                const Champ_Inc& inco )
{
  Cerr << "Op_Conv_ALE_VEF::associer" << finl;
  const Zone_VEF& zvef = ref_cast(Zone_VEF,zone_dis.valeur());
  const Zone_Cl_VEF& zclvef = ref_cast(Zone_Cl_VEF,zone_cl_dis.valeur());

  la_zone_vef = zvef;
  la_zcl_vef = zclvef;
  Op_Conv_ALE::associer(zone_dis,zone_cl_dis,inco);
}

DoubleTab& Op_Conv_ALE_VEF::ajouterALE(const DoubleTab& inco, DoubleTab& resu) const
{
  Cerr << "Op_Conv_ALE_VEF::ajouterALE" << finl;
  //   Cerr << "inco ******* :" << finl;
  //   Cerr << inco << finl;
  //   Cerr << "resu avant contribution ALE******* :" << finl;
  //   Cerr << resu << finl;

  int i,k,num_face;
  const Zone_Cl_VEF& zone_Cl_VEF = la_zcl_vef.valeur();
  const Zone_VEF& zone_VEF = la_zone_vef.valeur();
  const Zone_VF& zone_VF = ref_cast(Zone_VF, zone_VEF);

  const DoubleVect& volumes=zone_VF.volumes();
  const DoubleVect& volumes_entrelaces = zone_VEF.volumes_entrelaces();
  const DoubleVect& volumes_entrelaces_Cl = zone_Cl_VEF.volumes_entrelaces_Cl();
  //   Cerr << "ajouterALE::volumes_entrelaces_Cl : " << finl;
  //   Cerr << volumes_entrelaces_Cl << finl;
  //const Zone& zone = zone_VEF.zone();
  const int nb_faces = zone_VEF.nb_faces();
  //const int nb_elem = zone_VEF.nb_elem();
  const int nb_elem_tot = zone_VEF.nb_elem_tot();
  DoubleTab gradient_elem;//(nb_elem,dimension,dimension);
  //DoubleTab gradient_elem_scalaire(nb_elem,dimension);

  const DoubleTab& vitesse_faces = ref_cast(Domaine_ALE, dom.valeur()).vitesse_faces();
  const IntTab& face_voisins = zone_VF.face_voisins();

  //const int& nb_bord=zone.faces_bord().nb_bords();
  //const int& num_premiere_face_interne=zone_VF.premiere_face_int();
  int premiere_face_int = zone_VEF.premiere_face_int();
  //int nb_faces_bord;
  //const int& num_premiere_face_std=zone_VEF.premiere_face_std();

  double c_grad_scalaire_0;
  double c_grad_scalaire_1;
  DoubleVect c_grad_u_0(dimension);
  DoubleVect c_grad_u_1(dimension);

  int ncomp_ch_transporte=dimension;
  Motcle le_nom_eqn=equation().le_nom();
  if (inco.nb_dim() == 1)
    {
      ncomp_ch_transporte=1;
      if (le_nom_eqn=="pbNavier_Stokes_standard")
        {
          resu.resize(nb_faces);
        }
    }
  gradient_elem.resize(nb_elem_tot,ncomp_ch_transporte,dimension);
  //  Cerr << "resu avant contribution ALE******* :" << finl;
  //  Cerr << resu << finl;
  //   Cerr << "ncomp_ch_transporte : " << ncomp_ch_transporte << finl;
  //  gradient_elem.resize(nb_elem,ncomp_ch_transporte,dimension);
  // gradient_elem=0;
  calculer_gradientP1NC(inco,zone_VEF,zone_Cl_VEF ,gradient_elem);
  //   Cerr << "gradient_elem : " << finl;
  //   Cerr << gradient_elem << finl;
  //const DoubleTab& xv = zone_VEF.xv();
  if (ncomp_ch_transporte != 1)
    {
      // Traitement des bords

      const Conds_lim& les_cl = zone_Cl_VEF.les_conditions_limites();
      int nb_cl=les_cl.size();
      for (int num_cl=0; num_cl<nb_cl; num_cl++)
        {
          const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(num_cl);
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          //const IntTab& elem_faces = zone_VEF.elem_faces();
          int nb_faces_bord=le_bord.nb_faces();
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + nb_faces_bord;

          if (sub_type(Periodique,la_cl.valeur()))
            {
              const Periodique& la_cl_perio = ref_cast(Periodique,la_cl.valeur());
              // Patrick : pour eviter les mauvaises surprises
              IntVect fait(nb_faces_bord);
              fait = 0;
              for (num_face=num1; num_face<num2; num_face++)
                {
                  int face_associee=la_cl_perio.face_associee(num_face-num1);
                  if (fait(num_face-num1) == 0)
                    {
                      fait(num_face-num1) = 1;
                      fait(face_associee) = 1;
                      int elem_0 = face_voisins(num_face,0);
                      int elem_1 = face_voisins(num_face,1);
                      double vol_0=volumes(elem_0);
                      double vol_1=volumes(elem_1);
                      c_grad_u_0=0.;
                      c_grad_u_1=0.;

                      for (k=0; k<dimension; k++)
                        for (i=0; i<dimension; i++)
                          {
                            c_grad_u_0(k)+=vitesse_faces(num_face,i)*gradient_elem(elem_0,k,i);
                            c_grad_u_1(k)+=vitesse_faces(num_face,i)*gradient_elem(elem_1,k,i);
                          }

                      for (k=0; k<dimension; k++)
                        resu(num_face,k)+=(volumes_entrelaces_Cl(num_face)/(vol_0+vol_1))*(vol_0*c_grad_u_0(k) + vol_1*c_grad_u_1(k));

                      for (k=0; k<dimension; k++)
                        resu(num_face+num1,k)=resu(num_face,k);

                    }
                } // Fin de la boucle sur les faces
            } // Fin du cas periodique
          else
            {
              for (num_face=num1; num_face<num2; num_face++)
                {
                  int elem_0 = face_voisins(num_face,0);
                  //double vol_0=volumes(elem_0);
                  c_grad_u_0=0.;

                  for (k=0; k<dimension; k++)
                    for (i=0; i<dimension; i++)
                      c_grad_u_0(k)+=vitesse_faces(num_face,i)*gradient_elem(elem_0,k,i);

                  for (k=0; k<dimension; k++)
                    resu(num_face,k)+=volumes_entrelaces_Cl(num_face)*c_grad_u_0(k);
                }
            }

        }// fin traitement des bords


      // Traitement des faces internes

      for (num_face=premiere_face_int; num_face<nb_faces; num_face++)
        {
          int elem_0=face_voisins(num_face,0);
          int elem_1=face_voisins(num_face,1);
          double vol_0=volumes(elem_0);
          double vol_1=volumes(elem_1);
          c_grad_u_0=0.;
          c_grad_u_1=0.;

          for (k=0; k<dimension; k++)
            for(i=0; i<dimension; i++)
              {
                c_grad_u_0(k)+=vitesse_faces(num_face,i)*gradient_elem(elem_0,k,i);
                c_grad_u_1(k)+=vitesse_faces(num_face,i)*gradient_elem(elem_1,k,i);
              }

          for (k=0; k<dimension; k++)
            {
              resu(num_face,k)+=(volumes_entrelaces(num_face)/(vol_0+vol_1))*(vol_0*c_grad_u_0(k) + vol_1*c_grad_u_1(k));
            }
        }
    }
  else
    {
      const Conds_lim& les_cl = zone_Cl_VEF.les_conditions_limites();
      int nb_cl=les_cl.size();
      for (int num_cl=0; num_cl<nb_cl; num_cl++)
        {
          const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(num_cl);
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          //const IntTab& elem_faces = zone_VEF.elem_faces();
          int nb_faces_bord=le_bord.nb_faces();
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + nb_faces_bord;
          if (sub_type(Periodique,la_cl.valeur()))
            {
              const Periodique& la_cl_perio = ref_cast(Periodique,la_cl.valeur());
              // Patrick : pour eviter les mauvaises surprises
              IntVect fait(nb_faces_bord);
              fait = 0;
              for (num_face=num1; num_face<num2; num_face++)
                {
                  int face_associee=la_cl_perio.face_associee(num_face-num1);
                  if (fait(num_face-num1) == 0)
                    {
                      fait(num_face-num1) = 1;
                      fait(face_associee) = 1;
                      int elem_0 = face_voisins(num_face,0);
                      int elem_1 = face_voisins(num_face,1);
                      double vol_0=volumes(elem_0);
                      double vol_1=volumes(elem_1);
                      c_grad_scalaire_0=0;
                      c_grad_scalaire_1=0;
                      for (i=0; i<dimension; i++)
                        {
                          c_grad_scalaire_0+=vitesse_faces(num_face,i)*gradient_elem(elem_0,0,i);
                          c_grad_scalaire_1+=vitesse_faces(num_face,i)*gradient_elem(elem_1,0,i);
                        }

                      resu(num_face)+=(volumes_entrelaces_Cl(num_face)/(vol_0+vol_1))*(vol_0*c_grad_scalaire_0 + vol_1*c_grad_scalaire_1);
                      resu(num_face+num1)=resu(num_face);

                    }
                } // Fin de la boucle sur les faces
            } // Fin du cas periodique
          else
            {
              for (num_face=num1; num_face<num2; num_face++)
                {
                  int elem_0 = face_voisins(num_face,0);
                  //double vol_0=volumes(elem_0);
                  c_grad_scalaire_0=0;

                  for (i=0; i<dimension; i++)
                    c_grad_scalaire_0+=vitesse_faces(num_face,i)*gradient_elem(elem_0,0,i);

                  resu(num_face)+=volumes_entrelaces_Cl(num_face)*c_grad_scalaire_0;
                }
            }

        }// fin traitement des bords


      // Traitement des faces internes

      for (num_face=premiere_face_int; num_face<nb_faces; num_face++)
        {
          int elem_0=face_voisins(num_face,0);
          int elem_1=face_voisins(num_face,1);
          double vol_0=volumes(elem_0);
          double vol_1=volumes(elem_1);
          c_grad_scalaire_0=0;
          c_grad_scalaire_1=0;
          for(i=0; i<dimension; i++)
            {
              c_grad_scalaire_0+=vitesse_faces(num_face,i)*gradient_elem(elem_0,0,i);
              c_grad_scalaire_1+=vitesse_faces(num_face,i)*gradient_elem(elem_1,0,i);
            }
          resu(num_face)+=(volumes_entrelaces(num_face)/(vol_0+vol_1))*(vol_0*c_grad_scalaire_0 + vol_1*c_grad_scalaire_1);
        }
    }
  /*
  //Affichage variable pour postraitement manuel
  if (ncomp_ch_transporte != 1)
  {
  Cerr << "X" << "\tY" << "\tVITESSE U " << finl;
  //Affichage poiseuille
  //       Cerr << "Y" << "\tVITESSE U " << finl;
  for (num_face=0; num_face<nb_faces; num_face++)
  {
  Cerr << xv(num_face,0) << " \t" <<  xv(num_face,1) << "\t" <<  inco(num_face,0) << finl;
  // Affichage poiseuille
  //           if (xv(num_face,0)==xv(20,0)) Cerr << xv(num_face,1) << "\t" << inco(num_face,0) << finl;
  }
  }
  else
  {
  Cerr << " X " << "\tY" << "\tTEMPERATURE" << finl;
  for (num_face=0; num_face<nb_faces; num_face++)
  {
  // Cerr << " X " << xv(num_face,0) << "  Y " << xv(num_face,1) << " TEMPERATURE  " << inco(num_face) << " RESU " << resu(num_face) << " VOL_ENTR  " << volumes_entrelaces(num_face) << " VOL_ENTR_Cl  " << volumes_entrelaces_Cl(num_face) <<finl;
  Cerr << xv(num_face,0) << "\t" << xv(num_face,1) << "\t" << inco(num_face) << finl;
  }
  }
  */

  return resu;
}


DoubleTab& Op_Conv_ALE_VEF::supprimerALE(const DoubleTab& inco, DoubleTab& resu) const
{
  Cerr << "Op_Conv_ALE_VEF::ajouterALE" << finl;
  int i,k=-1,n,num_face;
  const Zone_Cl_VEF& zone_Cl_VEF = la_zcl_vef.valeur();
  const Zone_VEF& zone_VEF = la_zone_vef.valeur();
  const Zone_VF& zone_VF = ref_cast(Zone_VF, zone_VEF);

  const DoubleVect& volumes=zone_VF.volumes();
  const DoubleVect& volumes_entrelaces = zone_VEF.volumes_entrelaces();
  const DoubleVect& volumes_entrelaces_Cl = zone_Cl_VEF.volumes_entrelaces_Cl();

  const Zone& zone = zone_VEF.zone();
  const int nb_faces = zone_VEF.nb_faces();
  const int nb_elem = zone_VEF.nb_elem();

  DoubleTab gradient_elem;//(nb_elem,dimension,dimension);
  //DoubleTab gradient_elem_scalaire(nb_elem,dimension);

  const DoubleTab& vitesse_faces = ref_cast(Domaine_ALE, dom.valeur()).vitesse_faces();
  const IntTab& face_voisins = zone_VF.face_voisins();

  const int& nb_bord=zone.faces_bord().nb_bords();
  //const int& num_premiere_face_interne=zone_VF.premiere_face_int();
  int nb_faces_bord=0;
  const int& num_premiere_face_std=zone_VEF.premiere_face_std();

  double c_grad_scalaire_0;
  double c_grad_scalaire_1;
  DoubleVect c_grad_u_0(dimension);
  DoubleVect c_grad_u_1(dimension);

  int ncomp_ch_transporte=dimension;
  if (inco.nb_dim() == 1)
    ncomp_ch_transporte=1;
  gradient_elem.resize(nb_elem,ncomp_ch_transporte,dimension);
  gradient_elem=0;
  calculer_gradientP1NC(inco,zone_VEF,zone_Cl_VEF ,gradient_elem);
  if (ncomp_ch_transporte != 1)
    {
      for(n=0; n<nb_bord; n++)
        {
          const int nb_faces_bord_courant=zone.faces_bord()(n).nb_faces();
          const int num_premiere_face_bord_courant=zone.faces_bord()(n).num_premiere_face();
          const int num_derniere_face_bord_courant=num_premiere_face_bord_courant+ nb_faces_bord_courant-1;
          nb_faces_bord+=nb_faces_bord_courant;
          for (num_face=num_premiere_face_bord_courant; num_face<=num_derniere_face_bord_courant; num_face++)
            {
              int elem_non_std=face_voisins(num_face,0);
              c_grad_u_0=0;
              for (k=0; k<dimension; k++)
                {
                  for(i=0; i<dimension; i++)
                    {
                      c_grad_u_0(k)+=vitesse_faces(num_face,i)*gradient_elem(elem_non_std,k,i);
                    }
                }
              for (k=0; k<dimension; k++)
                {
                  resu(num_face,k)-=volumes_entrelaces_Cl(num_face)*c_grad_u_0(k);
                }
            }
        }
      Cout << finl;
      Cout << finl;
      Cout << finl;
      for (num_face=nb_faces_bord; num_face<num_premiere_face_std; num_face++)
        {
          int elem_0=face_voisins(num_face,0);
          double vol_0=volumes(elem_0);
          int elem_1=face_voisins(num_face,1);
          double vol_1=volumes(elem_1);
          c_grad_u_0=0;
          c_grad_u_1=0;
          for (k=0; k<dimension; k++)
            {
              for(i=0; i<dimension; i++)
                {
                  c_grad_u_0(k)+=vitesse_faces(num_face,i)*gradient_elem(elem_0,k,i);
                  c_grad_u_1(k)+=vitesse_faces(num_face,i)*gradient_elem(elem_1,k,i);
                }
            }
          for (k=0; k<dimension; k++)
            {
              resu(num_face,k)-=(volumes_entrelaces_Cl(num_face)/(vol_0+vol_1))*(vol_0*c_grad_u_0(k) + vol_1*c_grad_u_1(k));
            }
        }
      for (num_face=num_premiere_face_std; num_face<nb_faces; num_face++)
        {
          int elem_0=face_voisins(num_face,0);
          double vol_0=volumes(elem_0);
          int elem_1=face_voisins(num_face,1);
          double vol_1=volumes(elem_1);
          c_grad_u_0=0;
          c_grad_u_1=0;
          for (k=0; k<dimension; k++)
            {
              for(i=0; i<dimension; i++)
                {
                  c_grad_u_0(k)+=vitesse_faces(num_face,i)*gradient_elem(elem_0,k,i);
                  c_grad_u_1(k)+=vitesse_faces(num_face,i)*gradient_elem(elem_1,k,i);
                }
            }
          for (k=0; k<dimension; k++)
            {
              resu(num_face,k)-=(volumes_entrelaces(num_face)/(vol_0+vol_1))*(vol_0*c_grad_u_0(k) + vol_1*c_grad_u_1(k));
            }
        }
    }
  else
    {
      for(n=0; n<nb_bord; n++)
        {
          const int nb_faces_bord_courant=zone.faces_bord()(n).nb_faces();
          const int num_premiere_face_bord_courant=zone.faces_bord()(n).num_premiere_face();
          const int num_derniere_face_bord_courant=num_premiere_face_bord_courant+ nb_faces_bord_courant-1;
          nb_faces_bord+=nb_faces_bord_courant;
          for (num_face=num_premiere_face_bord_courant; num_face<=num_derniere_face_bord_courant; num_face++)
            {
              int elem_non_std=face_voisins(num_face,0);
              c_grad_scalaire_0=0;
              for(i=0; i<dimension; i++)
                {
                  c_grad_scalaire_0+=vitesse_faces(num_face,i)*gradient_elem(elem_non_std,0,i);
                }
              resu(num_face)-=volumes_entrelaces_Cl(num_face)*c_grad_scalaire_0;
            }
        }

      for (num_face=nb_faces_bord; num_face<num_premiere_face_std; num_face++)
        {
          int elem_0=face_voisins(num_face,0);
          double vol_0=volumes(elem_0);
          int elem_1=face_voisins(num_face,1);
          double vol_1=volumes(elem_1);
          c_grad_scalaire_0=0;
          c_grad_scalaire_1=0;
          for(i=0; i<dimension; i++)
            {
              c_grad_scalaire_0+=vitesse_faces(num_face,i)*gradient_elem(elem_0,0,i);
              c_grad_scalaire_1+=vitesse_faces(num_face,i)*gradient_elem(elem_1,0,i);
            }
          resu(num_face)-=(volumes_entrelaces_Cl(num_face)/(vol_0+vol_1))*(vol_0*c_grad_u_0(k) + vol_1*c_grad_u_1(k));
        }
      for (num_face=num_premiere_face_std; num_face<nb_faces; num_face++)
        {
          int elem_0=face_voisins(num_face,0);
          double vol_0=volumes(elem_0);
          int elem_1=face_voisins(num_face,1);
          double vol_1=volumes(elem_1);
          c_grad_scalaire_0=0;
          c_grad_scalaire_1=0;
          for(i=0; i<dimension; i++)
            {
              c_grad_scalaire_0+=vitesse_faces(num_face,i)*gradient_elem(elem_0,0,i);
              c_grad_scalaire_1+=vitesse_faces(num_face,i)*gradient_elem(elem_1,0,i);
            }
          resu(num_face)-=(volumes_entrelaces(num_face)/(vol_0+vol_1))*(vol_0*c_grad_scalaire_0 + vol_1*c_grad_scalaire_1);
        }
    }
  return resu;
}
