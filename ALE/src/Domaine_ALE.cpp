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
// File:        Domaine_ALE.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/ALE/src
// Version:     /main/21
//
//////////////////////////////////////////////////////////////////////////////

#include <Domaine_ALE.h>
#include <Probleme_base.h>
#include <Zone_VDF.h>
#include <Zone_VEF.h>
#include <Motcle.h>
#include <EFichier.h>
#include <Champ_front_ALE.h>

Implemente_instanciable(Domaine_ALE,"Domaine_ALE",Domaine);

Sortie& Domaine_ALE::printOn(Sortie& os) const
{
  Domaine::printOn(os);
  return os ;
}


Entree& Domaine_ALE::readOn(Entree& is)
{
  Domaine::readOn(is);
  return is ;
}

void Domaine_ALE::mettre_a_jour (double temps, Domaine_dis& le_domaine_dis, Probleme_base& pb)
{
  zone(0).invalide_octree();
  //Modification des coordonnees du maillage
  int N_som=nb_som();
  la_vitesse = calculer_vitesse(temps,le_domaine_dis,pb);
  for (int i=0; i<N_som; i++)
    {
      for (int k=0; k<dimension; k++)
        {
          coord(i,k)+=la_vitesse(i,k)*dt;
        }
    }

  //On recalcule les vitesses aux faces
  Zone_VF& la_zone_VF=ref_cast(Zone_VF,le_domaine_dis.zone_dis(0).valeur());
  int nb_faces=la_zone_VF.nb_faces();
  int nb_som_face=la_zone_VF.nb_som_face();
  IntTab& face_sommets=la_zone_VF.face_sommets();
  calculer_vitesse_faces(la_vitesse,nb_faces,nb_som_face,face_sommets);

  //On recalcule les metriques
  Zone la_zone=les_zones(0);
  la_zone_VF.volumes()=0;
  la_zone.calculer_volumes(la_zone_VF.volumes(),la_zone_VF.inverse_volumes());
  la_zone_VF.xp()=0;
  la_zone.calculer_centres_gravite(la_zone_VF.xp());

  DoubleTab& xv=la_zone_VF.xv();
  xv.set_smart_resize(1);
  xv.reset();
  Type_Face type_face=la_zone.type_elem().type_face();
  IntTab& elem_faces=la_zone_VF.elem_faces();
  IntTab& face_voisins=la_zone_VF.face_voisins();

  calculer_centres_gravite(xv, type_face,
                           sommets, face_sommets);
  if(sub_type(Zone_VDF, la_zone_VF))
    {
      Zone_VDF& la_zone_VDF=ref_cast(Zone_VDF,le_domaine_dis.zone_dis(0).valeur());
      la_zone_VF.volumes_entrelaces()=0;
      la_zone_VDF.calculer_volumes_entrelaces();

    }
  else if(sub_type(Zone_VEF, la_zone_VF))
    {
      Zone_VEF& la_zone_VEF=ref_cast(Zone_VEF,le_domaine_dis.zone_dis(0).valeur());
      DoubleTab& normales=la_zone_VEF.face_normales();
      DoubleTab& facette_normales_=la_zone_VEF.facette_normales();
      IntVect& rang_elem_non_standard=la_zone_VEF.rang_elem_non_std();
      la_zone_VF.volumes_entrelaces()=0;
      la_zone_VEF.calculer_volumes_entrelaces();
      // Recalcul des surfaces avec les normales:
      // PL: je trouve etonnant que le calcul des surfaces se fasse AVANT le calcul des normales
      // PL: il devrait se faire APRES
      int nb_faces_tot=face_sommets.dimension_tot(0);
      DoubleVect face_surfaces_(nb_faces_tot);
      for (int i=0; i<nb_faces_tot; i++)
        {
          double surf=0;
          for (int k=0; k<dimension; k++)
            surf += (la_zone_VF.face_normales(i,k)*la_zone_VF.face_normales(i,k));
          face_surfaces_(i) = sqrt(surf);
        }
      la_zone_VF.calculer_face_surfaces(face_surfaces_);

      la_zone_VEF.calculer_h_carre();
      const Elem_VEF& type_elem=la_zone_VEF.type_elem();
      //int nb_faces_tot=face_sommets.dimension_tot(0);
      // Recalcul des normales
      normales=0;
      for (int num_face=0; num_face<nb_faces_tot; num_face++)
        type_elem.normale(num_face,normales, face_sommets,
                          face_voisins,elem_faces,
                          la_zone) ;
      //       Cerr << "face_normale" << normales << finl;
      type_elem.creer_facette_normales(la_zone, facette_normales_, rang_elem_non_standard);
      Cerr << "carre_pas_du_maillage : " << la_zone_VEF.carre_pas_du_maillage() << finl;
      int nb_eqn=pb.nombre_d_equations();

      for(int num_eq=0; num_eq<nb_eqn; num_eq++)
        {
          Zone_Cl_dis& zcl_dis=pb.equation(num_eq).zone_Cl_dis();
          Zone_Cl_VEF& la_zcl_VEF=ref_cast(Zone_Cl_VEF, zcl_dis.valeur());
          la_zcl_VEF.remplir_volumes_entrelaces_Cl(la_zone_VEF);
          la_zcl_VEF.remplir_normales_facettes_Cl(la_zone_VEF );
        }
    }
  else
    {
      Cerr << "Discretisation non reconnue par ALE!" << finl;
      exit();
    }
}
void Domaine_ALE::initialiser (double temps, Domaine_dis& le_domaine_dis,Probleme_base& pb)
{
  deformable_=1;
  zone(0).invalide_octree();
  la_vitesse=calculer_vitesse(temps,le_domaine_dis,pb);
  //On initialise les vitesses aux faces
  Zone_VF& la_zone_VF=ref_cast(Zone_VF,le_domaine_dis.zone_dis(0).valeur());
  int nb_faces=la_zone_VF.nb_faces();
  int nb_som_face=la_zone_VF.nb_som_face();
  IntTab& face_sommets=la_zone_VF.face_sommets();
  vf.resize(nb_faces, dimension);
  calculer_vitesse_faces(la_vitesse,nb_faces,nb_som_face,face_sommets);
}

DoubleTab Domaine_ALE::calculer_vitesse(double temps, Domaine_dis& le_domaine_dis,Probleme_base& pb)
{

  /////////////////////////////////////////
  //Cas test de la conduction+hydraulique//
  /////////////////////////////////////////
  //   int i,j,k,n;
  //   double T_un=0.3;
  //   double T_deux=2*T_un;
  //   double w_un=(2*3.14)/T_un;
  //   double w_deux=(2*3.14)/T_deux;
  //   int N_som=nb_som();
  //   DoubleTab vit_maillage(N_som,dimension);
  //   Pave& le_pave=ref_cast(Pave,les_zones(0));
  //   int N_som_x=le_pave.Nx+1;
  //   int N_som_y=le_pave.Ny+1;
  //   double L_x=le_pave.Longueurs(0);
  //   double L_y=le_pave.Longueurs(1);
  //   double pas_x=L_x/N_som_x; //Pas de maillage regulier
  //   double pas_y=L_y/N_som_y;
  //   // WARNING
  //   //pas_x=0.1;
  //   //pas_y=0.1;

  //   DoubleVect w(N_som);
  //   double l=w_un;
  // //   for (i=0; i<N_som_x; i++)
  // //     {
  // //       l=(l==w_un) ? w_deux : w_un;
  // //       w(i)=l;
  // //     }

  //   for (i=0; i<N_som_x; i++)
  //     {
  //       l=(l==w_un) ? -w_un : w_un;
  //       w(i)=l;
  //     }
  //   for (n=1; n<N_som_y; n++)
  //     {
  //       for (i=n*N_som_x; i<(n+1)*N_som_x; i++)
  //         {
  //           w(i)=w(i-n*N_som_x);
  //         }
  //     }

  //   for (k=0; k<dimension; k++)
  //     {
  //       if (k!=0)
  //         {
  //           for (i=0; i<N_som; i++)
  //             {
  //               vit_maillage(i,k)=0;
  //             }
  //         }
  //       else
  //         {
  //         for (i=0; i<N_som; i++)
  //           {
  //             //vit_maillage(i,k)=(pas_x/6)*w(i)*cos(w(i)*temps);
  //             //vit_maillage (i,k)=2*coord(i,k)*sin(((2*3.14)/T_un)*temps); //Colonne d'eau
  //             vit_maillage(i,k)=pow(10,3)*pas_x;

  //           }
  //         }
  //     }

  //Zone_VEF& la_zone_VEF=ref_cast(Zone_VEF,le_domaine_dis.zone_dis(0).valeur());
  //const DoubleTab& xv = la_zone_VEF.xv();
  //Les noeuds des bords du domaine ne se deplacent pas

  int n; // A activer ou desactiver si on utilise le laplacien ou non
  int N_som=nb_som(); //A activer ou desactiver si on utilise le laplacien ou non
  //   int N_bord=les_zones(0).faces_bord().nb_bords();
  //   int N_tot_faces_bords=0;
  //   int N_som_par_face;
  //   if (dimension!=3)
  //     N_som_par_face=dimension;
  //   else
  //     N_som_par_face=dimension+1;

  //   for (n=0; n<N_bord; n++)
  //     {
  //       N_tot_faces_bords+=les_zones(0).faces_bord()(n).nb_faces();
  //     }
  //   som_faces_bords.resize(N_tot_faces_bords,N_som_par_face);
  //   N_tot_faces_bords=0;

  //   for (n=0; n<N_bord; n++)
  //     {
  //       N_tot_faces_bords+=les_zones(0).faces_bord()(n).nb_faces();
  //       int N_faces_bord_courant=les_zones(0).faces_bord()(n).nb_faces();
  //       IntTab som_bord_courant=les_zones(0).faces_bord()(n).les_sommets_des_faces();

  //       for (i=N_tot_faces_bords-N_faces_bord_courant; i<N_tot_faces_bords; i++)
  //         {
  //           for (j=0; j<N_som_par_face; j++)
  //             {
  //               som_faces_bords(i,j)=som_bord_courant(i-N_tot_faces_bords+N_faces_bord_courant,j);
  //             }
  //         }
  //     }

  //   for (i=0; i<N_tot_faces_bords; i++)
  //     {
  //       for (j=0;j<N_som_par_face; j++)
  //         {
  //           int t=som_faces_bords(i,j);
  //           for (k=0; k<dimension; k++)
  //             {
  //               vit_maillage(t,k)=0;
  //             }
  //         }
  //     }


  // Les noeuds des bords gauche et droit sont fixes.
  // Par contre ceux des bords haut et bas sont mobiles.
  // A utiliser pour cas test de convection (pic de temperature) , diffusion et constante
  //   int N_bord=les_zones(0).faces_bord().nb_bords();
  //   int N_tot_faces_bords=0;
  //   int N_som_par_face;
  //   if (dimension!=3)
  //     N_som_par_face=dimension;
  //   else
  //     N_som_par_face=dimension+1;

  //   N_tot_faces_bords=les_zones(0).faces_bord()(0).nb_faces()+les_zones(0).faces_bord()(1).nb_faces();
  //   Cerr << "N_tot_faces_bords : " << N_tot_faces_bords << finl;
  //   som_faces_bords.resize(N_tot_faces_bords,N_som_par_face);
  //   N_tot_faces_bords=0;

  //   for (n=0; n<2; n++)
  //     {
  //       N_tot_faces_bords+=les_zones(0).faces_bord()(3*n).nb_faces();
  //       int N_faces_bord_courant=les_zones(0).faces_bord()(3*n).nb_faces();
  //       Cerr << "faces du bord n= " << n << finl;
  //       for (i=les_zones(0).faces_bord()(3*n).num_premiere_face(); i<les_zones(0).faces_bord()(3*n).num_premiere_face()+N_faces_bord_courant; i++)
  //         Cerr << "X: " << xv(i,0) << "Y: " << xv(i,1) << finl;
  //       IntTab som_bord_courant=les_zones(0).faces_bord()(3*n).les_sommets_des_faces();

  //       for (i=N_tot_faces_bords-N_faces_bord_courant; i<N_tot_faces_bords; i++)
  //         {
  //           for (j=0; j<N_som_par_face; j++)
  //             {
  //               som_faces_bords(i,j)=som_bord_courant(i-N_tot_faces_bords+N_faces_bord_courant,j);
  //             }
  //         }
  //     }

  //   for (i=0; i<N_tot_faces_bords; i++)
  //     {
  //       for (j=0;j<N_som_par_face; j++)
  //         {
  //           int t=som_faces_bords(i,j);
  //           for (k=0; k<dimension; k++)
  //             {
  //               vit_maillage (t,k)=0;
  //             }
  //         }
  //     }

  // Les noeuds des bords haut et bas sont fixes.
  // Par contre ceux des bords gauche et droite sont mobiles.
  // A utiliser pour le cas test de poiseuille.

  //   int N_bord=les_zones(0).faces_bord().nb_bords();
  //   int N_tot_faces_bords=0;
  //   int N_som_par_face;
  //   if (dimension!=3)
  //     N_som_par_face=dimension;
  //   else
  //     N_som_par_face=dimension+1;

  //   N_tot_faces_bords=les_zones(0).faces_bord()(2).nb_faces()+les_zones(0).faces_bord()(3).nb_faces();
  //   Cerr << "N_tot_faces_bords : " << N_tot_faces_bords << finl;
  //   som_faces_bords.resize(N_tot_faces_bords,N_som_par_face);
  //   N_tot_faces_bords=0;

  //   for (n=1; n<3; n++)
  //     {
  //       N_tot_faces_bords+=les_zones(0).faces_bord()(n+1).nb_faces();
  //       int N_faces_bord_courant=les_zones(0).faces_bord()(n+1).nb_faces();
  // //       Cerr << "faces du bord n= " << n << finl;
  // //       for (i=les_zones(0).faces_bord()(n+1).num_premiere_face(); i<les_zones(0).faces_bord()(n+1).num_premiere_face()+N_faces_bord_courant; i++)
  // //         Cerr << "X: " << xv(i,0) << "Y: " << xv(i,1) << finl;
  //       IntTab som_bord_courant=les_zones(0).faces_bord()(n+1).les_sommets_des_faces();
  //       for (i=N_tot_faces_bords-N_faces_bord_courant; i<N_tot_faces_bords; i++)
  //         {
  //           for (j=0; j<N_som_par_face; j++)
  //             {
  //               som_faces_bords(i,j)=som_bord_courant(i-N_tot_faces_bords+N_faces_bord_courant,j);
  //             }
  //         }
  //     }
  //   for (i=0; i<N_tot_faces_bords; i++)
  //     {
  //       for (j=0;j<N_som_par_face; j++)
  //         {
  //           int t=som_faces_bords(i,j);
  //           for (k=0; k<dimension; k++)
  //             {
  //               vit_maillage (t,k)=0;
  //             }
  //         }
  //     }


  //   //////////////////////////////////////////////
  //   //Calcul de la vitesse avec le laplacien nul//
  //   //////////////////////////////////////////////
  //   Cerr << "Domaine_ALE::calculer vitesse" << finl;
  //   Cerr << "nombre bords ALE : " << nb_bords_ALE << finl;
  DoubleTab vit_maillage(N_som,dimension);
  DoubleTab vit_bords(N_som,dimension);
  vit_bords=0;
  DoubleTab tab_champ_front(N_som,dimension);
  for (n=0; n<nb_bords_ALE; n++)
    {
      const Nom& le_nom_bord_ALE=les_bords_ALE(n).le_nom();
      //       Cerr << "le_nom_bord_ALE : " << le_nom_bord_ALE << finl;
      int rang=les_zones(0).rang_frontiere(le_nom_bord_ALE);
      const Frontiere_dis_base& la_fr_dis=le_domaine_dis.zone_dis(0).frontiere_dis(rang);
      les_champs_front[n].valeur().associer_fr_dis_base(la_fr_dis);
      const Nom& le_nom_ch_front_courant=les_champs_front[n].valeur().que_suis_je();
      //       Cerr << le_nom_ch_front_courant << finl;
      if (le_nom_ch_front_courant == "Champ_front_ALE")
        {
          ref_cast(Champ_front_ALE, les_champs_front[n].valeur()).remplir_vit_som_bord_ALE(temps);
          vit_bords+=ref_cast(Champ_front_ALE, les_champs_front[n].valeur()).get_vit_som_bord_ALE();
        }
      else
        {
          Cerr << "Un champ front de type : "
               << les_champs_front[n].valeur().le_nom()
               << "ne peut etre utilise pour un probleme ALE pour le moment...."
               << finl;
          exit();
        }
    }
  //   Cerr << "vit_bords :" << vit_bords << finl;
  laplacien(le_domaine_dis, pb, vit_bords, vit_maillage);
  return vit_maillage;
}

DoubleTab& Domaine_ALE::calculer_vitesse_faces(DoubleTab& vit_maillage,int nb_faces,int nb_som_face, IntTab& face_sommets)
{
  int i,j,s;
  vf=0;
  for (j=0; j<dimension; j++)
    {
      for (i=0; i<nb_faces; i++)
        {
          for (s=0; s<nb_som_face; s++)
            {
              vf(i,j)+=vit_maillage(face_sommets(i,s),j);
            }
          vf(i,j)/=nb_som_face;
        }
    }
  return vf;

}

void Domaine_ALE::lecture_vit_bords_ALE(Entree& is)
{
  Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");
  Motcle motlu;
  Nom nomlu;
  is >> motlu;
  if (motlu != accolade_ouverte)
    {
      Cerr << "Erreur a la lecture des vitesses ALE aux bords\n";
      Cerr << "On attendait une " << accolade_ouverte << " a la place de \n"
           << motlu;
      exit();
    }
  is >> nb_bords_ALE;
  Cerr << "nombre de bords ALE : " <<  nb_bords_ALE << finl;
  les_champs_front.dimensionner(nb_bords_ALE);
  int compteur=0;
  while(1)
    {
      // lecture d'un nom de bord ou de }
      is >> nomlu;
      motlu=nomlu;
      if (motlu == accolade_fermee)
        break;
      Cerr << "Lecture des vitesses "
           << "imposees a " << nomlu << "....." << finl;
      int rang=les_zones(0).rang_frontiere(nomlu);
      //       le_bord_ALE=les_zones(0).faces_bord()(rang);
      //       is >> le_champ_front_uniforme;
      //       Cerr << le_champ_front_uniforme.valeurs() << finl;
      les_bords_ALE.add(les_zones(0).faces_bord()(rang));
      is >> les_champs_front[compteur];
      compteur++;
    }
  //   Cerr << "le_champ_front_uniforme : " << finl;
  //   Cerr << le_champ_front_uniforme.valeurs() << finl;
  //   for (int i=0; i<nb_bords_ALE; i++)
  //     {
  //       Cerr << "numero de bord ALE : " << i << finl;
  //       Cerr << les_champs_front[i].valeur().valeurs() << finl;
  //       Cerr << les_champs_front[i].valeur().valeurs()(0,1) << finl;
  //     }
  //Cerr << "les_bords_ALE : " << finl;
  //Cerr << les_bords_ALE << finl;
}

DoubleTab& Domaine_ALE::laplacien(Domaine_dis& le_domaine_dis,Probleme_base& pb, const DoubleTab& vit_bords, DoubleTab& ch_som)
{
  Cerr << "Domaine_ALE::laplacien" << finl;
  const Zone& lazone = les_zones(0);
  int nb_elem_tot = lazone.nb_elem_tot();
  int nbsom = lazone.nb_som();
  int nb_som_elem = lazone.nb_som_elem();
  //const DoubleTab& ch = cha.valeurs();
  const Zone_VEF& zone_VEF=ref_cast(Zone_VEF,le_domaine_dis.zone_dis(0).valeur());
  const DoubleTab& normales=zone_VEF.face_normales();
  const Zone_Cl_VEF& zone_Cl_VEF = ref_cast(Zone_Cl_VEF, pb.equation(0).zone_Cl_dis().valeur());
  const IntTab& elem_som = lazone.les_elems();
  const IntTab& elem_faces=zone_VEF.elem_faces();
  double mijK;
  int nb_comp=ch_som.dimension(1);
  //int dimension=Objet_U::dimension;

  if(dimension==2)
    {
      mijK=1./4.;
    }
  else
    {
      mijK=1./9.;
    }

  int elem;
  int n_bord;
  //  if ( cha.get_MatP1NC2P1().nb_lignes() <2 )
  //if (mat.nb_lignes() <2 )
  {
    EFichier fic("solveur.bar");
    fic >> solv;
    solv.nommer("ALE_solver");
    int rang;
    int nnz=0;
    IntLists voisins(nbsom);
    DoubleLists coeffs(nbsom);
    DoubleVect diag(nbsom);
    for (elem=0; elem<nb_elem_tot; elem++)
      {
        double volume=zone_VEF.volumes(elem);
        for (int isom=0; isom<nb_som_elem; isom++)
          {
            int facei=elem_faces(elem,isom);
            int ii=get_renum_som_perio(elem_som(elem,isom));
            for (int jsom=isom+1; jsom<nb_som_elem; jsom++)
              {
                int i=ii;
                int facej=elem_faces(elem,jsom);
                int j=get_renum_som_perio(elem_som(elem,jsom));

                if(i>j)
                  {
                    int tmp=i;
                    i=j;
                    j=tmp;
                  }
                double coeffij=0.0;
                for(int k=0; k<dimension; k++)
                  coeffij+=normales(facei,k)*normales(facej,k);

                coeffij*=zone_VEF.oriente_normale(facei,elem)*
                         zone_VEF.oriente_normale(facej,elem);
                coeffij*=mijK/volume;
                rang=voisins[i].rang(j);
                if(rang==-1)
                  {
                    voisins[i].add(j);
                    coeffs[i].add(coeffij);
                    nnz++;
                  }
                else
                  {
                    coeffs[i][rang]+=coeffij;
                  }
                diag[i]-=coeffij;
                diag[j]-=coeffij;
              }
          }
      }
    for (n_bord=0; n_bord<zone_VEF.nb_front_Cl(); n_bord++)
      {
        //for n_bord
        const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
        const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
        int num1 = le_bord.num_premiere_face();
        int num2 = num1 + le_bord.nb_faces();

        for (int face=num1; face<num2; face++)
          {
            elem=zone_VEF.face_voisins(face,0);
            for(int isom=0; isom<dimension; isom++)
              {
                int som=zone_VEF.face_sommets(face,isom);
                som=get_renum_som_perio(som);
                diag[som]=1.e6;
              }
          }
      }

    mat.dimensionner(nbsom, nnz+nbsom) ;
    mat.remplir(voisins, coeffs, diag) ;
    mat.compacte() ;
    mat.set_est_definie(1);
    Cerr << "Matrice de filtrage OK" << finl;
    Cerr << "Matrice de Laplacien P1 : " << finl;
    //mat.imprimer(Cerr) << finl;
  }


  //double coeff=-1./dimension;
  DoubleVect position(dimension);
  for(int comp=0; comp<nb_comp; comp++)
    {
      //for comp
      DoubleVect secmem(nbsom);
      DoubleVect solution(nbsom);
      //       for (elem=0; elem<nb_elem_tot; elem++)
      //         {//for elem
      //           double volume=zone_VEF.volumes(elem);
      //           for (int isom=0; isom<nb_som_elem; isom++)
      //             {
      //               int facei=elem_faces(elem,isom);
      //               int i=get_renum_som_perio(elem_som(elem,isom));
      //               for (int jsom=0; jsom<nb_som_elem; jsom++)
      //                 {
      //                   int facej=elem_faces(elem,jsom);
      //                   double coeffij=0;
      //                   for(k=0; k<dimension; k++)
      //                     coeffij+=normales(facei,k)*normales(facej,k);
      //                   coeffij*=zone_VEF.oriente_normale(facei,elem)*
      //                     zone_VEF.oriente_normale(facej,elem);
      //                   coeffij*=coeff/volume;
      //                   secmem(i)+=coeffij*ch(facej,comp);
      //                 }
      //             }
      //         }
      for (n_bord=0; n_bord<zone_VEF.nb_front_Cl(); n_bord++)
        {
          //for n_bord
          const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();

          for (int face=num1; face<num2; face++)
            {
              //elem=zone_VEF.face_voisins(face,0);
              for(int isom=0; isom<dimension; isom++)
                {
                  int som=zone_VEF.face_sommets(face,isom);
                  som=get_renum_som_perio(som);
                  //position(0)=coord(som, 0);
                  //position(1)=coord(som, 1);
                  //if(dimension == 3)
                  //position(2)=coord(som, 2);
                  //secmem(som)=1.e6*cha.valeur_a_elem_compo(position, elem, comp);
                  secmem(som)=1.e6*vit_bords(som,comp);
                }
            }
        }

      Cerr << "Resolution du systeme de filtrage" << finl;
      //solv.resoudre_systeme(cha.get_MatP1NC2P1(), secmem, solution);
      solv.resoudre_systeme(mat, secmem, solution);
      for(int som=0; som<nbsom; som++)
        ch_som(som,comp)=solution(get_renum_som_perio(som));

    }
  // Calcul de la norme du vecteur distribue
  {
    double x = mp_norme_vect(ch_som);
    Cerr << "norme(c) = " << x << finl;
  }
  return ch_som;


}

void Domaine_ALE::set_dt(double& dt_)
{
  dt=dt_;
}
