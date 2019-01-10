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
// File:        Op_Conv_ALE_VEF.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/ALE/src
// Version:     /main/14
//
//////////////////////////////////////////////////////////////////////////////

#include <Op_Conv_ALE_VEF.h>
#include <Op_Conv_VEF_base.h>
#include <Domaine_ALE.h>
#include <DoubleTrav.h>
#include <Op_Conv_VEF_Face.h>
#include <Probleme_base.h>
#include <Schema_Temps_base.h>
#include <Porosites_champ.h>
#include <Debog.h>


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

void Op_Conv_ALE_VEF::ajouterALE_Centre(const DoubleTab& inco,DoubleTab& resu) const
{
  // This in not exactly "center" scheme. It is center-like scheme, where inverse distance (power 3) interpolation is used.
  const Domaine_ALE& dom_ale=ref_cast(Domaine_ALE, dom.valeur());
  const DoubleTab& vitesse_faces_ALE = dom_ale.vitesse_faces();

  int i,k,num_face;
  const Zone_Cl_VEF& zone_Cl_VEF = la_zcl_vef.valeur();
  const Zone_VEF& zone_VEF = la_zone_vef.valeur();
  const Zone_VF& zone_VF = ref_cast(Zone_VF, zone_VEF);
  const DoubleVect& volumes=zone_VF.volumes();
  const DoubleVect& volumes_entrelaces = zone_VEF.volumes_entrelaces();
  const DoubleVect& volumes_entrelaces_Cl = zone_Cl_VEF.volumes_entrelaces_Cl();
  const int nb_faces = zone_VEF.nb_faces();
  const int nb_elem_tot = zone_VEF.nb_elem_tot();
  DoubleTab gradient_elem;//(nb_elem,dimension,dimension);
  const IntTab& face_voisins = zone_VF.face_voisins();
  int premiere_face_int = zone_VEF.premiere_face_int();
  double c_grad_scalaire_0;
  double c_grad_scalaire_1;
  DoubleVect c_grad_u_0(dimension);
  DoubleVect c_grad_u_1(dimension);
  int ncomp_ch_transporte=dimension;
  Motcle le_nom_eqn=equation().le_nom();

  /*if (inco.nb_dim() == 1)
    {
      ncomp_ch_transporte=1;
      if (le_nom_eqn=="pbNavier_Stokes_standard")
        {
          resu.resize(nb_faces);
        }
    }*/
  gradient_elem.resize(nb_elem_tot,ncomp_ch_transporte,dimension);
  calculer_gradientP1NC(inco,zone_VEF,zone_Cl_VEF ,gradient_elem);
  gradient_elem.echange_espace_virtuel();
  if (ncomp_ch_transporte != 1)
    {
      // Traitement des bords

      const Conds_lim& les_cl = zone_Cl_VEF.les_conditions_limites();
      int nb_cl=les_cl.size();
      for (int num_cl=0; num_cl<nb_cl; num_cl++)
        {
          const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(num_cl);
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
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
                            c_grad_u_0(k)+=vitesse_faces_ALE(num_face,i)*gradient_elem(elem_0,k,i);
                            c_grad_u_1(k)+=vitesse_faces_ALE(num_face,i)*gradient_elem(elem_1,k,i);
                          }

                      for (k=0; k<dimension; k++)
                        //Inverse distance weighting interpolation (p=3) is used for obtaining c_grad_u at the faces and it is multiplied by volumes_entrelaces_Cl
                        resu(num_face,k)+=((vol_1*c_grad_u_0(k) + vol_0*c_grad_u_1(k))/(vol_0+vol_1))*volumes_entrelaces_Cl(num_face);

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
                  c_grad_u_0=0.;

                  for (k=0; k<dimension; k++)
                    for (i=0; i<dimension; i++)
                      c_grad_u_0(k)+=vitesse_faces_ALE(num_face,i)*gradient_elem(elem_0,k,i);

                  for (k=0; k<dimension; k++)
                    resu(num_face,k)+=volumes_entrelaces_Cl(num_face)*c_grad_u_0(k);
                }
            }

        }// fin traitement des bords

      resu.echange_espace_virtuel();
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
                c_grad_u_0(k)+=vitesse_faces_ALE(num_face,i)*gradient_elem(elem_0,k,i);
                c_grad_u_1(k)+=vitesse_faces_ALE(num_face,i)*gradient_elem(elem_1,k,i);
              }

          for (k=0; k<dimension; k++)
            {
              //Inverse distance weighting interpolation (p=3) is used for obtaining c_grad_u at the faces and it is multiplied by volumes_entrelaces
              resu(num_face,k)+=((vol_1*c_grad_u_0(k) + vol_0*c_grad_u_1(k))/(vol_0+vol_1))*volumes_entrelaces(num_face);
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
                          c_grad_scalaire_0+=vitesse_faces_ALE(num_face,i)*gradient_elem(elem_0,0,i);
                          c_grad_scalaire_1+=vitesse_faces_ALE(num_face,i)*gradient_elem(elem_1,0,i);
                        }
                      //Inverse distance weighting interpolation (p=3) is used for obtaining c_grad_scalaire at the faces and it is multiplied by volumes_entrelaces_Cl
                      resu(num_face)+=((vol_1*c_grad_scalaire_0 + vol_0*c_grad_scalaire_1)/(vol_0+vol_1))*volumes_entrelaces_Cl(num_face);
                      resu(num_face+num1)=resu(num_face);

                    }
                } // Fin de la boucle sur les faces
            } // Fin du cas periodique
          else
            {
              for (num_face=num1; num_face<num2; num_face++)
                {
                  int elem_0 = face_voisins(num_face,0);
                  c_grad_scalaire_0=0;

                  for (i=0; i<dimension; i++)
                    c_grad_scalaire_0+=vitesse_faces_ALE(num_face,i)*gradient_elem(elem_0,0,i);

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
              c_grad_scalaire_0+=vitesse_faces_ALE(num_face,i)*gradient_elem(elem_0,0,i);
              c_grad_scalaire_1+=vitesse_faces_ALE(num_face,i)*gradient_elem(elem_1,0,i);
            }
          //Inverse distance weighting interpolation (p=3) is used for obtaining c_grad_scalaire at the faces and it is multiplied by volumes_entrelaces
          resu(num_face)+=((vol_1*c_grad_scalaire_0 + vol_0*c_grad_scalaire_1)/(vol_0+vol_1))*volumes_entrelaces(num_face);
        }
    }
}

void Op_Conv_ALE_VEF::ajouterALE_Amont(const DoubleTab& inco, DoubleTab& resu) const
{
  //This is based on Op_Conv_VEF_Face::ajouter(...) from where only "amont" related part is used. ALE convection term has the opposite sign compared to the normal convection term.
  const Domaine_ALE& dom_ale=ref_cast(Domaine_ALE, dom.valeur());
  const Zone_Cl_VEF& zone_Cl_VEF = la_zcl_vef.valeur();
  const Zone_VEF& zone_VEF = ref_cast(Zone_VEF, la_zone_vef.valeur());
  const Op_Conv_VEF_base& opConvVEFbase = ref_cast(Op_Conv_VEF_base, op_conv.valeur());

  const DoubleTab& vitesse_face_absolue_ALE=dom_ale.vitesse_faces();              //ALE mesh velocity.

  const DoubleVect& porosite_face = zone_VEF.porosite_face();
  int marq=opConvVEFbase.phi_u_transportant(equation());

  DoubleTab transporte_face_;
  DoubleTab vitesse_face_;

  const DoubleTab& transporte_face = modif_par_porosite_si_flag(inco,transporte_face_,!marq,porosite_face); //Fluid velocity = inco .
  const DoubleTab& vitesse_face_ALEmodif    = modif_par_porosite_si_flag(vitesse_face_absolue_ALE,vitesse_face_,marq,porosite_face);

  const IntTab& elem_faces = zone_VEF.elem_faces();
  const DoubleTab& facenormales = zone_VEF.face_normales();
  const DoubleTab& facette_normales = zone_VEF.facette_normales();
  const Zone& zone = zone_VEF.zone();
  const int nfa7 = zone_VEF.type_elem().nb_facette();
  const int nb_elem_tot = zone_VEF.nb_elem_tot();
  const IntVect& rang_elem_non_std = zone_VEF.rang_elem_non_std();
  const DoubleTab& normales_facettes_Cl = zone_Cl_VEF.normales_facettes_Cl();
  int nfac = zone.nb_faces_elem();
  int nsom = zone.nb_som_elem();
  int nb_bord = zone_VEF.nb_front_Cl();
  const IntTab& les_elems=zone.les_elems();

  // Permet d'avoir un flux_bord coherent avec les CLs (mais parfois diverge?)
  // Active uniquement pour ordre 3
  int option_appliquer_cl_dirichlet = 0;
  // Definition d'un tableau pour un traitement special des schemas pres des bords

  if (traitement_pres_bord_.size_array()!=nb_elem_tot)
    {
      traitement_pres_bord_.resize_array(nb_elem_tot);
      traitement_pres_bord_=0;


      // Pour le muscl/centre actuels on utilise un calcul de flux a l'ordre 1
      // aux mailles de bord ou aux mailles ayant un sommet de Dirichlet
      ArrOfInt est_un_sommet_de_bord_(zone_VEF.nb_som_tot());
      for (int n_bord=0; n_bord<nb_bord; n_bord++)
        {
          const Cond_lim_base& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord).valeur();
          if ( sub_type(Dirichlet,la_cl) || sub_type(Dirichlet_homogene,la_cl) )
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              int nb_faces_tot = le_bord.nb_faces_tot();
              int size = zone_VEF.face_sommets().dimension(1);
              for (int ind_face=0; ind_face<nb_faces_tot; ind_face++)
                for (int som=0; som<size; som++)
                  {
                    int face = le_bord.num_face(ind_face);
                    est_un_sommet_de_bord_(zone_VEF.face_sommets(face,som))=1;
                  }
            }
        }
      for (int elem=0; elem<nb_elem_tot; elem++)
        {
          if (rang_elem_non_std(elem)!=-1)
            traitement_pres_bord_(elem)=1;
          else
            {
              for (int n_som=0; n_som<nsom; n_som++)
                if (est_un_sommet_de_bord_(les_elems(elem,n_som)))
                  traitement_pres_bord_(elem)=1;
            }
        }

      // Construction du tableau est_une_face_de_dirichlet_
      est_une_face_de_dirichlet_.resize_array(zone_VEF.nb_faces_tot());
      est_une_face_de_dirichlet_=0;
      for (int n_bord=0; n_bord<nb_bord; n_bord++)
        {
          const Cond_lim_base& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord).valeur();
          if ( sub_type(Dirichlet,la_cl) || sub_type(Dirichlet_homogene,la_cl) )
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              int nb_faces_tot = le_bord.nb_faces_tot();
              for (int ind_face=0; ind_face<nb_faces_tot; ind_face++)
                {
                  int num_face = le_bord.num_face(ind_face);
                  est_une_face_de_dirichlet_(num_face) = 1;
                }
            }
        }
    }

  // Pour le traitement de la convection on distingue les polyedres
  // standard qui ne "voient" pas les conditions aux limites et les
  // polyedres non standard qui ont au moins une face sur le bord.
  // Un polyedre standard a n facettes sur lesquelles on applique le
  // schema de convection.
  // Pour un polyedre non standard qui porte des conditions aux limites
  // de Dirichlet, une partie des facettes sont portees par les faces.
  // En bref pour un polyedre le traitement de la convection depend
  // du type (triangle, tetraedre ...) et du nombre de faces de Dirichlet.

  const Elem_VEF_base& type_elemvef= zone_VEF.type_elem().valeur();
  int istetra=0;
  Nom nom_elem=type_elemvef.que_suis_je();
  if ((nom_elem=="Tetra_VEF")||(nom_elem=="Tri_VEF"))
    istetra=1;

  const DoubleVect& porosite_elem = zone_VEF.porosite_elem();
  double psc;
  int poly,face_adj,fa7,i,j,n_bord;
  int num_face, rang;
  int num10,num20;
  int ncomp_ch_transporte=(transporte_face.nb_dim() == 1?1:transporte_face.dimension(1));

  // Traitement particulier pour les faces de periodicite
  int nb_faces_perio = 0;
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          nb_faces_perio+=le_bord.nb_faces();
        }
    }

  DoubleTab tab;
  if (ncomp_ch_transporte == 1)
    tab.resize(nb_faces_perio);
  else
    tab.resize(nb_faces_perio,ncomp_ch_transporte);

  nb_faces_perio=0;
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (num_face=num1; num_face<num2; num_face++)
            {
              if (ncomp_ch_transporte == 1)
                tab(nb_faces_perio) = resu(num_face);
              else
                for (int comp=0; comp<ncomp_ch_transporte; comp++)
                  tab(nb_faces_perio,comp) = resu(num_face,comp);
              nb_faces_perio++;
            }
        }
    }

  int comp0;
  int nb_faces_ = zone_VEF.nb_faces();
  ArrOfInt face(nfac);

  // Tableau gradient base sur gradient_elem selon schema
  DoubleTab gradient_elem(nb_elem_tot,ncomp_ch_transporte,dimension);  // (du/dx du/dy dv/dx dv/dy) pour un poly
  DoubleTab gradient;

  ArrOfDouble vs(dimension);
  ArrOfDouble vc(dimension);
  ArrOfDouble cc(dimension);
  DoubleVect xc(dimension);
  DoubleTab vsom(nsom,dimension);
  DoubleTab xsom(nsom,dimension);

  const IntTab& KEL=type_elemvef.KEL();

  // Boucle ou non selon la valeur de alpha (uniquement a l'ordre 3 pour le moment)
  // Si alpha=1, la boucle se limite a une simple passe avec le schema choisi (muscl, amont, centre)
  // Si alpha<1, la boucle se compose de 2 passes:
  //                         -la premiere avec le schema choisi et une ponderation de alpha
  //                         -la seconde avec le schema centre et une ponderation de 1-alpha
  double alpha = 1;

  // Les polyedres non standard sont ranges en 2 groupes dans la Zone_VEF:
  //  - polyedres bords et joints
  //  - polyedres bords et non joints
  // On traite les polyedres en suivant l'ordre dans lequel ils figurent
  // dans la zone
  // boucle sur les polys

  for (poly=0; poly<nb_elem_tot; poly++)
    {
      int contrib = 0;
      // calcul des numeros des faces du polyedre
      for (face_adj=0; face_adj<nfac; face_adj++)
        {
          int face_ = elem_faces(poly,face_adj);
          face(face_adj)= face_;
          if (face_<nb_faces_) contrib=1; // Une face reelle sur l'element virtuel
        }
      //
      if (contrib)
        {
          for (j=0; j<dimension; j++)
            {
              vs(j) = vitesse_face_absolue_ALE(face(0),j)*porosite_face(face(0));
              for (i=1; i<nfac; i++)
                vs(j)+= vitesse_face_absolue_ALE(face(i),j)*porosite_face(face(i));
            }
          // calcul de la vitesse aux sommets des polyedres
          // On va utliser les fonctions de forme implementees dans la classe Champs_P1_impl ou Champs_Q1_impl
          if (istetra==1)
            {
              for (i=0; i<nsom; i++)
                for (j=0; j<dimension; j++)
                  vsom(i,j) = (vs(j) - dimension*vitesse_face_absolue_ALE(face(i),j)*porosite_face(face(i)));
            }
          else
            {
              Cout << "ALE is currently working only with Tetra_VEF (3D) or Tri_VEF (2D)." << finl;
              exit();
            }

          // Determination du type de CL selon le rang
          rang = rang_elem_non_std(poly);
          int itypcl = (rang==-1 ? 0 : zone_Cl_VEF.type_elem_Cl(rang));

          // calcul de vc (a l'intersection des 3 facettes) vc vs vsom proportionnelles a la porosite
          calcul_vc_ALE(face,vc,vs,vsom,vitesse_face_absolue_ALE,itypcl,porosite_face);

          // Gestion de la porosite
          if (marq==0)
            {
              double coeff=1./porosite_elem(poly);
              vsom*=coeff;
              vc*=coeff;
            }
          // Boucle sur les facettes du polyedre non standard:
          for (fa7=0; fa7<nfa7; fa7++)
            {
              num10 = face(KEL(0,fa7));
              num20 = face(KEL(1,fa7));

              // normales aux facettes
              if (rang==-1)
                for (i=0; i<dimension; i++)
                  cc[i] = facette_normales(poly, fa7, i);
              else
                for (i=0; i<dimension; i++)
                  cc[i] = normales_facettes_Cl(rang,fa7,i);

              // Calcul des vitesses en C,S,S2 les 3 extremites de la fa7 et M le centre de la fa7
              double psc_c=0,psc_s=0,psc_m,psc_s2=0;
              if (dimension==2)
                {
                  for (i=0; i<2; i++)
                    {
                      psc_c+=vc[i]*cc[i];
                      psc_s+=vsom(KEL(2,fa7),i)*cc[i];
                    }
                  psc_m=(psc_c+psc_s)/2.;
                }
              else
                {
                  for (i=0; i<3; i++)
                    {
                      psc_c+=vc[i]*cc[i];
                      psc_s+=vsom(KEL(2,fa7),i)*cc[i];
                      psc_s2+=vsom(KEL(3,fa7),i)*cc[i];
                    }
                  psc_m=(psc_c+psc_s+psc_s2)/3.;
                }
              // On applique les CL de Dirichlet si num1 ou num2 est une face avec CL de Dirichlet
              // auquel cas la fa7 coincide avec la face num1 ou num2 -> C est au centre de la face
              int appliquer_cl_dirichlet=0;
              if (option_appliquer_cl_dirichlet)
                if (est_une_face_de_dirichlet_(num10) || est_une_face_de_dirichlet_(num20))
                  {
                    appliquer_cl_dirichlet = 1;
                    psc_m = psc_c;
                  }

              // Determination de la face amont pour M
              int face_amont_m;
              if (psc_m >= 0)
                {
                  face_amont_m = num10;
                }
              else
                {
                  face_amont_m = num20;
                }

              for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                {
                  double flux;
                  double inco_m = (ncomp_ch_transporte==1?transporte_face(face_amont_m):transporte_face(face_amont_m,comp0));

                  if (true || appliquer_cl_dirichlet)
                    {
                      flux = inco_m*psc_m;
                    }

                  // Ponderation par coefficient alpha
                  flux*=alpha;

                  if (ncomp_ch_transporte == 1)
                    {
                      resu(num10) += flux;
                      resu(num20) -= flux;
                    }
                  else
                    {
                      resu(num10,comp0) += flux;
                      resu(num20,comp0) -= flux;
                    }

                }// boucle sur comp
            } // fin de la boucle sur les facettes
        }
    } // fin de la boucle
  alpha = 1 - alpha;
  resu.echange_espace_virtuel();
  int voisine;
  nb_faces_perio = 0;
  double diff1,diff2;

  // Boucle sur les bords pour traiter les conditions aux limites
  // il y a prise en compte d'un terme de convection pour les
  // conditions aux limites de Neumann_sortie_libre seulement
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);

      if (sub_type(Neumann_sortie_libre,la_cl.valeur()))
        {
          const Neumann_sortie_libre& la_sortie_libre = ref_cast(Neumann_sortie_libre, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (num_face=num1; num_face<num2; num_face++)
            {
              psc =0;
              for (i=0; i<dimension; i++)
                psc += vitesse_face_ALEmodif(num_face,i)*facenormales(num_face,i);
              if (psc>0)
                if (ncomp_ch_transporte == 1)
                  {
                    resu(num_face) += psc*transporte_face(num_face);
                  }
                else
                  for (i=0; i<ncomp_ch_transporte; i++)
                    {
                      resu(num_face,i) += psc*transporte_face(num_face,i);
                    }
              else
                {
                  if (ncomp_ch_transporte == 1)
                    {
                      resu(num_face) += psc*la_sortie_libre.val_ext(num_face-num1);
                    }
                  else
                    for (i=0; i<ncomp_ch_transporte; i++)
                      {
                        resu(num_face,i) += psc*la_sortie_libre.val_ext(num_face-num1,i);
                      }
                }
            }
        }
      else if (sub_type(Periodique,la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          ArrOfInt fait(le_bord.nb_faces());
          fait = 0;
          for (num_face=num1; num_face<num2; num_face++)
            {
              if (fait[num_face-num1] == 0)
                {
                  voisine = la_cl_perio.face_associee(num_face-num1) + num1;

                  if (ncomp_ch_transporte == 1)
                    {
                      diff1 = resu(num_face)-tab(nb_faces_perio);
                      diff2 = resu(voisine)-tab(nb_faces_perio+voisine-num_face);
                      resu(voisine)  -= diff1;
                      resu(num_face) -= diff2;
                    }
                  else
                    for (int comp=0; comp<ncomp_ch_transporte; comp++)
                      {
                        diff1 = resu(num_face,comp)-tab(nb_faces_perio,comp);
                        diff2 = resu(voisine,comp)-tab(nb_faces_perio+voisine-num_face,comp);
                        resu(voisine,comp)  -= diff1;
                        resu(num_face,comp) -= diff2;
                      }

                  fait[num_face-num1]= 1;
                  fait[voisine-num1] = 1;
                }
              nb_faces_perio++;
            }
        }
    }
}

void Op_Conv_ALE_VEF::ajouterALE_Muscl(const DoubleTab& inco,DoubleTab& resu) const
{

  //This is based on Op_Conv_VEF_Face::ajouter(...) from where only "amont" related part is used. ALE convection term has the opposite sign compared to the normal convection term.

  const Domaine_ALE& dom_ale=ref_cast(Domaine_ALE, dom.valeur());
  const Op_Conv_VEF_base& opConvVEFbase = ref_cast(Op_Conv_VEF_base, op_conv.valeur());
  const Op_Conv_VEF_Face& opConvVEFFace = ref_cast(Op_Conv_VEF_Face, op_conv.valeur());
  int ordre;
  opConvVEFFace.get_ordre(ordre);

  double alpha_;
  opConvVEFFace.get_alpha(alpha_);

  Motcle type_limit;
  if(op_conv.type()=="MUSCL") //In case of: "convection {  ALE { muscl } }" , type_lumit has default value set in Op_Conv_Muscl_VEF_Face.cpp .
    {
      type_limit="vanleer";
    }
  else // In case of: "convection {  ALE { generic muscl [limiter] [order of accuracy] [alpha] } }". op_conv.type()=="GENERIC" .
    {
      opConvVEFFace.get_type_lim(type_limit);
    }

  const Zone_Cl_VEF& zone_Cl_VEF = la_zcl_vef.valeur();
  const Zone_VEF& zone_VEF = ref_cast(Zone_VEF, la_zone_vef.valeur());

  const DoubleTab& vitesse_face_absolue_ALE=dom_ale.vitesse_faces();              //ALE mesh velocity.

  const DoubleVect& porosite_face = zone_VEF.porosite_face();

  int marq=opConvVEFbase.phi_u_transportant(equation());
  DoubleTab transporte_face_;
  DoubleTab vitesse_face_;
  // soit on a transporte_face=phi*transporte et vitesse_face=vitesse
  // soit on a transporte_face=transporte et vitesse_face=phi*vitesse
  // cela depend si on transporte avec phi*u ou avec u.
  const DoubleTab& transporte_face = modif_par_porosite_si_flag(inco,transporte_face_,!marq,porosite_face); //Fluid velocity = inco .
  const DoubleTab& vitesse_face_ALEmodif    = modif_par_porosite_si_flag(vitesse_face_absolue_ALE,vitesse_face_,marq,porosite_face);

  const IntTab& elem_faces = zone_VEF.elem_faces();
  const DoubleTab& facenormales = zone_VEF.face_normales();
  const DoubleTab& facette_normales = zone_VEF.facette_normales();
  const Zone& zone = zone_VEF.zone();
  const int nfa7 = zone_VEF.type_elem().nb_facette();
  const int nb_elem_tot = zone_VEF.nb_elem_tot();
  const IntVect& rang_elem_non_std = zone_VEF.rang_elem_non_std();
  const IntTab& face_voisins = zone_VEF.face_voisins();
  const DoubleTab& normales_facettes_Cl = zone_Cl_VEF.normales_facettes_Cl();
  int premiere_face_int = zone_VEF.premiere_face_int();
  int nfac = zone.nb_faces_elem();
  int nsom = zone.nb_som_elem();
  const IntTab& sommet_elem = zone.les_elems();
  const DoubleTab& vecteur_face_facette = ref_cast_non_const(Zone_VEF,zone_VEF).vecteur_face_facette();
  const  DoubleTab& vecteur_face_facette_Cl = zone_Cl_VEF.vecteur_face_facette_Cl();
  int nb_bord = zone_VEF.nb_front_Cl();
  const IntTab& les_elems=zone.les_elems();

  // Permet d'avoir un flux_bord coherent avec les CLs (mais parfois diverge?)
  // Active uniquement pour ordre 3
  int option_appliquer_cl_dirichlet = (ordre==3 ? 1 : 0);
  int option_calcul_flux_en_un_point = 0;//(ordre==3 ? 1 : 0);
  // Definition d'un tableau pour un traitement special des schemas pres des bords
  if (traitement_pres_bord_.size_array()!=nb_elem_tot)
    {
      traitement_pres_bord_.resize_array(nb_elem_tot);
      traitement_pres_bord_=0;
      // Pour muscl3 on applique le minmod sur les elements ayant une face de Dirichlet
      if (ordre==3)
        {
          for (int n_bord=0; n_bord<nb_bord; n_bord++)
            {
              const Cond_lim_base& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord).valeur();
              if ( sub_type(Dirichlet,la_cl) || sub_type(Dirichlet_homogene,la_cl) )
                {
                  const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
                  int nb_faces_tot = le_bord.nb_faces_tot();
                  for (int ind_face=0; ind_face<nb_faces_tot; ind_face++)
                    {
                      int num_face = le_bord.num_face(ind_face);
                      int elem = face_voisins(num_face,0);
                      traitement_pres_bord_(elem)=1;
                    }
                }
            }
        }
      else
        {
          // Pour le muscl/centre actuels on utilise un calcul de flux a l'ordre 1
          // aux mailles de bord ou aux mailles ayant un sommet de Dirichlet
          ArrOfInt est_un_sommet_de_bord_(zone_VEF.nb_som_tot());
          for (int n_bord=0; n_bord<nb_bord; n_bord++)
            {
              const Cond_lim_base& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord).valeur();
              if ( sub_type(Dirichlet,la_cl) || sub_type(Dirichlet_homogene,la_cl) )
                {
                  const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
                  int nb_faces_tot = le_bord.nb_faces_tot();
                  int size = zone_VEF.face_sommets().dimension(1);
                  for (int ind_face=0; ind_face<nb_faces_tot; ind_face++)
                    for (int som=0; som<size; som++)
                      {
                        int face = le_bord.num_face(ind_face);
                        est_un_sommet_de_bord_(zone_VEF.face_sommets(face,som))=1;
                      }
                }
            }
          for (int elem=0; elem<nb_elem_tot; elem++)
            {
              if (rang_elem_non_std(elem)!=-1)
                traitement_pres_bord_(elem)=1;
              else
                {
                  for (int n_som=0; n_som<nsom; n_som++)
                    if (est_un_sommet_de_bord_(les_elems(elem,n_som)))
                      traitement_pres_bord_(elem)=1;
                }
            }
        }
      // Construction du tableau est_une_face_de_dirichlet_
      est_une_face_de_dirichlet_.resize_array(zone_VEF.nb_faces_tot());
      est_une_face_de_dirichlet_=0;
      for (int n_bord=0; n_bord<nb_bord; n_bord++)
        {
          const Cond_lim_base& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord).valeur();
          if ( sub_type(Dirichlet,la_cl) || sub_type(Dirichlet_homogene,la_cl) )
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              int nb_faces_tot = le_bord.nb_faces_tot();
              for (int ind_face=0; ind_face<nb_faces_tot; ind_face++)
                {
                  int num_face = le_bord.num_face(ind_face);
                  est_une_face_de_dirichlet_(num_face) = 1;
                }
            }
        }
    }

  // Pour le traitement de la convection on distingue les polyedres
  // standard qui ne "voient" pas les conditions aux limites et les
  // polyedres non standard qui ont au moins une face sur le bord.
  // Un polyedre standard a n facettes sur lesquelles on applique le
  // schema de convection.
  // Pour un polyedre non standard qui porte des conditions aux limites
  // de Dirichlet, une partie des facettes sont portees par les faces.
  // En bref pour un polyedre le traitement de la convection depend
  // du type (triangle, tetraedre ...) et du nombre de faces de Dirichlet.

  const Elem_VEF_base& type_elemvef= zone_VEF.type_elem().valeur();
  int istetra=0;
  Nom nom_elem=type_elemvef.que_suis_je();
  if ((nom_elem=="Tetra_VEF")||(nom_elem=="Tri_VEF"))
    istetra=1;

  const DoubleVect& porosite_elem = zone_VEF.porosite_elem();
  double psc;
  int poly,face_adj,fa7,i,j,n_bord;
  int num_face, rang;
  int num10,num20;
  int ncomp_ch_transporte=(transporte_face.nb_dim() == 1?1:transporte_face.dimension(1));

  // Traitement particulier pour les faces de periodicite
  int nb_faces_perio = 0;
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          nb_faces_perio+=le_bord.nb_faces();
        }
    }

  DoubleTab tab;
  if (ncomp_ch_transporte == 1)
    tab.resize(nb_faces_perio);
  else
    tab.resize(nb_faces_perio,ncomp_ch_transporte);

  nb_faces_perio=0;
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (num_face=num1; num_face<num2; num_face++)
            {
              if (ncomp_ch_transporte == 1)
                tab(nb_faces_perio) = resu(num_face);
              else
                for (int comp=0; comp<ncomp_ch_transporte; comp++)
                  tab(nb_faces_perio,comp) = resu(num_face,comp);
              nb_faces_perio++;
            }
        }
    }

  int fac=0,elem1,elem2,comp0;
  int nb_faces_ = zone_VEF.nb_faces();
  ArrOfInt face(nfac);
  //statistiques().end_count(m1);
  //statistiques().begin_count(m2);

  // Tableau gradient base sur gradient_elem selon schema
  DoubleTab gradient_elem(nb_elem_tot,ncomp_ch_transporte,dimension);  // (du/dx du/dy dv/dx dv/dy) pour un poly

  Champ_P1NC::calcul_gradient(transporte_face,gradient_elem,zone_Cl_VEF);

  DoubleTab gradient;

  //  application du limiteur
  gradient.resize(0, ncomp_ch_transporte, dimension);     // (du/dx du/dy dv/dx dv/dy) pour une face
  zone_VEF.creer_tableau_faces(gradient);
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
      int num1 = le_bord.num_premiere_face();
      int num2 = num1 + le_bord.nb_faces();
      if (sub_type(Periodique,la_cl.valeur()))
        {
          for (fac=num1; fac<num2; fac++)
            {
              elem1=face_voisins(fac,0);
              elem2=face_voisins(fac,1);
              for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                for (i=0; i<dimension; i++)
                  {
                    double grad1=gradient_elem(elem1, comp0, i);
                    double grad2=gradient_elem(elem2, comp0, i);
                    //gradient(fac, comp0, i) =(*LIMITEUR)(grad1, grad2);
                    gradient(fac, comp0, i) =application_LIMITEUR(grad1, grad2, type_limit);
                  }
            }
        }
      else if (sub_type(Symetrie,la_cl.valeur()))
        {
          for (fac=num1; fac<num2; fac++)
            {
              elem1=face_voisins(fac,0);
              for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                for (i=0; i<dimension; i++)
                  gradient(fac, comp0, i) = gradient_elem(elem1, comp0, i);

              if (ordre==3)
                {
                  // On enleve la composante normale (on pourrait le faire pour les autres schemas...)
                  // mais pour le moment, on ne veut pas changer le comportement par defaut du muscl...
                  //const DoubleTab& facenormales = zone_VEF.face_normales();
                  for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                    for (i=0; i<dimension; i++)
                      {
                        double carre_surface=0;
                        double tmp=0;
                        for (j=0; j<dimension; j++)
                          {
                            double ndS=facenormales(fac,j);
                            carre_surface += ndS*ndS;
                            tmp += gradient(fac, comp0, j)*ndS;
                          }
                        gradient(fac, comp0, i) -= tmp*facenormales(fac,i)/carre_surface;
                      }
                }
            }
        }
    }

  for (fac=premiere_face_int; fac<nb_faces_; fac++)
    {
      elem1=face_voisins(fac,0);
      elem2=face_voisins(fac,1);
      int minmod_pres_du_bord = 0;
      if (ordre==3 && (traitement_pres_bord_(elem1) || traitement_pres_bord_(elem2))) minmod_pres_du_bord = 1;
      for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
        for (i=0; i<dimension; i++)
          {
            double grad1=gradient_elem(elem1, comp0, i);
            double grad2=gradient_elem(elem2, comp0, i);
            if (minmod_pres_du_bord)
              gradient(fac, comp0, i) = minmod(grad1, grad2);
            else
              //gradient(fac, comp0, i) = (*LIMITEUR)(grad1, grad2);
              gradient(fac, comp0, i) =application_LIMITEUR(grad1, grad2, type_limit);
          }
    } // fin du for faces
  gradient.echange_espace_virtuel();

  ArrOfDouble vs(dimension);
  ArrOfDouble vc(dimension);
  ArrOfDouble cc(dimension);
  DoubleVect xc(dimension);
  DoubleTab vsom(nsom,dimension);
  DoubleTab xsom(nsom,dimension);

  const IntTab& KEL=type_elemvef.KEL();
  const DoubleTab& xv=zone_VEF.xv();
  const Domaine& domaine=zone.domaine();
  const DoubleTab& coord_sommets=domaine.coord_sommets();

  // Boucle ou non selon la valeur de alpha (uniquement a l'ordre 3 pour le moment)
  // Si alpha=1, la boucle se limite a une simple passe avec le schema choisi (muscl, amont, centre)
  // Si alpha<1, la boucle se compose de 2 passes:
  //                         -la premiere avec le schema choisi et une ponderation de alpha
  //                         -la seconde avec le schema centre et une ponderation de 1-alpha
  double alpha = alpha_;
  int nombre_passes = (alpha==1 ? 1 : 2);
  for (int passe=1; passe<=nombre_passes; passe++)
    {
      if (passe==2)
        {
          gradient.ref(gradient_elem);
        }
      // Les polyedres non standard sont ranges en 2 groupes dans la Zone_VEF:
      //  - polyedres bords et joints
      //  - polyedres bords et non joints
      // On traite les polyedres en suivant l'ordre dans lequel ils figurent
      // dans la zone
      // boucle sur les polys
      for (poly=0; poly<nb_elem_tot; poly++)
        {
          int contrib = 0;
          // calcul des numeros des faces du polyedre
          for (face_adj=0; face_adj<nfac; face_adj++)
            {
              int face_ = elem_faces(poly,face_adj);
              face(face_adj)= face_;
              if (face_<nb_faces_) contrib=1; // Une face reelle sur l'element virtuel
            }
          //
          if (contrib)
            {
              int calcul_flux_en_un_point = (ordre != 3) && (ordre==1 || traitement_pres_bord_(poly));
              for (j=0; j<dimension; j++)
                {
                  vs(j) = vitesse_face_absolue_ALE(face(0),j)*porosite_face(face(0));
                  for (i=1; i<nfac; i++)
                    vs(j)+= vitesse_face_absolue_ALE(face(i),j)*porosite_face(face(i));
                }
              // calcul de la vitesse aux sommets des polyedres
              // On va utliser les fonctions de forme implementees dans la classe Champs_P1_impl ou Champs_Q1_impl
              if (istetra==1)
                {
                  for (i=0; i<nsom; i++)
                    for (j=0; j<dimension; j++)
                      vsom(i,j) = (vs(j) - dimension*vitesse_face_absolue_ALE(face(i),j)*porosite_face(face(i)));
                }
              else
                {
                  Cout << "ALE is currently working only with Tetra_VEF (3D) or Tri_VEF (2D)." << finl;
                  exit();
                }

              // Determination du type de CL selon le rang
              rang = rang_elem_non_std(poly);
              int itypcl = (rang==-1 ? 0 : zone_Cl_VEF.type_elem_Cl(rang));

              // calcul de vc (a l'intersection des 3 facettes) vc vs vsom proportionnelles a la porosite
              calcul_vc_ALE(face,vc,vs,vsom,vitesse_face_absolue_ALE,itypcl,porosite_face);

              // calcul de xc (a l'intersection des 3 facettes) necessaire pour muscl3
              if (ordre==3)
                {
                  int idirichlet;
                  int n1,n2,n3;
                  for (i=0; i<nsom; i++)
                    for (j=0; j<dimension; j++)
                      xsom(i,j) = coord_sommets(les_elems(poly,i),j);
                  type_elemvef.calcul_xg(xc,xsom,itypcl,idirichlet,n1,n2,n3);
                }

              // Gestion de la porosite
              if (marq==0)
                {
                  double coeff=1./porosite_elem(poly);
                  vsom*=coeff;
                  vc*=coeff;
                }
              // Boucle sur les facettes du polyedre non standard:
              for (fa7=0; fa7<nfa7; fa7++)
                {
                  num10 = face(KEL(0,fa7));
                  num20 = face(KEL(1,fa7));
                  // normales aux facettes
                  if (rang==-1)
                    for (i=0; i<dimension; i++)
                      cc[i] = facette_normales(poly, fa7, i);
                  else
                    for (i=0; i<dimension; i++)
                      cc[i] = normales_facettes_Cl(rang,fa7,i);

                  // Calcul des vitesses en C,S,S2 les 3 extremites de la fa7 et M le centre de la fa7
                  double psc_c=0,psc_s=0,psc_m,psc_s2=0;
                  if (dimension==2)
                    {
                      for (i=0; i<2; i++)
                        {
                          psc_c+=vc[i]*cc[i];
                          psc_s+=vsom(KEL(2,fa7),i)*cc[i];
                        }
                      psc_m=(psc_c+psc_s)/2.;
                    }
                  else
                    {
                      for (i=0; i<3; i++)
                        {
                          psc_c+=vc[i]*cc[i];
                          psc_s+=vsom(KEL(2,fa7),i)*cc[i];
                          psc_s2+=vsom(KEL(3,fa7),i)*cc[i];
                        }
                      psc_m=(psc_c+psc_s+psc_s2)/3.;
                    }
                  // On applique les CL de Dirichlet si num1 ou num2 est une face avec CL de Dirichlet
                  // auquel cas la fa7 coincide avec la face num1 ou num2 -> C est au centre de la face
                  int appliquer_cl_dirichlet=0;
                  if (option_appliquer_cl_dirichlet)
                    if (est_une_face_de_dirichlet_(num10) || est_une_face_de_dirichlet_(num20))
                      {
                        appliquer_cl_dirichlet = 1;
                        psc_m = psc_c;
                      }

                  // Determination de la face amont pour M
                  int face_amont_m,dir;
                  if (psc_m >= 0)
                    {
                      face_amont_m = num10;
                      dir=0;
                    }
                  else
                    {
                      face_amont_m = num20;
                      dir=1;
                    }
                  // Determination des faces amont pour les points C,S,S2
                  int face_amont_c=face_amont_m;
                  int face_amont_s=face_amont_m;
                  int face_amont_s2=face_amont_m;
                  if (passe==1 && ordre==3)
                    {
                      face_amont_c  = (psc_c >= 0)  ? num10 : num20;
                      face_amont_s  = (psc_s >= 0)  ? num10 : num20;
                      face_amont_s2 = (psc_s2 >= 0) ? num10 : num20;
                    }
                  // gradient aux items element (schema centre) ou aux items face (schemas muscl)
                  int item_m=poly;
                  int item_c=poly;
                  int item_s=poly;
                  int item_s2=poly;
                  if (passe==1)
                    {
                      item_m = face_amont_m;
                      item_c = face_amont_c;
                      item_s = face_amont_s;
                      item_s2 = face_amont_s2;
                    }

                  for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                    {
                      double flux;
                      double inco_m = (ncomp_ch_transporte==1?transporte_face(face_amont_m):transporte_face(face_amont_m,comp0));
                      if (appliquer_cl_dirichlet)
                        {
                          flux = inco_m*psc_m;
                        }

                      else // muscl ou centre
                        {
                          // Calcul de l'inconnue au centre M de la fa7
                          if (rang==-1)
                            for (j=0; j<dimension; j++)
                              inco_m+= gradient(item_m,comp0,j)*vecteur_face_facette(poly,fa7,j,dir);
                          else
                            for (j=0; j<dimension; j++)
                              inco_m+= gradient(item_m,comp0,j)*vecteur_face_facette_Cl(rang,fa7,j,dir);

                          // Calcul de l'inconnue au sommet S, une premiere extremite de la fa7
                          double inco_s = (ncomp_ch_transporte==1?transporte_face(face_amont_s):transporte_face(face_amont_s,comp0));
                          int sommet_s = sommet_elem(poly,KEL(2,fa7));
                          for (j=0; j<dimension; j++)
                            inco_s+= gradient(item_s,comp0,j)*(-xv(face_amont_s,j)+coord_sommets(sommet_s,j));

                          // Calcul de l'inconnue au sommet S2, la derniere extremite de la fa7 en 3D
                          double inco_s2=0;
                          if (dimension==3)
                            {
                              inco_s2 = (ncomp_ch_transporte==1?transporte_face(face_amont_s2):transporte_face(face_amont_s2,comp0));
                              int sommet_s2 = sommet_elem(poly,KEL(3,fa7));
                              for (j=0; j<dimension; j++)
                                inco_s2+= gradient(item_s2,comp0,j)*(-xv(face_amont_s2,j)+coord_sommets(sommet_s2,j));
                            }

                          // Calcul de l'inconnue a C, une autre extremite de la fa7, intersection avec les autres fa7
                          // du polyedre. C=G centre du polyedre si volume non etendu
                          // xc donne par elemvef.calcul_xg()
                          double inco_c;
                          if (ordre==3)
                            {
                              inco_c = (ncomp_ch_transporte==1?transporte_face(face_amont_c):transporte_face(face_amont_c,comp0));
                              for (j=0; j<dimension; j++)
                                inco_c+= gradient(item_c,comp0,j)*(-xv(face_amont_c,j)+xc(j));
                            }
                          else
                            {
                              inco_c = dimension*inco_m-inco_s-inco_s2;
                            }

                          // Calcul du flux sur 1 point
                          if (calcul_flux_en_un_point || option_calcul_flux_en_un_point)
                            {
                              flux = inco_m*psc_m;
                            }
                          else
                            {
                              // Calcul du flux sur 3 points
                              flux = (dimension==2) ? (inco_c*psc_c + inco_s*psc_s + 4*inco_m*psc_m)/6
                                     : (inco_c*psc_c + inco_s*psc_s + inco_s2*psc_s2 + 9*inco_m*psc_m)/12;
                            }

                        }
                      // Ponderation par coefficient alpha
                      flux*=alpha;

                      if (ncomp_ch_transporte == 1)
                        {
                          resu(num10) += flux;
                          resu(num20) -= flux;
                        }
                      else
                        {
                          resu(num10,comp0) += flux;
                          resu(num20,comp0) -= flux;
                        }

                    }// boucle sur comp
                } // fin de la boucle sur les facettes
            }
        } // fin de la boucle
      alpha = 1 - alpha;
    } // fin de la boucle
  //statistiques().end_count(m2);
  //statistiques().begin_count(m3);
  int voisine;
  nb_faces_perio = 0;
  double diff1,diff2;

  // Boucle sur les bords pour traiter les conditions aux limites
  // il y a prise en compte d'un terme de convection pour les
  // conditions aux limites de Neumann_sortie_libre seulement
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);

      if (sub_type(Neumann_sortie_libre,la_cl.valeur()))
        {
          const Neumann_sortie_libre& la_sortie_libre = ref_cast(Neumann_sortie_libre, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (num_face=num1; num_face<num2; num_face++)
            {
              psc =0;
              for (i=0; i<dimension; i++)
                psc += vitesse_face_ALEmodif(num_face,i)*facenormales(num_face,i);
              if (psc>0)
                if (ncomp_ch_transporte == 1)
                  {
                    resu(num_face) += psc*transporte_face(num_face);
                  }
                else
                  for (i=0; i<ncomp_ch_transporte; i++)
                    {
                      resu(num_face,i) += psc*transporte_face(num_face,i);
                    }
              else
                {
                  if (ncomp_ch_transporte == 1)
                    {
                      resu(num_face) += psc*la_sortie_libre.val_ext(num_face-num1);
                    }
                  else
                    for (i=0; i<ncomp_ch_transporte; i++)
                      {
                        resu(num_face,i) += psc*la_sortie_libre.val_ext(num_face-num1,i);
                      }
                }
            }
        }
      else if (sub_type(Periodique,la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          ArrOfInt fait(le_bord.nb_faces());
          fait = 0;
          for (num_face=num1; num_face<num2; num_face++)
            {
              if (fait[num_face-num1] == 0)
                {
                  voisine = la_cl_perio.face_associee(num_face-num1) + num1;

                  if (ncomp_ch_transporte == 1)
                    {
                      diff1 = resu(num_face)-tab(nb_faces_perio);
                      diff2 = resu(voisine)-tab(nb_faces_perio+voisine-num_face);
                      resu(voisine)  -= diff1;
                      resu(num_face) -= diff2;
                    }
                  else
                    for (int comp=0; comp<ncomp_ch_transporte; comp++)
                      {
                        diff1 = resu(num_face,comp)-tab(nb_faces_perio,comp);
                        diff2 = resu(voisine,comp)-tab(nb_faces_perio+voisine-num_face,comp);
                        resu(voisine,comp)  -= diff1;
                        resu(num_face,comp) -= diff2;
                      }

                  fait[num_face-num1]= 1;
                  fait[voisine-num1] = 1;
                }
              nb_faces_perio++;
            }
        }
    }

}


DoubleTab& Op_Conv_ALE_VEF::ajouterALE(const DoubleTab& inco, DoubleTab& resu) const
{
  Cerr << "Op_Conv_ALE_VEF::ajouterALE" << finl;

  const Domaine_ALE& dom_ale=ref_cast(Domaine_ALE, dom.valeur());
  const Op_Conv_VEF_base& opConvVEFbase = ref_cast(Op_Conv_VEF_base, op_conv.valeur());

  //If ALE mesh velocity is zero, then ajouterALE is directly zero (resu is initialized with zero values in Navier_Stokes_std). Otherwise it is calculated.
  if(dom_ale.vitesse_faces().mp_max_abs_vect() >=1.e-12)
    {
      //Check convection scheme discretization type:
      const Op_Conv_VEF_Face& opConvVEFFace = ref_cast(Op_Conv_VEF_Face, op_conv.valeur());
      int convectionSchemeDiscrType; // amont=0, muscl=1, centre=2
      opConvVEFFace.get_type_op(convectionSchemeDiscrType);

      if(convectionSchemeDiscrType==2) // Convection scheme discr "centre".
        {
          ajouterALE_Centre(inco, resu);
          opConvVEFbase.modifier_flux(*this); // Multiplication by density in case of incompressible Navier Stokes
        }// end centre

      else if (convectionSchemeDiscrType==0) // Convection scheme discr "amont".
        {
          ajouterALE_Amont(inco, resu);
          opConvVEFbase.modifier_flux(*this);

        }// end amont

      else if (convectionSchemeDiscrType==1) // Convection scheme discr "muscl".
        {

          ajouterALE_Muscl(inco, resu);
          opConvVEFbase.modifier_flux(*this);

        }// end muscl

      else
        {
          Cout << "ALE cannot use currently convection scheme " << op_conv.type() <<". Try available ALE convection schemes: centre , amont , muscl." << finl;
          exit();
        }

    }
  resu.echange_espace_virtuel();
  Debog::verifier("resu dans Op_Conv_ALE_VEF::ajouterALE", resu);
  return resu;
}


double Op_Conv_ALE_VEF::calculer_dt_stab() const
{
  const Zone_Cl_VEF& zone_Cl_VEF = la_zcl_vef.valeur();
  const Zone_VEF& zone_VEF = la_zone_vef.valeur();
  const DoubleVect& volumes_entrelaces = zone_VEF.volumes_entrelaces();
  const DoubleVect& volumes_entrelaces_Cl = zone_Cl_VEF.volumes_entrelaces_Cl();
  DoubleTrav fluent(volumes_entrelaces);

  remplir_fluent_ALEincluded(fluent);

  double dt_face,dt_stab =1.e30;

  // On traite les conditions aux limites
  // Si une face porte une condition de Dirichlet on n'en tient pas compte
  // dans le calcul de dt_stab
  for (int n_bord=0; n_bord<zone_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
      if ( (sub_type(Dirichlet,la_cl.valeur()))
           ||
           (sub_type(Dirichlet_homogene,la_cl.valeur()))
         )
        ;
      else
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();
          for (int num_face=ndeb; num_face<nfin; num_face++)
            {
              dt_face = volumes_entrelaces_Cl(num_face)/(fluent[num_face]+1.e-30);
              dt_stab = (dt_face < dt_stab) ? dt_face : dt_stab;
            }
        }
    }

  // On traite les faces internes non standard
  int ndeb = zone_VEF.premiere_face_int();
  int nfin = zone_VEF.premiere_face_std();

  for (int num_face=ndeb; num_face<nfin; num_face++)
    {
      dt_face = volumes_entrelaces_Cl(num_face)/(fluent[num_face]+DMINFLOAT);
      dt_stab =(dt_face < dt_stab) ? dt_face : dt_stab;
    }

  // On traite les faces internes standard
  ndeb = nfin;
  nfin = zone_VEF.nb_faces();
  for (int num_face=ndeb; num_face<nfin; num_face++)
    {
      dt_face = volumes_entrelaces(num_face)/(fluent[num_face]+DMINFLOAT);
      dt_stab =(dt_face < dt_stab) ? dt_face : dt_stab;
    }

  dt_stab = Process::mp_min(dt_stab);
  // astuce pour contourner le type const de la methode
  Op_Conv_ALE_VEF& op = ref_cast_non_const(Op_Conv_ALE_VEF,*this);
  op.fixer_dt_stab_conv(dt_stab);

  return dt_stab;
}

void Op_Conv_ALE_VEF::remplir_fluent_ALEincluded(DoubleVect& tab_fluent) const
{
// Taken from Op_Conv_EF_VEF_P1NC_Stab::remplir_fluent(DoubleVect& fluent_)
//
  const Zone_VEF& zone_VEF = la_zone_vef.valeur();
//old const Champ_Inc_base& la_vitesse=vitesse_.valeur();
//old const DoubleTab& tab_vitesse=la_vitesse.valeurs();
//new DoubleTab velocity= equation().inconnue().valeurs();
  const DoubleTab& tab_vitesse= equation().inconnue().valeurs();
  const DoubleTab& tab_vitesse_faces_ALE = ref_cast(Domaine_ALE, dom.valeur()).vitesse_faces();

  const DoubleTab& face_normales=zone_VEF.face_normales();
  const int nb_faces = zone_VEF.nb_faces();
  // calcul de la CFL.
  double psc;
  // On remet a zero le tableau qui sert pour
  // le calcul du pas de temps de stabilite


  //Initialization step before projection:
  if(tab_vitesse_faces_ALE.nb_dim()==1 && tab_vitesse_faces_ALE.mp_norme_vect()<=1.e-12)
    {
      for(int num_face=0; num_face<nb_faces; num_face++)
        {
          psc=0.;
          for (int i=0; i<dimension; i++)
            psc+=tab_vitesse(num_face,i)*face_normales(num_face,i);
          tab_fluent(num_face)=dabs(psc);
        }
    }

  //Projection step and solving problem:
  else
    {
      for(int num_face=0; num_face<nb_faces; num_face++)
        {
          psc=0.;
          for (int i=0; i<dimension; i++)
            psc+=(tab_vitesse(num_face,i)-tab_vitesse_faces_ALE(num_face,i))*face_normales(num_face,i);
          tab_fluent(num_face)=dabs(psc);
        }
    }

}

DoubleTab& Op_Conv_ALE_VEF::calculateALEMeshVelocityGradientOnFaces( DoubleTab& ALE_Mesh_Velocity_Gradient) const
{

  const DoubleTab& vitesse_faces_ALE = ref_cast(Domaine_ALE, dom.valeur()).vitesse_faces();
  int i,k,num_face;
  const Zone_Cl_VEF& zone_Cl_VEF = la_zcl_vef.valeur();
  const Zone_VEF& zone_VEF = la_zone_vef.valeur();
  const Zone_VF& zone_VF = ref_cast(Zone_VF, zone_VEF);
  const DoubleVect& volumes=zone_VF.volumes();
  const int nb_faces = zone_VEF.nb_faces();
  const int nb_elem_tot = zone_VEF.nb_elem_tot();
  DoubleTab gradient_elem;//(nb_elem,dimension,dimension);
  const IntTab& face_voisins = zone_VF.face_voisins();
  int premiere_face_int = zone_VEF.premiere_face_int();
  int ncomp_ch_transporte=dimension;
  Motcle le_nom_eqn=equation().le_nom();

  if (dimension == 1)
    {
      ncomp_ch_transporte=1;
      if (le_nom_eqn=="pbNavier_Stokes_standard")
        {
          ALE_Mesh_Velocity_Gradient.resize(nb_faces,dimension);
        }
    }
  gradient_elem.resize(nb_elem_tot,ncomp_ch_transporte,dimension);
  calculer_gradientP1NC(vitesse_faces_ALE,zone_VEF,zone_Cl_VEF ,gradient_elem);

  if (ncomp_ch_transporte != 1)
    {
      // Traitement des bords

      const Conds_lim& les_cl = zone_Cl_VEF.les_conditions_limites();
      int nb_cl=les_cl.size();
      for (int num_cl=0; num_cl<nb_cl; num_cl++)
        {
          const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(num_cl);
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int nb_faces_bord=le_bord.nb_faces();
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + nb_faces_bord;

          if (sub_type(Periodique,la_cl.valeur()))
            {
              const Periodique& la_cl_perio = ref_cast(Periodique,la_cl.valeur());
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

                      for (k=0; k<dimension; k++)
                        for (i=0; i<dimension; i++)
                          {
                            //Inverse distance weighting interpolation (p=3) is used for obtaining ALE mesh velocity gradient at the faces
                            ALE_Mesh_Velocity_Gradient(num_face,k,i)=(vol_1*gradient_elem(elem_0,k,i)+vol_0*gradient_elem(elem_1,k,i))/(vol_0+vol_1);
                          }
                      for (k=0; k<dimension; k++)
                        for (i=0; i<dimension; i++)
                          {
                            ALE_Mesh_Velocity_Gradient(num_face+num1,k,i)=ALE_Mesh_Velocity_Gradient(num_face,k,i);
                          }

                    }
                } // Fin de la boucle sur les faces
            } // Fin du cas periodique
          else
            {
              for (num_face=num1; num_face<num2; num_face++)
                {
                  int elem_0 = face_voisins(num_face,0);
                  for (k=0; k<dimension; k++)
                    for (i=0; i<dimension; i++)
                      {
                        ALE_Mesh_Velocity_Gradient(num_face,k,i)=gradient_elem(elem_0,k,i);
                      }
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

          for (k=0; k<dimension; k++)
            for (i=0; i<dimension; i++)
              {
                //Inverse distance weighting interpolation (p=3) is used for obtaining ALE mesh velocity gradient at the faces
                ALE_Mesh_Velocity_Gradient(num_face,k,i)=(vol_1*gradient_elem(elem_0,k,i)+vol_0*gradient_elem(elem_1,k,i))/(vol_0+vol_1);
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
          int nb_faces_bord=le_bord.nb_faces();
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + nb_faces_bord;
          if (sub_type(Periodique,la_cl.valeur()))
            {
              const Periodique& la_cl_perio = ref_cast(Periodique,la_cl.valeur());
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

                      for (i=0; i<dimension; i++)
                        {
                          //Inverse distance weighting interpolation (p=3) is used for obtaining ALE mesh velocity gradient at the faces
                          ALE_Mesh_Velocity_Gradient(num_face,i)=(vol_1*gradient_elem(elem_0,0,i) + vol_0*gradient_elem(elem_1,0,i))/(vol_0+vol_1);
                          ALE_Mesh_Velocity_Gradient(num_face+num1,i)=ALE_Mesh_Velocity_Gradient(num_face,i);
                        }

                    }
                } // Fin de la boucle sur les faces
            } // Fin du cas periodique

          else
            {
              for (num_face=num1; num_face<num2; num_face++)
                {
                  int elem_0 = face_voisins(num_face,0);

                  for (i=0; i<dimension; i++)
                    ALE_Mesh_Velocity_Gradient(num_face,i)=gradient_elem(elem_0,0,i);
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
          for(i=0; i<dimension; i++)
            {
              //Inverse distance weighting interpolation (p=3) is used for obtaining ALE mesh velocity gradient at the faces
              ALE_Mesh_Velocity_Gradient(num_face,i)=(vol_1*gradient_elem(elem_0,0,i) + vol_0*gradient_elem(elem_1,0,i))/(vol_0+vol_1);
            }

        }
    }

  return ALE_Mesh_Velocity_Gradient;
}

DoubleTab& Op_Conv_ALE_VEF::calculateALEjacobian(DoubleTab& jacobianALE) const
{
//////////////////////////////////////////////////////////////////////////////////
//Jacobian determinant in the case when the mesh is moved as x_new=x_old+dt*v_ALE
//////////////////////////////////////////////////////////////////////////////////
//Case with periodic boundary conditions are not implemented yet!

  int num_face;
  double timestep=equation().probleme().schema_temps().pas_de_temps();
  const Zone_VEF& zone_VEF = la_zone_vef.valeur();
  const int nb_faces_tot = zone_VEF.nb_faces_tot();
  DoubleTab ALEmeshVelocityGradient(nb_faces_tot,dimension,dimension);
  jacobianALE.resize(nb_faces_tot,dimension);
  calculateALEMeshVelocityGradientOnFaces(ALEmeshVelocityGradient);

  //Multiply ALE mesh velocity gradient values with time step and used the result to calculate Jacobian determinants on faces.
  ALEmeshVelocityGradient*=timestep;

  if (dimension == 2)
    {
      for (num_face=0; num_face<nb_faces_tot; num_face++)
        {
          jacobianALE(num_face,0)=(1+ALEmeshVelocityGradient(num_face,0,0))*(1+ALEmeshVelocityGradient(num_face,1,1))-ALEmeshVelocityGradient(num_face,0,1)*ALEmeshVelocityGradient(num_face,1,0);
          jacobianALE(num_face,1)=jacobianALE(num_face,0);
        }
    }
  else if (dimension == 3)
    {
      for (num_face=0; num_face<nb_faces_tot; num_face++)
        {
          jacobianALE(num_face,0)=(1+ALEmeshVelocityGradient(num_face,0,0))*(1+ALEmeshVelocityGradient(num_face,1,1))*(1+ALEmeshVelocityGradient(num_face,2,2))+ALEmeshVelocityGradient(num_face,0,1)*ALEmeshVelocityGradient(num_face,1,2)*ALEmeshVelocityGradient(num_face,2,0)+ALEmeshVelocityGradient(num_face,0,2)*ALEmeshVelocityGradient(num_face,1,0)*ALEmeshVelocityGradient(num_face,2,1)-ALEmeshVelocityGradient(num_face,0,2)*(1+ALEmeshVelocityGradient(num_face,1,1))*ALEmeshVelocityGradient(num_face,2,0)-ALEmeshVelocityGradient(num_face,0,1)*ALEmeshVelocityGradient(num_face,1,0)*(1*ALEmeshVelocityGradient(num_face,2,2))-(1+ALEmeshVelocityGradient(num_face,0,0))*ALEmeshVelocityGradient(num_face,1,2)*ALEmeshVelocityGradient(num_face,2,1);
          jacobianALE(num_face,1)=jacobianALE(num_face,0);
          jacobianALE(num_face,2)=jacobianALE(num_face,0);
        }
    }
  else
    {
      for (num_face=0; num_face<nb_faces_tot; num_face++)
        {
          jacobianALE(num_face,0)=(1+ALEmeshVelocityGradient(num_face,0,0));
        }
    }

  return jacobianALE;
}

void Op_Conv_ALE_VEF::calcul_vc_ALE(const ArrOfInt& Face, ArrOfDouble& vc, const ArrOfDouble& vs, const DoubleTab& vsom, const DoubleTab& vitesse_face_absolue_ALE,int type_cl, const DoubleVect& porosite_face) const
{
// Op_Conv_ALE_VEF::calcul_vc_ALE(...) is based on Tetra_VEF::calcul_vc(...) and on Tri_VEF::calcul_vc(...). It was created in order to have vitesse_face_absolue_ALE input as DoubleTab& compared to the original Champ_Inc_base& .
  if (dimension == 2) //Tri_VEF
    {

      switch(type_cl)
        {
        case 0: // le triangle n'a pas de Face de Dirichlet
          {
            vc[0] = vs[0]/3;
            vc[1] = vs[1]/3;
            break;
          }

        case 1: // le triangle a une Face de Dirichlet :la Face 2
          {
            vc[0]= vitesse_face_absolue_ALE(Face[2],0)*porosite_face(Face[2]);
            vc[1]= vitesse_face_absolue_ALE(Face[2],1)*porosite_face(Face[2]);
            break;
          }

        case 2: // le triangle a une Face de Dirichlet :la Face 1
          {
            vc[0]= vitesse_face_absolue_ALE(Face[1],0)*porosite_face(Face[1]);
            vc[1]= vitesse_face_absolue_ALE(Face[1],1)*porosite_face(Face[1]);;
            break;
          }

        case 4: // le triangle a une Face de Dirichlet :la Face 0
          {
            vc[0]= vitesse_face_absolue_ALE(Face[0],0)*porosite_face(Face[0]);
            vc[1]= vitesse_face_absolue_ALE(Face[0],1)*porosite_face(Face[0]);
            break;
          }

        case 3: // le triangle a deux faces de Dirichlet :les faces 1 et 2
          {
            vc[0]= vsom(0,0);
            vc[1]= vsom(0,1);
            break;
          }

        case 5: // le triangle a deux faces de Dirichlet :les faces 0 et 2
          {
            vc[0]= vsom(1,0);
            vc[1]= vsom(1,1);
            break;
          }

        case 6: // le triangle a deux faces de Dirichlet :les faces 0 et 1
          {
            vc[0]= vsom(2,0);
            vc[1]= vsom(2,1);
            break;
          }

        } // switch end

    } // 2D end

  else if (dimension == 3) //Tetra_VEF
    {

      int comp;
      switch(type_cl)
        {

        case 0: // le tetraedre n'a pas de Face de Dirichlet
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = 0.25*vs[comp];
            break;
          }

        case 1: // le tetraedre a une Face de Dirichlet : KEL3
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = vitesse_face_absolue_ALE(Face[3],comp)*porosite_face(Face[3]);
            break;
          }

        case 2: // le tetraedre a une Face de Dirichlet : KEL2
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = vitesse_face_absolue_ALE(Face[2],comp)*porosite_face(Face[2]);
            break;
          }

        case 4: // le tetraedre a une Face de Dirichlet : KEL1
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = vitesse_face_absolue_ALE(Face[1],comp)*porosite_face(Face[1]);
            break;
          }

        case 8: // le tetraedre a une Face de Dirichlet : KEL0
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = vitesse_face_absolue_ALE(Face[0],comp)*porosite_face(Face[0]);
            break;
          }

        case 3: // le tetraedre a deux faces de Dirichlet : KEL3 et KEL2
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = 0.5* (vsom(0,comp) + vsom(1,comp));
            break;
          }

        case 5: // le tetraedre a deux faces de Dirichlet : KEL3 et KEL1
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = 0.5* (vsom(0,comp) + vsom(2,comp));
            break;
          }

        case 6: // le tetraedre a deux faces de Dirichlet : KEL1 et KEL2
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = 0.5* (vsom(0,comp) + vsom(3,comp));
            break;
          }

        case 9: // le tetraedre a deux faces de Dirichlet : KEL0 et KEL3
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = 0.5* (vsom(1,comp) + vsom(2,comp));
            break;
          }

        case 10: // le tetraedre a deux faces de Dirichlet : KEL0 et KEL2
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = 0.5* (vsom(1,comp) + vsom(3,comp));
            break;
          }

        case 12: // le tetraedre a deux faces de Dirichlet : KEL0 et KEL1
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = 0.5*(vsom(2,comp) + vsom(3,comp));
            break;
          }

        case 7: // le tetraedre a trois faces de Dirichlet : KEL1, KEL2 et KEL3
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = vsom(0,comp);
            break;
          }

        case 11: // le tetraedre a trois faces de Dirichlet : KEL0,KEL2 et KEL3
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = vsom(1,comp);
            break;
          }

        case 13: // le tetraedre a trois faces de Dirichlet : KEL0, KEL1 et KEL3
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = vsom(2,comp);
            break;
          }

        case 14: // le tetraedre a trois faces de Dirichlet : KEL0, KEL1 et KEL2
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = vsom(3,comp);
            break;
          }

        } // switch end

    } // 3D end
}


double Op_Conv_ALE_VEF::application_LIMITEUR(double grad1, double grad2, Motcle& type_limit) const
{
  //////////////////////////////////////////////////////////////
  //   Fonctions limiteurs de MUSCL
  ////////////////////////////////////////////////////////////////
  // Inlining inutile car ces fonctions sont appelees dynamiquement
  // Pour optimiser, passer par une macro ?
  if ( !(type_limit=="minmod") && !(type_limit=="vanleer") && !(type_limit=="vanalbada")
       &&  !(type_limit=="chakravarthy") && !(type_limit=="superbee") )
    {
      Cerr << type_limit << " is not implemented. " << finl;
      Cerr << " Choose from: minmod - vanleer - vanalbada - chakravarthy - superbee " << finl;
      exit();
    }
  double gradlim=0.;
  if (type_limit=="minmod")
    {
      if(grad1*grad2>0.) (dabs(grad1)<dabs(grad2)) ? gradlim=grad1 : gradlim=grad2 ;
    }
  if (type_limit=="vanleer")
    {
      if(grad1*grad2>0.) gradlim=2.*grad1*grad2/(grad1+grad2) ;
    }
  if (type_limit=="vanalbada")
    {
      if(grad1*grad2>0.) gradlim=grad1*grad2*(grad1+grad2)/(grad1*grad1+grad2*grad2) ;
    }
  if (type_limit=="chakravarthy")
    {
      /*
        Cerr << " limiteur chakavarthy non preconise (non symetrique) " << finl;
        exit();
        return 0;
      */
      if ((grad1*grad2)>0)
        {
          gradlim=dmin(grad1/grad2,1.8); // 1<<beta<<2
          gradlim=dmax(gradlim,0.);
          gradlim*=grad2;
        }
    }
  if (type_limit=="superbee")
    {
      /*
        Cerr << " limiteur superbee non preconise (source d'instabilites) " << finl;
        exit();
        return 0;
      */
      if ((grad1*grad2)>0)
        {
          double gradlim1,gradlim2;
          gradlim1=dmin(2*(grad1/grad2),1);
          gradlim2=dmin(grad1/grad2,2);
          gradlim=dmax(gradlim1,gradlim2);
          gradlim=dmax(gradlim,0.);
          gradlim*=grad2;
        }
    }
  return gradlim;
}


DoubleTab& Op_Conv_ALE_VEF::supprimerALE(const DoubleTab& inco, DoubleTab& resu) const
{
  const Op_Conv_VEF_base& opConvVEFbase = ref_cast(Op_Conv_VEF_base, op_conv.valeur());
  //Check convection scheme discretization type:
  const Op_Conv_VEF_Face& opConvVEFFace = ref_cast(Op_Conv_VEF_Face, op_conv.valeur());
  int convectionSchemeDiscrType; // amont=0, muscl=1, centre=2
  opConvVEFFace.get_type_op(convectionSchemeDiscrType);

  if(convectionSchemeDiscrType==2) // Convection scheme discr "centre".
    {

      supprimerALE_Centre(inco, resu);
      opConvVEFbase.modifier_flux(*this);

    }// end centre

  else if (convectionSchemeDiscrType==0) // Convection scheme discr "amont".
    {
      supprimerALE_Amont(inco, resu);
      opConvVEFbase.modifier_flux(*this);

    }// end amont

  else if (convectionSchemeDiscrType==1) // Convection scheme discr "muscl".
    {
      supprimerALE_Muscl(inco, resu);
      opConvVEFbase.modifier_flux(*this);

    }// end muscl


  else
    {
      Cout << "ALE cannot use currently convection scheme discretization type " << op_conv.type() <<". Try available ALE convection scheme discretization types : centre , amont , muscl." << finl;
      exit();
    }

  return resu;
}

void Op_Conv_ALE_VEF::supprimerALE_Centre(const DoubleTab& inco, DoubleTab& resu) const
{

  Cerr << "Op_Conv_ALE_VEF::supprimerALE" << finl;
  const Domaine_ALE& dom_ale=ref_cast(Domaine_ALE, dom.valeur());
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
  const DoubleTab& vitesse_faces_ALE = dom_ale.vitesse_faces();
  const IntTab& face_voisins = zone_VF.face_voisins();
  const int& nb_bord=zone.faces_bord().nb_bords();
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
                      c_grad_u_0(k)+=vitesse_faces_ALE(num_face,i)*gradient_elem(elem_non_std,k,i);
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
                  c_grad_u_0(k)+=vitesse_faces_ALE(num_face,i)*gradient_elem(elem_0,k,i);
                  c_grad_u_1(k)+=vitesse_faces_ALE(num_face,i)*gradient_elem(elem_1,k,i);
                }
            }
          for (k=0; k<dimension; k++)
            {
              resu(num_face,k)-=(volumes_entrelaces_Cl(num_face)/(vol_0+vol_1))*(vol_1*c_grad_u_0(k) + vol_0*c_grad_u_1(k));
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
                  c_grad_u_0(k)+=vitesse_faces_ALE(num_face,i)*gradient_elem(elem_0,k,i);
                  c_grad_u_1(k)+=vitesse_faces_ALE(num_face,i)*gradient_elem(elem_1,k,i);
                }
            }
          for (k=0; k<dimension; k++)
            {
              resu(num_face,k)-=(volumes_entrelaces(num_face)/(vol_0+vol_1))*(vol_1*c_grad_u_0(k) + vol_0*c_grad_u_1(k));
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
                  c_grad_scalaire_0+=vitesse_faces_ALE(num_face,i)*gradient_elem(elem_non_std,0,i);
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
              c_grad_scalaire_0+=vitesse_faces_ALE(num_face,i)*gradient_elem(elem_0,0,i);
              c_grad_scalaire_1+=vitesse_faces_ALE(num_face,i)*gradient_elem(elem_1,0,i);
            }
          resu(num_face)-=(volumes_entrelaces_Cl(num_face)/(vol_0+vol_1))*(vol_1*c_grad_u_0(k) + vol_0*c_grad_u_1(k));
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
              c_grad_scalaire_0+=vitesse_faces_ALE(num_face,i)*gradient_elem(elem_0,0,i);
              c_grad_scalaire_1+=vitesse_faces_ALE(num_face,i)*gradient_elem(elem_1,0,i);
            }
          resu(num_face)-=(volumes_entrelaces(num_face)/(vol_0+vol_1))*(vol_1*c_grad_scalaire_0 + vol_0*c_grad_scalaire_1);
        }
    }


}
void Op_Conv_ALE_VEF::supprimerALE_Amont(const DoubleTab& inco, DoubleTab& resu) const
{

  const Domaine_ALE& dom_ale=ref_cast(Domaine_ALE, dom.valeur());
  const Zone_Cl_VEF& zone_Cl_VEF = la_zcl_vef.valeur();
  const Op_Conv_VEF_base& opConvVEFbase = ref_cast(Op_Conv_VEF_base, op_conv.valeur());

  const Zone_VEF& zone_VEF = ref_cast(Zone_VEF, la_zone_vef.valeur());

  const DoubleTab& vitesse_face_absolue_ALE=dom_ale.vitesse_faces();              //ALE mesh velocity.

  const DoubleVect& porosite_face = zone_VEF.porosite_face();
  int marq=opConvVEFbase.phi_u_transportant(equation());

  DoubleTab transporte_face_;
  DoubleTab vitesse_face_;

  const DoubleTab& transporte_face = modif_par_porosite_si_flag(inco,transporte_face_,!marq,porosite_face); //Fluid velocity = inco .
  const DoubleTab& vitesse_face_ALEmodif    = modif_par_porosite_si_flag(vitesse_face_absolue_ALE,vitesse_face_,marq,porosite_face);

  const IntTab& elem_faces = zone_VEF.elem_faces();
  const DoubleTab& facenormales = zone_VEF.face_normales();
  const DoubleTab& facette_normales = zone_VEF.facette_normales();
  const Zone& zone = zone_VEF.zone();
  const int nfa7 = zone_VEF.type_elem().nb_facette();
  const int nb_elem_tot = zone_VEF.nb_elem_tot();
  const IntVect& rang_elem_non_std = zone_VEF.rang_elem_non_std();
  const DoubleTab& normales_facettes_Cl = zone_Cl_VEF.normales_facettes_Cl();
  int nfac = zone.nb_faces_elem();
  int nsom = zone.nb_som_elem();
  int nb_bord = zone_VEF.nb_front_Cl();
  const IntTab& les_elems=zone.les_elems();

  // Permet d'avoir un flux_bord coherent avec les CLs (mais parfois diverge?)
  // Active uniquement pour ordre 3
  int option_appliquer_cl_dirichlet = 0;
  // Definition d'un tableau pour un traitement special des schemas pres des bords

  if (traitement_pres_bord_.size_array()!=nb_elem_tot)
    {
      traitement_pres_bord_.resize_array(nb_elem_tot);
      traitement_pres_bord_=0;


      // Pour le muscl/centre actuels on utilise un calcul de flux a l'ordre 1
      // aux mailles de bord ou aux mailles ayant un sommet de Dirichlet
      ArrOfInt est_un_sommet_de_bord_(zone_VEF.nb_som_tot());
      for (int n_bord=0; n_bord<nb_bord; n_bord++)
        {
          const Cond_lim_base& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord).valeur();
          if ( sub_type(Dirichlet,la_cl) || sub_type(Dirichlet_homogene,la_cl) )
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              int nb_faces_tot = le_bord.nb_faces_tot();
              int size = zone_VEF.face_sommets().dimension(1);
              for (int ind_face=0; ind_face<nb_faces_tot; ind_face++)
                for (int som=0; som<size; som++)
                  {
                    int face = le_bord.num_face(ind_face);
                    est_un_sommet_de_bord_(zone_VEF.face_sommets(face,som))=1;
                  }
            }
        }
      for (int elem=0; elem<nb_elem_tot; elem++)
        {
          if (rang_elem_non_std(elem)!=-1)
            traitement_pres_bord_(elem)=1;
          else
            {
              for (int n_som=0; n_som<nsom; n_som++)
                if (est_un_sommet_de_bord_(les_elems(elem,n_som)))
                  traitement_pres_bord_(elem)=1;
            }
        }

      // Construction du tableau est_une_face_de_dirichlet_
      est_une_face_de_dirichlet_.resize_array(zone_VEF.nb_faces_tot());
      est_une_face_de_dirichlet_=0;
      for (int n_bord=0; n_bord<nb_bord; n_bord++)
        {
          const Cond_lim_base& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord).valeur();
          if ( sub_type(Dirichlet,la_cl) || sub_type(Dirichlet_homogene,la_cl) )
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              int nb_faces_tot = le_bord.nb_faces_tot();
              for (int ind_face=0; ind_face<nb_faces_tot; ind_face++)
                {
                  int num_face = le_bord.num_face(ind_face);
                  est_une_face_de_dirichlet_(num_face) = 1;
                }
            }
        }
    }

  // Pour le traitement de la convection on distingue les polyedres
  // standard qui ne "voient" pas les conditions aux limites et les
  // polyedres non standard qui ont au moins une face sur le bord.
  // Un polyedre standard a n facettes sur lesquelles on applique le
  // schema de convection.
  // Pour un polyedre non standard qui porte des conditions aux limites
  // de Dirichlet, une partie des facettes sont portees par les faces.
  // En bref pour un polyedre le traitement de la convection depend
  // du type (triangle, tetraedre ...) et du nombre de faces de Dirichlet.

  const Elem_VEF_base& type_elemvef= zone_VEF.type_elem().valeur();
  int istetra=0;
  Nom nom_elem=type_elemvef.que_suis_je();
  if ((nom_elem=="Tetra_VEF")||(nom_elem=="Tri_VEF"))
    istetra=1;

  const DoubleVect& porosite_elem = zone_VEF.porosite_elem();
  double psc;
  int poly,face_adj,fa7,i,j,n_bord;
  int num_face, rang;
  int num10,num20;
  int ncomp_ch_transporte=(transporte_face.nb_dim() == 1?1:transporte_face.dimension(1));

  // Traitement particulier pour les faces de periodicite
  int nb_faces_perio = 0;
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          nb_faces_perio+=le_bord.nb_faces();
        }
    }

  DoubleTab tab;
  if (ncomp_ch_transporte == 1)
    tab.resize(nb_faces_perio);
  else
    tab.resize(nb_faces_perio,ncomp_ch_transporte);

  nb_faces_perio=0;
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (num_face=num1; num_face<num2; num_face++)
            {
              if (ncomp_ch_transporte == 1)
                tab(nb_faces_perio) = resu(num_face);
              else
                for (int comp=0; comp<ncomp_ch_transporte; comp++)
                  tab(nb_faces_perio,comp) = resu(num_face,comp);
              nb_faces_perio++;
            }
        }
    }

  int comp0;
  int nb_faces_ = zone_VEF.nb_faces();
  ArrOfInt face(nfac);

  // Tableau gradient base sur gradient_elem selon schema
  DoubleTab gradient_elem(nb_elem_tot,ncomp_ch_transporte,dimension);  // (du/dx du/dy dv/dx dv/dy) pour un poly
  DoubleTab gradient;

  ArrOfDouble vs(dimension);
  ArrOfDouble vc(dimension);
  ArrOfDouble cc(dimension);
  DoubleVect xc(dimension);
  DoubleTab vsom(nsom,dimension);
  DoubleTab xsom(nsom,dimension);

  const IntTab& KEL=type_elemvef.KEL();

  // Boucle ou non selon la valeur de alpha (uniquement a l'ordre 3 pour le moment)
  // Si alpha=1, la boucle se limite a une simple passe avec le schema choisi (muscl, amont, centre)
  // Si alpha<1, la boucle se compose de 2 passes:
  //                         -la premiere avec le schema choisi et une ponderation de alpha
  //                         -la seconde avec le schema centre et une ponderation de 1-alpha
  double alpha = 1;

  // Les polyedres non standard sont ranges en 2 groupes dans la Zone_VEF:
  //  - polyedres bords et joints
  //  - polyedres bords et non joints
  // On traite les polyedres en suivant l'ordre dans lequel ils figurent
  // dans la zone
  // boucle sur les polys
  for (poly=0; poly<nb_elem_tot; poly++)
    {
      int contrib = 0;
      // calcul des numeros des faces du polyedre
      for (face_adj=0; face_adj<nfac; face_adj++)
        {
          int face_ = elem_faces(poly,face_adj);
          face(face_adj)= face_;
          if (face_<nb_faces_) contrib=1; // Une face reelle sur l'element virtuel
        }
      //
      if (contrib)
        {
          for (j=0; j<dimension; j++)
            {
              vs(j) = vitesse_face_absolue_ALE(face(0),j)*porosite_face(face(0));
              for (i=1; i<nfac; i++)
                vs(j)+= vitesse_face_absolue_ALE(face(i),j)*porosite_face(face(i));
            }
          // calcul de la vitesse aux sommets des polyedres
          // On va utliser les fonctions de forme implementees dans la classe Champs_P1_impl ou Champs_Q1_impl
          if (istetra==1)
            {
              for (i=0; i<nsom; i++)
                for (j=0; j<dimension; j++)
                  vsom(i,j) = (vs(j) - dimension*vitesse_face_absolue_ALE(face(i),j)*porosite_face(face(i)));
            }
          else
            {
              Cout << "ALE is currently working only with Tetra_VEF (3D) or Tri_VEF (2D)." << finl;
              exit();
            }

          // Determination du type de CL selon le rang
          rang = rang_elem_non_std(poly);
          int itypcl = (rang==-1 ? 0 : zone_Cl_VEF.type_elem_Cl(rang));

          // calcul de vc (a l'intersection des 3 facettes) vc vs vsom proportionnelles a la porosite
          calcul_vc_ALE(face,vc,vs,vsom,vitesse_face_absolue_ALE,itypcl,porosite_face);

          // Gestion de la porosite
          if (marq==0)
            {
              double coeff=1./porosite_elem(poly);
              vsom*=coeff;
              vc*=coeff;
            }
          // Boucle sur les facettes du polyedre non standard:
          for (fa7=0; fa7<nfa7; fa7++)
            {
              num10 = face(KEL(0,fa7));
              num20 = face(KEL(1,fa7));

              // normales aux facettes
              if (rang==-1)
                for (i=0; i<dimension; i++)
                  cc[i] = facette_normales(poly, fa7, i);
              else
                for (i=0; i<dimension; i++)
                  cc[i] = normales_facettes_Cl(rang,fa7,i);

              // Calcul des vitesses en C,S,S2 les 3 extremites de la fa7 et M le centre de la fa7
              double psc_c=0,psc_s=0,psc_m,psc_s2=0;
              if (dimension==2)
                {
                  for (i=0; i<2; i++)
                    {
                      psc_c+=vc[i]*cc[i];
                      psc_s+=vsom(KEL(2,fa7),i)*cc[i];
                    }
                  psc_m=(psc_c+psc_s)/2.;
                }
              else
                {
                  for (i=0; i<3; i++)
                    {
                      psc_c+=vc[i]*cc[i];
                      psc_s+=vsom(KEL(2,fa7),i)*cc[i];
                      psc_s2+=vsom(KEL(3,fa7),i)*cc[i];
                    }
                  psc_m=(psc_c+psc_s+psc_s2)/3.;
                }
              // On applique les CL de Dirichlet si num1 ou num2 est une face avec CL de Dirichlet
              // auquel cas la fa7 coincide avec la face num1 ou num2 -> C est au centre de la face
              int appliquer_cl_dirichlet=0;
              if (option_appliquer_cl_dirichlet)
                if (est_une_face_de_dirichlet_(num10) || est_une_face_de_dirichlet_(num20))
                  {
                    appliquer_cl_dirichlet = 1;
                    psc_m = psc_c;
                  }

              // Determination de la face amont pour M
              int face_amont_m;
              if (psc_m >= 0)
                {
                  face_amont_m = num10;
                }
              else
                {
                  face_amont_m = num20;
                }

              for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                {
                  double flux;
                  double inco_m = (ncomp_ch_transporte==1?transporte_face(face_amont_m):transporte_face(face_amont_m,comp0));

                  if (true || appliquer_cl_dirichlet)
                    {
                      flux = inco_m*psc_m;
                    }

                  // Ponderation par coefficient alpha
                  flux*=alpha;

                  if (ncomp_ch_transporte == 1)
                    {
                      resu(num10) -= flux;
                      resu(num20) += flux;
                    }
                  else
                    {
                      resu(num10,comp0) -= flux;
                      resu(num20,comp0) += flux;
                    }

                }// boucle sur comp
            } // fin de la boucle sur les facettes
        }
    } // fin de la boucle
  alpha = 1 - alpha;

  int voisine;
  nb_faces_perio = 0;
  double diff1,diff2;

  // Boucle sur les bords pour traiter les conditions aux limites
  // il y a prise en compte d'un terme de convection pour les
  // conditions aux limites de Neumann_sortie_libre seulement
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);

      if (sub_type(Neumann_sortie_libre,la_cl.valeur()))
        {
          const Neumann_sortie_libre& la_sortie_libre = ref_cast(Neumann_sortie_libre, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (num_face=num1; num_face<num2; num_face++)
            {
              psc =0;
              for (i=0; i<dimension; i++)
                psc += vitesse_face_ALEmodif(num_face,i)*facenormales(num_face,i);
              if (psc>0)
                if (ncomp_ch_transporte == 1)
                  {
                    resu(num_face) -= psc*transporte_face(num_face);
                  }
                else
                  for (i=0; i<ncomp_ch_transporte; i++)
                    {
                      resu(num_face,i) -= psc*transporte_face(num_face,i);
                    }
              else
                {
                  if (ncomp_ch_transporte == 1)
                    {
                      resu(num_face) -= psc*la_sortie_libre.val_ext(num_face-num1);
                    }
                  else
                    for (i=0; i<ncomp_ch_transporte; i++)
                      {
                        resu(num_face,i) -= psc*la_sortie_libre.val_ext(num_face-num1,i);
                      }
                }
            }
        }
      else if (sub_type(Periodique,la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          ArrOfInt fait(le_bord.nb_faces());
          fait = 0;
          for (num_face=num1; num_face<num2; num_face++)
            {
              if (fait[num_face-num1] == 0)
                {
                  voisine = la_cl_perio.face_associee(num_face-num1) + num1;

                  if (ncomp_ch_transporte == 1)
                    {
                      diff1 = resu(num_face)-tab(nb_faces_perio);
                      diff2 = resu(voisine)-tab(nb_faces_perio+voisine-num_face);
                      resu(voisine)  += diff1;
                      resu(num_face) += diff2;
                    }
                  else
                    for (int comp=0; comp<ncomp_ch_transporte; comp++)
                      {
                        diff1 = resu(num_face,comp)-tab(nb_faces_perio,comp);
                        diff2 = resu(voisine,comp)-tab(nb_faces_perio+voisine-num_face,comp);
                        resu(voisine,comp)  += diff1;
                        resu(num_face,comp) += diff2;
                      }

                  fait[num_face-num1]= 1;
                  fait[voisine-num1] = 1;
                }
              nb_faces_perio++;
            }
        }
    }

}
void Op_Conv_ALE_VEF::supprimerALE_Muscl(const DoubleTab& inco, DoubleTab& resu) const
{

  const Domaine_ALE& dom_ale=ref_cast(Domaine_ALE, dom.valeur());
  const Op_Conv_VEF_base& opConvVEFbase = ref_cast(Op_Conv_VEF_base, op_conv.valeur());
  const Op_Conv_VEF_Face& opConvVEFFace = ref_cast(Op_Conv_VEF_Face, op_conv.valeur());

  //This is based on Op_Conv_VEF_Face::ajouter(...) from where only "muscl" related part is used. ALE convection term has the opposite sign compared to the normal convection term.


  int ordre;
  opConvVEFFace.get_ordre(ordre);

  double alpha_;
  opConvVEFFace.get_alpha(alpha_);

  Motcle type_limit;
  if(op_conv.type()=="MUSCL") //In case of: "convection {  ALE { muscl } }" , type_lumit has default value set in Op_Conv_Muscl_VEF_Face.cpp .
    {
      type_limit="vanleer";
    }
  else // In case of: "convection {  ALE { generic muscl [limiter] [order of accuracy] [alpha] } }". op_conv.type()=="GENERIC" .
    {
      opConvVEFFace.get_type_lim(type_limit);
    }

  const Zone_Cl_VEF& zone_Cl_VEF = la_zcl_vef.valeur();
  const Zone_VEF& zone_VEF = ref_cast(Zone_VEF, la_zone_vef.valeur());

  const DoubleTab& vitesse_face_absolue_ALE=dom_ale.vitesse_faces();              //ALE mesh velocity.

  const DoubleVect& porosite_face = zone_VEF.porosite_face();

  int marq=opConvVEFbase.phi_u_transportant(equation());
  DoubleTab transporte_face_;
  DoubleTab vitesse_face_;
  // soit on a transporte_face=phi*transporte et vitesse_face=vitesse
  // soit on a transporte_face=transporte et vitesse_face=phi*vitesse
  // cela depend si on transporte avec phi*u ou avec u.
  const DoubleTab& transporte_face = modif_par_porosite_si_flag(inco,transporte_face_,!marq,porosite_face); //Fluid velocity = inco .
  const DoubleTab& vitesse_face_ALEmodif    = modif_par_porosite_si_flag(vitesse_face_absolue_ALE,vitesse_face_,marq,porosite_face);

  const IntTab& elem_faces = zone_VEF.elem_faces();
  const DoubleTab& facenormales = zone_VEF.face_normales();
  const DoubleTab& facette_normales = zone_VEF.facette_normales();
  const Zone& zone = zone_VEF.zone();
  const int nfa7 = zone_VEF.type_elem().nb_facette();
  const int nb_elem_tot = zone_VEF.nb_elem_tot();
  const IntVect& rang_elem_non_std = zone_VEF.rang_elem_non_std();
  const IntTab& face_voisins = zone_VEF.face_voisins();
  const DoubleTab& normales_facettes_Cl = zone_Cl_VEF.normales_facettes_Cl();
  int premiere_face_int = zone_VEF.premiere_face_int();
  int nfac = zone.nb_faces_elem();
  int nsom = zone.nb_som_elem();
  const IntTab& sommet_elem = zone.les_elems();
  const DoubleTab& vecteur_face_facette = ref_cast_non_const(Zone_VEF,zone_VEF).vecteur_face_facette();
  const  DoubleTab& vecteur_face_facette_Cl = zone_Cl_VEF.vecteur_face_facette_Cl();
  int nb_bord = zone_VEF.nb_front_Cl();
  const IntTab& les_elems=zone.les_elems();

  // Permet d'avoir un flux_bord coherent avec les CLs (mais parfois diverge?)
  // Active uniquement pour ordre 3
  int option_appliquer_cl_dirichlet = (ordre==3 ? 1 : 0);
  int option_calcul_flux_en_un_point = 0;//(ordre==3 ? 1 : 0);
  // Definition d'un tableau pour un traitement special des schemas pres des bords
  if (traitement_pres_bord_.size_array()!=nb_elem_tot)
    {
      traitement_pres_bord_.resize_array(nb_elem_tot);
      traitement_pres_bord_=0;
      // Pour muscl3 on applique le minmod sur les elements ayant une face de Dirichlet
      if (ordre==3)
        {
          for (int n_bord=0; n_bord<nb_bord; n_bord++)
            {
              const Cond_lim_base& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord).valeur();
              if ( sub_type(Dirichlet,la_cl) || sub_type(Dirichlet_homogene,la_cl) )
                {
                  const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
                  int nb_faces_tot = le_bord.nb_faces_tot();
                  for (int ind_face=0; ind_face<nb_faces_tot; ind_face++)
                    {
                      int num_face = le_bord.num_face(ind_face);
                      int elem = face_voisins(num_face,0);
                      traitement_pres_bord_(elem)=1;
                    }
                }
            }
        }
      else
        {
          // Pour le muscl/centre actuels on utilise un calcul de flux a l'ordre 1
          // aux mailles de bord ou aux mailles ayant un sommet de Dirichlet
          ArrOfInt est_un_sommet_de_bord_(zone_VEF.nb_som_tot());
          for (int n_bord=0; n_bord<nb_bord; n_bord++)
            {
              const Cond_lim_base& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord).valeur();
              if ( sub_type(Dirichlet,la_cl) || sub_type(Dirichlet_homogene,la_cl) )
                {
                  const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
                  int nb_faces_tot = le_bord.nb_faces_tot();
                  int size = zone_VEF.face_sommets().dimension(1);
                  for (int ind_face=0; ind_face<nb_faces_tot; ind_face++)
                    for (int som=0; som<size; som++)
                      {
                        int face = le_bord.num_face(ind_face);
                        est_un_sommet_de_bord_(zone_VEF.face_sommets(face,som))=1;
                      }
                }
            }
          for (int elem=0; elem<nb_elem_tot; elem++)
            {
              if (rang_elem_non_std(elem)!=-1)
                traitement_pres_bord_(elem)=1;
              else
                {
                  for (int n_som=0; n_som<nsom; n_som++)
                    if (est_un_sommet_de_bord_(les_elems(elem,n_som)))
                      traitement_pres_bord_(elem)=1;
                }
            }
        }
      // Construction du tableau est_une_face_de_dirichlet_
      est_une_face_de_dirichlet_.resize_array(zone_VEF.nb_faces_tot());
      est_une_face_de_dirichlet_=0;
      for (int n_bord=0; n_bord<nb_bord; n_bord++)
        {
          const Cond_lim_base& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord).valeur();
          if ( sub_type(Dirichlet,la_cl) || sub_type(Dirichlet_homogene,la_cl) )
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              int nb_faces_tot = le_bord.nb_faces_tot();
              for (int ind_face=0; ind_face<nb_faces_tot; ind_face++)
                {
                  int num_face = le_bord.num_face(ind_face);
                  est_une_face_de_dirichlet_(num_face) = 1;
                }
            }
        }
    }

  // Pour le traitement de la convection on distingue les polyedres
  // standard qui ne "voient" pas les conditions aux limites et les
  // polyedres non standard qui ont au moins une face sur le bord.
  // Un polyedre standard a n facettes sur lesquelles on applique le
  // schema de convection.
  // Pour un polyedre non standard qui porte des conditions aux limites
  // de Dirichlet, une partie des facettes sont portees par les faces.
  // En bref pour un polyedre le traitement de la convection depend
  // du type (triangle, tetraedre ...) et du nombre de faces de Dirichlet.

  const Elem_VEF_base& type_elemvef= zone_VEF.type_elem().valeur();
  int istetra=0;
  Nom nom_elem=type_elemvef.que_suis_je();
  if ((nom_elem=="Tetra_VEF")||(nom_elem=="Tri_VEF"))
    istetra=1;

  const DoubleVect& porosite_elem = zone_VEF.porosite_elem();
  double psc;
  int poly,face_adj,fa7,i,j,n_bord;
  int num_face, rang;
  int num10,num20;
  int ncomp_ch_transporte=(transporte_face.nb_dim() == 1?1:transporte_face.dimension(1));

  // Traitement particulier pour les faces de periodicite
  int nb_faces_perio = 0;
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          nb_faces_perio+=le_bord.nb_faces();
        }
    }

  DoubleTab tab;
  if (ncomp_ch_transporte == 1)
    tab.resize(nb_faces_perio);
  else
    tab.resize(nb_faces_perio,ncomp_ch_transporte);

  nb_faces_perio=0;
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (num_face=num1; num_face<num2; num_face++)
            {
              if (ncomp_ch_transporte == 1)
                tab(nb_faces_perio) = resu(num_face);
              else
                for (int comp=0; comp<ncomp_ch_transporte; comp++)
                  tab(nb_faces_perio,comp) = resu(num_face,comp);
              nb_faces_perio++;
            }
        }
    }

  int fac=0,elem1,elem2,comp0;
  int nb_faces_ = zone_VEF.nb_faces();
  ArrOfInt face(nfac);
  //statistiques().end_count(m1);
  //statistiques().begin_count(m2);

  // Tableau gradient base sur gradient_elem selon schema
  DoubleTab gradient_elem(nb_elem_tot,ncomp_ch_transporte,dimension);  // (du/dx du/dy dv/dx dv/dy) pour un poly

  Champ_P1NC::calcul_gradient(transporte_face,gradient_elem,zone_Cl_VEF);

  DoubleTab gradient;

  //  application du limiteur
  gradient.resize(0, ncomp_ch_transporte, dimension);     // (du/dx du/dy dv/dx dv/dy) pour une face
  zone_VEF.creer_tableau_faces(gradient);
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
      int num1 = le_bord.num_premiere_face();
      int num2 = num1 + le_bord.nb_faces();
      if (sub_type(Periodique,la_cl.valeur()))
        {
          for (fac=num1; fac<num2; fac++)
            {
              elem1=face_voisins(fac,0);
              elem2=face_voisins(fac,1);
              for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                for (i=0; i<dimension; i++)
                  {
                    double grad1=gradient_elem(elem1, comp0, i);
                    double grad2=gradient_elem(elem2, comp0, i);
                    //gradient(fac, comp0, i) =(*LIMITEUR)(grad1, grad2);
                    gradient(fac, comp0, i) =application_LIMITEUR(grad1, grad2, type_limit);
                  }
            }
        }
      else if (sub_type(Symetrie,la_cl.valeur()))
        {
          for (fac=num1; fac<num2; fac++)
            {
              elem1=face_voisins(fac,0);
              for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                for (i=0; i<dimension; i++)
                  gradient(fac, comp0, i) = gradient_elem(elem1, comp0, i);

              if (ordre==3)
                {
                  // On enleve la composante normale (on pourrait le faire pour les autres schemas...)
                  // mais pour le moment, on ne veut pas changer le comportement par defaut du muscl...
                  //const DoubleTab& facenormales = zone_VEF.face_normales();
                  for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                    for (i=0; i<dimension; i++)
                      {
                        double carre_surface=0;
                        double tmp=0;
                        for (j=0; j<dimension; j++)
                          {
                            double ndS=facenormales(fac,j);
                            carre_surface += ndS*ndS;
                            tmp += gradient(fac, comp0, j)*ndS;
                          }
                        gradient(fac, comp0, i) -= tmp*facenormales(fac,i)/carre_surface;
                      }
                }
            }
        }
    }

  for (fac=premiere_face_int; fac<nb_faces_; fac++)
    {
      elem1=face_voisins(fac,0);
      elem2=face_voisins(fac,1);
      int minmod_pres_du_bord = 0;
      if (ordre==3 && (traitement_pres_bord_(elem1) || traitement_pres_bord_(elem2))) minmod_pres_du_bord = 1;
      for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
        for (i=0; i<dimension; i++)
          {
            double grad1=gradient_elem(elem1, comp0, i);
            double grad2=gradient_elem(elem2, comp0, i);
            if (minmod_pres_du_bord)
              gradient(fac, comp0, i) = minmod(grad1, grad2);
            else
              //gradient(fac, comp0, i) = (*LIMITEUR)(grad1, grad2);
              gradient(fac, comp0, i) =application_LIMITEUR(grad1, grad2, type_limit);
          }
    } // fin du for faces
  gradient.echange_espace_virtuel();

  ArrOfDouble vs(dimension);
  ArrOfDouble vc(dimension);
  ArrOfDouble cc(dimension);
  DoubleVect xc(dimension);
  DoubleTab vsom(nsom,dimension);
  DoubleTab xsom(nsom,dimension);

  const IntTab& KEL=type_elemvef.KEL();
  const DoubleTab& xv=zone_VEF.xv();
  const Domaine& domaine=zone.domaine();
  const DoubleTab& coord_sommets=domaine.coord_sommets();

  // Boucle ou non selon la valeur de alpha (uniquement a l'ordre 3 pour le moment)
  // Si alpha=1, la boucle se limite a une simple passe avec le schema choisi (muscl, amont, centre)
  // Si alpha<1, la boucle se compose de 2 passes:
  //                         -la premiere avec le schema choisi et une ponderation de alpha
  //                         -la seconde avec le schema centre et une ponderation de 1-alpha
  double alpha = alpha_;
  int nombre_passes = (alpha==1 ? 1 : 2);
  for (int passe=1; passe<=nombre_passes; passe++)
    {
      if (passe==2)
        {
          gradient.ref(gradient_elem);
        }
      // Les polyedres non standard sont ranges en 2 groupes dans la Zone_VEF:
      //  - polyedres bords et joints
      //  - polyedres bords et non joints
      // On traite les polyedres en suivant l'ordre dans lequel ils figurent
      // dans la zone
      // boucle sur les polys
      for (poly=0; poly<nb_elem_tot; poly++)
        {
          int contrib = 0;
          // calcul des numeros des faces du polyedre
          for (face_adj=0; face_adj<nfac; face_adj++)
            {
              int face_ = elem_faces(poly,face_adj);
              face(face_adj)= face_;
              if (face_<nb_faces_) contrib=1; // Une face reelle sur l'element virtuel
            }
          //
          if (contrib)
            {
              int calcul_flux_en_un_point = (ordre != 3) && (ordre==1 || traitement_pres_bord_(poly));
              for (j=0; j<dimension; j++)
                {
                  vs(j) = vitesse_face_absolue_ALE(face(0),j)*porosite_face(face(0));
                  for (i=1; i<nfac; i++)
                    vs(j)+= vitesse_face_absolue_ALE(face(i),j)*porosite_face(face(i));
                }
              // calcul de la vitesse aux sommets des polyedres
              // On va utliser les fonctions de forme implementees dans la classe Champs_P1_impl ou Champs_Q1_impl
              if (istetra==1)
                {
                  for (i=0; i<nsom; i++)
                    for (j=0; j<dimension; j++)
                      vsom(i,j) = (vs(j) - dimension*vitesse_face_absolue_ALE(face(i),j)*porosite_face(face(i)));
                }
              else
                {
                  Cout << "ALE is currently working only with Tetra_VEF (3D) or Tri_VEF (2D)." << finl;
                  exit();
                }

              // Determination du type de CL selon le rang
              rang = rang_elem_non_std(poly);
              int itypcl = (rang==-1 ? 0 : zone_Cl_VEF.type_elem_Cl(rang));

              // calcul de vc (a l'intersection des 3 facettes) vc vs vsom proportionnelles a la porosite
              calcul_vc_ALE(face,vc,vs,vsom,vitesse_face_absolue_ALE,itypcl,porosite_face);

              // calcul de xc (a l'intersection des 3 facettes) necessaire pour muscl3
              if (ordre==3)
                {
                  int idirichlet;
                  int n1,n2,n3;
                  for (i=0; i<nsom; i++)
                    for (j=0; j<dimension; j++)
                      xsom(i,j) = coord_sommets(les_elems(poly,i),j);
                  type_elemvef.calcul_xg(xc,xsom,itypcl,idirichlet,n1,n2,n3);
                }

              // Gestion de la porosite
              if (marq==0)
                {
                  double coeff=1./porosite_elem(poly);
                  vsom*=coeff;
                  vc*=coeff;
                }
              // Boucle sur les facettes du polyedre non standard:
              for (fa7=0; fa7<nfa7; fa7++)
                {
                  num10 = face(KEL(0,fa7));
                  num20 = face(KEL(1,fa7));
                  // normales aux facettes
                  if (rang==-1)
                    for (i=0; i<dimension; i++)
                      cc[i] = facette_normales(poly, fa7, i);
                  else
                    for (i=0; i<dimension; i++)
                      cc[i] = normales_facettes_Cl(rang,fa7,i);

                  // Calcul des vitesses en C,S,S2 les 3 extremites de la fa7 et M le centre de la fa7
                  double psc_c=0,psc_s=0,psc_m,psc_s2=0;
                  if (dimension==2)
                    {
                      for (i=0; i<2; i++)
                        {
                          psc_c+=vc[i]*cc[i];
                          psc_s+=vsom(KEL(2,fa7),i)*cc[i];
                        }
                      psc_m=(psc_c+psc_s)/2.;
                    }
                  else
                    {
                      for (i=0; i<3; i++)
                        {
                          psc_c+=vc[i]*cc[i];
                          psc_s+=vsom(KEL(2,fa7),i)*cc[i];
                          psc_s2+=vsom(KEL(3,fa7),i)*cc[i];
                        }
                      psc_m=(psc_c+psc_s+psc_s2)/3.;
                    }
                  // On applique les CL de Dirichlet si num1 ou num2 est une face avec CL de Dirichlet
                  // auquel cas la fa7 coincide avec la face num1 ou num2 -> C est au centre de la face
                  int appliquer_cl_dirichlet=0;
                  if (option_appliquer_cl_dirichlet)
                    if (est_une_face_de_dirichlet_(num10) || est_une_face_de_dirichlet_(num20))
                      {
                        appliquer_cl_dirichlet = 1;
                        psc_m = psc_c;
                      }

                  // Determination de la face amont pour M
                  int face_amont_m,dir;
                  if (psc_m >= 0)
                    {
                      face_amont_m = num10;
                      dir=0;
                    }
                  else
                    {
                      face_amont_m = num20;
                      dir=1;
                    }
                  // Determination des faces amont pour les points C,S,S2
                  int face_amont_c=face_amont_m;
                  int face_amont_s=face_amont_m;
                  int face_amont_s2=face_amont_m;
                  if (passe==1 && ordre==3)
                    {
                      face_amont_c  = (psc_c >= 0)  ? num10 : num20;
                      face_amont_s  = (psc_s >= 0)  ? num10 : num20;
                      face_amont_s2 = (psc_s2 >= 0) ? num10 : num20;
                    }
                  // gradient aux items element (schema centre) ou aux items face (schemas muscl)
                  int item_m=poly;
                  int item_c=poly;
                  int item_s=poly;
                  int item_s2=poly;
                  if (passe==1)
                    {
                      item_m = face_amont_m;
                      item_c = face_amont_c;
                      item_s = face_amont_s;
                      item_s2 = face_amont_s2;
                    }

                  for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
                    {
                      double flux;
                      double inco_m = (ncomp_ch_transporte==1?transporte_face(face_amont_m):transporte_face(face_amont_m,comp0));
                      if (appliquer_cl_dirichlet)
                        {
                          flux = inco_m*psc_m;
                        }

                      else // muscl ou centre
                        {
                          // Calcul de l'inconnue au centre M de la fa7
                          if (rang==-1)
                            for (j=0; j<dimension; j++)
                              inco_m+= gradient(item_m,comp0,j)*vecteur_face_facette(poly,fa7,j,dir);
                          else
                            for (j=0; j<dimension; j++)
                              inco_m+= gradient(item_m,comp0,j)*vecteur_face_facette_Cl(rang,fa7,j,dir);

                          // Calcul de l'inconnue au sommet S, une premiere extremite de la fa7
                          double inco_s = (ncomp_ch_transporte==1?transporte_face(face_amont_s):transporte_face(face_amont_s,comp0));
                          int sommet_s = sommet_elem(poly,KEL(2,fa7));
                          for (j=0; j<dimension; j++)
                            inco_s+= gradient(item_s,comp0,j)*(-xv(face_amont_s,j)+coord_sommets(sommet_s,j));

                          // Calcul de l'inconnue au sommet S2, la derniere extremite de la fa7 en 3D
                          double inco_s2=0;
                          if (dimension==3)
                            {
                              inco_s2 = (ncomp_ch_transporte==1?transporte_face(face_amont_s2):transporte_face(face_amont_s2,comp0));
                              int sommet_s2 = sommet_elem(poly,KEL(3,fa7));
                              for (j=0; j<dimension; j++)
                                inco_s2+= gradient(item_s2,comp0,j)*(-xv(face_amont_s2,j)+coord_sommets(sommet_s2,j));
                            }

                          // Calcul de l'inconnue a C, une autre extremite de la fa7, intersection avec les autres fa7
                          // du polyedre. C=G centre du polyedre si volume non etendu
                          // xc donne par elemvef.calcul_xg()
                          double inco_c;
                          if (ordre==3)
                            {
                              inco_c = (ncomp_ch_transporte==1?transporte_face(face_amont_c):transporte_face(face_amont_c,comp0));
                              for (j=0; j<dimension; j++)
                                inco_c+= gradient(item_c,comp0,j)*(-xv(face_amont_c,j)+xc(j));
                            }
                          else
                            {
                              inco_c = dimension*inco_m-inco_s-inco_s2;
                            }

                          // Calcul du flux sur 1 point
                          if (calcul_flux_en_un_point || option_calcul_flux_en_un_point)
                            {
                              flux = inco_m*psc_m;
                            }
                          else
                            {
                              // Calcul du flux sur 3 points
                              flux = (dimension==2) ? (inco_c*psc_c + inco_s*psc_s + 4*inco_m*psc_m)/6
                                     : (inco_c*psc_c + inco_s*psc_s + inco_s2*psc_s2 + 9*inco_m*psc_m)/12;
                            }

                        }
                      // Ponderation par coefficient alpha
                      flux*=alpha;

                      if (ncomp_ch_transporte == 1)
                        {
                          resu(num10) -= flux;
                          resu(num20) += flux;
                        }
                      else
                        {
                          resu(num10,comp0) -= flux;
                          resu(num20,comp0) += flux;
                        }

                    }// boucle sur comp
                } // fin de la boucle sur les facettes
            }
        } // fin de la boucle
      alpha = 1 - alpha;
    } // fin de la boucle
  //statistiques().end_count(m2);
  //statistiques().begin_count(m3);
  int voisine;
  nb_faces_perio = 0;
  double diff1,diff2;

  // Boucle sur les bords pour traiter les conditions aux limites
  // il y a prise en compte d'un terme de convection pour les
  // conditions aux limites de Neumann_sortie_libre seulement
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);

      if (sub_type(Neumann_sortie_libre,la_cl.valeur()))
        {
          const Neumann_sortie_libre& la_sortie_libre = ref_cast(Neumann_sortie_libre, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (num_face=num1; num_face<num2; num_face++)
            {
              psc =0;
              for (i=0; i<dimension; i++)
                psc += vitesse_face_ALEmodif(num_face,i)*facenormales(num_face,i);
              if (psc>0)
                if (ncomp_ch_transporte == 1)
                  {
                    resu(num_face) -= psc*transporte_face(num_face);
                  }
                else
                  for (i=0; i<ncomp_ch_transporte; i++)
                    {
                      resu(num_face,i) -= psc*transporte_face(num_face,i);
                    }
              else
                {
                  if (ncomp_ch_transporte == 1)
                    {
                      resu(num_face) -= psc*la_sortie_libre.val_ext(num_face-num1);
                    }
                  else
                    for (i=0; i<ncomp_ch_transporte; i++)
                      {
                        resu(num_face,i) -= psc*la_sortie_libre.val_ext(num_face-num1,i);
                      }
                }
            }
        }
      else if (sub_type(Periodique,la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          ArrOfInt fait(le_bord.nb_faces());
          fait = 0;
          for (num_face=num1; num_face<num2; num_face++)
            {
              if (fait[num_face-num1] == 0)
                {
                  voisine = la_cl_perio.face_associee(num_face-num1) + num1;

                  if (ncomp_ch_transporte == 1)
                    {
                      diff1 = resu(num_face)-tab(nb_faces_perio);
                      diff2 = resu(voisine)-tab(nb_faces_perio+voisine-num_face);
                      resu(voisine)  += diff1;
                      resu(num_face) += diff2;
                    }
                  else
                    for (int comp=0; comp<ncomp_ch_transporte; comp++)
                      {
                        diff1 = resu(num_face,comp)-tab(nb_faces_perio,comp);
                        diff2 = resu(voisine,comp)-tab(nb_faces_perio+voisine-num_face,comp);
                        resu(voisine,comp)  += diff1;
                        resu(num_face,comp) += diff2;
                      }

                  fait[num_face-num1]= 1;
                  fait[voisine-num1] = 1;
                }
              nb_faces_perio++;
            }
        }
    }


}

// contribuer_a_avec() and modifier_pour_Cl() are needed when using Scheme_euler_implicit time scheme.

//Description:
//On assemble la matrice des inconnues implicite.
void Op_Conv_ALE_VEF::contribuer_a_avec(const DoubleTab& inco, Matrice_Morse& matrice) const
{
  const Op_Conv_VEF_Face& opConvVEFFace = ref_cast(Op_Conv_VEF_Face, op_conv.valeur());
  opConvVEFFace.ajouter_contribution(inco, matrice);
}

// Description:
// On modifie le second membre et la matrice dans le cas des
// conditions de dirichlet.
void Op_Conv_ALE_VEF::modifier_pour_Cl(Matrice_Morse& matrice, DoubleTab& secmem) const
{
  const Op_Conv_VEF_Face& opConvVEFFace = ref_cast(Op_Conv_VEF_Face, op_conv.valeur());
  opConvVEFFace.modifier_pour_Cl(matrice, secmem);
}
