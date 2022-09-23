/****************************************************************************
* Copyright (c) 2022, CEA
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

#include <Op_Conv_RT_VEF_Face.h>
#include <Champ_P1NC.h>
#include <Porosites_champ.h>
#include <stat_counters.h>
#include <Zone_VEF_PreP1b.h>
#include <CL_Types_include.h>

extern double calculer_coef_som(int elem, int& nb_face_diri, ArrOfInt& indice_diri, const Zone_Cl_VEF& zcl, const Zone_VEF& zone_VEF);
Implemente_instanciable_sans_constructeur(Op_Conv_RT_VEF_Face,"Op_Conv_RT_VEF_P1NC",Op_Conv_VEF_Face);
// XD convection_RT convection_deriv RT 0 Keyword to use RT projection for P1NCP0RT discretization

//// printOn
//
Sortie& Op_Conv_RT_VEF_Face::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}

//// readOn
//
Entree& Op_Conv_RT_VEF_Face::readOn(Entree& s )
{
  ordre=1;
  return s ;
}

//
//   Fonctions de la classe Op_Conv_VEF_Face
//
////////////////////////////////////////////////////////////////////
//
//                      Implementation des fonctions
//
//                   de la classe Op_Conv_VEF_Face
//
////////////////////////////////////////////////////////////////////




DoubleTab& Op_Conv_RT_VEF_Face::ajouter(const DoubleTab& transporte,
                                        DoubleTab& resu) const
{
  //statistiques().begin_count(m1);
  const Zone_Cl_VEF& zone_Cl_VEF = la_zcl_vef.valeur();
  const Zone_VEF& zone_VEF = ref_cast(Zone_VEF, la_zone_vef.valeur());
  const Champ_Inc_base& la_vitesse=vitesse();
  const DoubleTab& vitesse_face_absolue=la_vitesse.valeurs();
  const DoubleVect& porosite_face = equation().milieu().porosite_face();

  int marq=phi_u_transportant(equation());
  DoubleTab transporte_face_;
  DoubleTab vitesse_face_;
  // soit on a transporte_face=phi*transporte et vitesse_face=vitesse
  // soit on a transporte_face=transporte et vitesse_face=phi*vitesse
  // cela depend si on transporte avec phi*u ou avec u.
  const DoubleTab& transporte_face = modif_par_porosite_si_flag(transporte,transporte_face_,!marq,porosite_face);
  const DoubleTab& vitesse_face    = modif_par_porosite_si_flag(vitesse_face_absolue,vitesse_face_,marq,porosite_face);

  const IntTab& elem_faces = zone_VEF.elem_faces();
  const DoubleTab& facenormales = zone_VEF.face_normales();
  const DoubleTab& facette_normales = zone_VEF.facette_normales();
  const Zone& zone = zone_VEF.zone();
  const int nb_elem_tot = zone_VEF.nb_elem_tot();
  const IntVect& rang_elem_non_std = zone_VEF.rang_elem_non_std();
  const IntTab& face_voisins = zone_VEF.face_voisins();
  const DoubleTab& normales_facettes_Cl = zone_Cl_VEF.normales_facettes_Cl();
  const IntTab& elem_sommets = zone.les_elems();
  const DoubleTab& coord_sommets=zone.domaine().les_sommets();
  const DoubleTab& vecteur_face_facette = ref_cast_non_const(Zone_VEF,zone_VEF).vecteur_face_facette();
  const  DoubleTab& vecteur_face_facette_Cl = zone_Cl_VEF.vecteur_face_facette_Cl();
  const IntTab& les_elems=zone.les_elems();

  // Permet d'avoir un flux_bord coherent avec les CLs (mais parfois diverge?)
  // Active uniquement pour ordre 3
  //

  {
    const DoubleVect& volumes = zone_VEF.volumes();
    const DoubleTab& face_normales = zone_VEF.face_normales();
    double volume;
    // Cout<<"Op_Conv_RT_VEF_Face::ajouter RT\n";
    // Traitement des CL de Dirichlet
    int nb_face_diri=0;
    ArrOfInt indice_diri(dimension+1);
    int modif_traitement_diri=0;
    if (sub_type(Zone_VEF_PreP1b,zone_VEF))
      modif_traitement_diri=ref_cast(Zone_VEF_PreP1b,zone_VEF).get_modif_div_face_dirichlet();
    int elem,i,j,alfa,dim;
    i=0;
    if (dimension==2)
      {
        for (elem=0; elem<nb_elem_tot; elem++)
          {
            DoubleTab coordSommet(dimension+1,dimension);
            for (i=0; i<=dimension; i++)
              {
                int numSom=elem_sommets(elem, i);
                for (dim=0; dim<dimension; dim++) coordSommet(i,dim)=coord_sommets(numSom,dim);
              }

            if (modif_traitement_diri)
              calculer_coef_som(elem,nb_face_diri,indice_diri,zone_Cl_VEF,zone_VEF);
            volume=volumes(elem);
            double invVol  = 1./(12*volume);
            double invVol2 = invVol/volume;
            // somme (vitesse_dir) fois (face normale_dir)
            DoubleTab FacesNormales(dimension+1,dimension);
            DoubleVect vitFaceNormale(dimension+1);
            for (j=0; j<=dimension; j++)
              {
                int num_face=elem_faces(elem, j);
                for (dim=0; dim<dimension; dim++) FacesNormales(j,dim)=face_normales(num_face,dim);
                if (elem!=face_voisins(num_face,0))
                  {
                    for (dim=0; dim<dimension; dim++) FacesNormales(j,dim)=-FacesNormales(j,dim);
                  }
                vitFaceNormale(j)=0.;
                for (dim=0; dim<dimension; dim++)
                  vitFaceNormale(j)+=vitesse_face_absolue(num_face,dim)*FacesNormales(j,dim);
              }
            // Calculer le rot de la vitesse dans T
            double rotVit=0.;
            for (j=0; j<=dimension; j++)
              {
                int num_face=elem_faces(elem, j);
                rotVit+=vitesse_face_absolue(num_face,1)*FacesNormales(j,0)-vitesse_face_absolue(num_face,0)*FacesNormales(j,1);
              }
            rotVit*=invVol2;
            DoubleTab FiFj(dimension+1,dimension+1);
            for (i=0; i<=dimension; i++)
              {
                for (j=i+1; j<=dimension; j++)
                  {
                    FiFj(i,j) = FacesNormales(j,0)*FacesNormales(i,1)-FacesNormales(j,1)*FacesNormales(i,0);
                    FiFj(j,i) =-FiFj(i,j);
                  }
                FiFj(i,i)=0.;
              }
            DoubleVect resu_face(dimension+1);
            for (i=0; i<=dimension; i++)
              {
                resu_face(i)=0;
                for (j=0; j<=dimension; j++)
                  {
                    if (i!=j) resu_face(i) += vitFaceNormale(j)*FiFj(j,i);
                  }
              }
            for (i=0; i<=dimension; i++)
              {
                int num_face=elem_faces(elem, i);
                // Code modifié pour prise en compte CL de Dirichlet
                for (dim=0; dim<dimension; dim++)
                  {
                    double contrib_resu=rotVit*FacesNormales(i,dim)*resu_face(i);
                    resu(num_face, dim)+=contrib_resu;
                    // Traitement CL de Dirichlet
                    for (int fdiri=0; fdiri<nb_face_diri; fdiri++)
                      {
                        int indice=indice_diri[fdiri];
                        int facel = elem_faces(elem,indice);
                        if (num_face==facel)
                          {
                            resu(facel,dim)-=contrib_resu;
                            double contrib_resu2=contrib_resu/(dimension+1-nb_face_diri);
                            for (int f2=0; f2<dimension+1; f2++)
                              {
                                // Cerr<<num_face_i<<" "<<elem<<" la "<< f2<<" "<<fdiri<<" "<<nb_face_diri<<" dim "<< dim <<finl;
                                int face2=elem_faces(elem,f2);
                                resu(face2,dim)+=contrib_resu2;
                              }
                          }
                      }
                  }

              }
          }
      }
    else if (dimension==3)
      {
        for ( elem=0; elem<nb_elem_tot; elem++)
          {
            DoubleTab coordSommet(dimension+1,dimension);
            for (i=0; i<=dimension; i++)
              {
                int numSom=elem_sommets(elem, i);
                for (dim=0; dim<dimension; dim++) coordSommet(i,dim)=coord_sommets(numSom,dim);
              }
            // Calculer le rot de la vitesse dans elem
            DoubleVect rotVit(dimension);
            // Calculer rotVit.FiFj
            DoubleTab  rotVitFiFj(dimension+1,dimension+1);
            DoubleTab FacesNormales(dimension+1,dimension);
            DoubleVect vitFaceNormale(dimension+1);
            if (modif_traitement_diri)
              calculer_coef_som(elem,nb_face_diri,indice_diri,zone_Cl_VEF,zone_VEF);
            volume=volumes(elem);
            for (dim=0; dim<dimension; dim++) rotVit(dim)=0.;
            for (i=0; i<=dimension; i++)
              {
                int num_face=elem_faces(elem, i);
                for (dim=0; dim<dimension; dim++)
                  FacesNormales(i,dim)=face_normales(num_face,dim);
                if (elem==face_voisins(num_face,0))
                  {
                    for (dim=0; dim<dimension; dim++)
                      FacesNormales(i,dim)=-FacesNormales(i,dim);
                  }
                //
                rotVit(0)+=vitesse_face_absolue(num_face,2)*FacesNormales(i,1)-vitesse_face_absolue(num_face,1)*FacesNormales(i,2);
                rotVit(1)+=vitesse_face_absolue(num_face,0)*FacesNormales(i,2)-vitesse_face_absolue(num_face,2)*FacesNormales(i,0);
                rotVit(2)+=vitesse_face_absolue(num_face,1)*FacesNormales(i,0)-vitesse_face_absolue(num_face,0)*FacesNormales(i,1);
                for (j=0; j<=dimension; j++) rotVitFiFj(i,j)=0.;
                vitFaceNormale(i)=0.;
                for (dim=0; dim<dimension; dim++) vitFaceNormale(i)+=vitesse_face_absolue(num_face,dim)*FacesNormales(i,dim);

              } // end_for (i=0; i<=dimension; i++)
            double rotVol=1./(18.*volume*volume);
            for (dim=0; dim<dimension; dim++) rotVit(dim)*=rotVol;
            // remplissage de rot u.SimoinsSj

            for (dim=0; dim<dimension; dim++)
              {
                rotVitFiFj(0,1)+=rotVit(dim)*(FacesNormales(2,dim)-FacesNormales(3,dim));
                rotVitFiFj(0,2)+=rotVit(dim)*(FacesNormales(3,dim)-FacesNormales(1,dim));
                rotVitFiFj(0,3)+=rotVit(dim)*(FacesNormales(1,dim)-FacesNormales(2,dim));
                rotVitFiFj(1,2)+=rotVit(dim)*(FacesNormales(0,dim)-FacesNormales(3,dim));
                rotVitFiFj(1,3)+=rotVit(dim)*(FacesNormales(2,dim)-FacesNormales(0,dim));
                rotVitFiFj(2,3)+=rotVit(dim)*(FacesNormales(0,dim)-FacesNormales(1,dim));
              }
            for (i=0; i<=dimension; i++)
              {
                for (j=i+1; j<=dimension; j++) rotVitFiFj(j,i)=-rotVitFiFj(i,j);
              }
            DoubleVect resu_face(dimension+1);
            for (i=0; i<=dimension; i++)
              {
                resu_face(i)=0.;
                for (j=0; j<=dimension; j++)
                  {
                    if (i!=j) resu_face(i)+=vitFaceNormale(j)*rotVitFiFj(i,j);
                  }
              }
            for (i=0; i<=dimension; i++)
              {
                int num_face=elem_faces(elem,i);
                for (alfa=0; alfa<dimension; alfa++)
                  {
                    double contrib=resu_face(i)*FacesNormales(i,alfa);
                    resu(num_face,alfa)+=contrib;
                    // Traitement CL de Dirichlet
                    for (int fdiri=0; fdiri<nb_face_diri; fdiri++)
                      {
                        int indice=indice_diri[fdiri];
                        int facel = elem_faces(elem,indice);
                        if (num_face==facel)
                          {
                            resu(num_face,alfa)-=contrib;
                            double contrib2=contrib/(dimension+1-nb_face_diri);
                            for (int f2=0; f2<dimension+1; f2++)
                              {
                                // Cerr<<num_face_i<<" "<<elem<<" la "<< f2<<" "<<fdiri<<" "<<nb_face_diri<<" dim "<< dim <<finl;
                                int face2=elem_faces(elem,f2);
                                resu(face2,alfa)+=contrib2;
                              } // end_for (int f2=0; f2<dimension+1; f2++)
                          } // end_if (num_face_i==facel)
                      } // end_for (int fdiri=0; fdiri<nb_face_diri; fdiri++)

                  } // end_for (dim=0; dim<dimension;dim++)
              } // end_for (i=0;i<=dimension;i++)
          } // end_for (elem=0; elem<nb_elem_tot; elem++)
      } // end_if (dim==3)
  } // end_if (type_op==RT)
  modifier_flux(*this);
  //statistiques().end_count(m3);
  return resu;
}

void Op_Conv_RT_VEF_Face::ajouter_contribution(const DoubleTab& transporte, Matrice_Morse& matrice ) const
{
  modifier_matrice_pour_periodique_avant_contribuer(matrice,equation());
  const Zone_Cl_VEF& zone_Cl_VEF = la_zcl_vef.valeur();
  const Zone_VEF& zone_VEF = la_zone_vef.valeur();
  const Champ_Inc_base& la_vitesse=vitesse();
  const DoubleTab& vitesse_face_absolue=la_vitesse.valeurs();
  const IntTab& elem_faces = zone_VEF.elem_faces();
  const DoubleTab& face_normales = zone_VEF.face_normales();
  const DoubleTab& facette_normales = zone_VEF.facette_normales();
  const Zone& zone = zone_VEF.zone();
  const Elem_VEF& type_elem = zone_VEF.type_elem();
  const int nb_elem_tot = zone_VEF.nb_elem_tot();
  const IntVect& rang_elem_non_std = zone_VEF.rang_elem_non_std();


  const DoubleVect& porosite_face = equation().milieu().porosite_face();
  const DoubleVect& porosite_elem = equation().milieu().porosite_elem();
  const DoubleTab& normales_facettes_Cl = zone_Cl_VEF.normales_facettes_Cl();

  int nfac = zone.nb_faces_elem();
  int nsom = zone.nb_som_elem();
  const IntTab& elem_sommets = zone.les_elems();
  //const DoubleTab& coord_sommets=zone.domaine().les_sommets();
  // Pour le traitement de la convection on distingue les polyedres
  // standard qui ne "voient" pas les conditions aux limites et les
  // polyedres non standard qui ont au moins une face sur le bord.
  // Un polyedre standard a n facettes sur lesquelles on applique le
  // schema de convection.
  // Pour un polyedre non standard qui porte des conditions aux limites
  // de Dirichlet, une partie des facettes sont portees par les faces.
  // En bref pour un polyedre le traitement de la convection depend
  // du type (triangle, tetraedre ...) et du nombre de faces de Dirichlet.

  double psc;
  DoubleTab pscl=0;
  int poly,face_adj,fa7,i,j,n_bord;
  int num_face, rang ,itypcl;
  int num10,num20,num_som;
  int ncomp_ch_transporte;
  if (transporte.nb_dim() == 1)
    ncomp_ch_transporte=1;
  else
    ncomp_ch_transporte= transporte.dimension(1);


  int marq=phi_u_transportant(equation());

  DoubleTab vitesse_face_;
  // soit on a transporte=phi*transporte_ et vitesse_face=vitesse_
  // soit transporte=transporte_ et vitesse_face=phi*vitesse_
  // cela depend si on transporte avec phi u ou avec u.
  const DoubleTab& vitesse_face=modif_par_porosite_si_flag(vitesse_face_absolue,vitesse_face_,marq,porosite_face);
  ArrOfInt face(nfac);
  ArrOfDouble vs(dimension);
  ArrOfDouble vc(dimension);
  DoubleTab vsom(nsom,dimension);
  ArrOfDouble cc(dimension);
  const Elem_VEF_base& type_elemvef= zone_VEF.type_elem().valeur();

  Nom nom_elem=type_elemvef.que_suis_je();



  // Les polyedres non standard sont ranges en 2 groupes dans la Zone_VEF:
  //  - polyedres bords et joints
  //  - polyedres bords et non joints
  // On traite les polyedres en suivant l'ordre dans lequel ils figurent
  // dans la zone

  // boucle sur les polys
  const IntTab& KEL=zone_VEF.type_elem().valeur().KEL();
  int phi_u_transportant_yes=phi_u_transportant(equation());

  {
    Cout<<"Op_Conv_RT_VEF_Face::ajouter_contribution RT\n";
    const DoubleVect& volumes = zone_VEF.volumes();
    const IntTab& face_voisins  = zone_VEF.face_voisins();
    const DoubleTab& coord_sommets=zone_VEF.zone().domaine().les_sommets();
    double volume=0.;
    int dim,elem,alfa,beta,k;
    if (dimension==2)
      {
        for (elem=0; elem<nb_elem_tot; elem++)
          {
            // Orientation des triangles
            DoubleTab coordSommet(dimension+1,dimension);
            for (i=0; i<=dimension; i++)
              {
                int numSom=elem_sommets(elem, i);
                for (dim=0; dim<dimension; dim++) coordSommet(i,dim)=coord_sommets(numSom,dim);
              }
#ifdef DEBUG
            double det=coordSommet(0,0)*(coordSommet(1,1)-coordSommet(2,1));
            det-=coordSommet(1,0)*(coordSommet(0,1)-coordSommet(2,1));
            det+=coordSommet(2,0)*(coordSommet(0,1)-coordSommet(1,1));
            if (det<0)
              {
                // Changer le sens des sommets du triangle
                Cout<<"Changer sens triangle " << elem << "\n";
                for (dim=0; dim<dimension; dim++)
                  {
                    double temp=coordSommet(1,dim);
                    coordSommet(1,dim)=coordSommet(2,dim);
                    coordSommet(2,dim)=temp;
                  }
              }
#endif
            DoubleTab FacesNormales(dimension+1,dimension);
            DoubleVect vitFaceNormale(dimension+1);
            for (j=0; j<=dimension; j++)
              {
                num_face=elem_faces(elem, j);
                for (dim=0; dim<dimension; dim++) FacesNormales(j,dim)=face_normales(num_face,dim);
                if (elem!=face_voisins(num_face,0))
                  {
                    for (dim=0; dim<dimension; dim++) FacesNormales(j,dim)=-FacesNormales(j,dim);
                  }
                vitFaceNormale(j)=0.;
                for (dim=0; dim<dimension; dim++) vitFaceNormale(j)+=vitesse_face_absolue(num_face,dim)*FacesNormales(j,dim);
              }
            DoubleTab FiaFib(dimension+1,dimension,dimension+1,dimension);

            for (i=0; i<=dimension; i++)
              {
                for (alfa=0; alfa<dimension; alfa++)
                  {
                    for (j=i+1; j<=dimension; j++)
                      {
                        for (beta=0; beta<dimension; beta++)
                          {
                            FiaFib(i,alfa,j,beta)=FacesNormales(i,alfa)*FacesNormales(j,beta);
                            FiaFib(j,beta,i,beta)=FiaFib(i,alfa,j,beta);
                          }
                      }
                  }
              }
            DoubleTab Fij(dimension+1,dimension+1);
            volume=volumes(elem);
            double invVol=1./(12*volume);
            double invVolvol=invVol/volume;
            for (i=0; i<=dimension; i++)
              {
                for (j=i+1; j<=dimension; j++)
                  {
                    Fij(i,j)= FiaFib(i,0,j,1)-FiaFib(i,1,j,0);
                    Fij(j,i)=-Fij(i,j);
                  } // end_for (j=i+1; j<=dimension; j++)
                Fij(i,i)=0.;
              } // end_for (i=0; i<=dimension; i++)
            // calcul qui ne depend que de i
            DoubleVect vitFaceNormaleFij(dimension+1);
            for (i=0; i<=dimension; i++)
              {
                vitFaceNormaleFij(i)=0.;
                for (k=0; k<=dimension; k++)
                  {
                    if (k!=i) vitFaceNormaleFij(i)+=vitFaceNormale(k)*Fij(k,i);
                  }
                vitFaceNormaleFij(i)*=invVolvol;
              }

            for (i=0; i<=dimension; i++)
              {
                int num_face_i=elem_faces(elem, i);
                int num_i = num_face_i*dimension;
                for (j=i; j<=dimension; j++)
                  {
                    int num_face_j=elem_faces(elem, j);
                    int num_j = num_face_j*dimension;

                    matrice(num_i,num_j)     -= FiaFib(i,0,j,0)*vitFaceNormaleFij(i);
                    matrice(num_i,num_j+1)   += FiaFib(i,0,j,1)*vitFaceNormaleFij(i);
                    matrice(num_i+1,num_j)   += FiaFib(i,1,j,0)*vitFaceNormaleFij(i);
                    matrice(num_i+1,num_j+1) -= FiaFib(i,1,j,1)*vitFaceNormaleFij(i);
                    if (i!=j)
                      {
                        matrice(num_j,num_i)     -= FiaFib(j,0,i,0)*vitFaceNormaleFij(j);
                        matrice(num_j,num_i+1)   += FiaFib(j,0,i,1)*vitFaceNormaleFij(j);
                        matrice(num_j+1,num_i)   += FiaFib(j,1,i,0)*vitFaceNormaleFij(j);
                        matrice(num_j+1,num_i+1) -= FiaFib(j,1,i,1)*vitFaceNormaleFij(j);
                      }
                  }
              }
            {
              DoubleTab resu_j(dimension+1,dimension);
              // NB: partie de code dupliquee ...
              DoubleTab SiMj(dimension+1,dimension+1,dimension);
              DoubleTab demi_coordSommet(dimension+1,dimension);
              for (i=0; i<=dimension; i++)
                {
                  for (dim=0; dim<dimension; dim++)
                    {
                      demi_coordSommet(i,dim)=0.5*coordSommet(i,dim);
                      SiMj(i,i,dim) = -coordSommet(i,dim);
                    }
                }
              for (dim=0; dim<dimension; dim++)
                {
                  SiMj(0,1,dim) =  demi_coordSommet(2,dim)-demi_coordSommet(0,dim);
                  SiMj(2,1,dim) = -SiMj(0,1,dim);
                  SiMj(0,2,dim) =  demi_coordSommet(1,dim)-demi_coordSommet(0,dim);
                  SiMj(1,2,dim) = -SiMj(0,2,dim);
                  SiMj(1,0,dim) =  demi_coordSommet(2,dim)-demi_coordSommet(1,dim);
                  SiMj(2,0,dim) = -SiMj(1,0,dim);
                  SiMj(0,0,dim) +=  demi_coordSommet(2,dim)+demi_coordSommet(1,dim);
                  SiMj(1,1,dim) +=  demi_coordSommet(2,dim)+demi_coordSommet(0,dim);
                  SiMj(2,2,dim) +=  demi_coordSommet(0,dim)+demi_coordSommet(1,dim);

                }
              //
              for (j=0; j<=dimension; j++)
                {
                  for (beta=0; beta<dimension; beta++)
                    resu_j(j,beta)=0.;
                  for (k=0; k<=dimension; k++)
                    {
                      for (dim=0; dim<dimension; dim++) resu_j(j,dim)+=vitFaceNormale(k)*SiMj(k,j,dim);
                    }
                  for (dim=0; dim<dimension; dim++) resu_j(j,dim)*=invVol;
                }
              for (i=0; i<=dimension; i++)
                {
                  int num_face_i=elem_faces(elem, i);
                  int num_i = num_face_i*dimension;
                  for (alfa=0; alfa<dimension; alfa++)
                    {
                      for (j=i+1; j<=dimension; j++)
                        {
                          int num_face_j=elem_faces(elem, j);
                          int num_j = num_face_j*dimension;
                          for (beta=0; beta<dimension; beta++)
                            {
                              matrice(num_i,num_j)-= FacesNormales(i,alfa)*resu_j(j,beta);
                              matrice(num_j,num_i)-= FacesNormales(j,beta)*resu_j(i,alfa);
                              num_j++;
                            }
                        }
                      num_i++;
                    }
                }
            }
          } // end_for (elem=0; elem<nb_elem_tot; elem++)
      }
    else if (dimension==3)
      {
        for (elem=0; elem<nb_elem_tot; elem++)
          {
            // Orientation des tetra
            DoubleTab coordSommet(dimension+1,dimension);
            for (i=0; i<=dimension; i++)
              {
                int numSom=elem_sommets(elem, i);
                for (dim=0; dim<dimension; dim++) coordSommet(i,dim)=coord_sommets(numSom,dim);
              }

#ifdef DEBUG
            DoubleVect S1S2vectS1S3(dimension);
            DoubleTab S1Si(dimension,dimension);
            for (i=0; i<dimension; i++)
              {
                for (dim=0; dim<dimension; dim++) S1Si(i,dim)=coordSommet(i+1,dim)-coordSommet(0,dim);
              }
            S1S2vectS1S3(0)=S1Si(0,1)*S1Si(1,2)-S1Si(0,2)*S1Si(1,1);
            S1S2vectS1S3(1)=S1Si(0,2)*S1Si(1,0)-S1Si(0,0)*S1Si(1,2);
            S1S2vectS1S3(2)=S1Si(0,0)*S1Si(1,1)-S1Si(0,1)*S1Si(1,0);
            double det=0.;
            for (dim=0; dim<dimension; dim++) det+=S1S2vectS1S3(dim)*S1Si(2,dim);

            if (det<0)
              {
                // Changer le sens des sommets du triangle
                Cout<<"Changer sens triangle " << elem << "\n";
                for (dim=0; dim<dimension; dim++)
                  {
                    double temp=coordSommet(2,dim);
                    coordSommet(2,dim)=coordSommet(3,dim);
                    coordSommet(3,dim)=temp;
                  }
              }
#endif
            DoubleTab FacesNormales(dimension+1,dimension);
            DoubleVect vitFaceNormale(dimension+1);
            DoubleTab FiFj(dimension+1,dimension+1,dimension);
            IntTab ijkl(dimension+1,dimension+1,2);
            for (i=0; i<=dimension; i++)
              {
                num_face=elem_faces(elem, i);
                for (dim=0; dim<dimension; dim++) FacesNormales(i,dim)=face_normales(num_face,dim);
                if (elem!=face_voisins(num_face,0))
                  {
                    for (dim=0; dim<dimension; dim++) FacesNormales(i,dim)=-FacesNormales(i,dim);
                  }
                vitFaceNormale(i)=0.;
                for (dim=0; dim<dimension; dim++) vitFaceNormale(i) += vitesse_face_absolue(num_face,dim)*FacesNormales(i,dim);

                for (j=i+1; j<=dimension; j++)
                  {
                    FiFj(i,j,0)=FacesNormales(i,1)*FacesNormales(j,2)-FacesNormales(i,2)*FacesNormales(j,1);
                    FiFj(i,j,1)=FacesNormales(i,2)*FacesNormales(j,0)-FacesNormales(i,0)*FacesNormales(j,2);
                    FiFj(i,j,2)=FacesNormales(i,0)*FacesNormales(j,1)-FacesNormales(i,1)*FacesNormales(j,0);
                    for (dim=0; dim<dimension; dim++) FiFj(j,i,dim)=-FiFj(i,j,dim);
                  }
              }
            ijkl(0,1,0)=2;
            ijkl(0,1,1)=3;
            ijkl(0,2,0)=3;
            ijkl(0,2,1)=1;
            ijkl(0,3,0)=1;
            ijkl(0,3,1)=2;
            ijkl(1,2,0)=0;
            ijkl(1,2,1)=3;
            ijkl(1,3,0)=2;
            ijkl(1,3,1)=0;
            ijkl(2,3,0)=0;
            ijkl(2,3,1)=1;
            /*
            for (i=0; i<=dimension; i++) {
              ijkl(i,i,0)=-1; ijkl(i,i,1)=-1;
              for(j=i+1;j<=dimension;j++) {
            	  ijkl(j,i,0)=ijkl(i,j,1); ijkl(j,i,1)=ijkl(i,j,0);
              }
            }*/
            DoubleTab FkiFj(dimension+1,dimension+1,dimension+1,dimension);
            for (k=0; k<=dimension; k++)
              {
                for (i=k+1; i<=dimension; i++)
                  {
                    int kk=ijkl(k,i,0);
                    int ii=ijkl(k,i,1);
                    for (dim=0; dim<dimension; dim++)
                      {
                        FkiFj(k,i,k,dim)=FiFj(kk,k,dim)-FiFj(ii,k,dim);
                        FkiFj(i,k,k,dim)=-FkiFj(k,i,k,dim);
                        FkiFj(k,i,i,dim)=FiFj(kk,i,dim)-FiFj(ii,i,dim);
                        FkiFj(i,k,i,dim)=-FkiFj(k,i,i,dim);
                      }
                  }
              }

            volume=volumes(elem);
            double invVol=1./(18.*volume*volume);
            for (i=0; i<=dimension; i++)
              {
                int num_face_i=elem_faces(elem, i);
                int num_i=num_face_i*dimension;
                IntTab num_ia(dimension);
                for (dim=0; dim<dimension; dim++) num_ia(dim)=num_i+dim;
                for (j=0; j<=dimension; j++)
                  {
                    int num_face_j=elem_faces(elem, j);
                    int num_j = num_face_j*dimension;
                    for (beta=0; beta<dimension; beta++)
                      {
                        double resu_loc=0.;
                        for (k=0; k<=dimension; k++)
                          {
                            if (k!=i) resu_loc+=vitFaceNormale(k)*FkiFj(k,i,j,beta);
                          }
                        for (alfa=0; alfa<dimension; alfa++) matrice(num_ia(alfa),num_j)+=resu_loc*FacesNormales(i,alfa);
                        num_j++;
                      }
                  }
              }


            {
              // Schema d'Hammer-Stroud pour integration de polynomes d'ordre 2
              int nbpt=4;
              double c0= 0.13819660112501051518; // (5-sqrt(5))/20
              double c1= 0.58541019662496845446; // 1-3*c0 ou (5+3*sqrt(5))/20
              double c2=-0.75623058987490536338; // 1-3*c1 ou (5-9*sqrt(5))/20
              // lambda des points (coordonnees barycentriques)
              DoubleTab lambda(nbpt,dimension+1);
              lambda(0,0)=lambda(0,1)=lambda(0,2)=c0;
              lambda(0,3)=c1;
              lambda(1,0)=lambda(1,1)=lambda(1,3)=c0;
              lambda(1,2)=c1;
              lambda(2,0)=lambda(2,2)=lambda(2,3)=c0;
              lambda(2,1)=c1;
              lambda(3,2)=lambda(3,1)=lambda(3,3)=c0;
              lambda(3,0)=c1;
              // psi(p,k)=1-3*lambda(p,k)   // valeur fonction de base CR en OXp
              DoubleTab psi(nbpt,dimension+1);
              psi(0,0)=psi(0,1)=psi(0,2)=c1;
              psi(0,3)=c2;
              psi(1,0)=psi(1,1)=psi(1,3)=c1;
              psi(1,2)=c2;
              psi(2,0)=psi(2,2)=psi(2,3)=c1;
              psi(2,1)=c2;
              psi(3,2)=psi(3,1)=psi(3,3)=c1;
              psi(3,0)=c2;
              int pt,gamma;
              // les points d'integration
              DoubleTab OXp(nbpt,dimension);
              for (pt=0; pt<nbpt; pt++)
                {
                  for (beta=0; beta<dimension; beta++)
                    {
                      OXp(pt,beta)=0;
                      for (j=0; j<=dimension; j++) OXp(pt,beta)+=lambda(pt,j)*coordSommet(j,beta);
                    }
                }
              invVol=1./(24.*volume);
              // les vecteurs SkXp
              DoubleTab SkXp(dimension+1,nbpt,dimension);
              for (k=0; k<=dimension; k++)
                {
                  for (gamma=0; gamma<dimension; gamma++)
                    {
                      for (pt=0; pt<nbpt; pt++) SkXp(k,pt,gamma)=OXp(pt,gamma)-coordSommet(k,gamma);
                    }
                }
              for (i=0; i<=dimension; i++)
                {
                  int num_face_i=elem_faces(elem, i);
                  int num_i=num_face_i*dimension;
                  for (alfa=0; alfa<dimension; alfa++)
                    {
                      for (j=0; j<=dimension; j++)
                        {
                          int num_face_j=elem_faces(elem,j);
                          int num_j=num_face_j*dimension;
                          for (beta=0; beta<dimension; beta++)
                            {
                              double resu_j_beta=0.;
                              for (k=0; k<=dimension; k++)
                                {
                                  double resu_k=0.;
                                  // boucle sur les points d'integration
                                  for (pt=0; pt<nbpt; pt++) resu_k+=SkXp(k,pt,beta)*psi(pt,j);
                                  resu_k*=vitFaceNormale(k);
                                  resu_j_beta+=resu_k;
                                } // end_for (k=0; k<=dimension; k++)
                              matrice(num_i,num_j)-=invVol*FacesNormales(i,alfa)*resu_j_beta ;
                              num_j++;
                            } // end_for (beta=0; beta<dimension; beta++)
                        } // end_for (j=0; j<=dimension; j++)
                      num_i++;
                    } // end_for (alfa=0; alfa<dimension;alfa++)
                } // end_for (i=0; i<=dimension;i++)
            } // end_if (type_op==RT)
          } // end_for (int elem=0; elem<nb_elem_tot; elem++)
      }
  }

  // Boucle sur les bords pour traiter les conditions aux limites
  // il y a prise en compte d'un terme de convection pour les
  // conditions aux limites de Neumann_sortie_libre seulement
  int nb_bord = zone_VEF.nb_front_Cl();
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {

      const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);

      if (sub_type(Neumann_sortie_libre,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (num_face=num1; num_face<num2; num_face++)
            {
              psc =0;
              for (i=0; i<dimension; i++)
                psc += vitesse_face(num_face,i)*face_normales(num_face,i);
              if (!phi_u_transportant_yes)
                psc*=porosite_face(num_face);
              if (psc>0)
                {
                  for (j=0; j<ncomp_ch_transporte; j++)
                    {
                      int n0=num_face*ncomp_ch_transporte+j;
                      matrice(n0,n0)+=psc;
                    }
                }
            }
        }
    }
  modifier_matrice_pour_periodique_apres_contribuer(matrice,equation());
}

