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
// File:        Modele_Launder_Sharma_VEF.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Modeles_Turbulence/RANS/Fonc
//
//////////////////////////////////////////////////////////////////////////////

#include <Modele_Launder_Sharma_VEF.h>
#include <Domaine_VEF.h>
#include <Champ_Uniforme.h>
#include <Champ_P1NC.h>

Implemente_instanciable(Modele_Launder_Sharma_VEF,"Modele_Launder_Sharma_VEF",Modele_Lam_Bremhorst_VEF);



///////////////////////////////////////////////////////////////
//   Implementation des fonctions de la classe
///////////////////////////////////////////////////////////////
// printOn et readOn

Sortie& Modele_Launder_Sharma_VEF::printOn(Sortie& s ) const
{
  return s;
}

Entree& Modele_Launder_Sharma_VEF::readOn(Entree& is )
{
  Param param(que_suis_je());
  set_param(param);
  param.lire_avec_accolades_depuis(is);
  if (is_Reynolds_stress_isotrope_ == 0)
    is_Cmu_constant_ = 0;
  return is;
}

void Modele_Launder_Sharma_VEF::lire_distance_paroi( )
{

// Pas besoin de la distance a la proi dans ce modele

}


DoubleTab& Modele_Launder_Sharma_VEF::Calcul_D(DoubleTab& D,const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis,
                                               const DoubleTab& vitesse,const DoubleTab& K_eps_Bas_Re, const Champ_Don& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco->valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco.valeur());
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,domaine_dis);
  //  const Domaine_Cl_VEF& le_dom_Cl = ref_cast(Domaine_Cl_VEF,domaine_Cl_dis);
  const DoubleVect& volumes = le_dom.volumes();
  //  int nb_faces = le_dom.nb_faces();
  int nb_faces_tot = le_dom.nb_faces();
  //  int nb_elem = le_dom.nb_elem();
  int nb_elem_tot = le_dom.nb_elem_tot();
  //  int nb_elem_tot = le_dom.nb_elem_tot();
  const Domaine& domaine=le_dom.domaine();
  const IntTab& elem_faces = le_dom.elem_faces();
  const int nb_faces_elem = domaine.nb_faces_elem();
  const DoubleTab& face_normales = le_dom.face_normales();
  const DoubleVect& vol_ent = le_dom.volumes_entrelaces();

  D = 0;
  //return D;
  int num_elem,i,face_loc,face_glob,num_face;
  double deriv,Vol;

  // Algo :
  //   * Boucle sur les elements
  //      * boucle locale dans l'element sur les faces -> calcul du gradient de racine de k
  //      * distribution de la valeur de l'integrale sur les faces de l'element
  //  Cela doit etre OK!!!
  // Rq : k et epsilon sont definis comme la vitesse P1NC.

  // ET le // ????? nb_elem ou nb_elem_tot
  // Pbls des C.L. en paroi???? non car on impose epsilon-D = 0!!
  for (num_elem=0; num_elem<nb_elem_tot; num_elem++)
    {
      Vol = volumes(num_elem);
      for(i=0; i<dimension; i++)
        {
          deriv = 0;
          for(face_loc=0; face_loc<nb_faces_elem; face_loc++)
            {
              face_glob = elem_faces(num_elem,face_loc);
              deriv += sqrt(K_eps_Bas_Re(face_glob,0))*face_normales(face_glob,i)*le_dom.oriente_normale(face_glob,num_elem);
            }
          deriv /= Vol;
          for(face_loc=0; face_loc<nb_faces_elem; face_loc++)
            {
              face_glob = elem_faces(num_elem,face_loc);
              /*              Cerr<<"face_glob = "<<face_glob<<finl;
                              Cerr<<"K_eps_Bas_Re(face_glob,0) = "<<K_eps_Bas_Re(face_glob,0)<<finl;
                              Cerr<<"deriv = "<<deriv<<", Vol = "<<Vol<<", nb_faces_elem = "<<nb_faces_elem<<finl;*/
              if (!is_visco_const)
                visco=tab_visco(num_elem);
              // on divise par nb_faces_elem pour dire que chaque elem
              // contribue pour 1/3 de son volume en 2D.
              D(face_glob) += visco*deriv*deriv*Vol/nb_faces_elem;
            }
        }
    }
  // on divise par vol_ent car en 2d c 2/3 du volume. Donc en tout en 2D
  // chaque face recoit en tout 2 * (1/3) * vol / (vol * (2/3))=1 !
  for (num_face=0; num_face<nb_faces_tot; num_face++)
    D(num_face) *= 2./vol_ent(num_face);

  // RQ : il faut diviser par le volume entrelace car ce que nous calculons c est l integrale en V
  // Par contre dans les termes sources on multiplie par le volume entrelace donc pour rester coherent avec le reste....
  //D=0;
  //  D*=0.5;
  return D;
}


// Fonction utile visc
// Pour le calcul de la derivee seconde
// 1/|v|*(Si,elem)ind_der*(Sj,elem)ind_der
/* pas appele
   static double viscA(const Domaine_VEF& le_dom,int num_face,int num2,int dimension,
   int num_elem, int ind_der)
   {
   double pscal;

   pscal = le_dom.face_normales(num_face,ind_der)*le_dom.face_normales(num2,ind_der);

   if ( (le_dom.face_voisins(num_face,0) == le_dom.face_voisins(num2,0)) ||
   (le_dom.face_voisins(num_face,1) == le_dom.face_voisins(num2,1)))
   return -pscal/le_dom.volumes(num_elem);
   else
   return pscal/le_dom.volumes(num_elem);
   }
*/

DoubleTab& Modele_Launder_Sharma_VEF::Calcul_E(DoubleTab& E,const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis, const DoubleTab& transporte,const DoubleTab& K_eps_Bas_Re,const Champ_Don& ch_visco, const DoubleTab& visco_turb ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco->valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco.valeur());
  if (is_visco_const)
    visco=tab_visco(0,0);
  //const Domaine_Cl_VEF& domaine_Cl_VEF = ref_cast(Domaine_Cl_VEF,domaine_Cl_dis);
  const Domaine_VEF& domaine_VEF =  ref_cast(Domaine_VEF,domaine_dis);
  E = 0;
  //return E;
  //  const IntTab& elem_faces = domaine_VEF.elem_faces();
  const IntTab& face_voisins = domaine_VEF.face_voisins();
  const DoubleVect& volumes = domaine_VEF.volumes();
  //const DoubleTab& face_normales = domaine_VEF.face_normales();
  //  const DoubleTab& facette_normales = domaine_VEF.facette_normales();
  //  const DoubleVect& volumes_entrelaces = domaine_VEF.volumes_entrelaces();
  //  const Domaine& domaine = domaine_VEF.domaine();
  const int nb_faces = domaine_VEF.nb_faces();
  //  int nb_faces_tot = domaine_VEF.nb_faces_tot();
  //  const int nfa7 = domaine_VEF.type_elem().nb_facette();
  const int nb_elem = domaine_VEF.nb_elem();
  //  const IntVect& rang_elem_non_std = domaine_VEF.rang_elem_non_std();
  //int premiere_face_int = domaine_VEF.premiere_face_int();
  int premiere_face_int = domaine_VEF.premiere_face_int();

  int i,j;
  int elem,elem0,elem2;
  //int n_bord;//alpha,beta,elem0;
  //  int nb_faces_elem = domaine_VEF.domaine().nb_faces_elem();

  //  double val;

  const int ncomp_ch_transporte = transporte.line_size();

  DoubleTab gradient(0, ncomp_ch_transporte, dimension);
  domaine_VEF.creer_tableau_faces(gradient);

  DoubleTab deriv_seconde_de_gradient_elem(0, ncomp_ch_transporte, dimension, dimension);
  domaine_VEF.domaine().creer_tableau_elements(deriv_seconde_de_gradient_elem);
  deriv_seconde_de_gradient_elem=0.;

  DoubleTab champ_face(0,dimension);
  domaine_VEF.creer_tableau_faces(champ_face);

  // Calcul de la derivee premiere de U
  // stokage d'un gradient_elem par element.
  DoubleTab gradient_elem(0, ncomp_ch_transporte, dimension);
  domaine_VEF.domaine().creer_tableau_elements(gradient_elem);

  DoubleTab champ_elem(0, ncomp_ch_transporte, dimension);
  domaine_VEF.domaine().creer_tableau_elements(champ_elem);

  DoubleTab deriv_seconde_de_gradient(0, ncomp_ch_transporte, dimension, dimension);
  domaine_VEF.creer_tableau_faces(deriv_seconde_de_gradient);
  deriv_seconde_de_gradient=0.;

  int fac=-1,elem1,comp;
  // Rque methode non const Pourquoi ?

  const Champ_P1NC& vitesse = ref_cast(Champ_P1NC,eq_hydraulique->inconnue().valeur());
  ref_cast_non_const(Champ_P1NC,vitesse).calcul_gradient(transporte,gradient_elem,ref_cast(Domaine_Cl_VEF,eq_hydraulique->domaine_Cl_dis()));

  gradient_elem.echange_espace_virtuel();
  // On a les gradient_elem par elements
  // Boucle sur les faces
  //
  for (fac=0; fac< premiere_face_int; fac++)
    for (comp=0; comp<ncomp_ch_transporte; comp++)
      for (i=0; i<dimension; i++)
        {
          elem1=face_voisins(fac,0);
          if (elem1 != -1)
            gradient(fac, comp, i) = gradient_elem(elem1, comp, i);
          else
            {
              elem1=face_voisins(fac,1);
              gradient(fac, comp, i) = gradient_elem(elem1, comp, i);
            }
        }
// fin du for faces

  for (; fac<nb_faces; fac++)
    {
      elem1=face_voisins(fac,0);
      elem2=face_voisins(fac,1);
      double vol1=volumes(elem1);
      double vol2=volumes(elem2);
      double voltot=vol1+vol2;
      for (comp=0; comp<ncomp_ch_transporte; comp++)
        for (i=0; i<dimension; i++)
          {
            double grad1=gradient_elem(elem1, comp, i);
            double grad2=gradient_elem(elem2, comp, i);
            gradient(fac, comp, i) =  (vol1*grad1 + vol2*grad2)/voltot;
          }
    } // fin du for faces

  gradient.echange_espace_virtuel();
//on a les gradients sur les faces.


  // maintenant on calcul la derivee seconde
  // stokage de la derivee seconde par element.
  for (i=0; i<dimension; i++)
    {
      for (comp=0; comp<ncomp_ch_transporte; comp++)
        for (fac=0; fac<nb_faces; fac++)
          champ_face(fac,comp)=gradient(fac, comp, i);

      ref_cast_non_const(Champ_P1NC,vitesse).calcul_gradient(champ_face,champ_elem,ref_cast(Domaine_Cl_VEF,eq_hydraulique->domaine_Cl_dis()));

      for (comp=0; comp<ncomp_ch_transporte; comp++)
        for ( elem=0; elem<nb_elem; elem++)
          for (j=0; j<dimension; j++)
            deriv_seconde_de_gradient_elem(elem,comp,i,j) = champ_elem(elem,comp,j);

    }

  deriv_seconde_de_gradient_elem.echange_espace_virtuel();
  // On a les gradient_elem par elements

  // Boucle sur les faces
  //

  for (fac=0; fac< premiere_face_int; fac++)
    {
      for (comp=0; comp<ncomp_ch_transporte; comp++)
        for (j=0; j<dimension; j++)
          for (i=0; i<dimension; i++)
            {
              elem1=face_voisins(fac,0);
              if (elem1 != -1)
                deriv_seconde_de_gradient(fac, comp, i, j) = deriv_seconde_de_gradient_elem(elem1, comp, i, j);
              else
                {
                  elem1=face_voisins(fac,1);
                  deriv_seconde_de_gradient(fac, comp, i, j) = deriv_seconde_de_gradient_elem(elem1, comp, i, j);
                }
            }

    } // fin du for faces

  for (; fac<nb_faces; fac++)
    {
      elem1=face_voisins(fac,0);
      elem2=face_voisins(fac,1);
      double vol1=volumes(elem1);
      double vol2=volumes(elem2);
      double voltot=vol1+vol2;
      for (comp=0; comp<ncomp_ch_transporte; comp++)
        for (j=0; j<dimension; j++)
          for (i=0; i<dimension; i++)
            {
              double grad1=deriv_seconde_de_gradient_elem(elem1, comp, i, j);
              double grad2=deriv_seconde_de_gradient_elem(elem2, comp, i, j);
              deriv_seconde_de_gradient(fac, comp, i, j) = (vol1*grad1 + vol2*grad2)/voltot;
            }
    } // fin du for faces

  // Calcul de E
  double deriv,nuturb;
  for (fac=0; fac<nb_faces; fac++)
    {
      deriv = 0;
      if (dimension == 2)
        {
          double val2 = deriv_seconde_de_gradient(fac,0,1,1)*deriv_seconde_de_gradient(fac,0,1,1);
          double val3 = deriv_seconde_de_gradient(fac,1,0,0)*deriv_seconde_de_gradient(fac,1,0,0);
          deriv = val2 + val3;
          //deriv = gradient(fac, 0, 1);
        }
      else
        {
          double val2 = deriv_seconde_de_gradient(fac,1,0,0)+deriv_seconde_de_gradient(fac,2,0,0); //d2v_dx2(fac)+d2w_dx2(fac);
          val2 *= val2;
          double val3 = deriv_seconde_de_gradient(fac,0,1,1)+deriv_seconde_de_gradient(fac,2,1,1); //d2u_dy2(fac)+d2w_dy2(fac);
          val3 *= val3;
          double val4 = deriv_seconde_de_gradient(fac,0,2,2)+deriv_seconde_de_gradient(fac,1,2,2); //d2u_dz2(fac)+d2v_dz2(fac);
          val4 *= val4;
          deriv = val2 + val3 + val4;

        }
      elem0 = face_voisins(fac,0);
      elem1 = face_voisins(fac,1);

      if (elem1!=-1)
        {
          nuturb = visco_turb(elem0)*volumes(elem0)+visco_turb(elem1)*volumes(elem1);
          nuturb /= volumes(elem0) + volumes(elem1);
          if (!is_visco_const)
            visco=(tab_visco(elem0)+tab_visco(elem1))/2.;
        }
      else
        {
          nuturb =  visco_turb(elem0);
          if (!is_visco_const)
            visco=(tab_visco(elem0));
        }

      E(fac) = 2.*visco*deriv*nuturb;

    }
  //E=0;
//  Cerr<<E.mp_min_vect()<<" EEEEEEEEEEEEEEEEEEEEEEEE "<<E.mp_max_vect()<<finl;

  return E;
}

DoubleTab& Modele_Launder_Sharma_VEF::Calcul_F1( DoubleTab& F1, const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis, const DoubleTab& P, const DoubleTab& K_eps_Bas_Re,const Champ_base& ch_visco) const
{
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,domaine_dis);
  int nb_faces = le_dom.nb_faces();
  for (int num_face=0; num_face <nb_faces; num_face ++ )
    F1[num_face] = 1.;
  return F1;
}

DoubleTab& Modele_Launder_Sharma_VEF::Calcul_F2( DoubleTab& F2, DoubleTab& Deb, const Domaine_dis_base& domaine_dis,const DoubleTab& K_eps_Bas_Re,const Champ_base& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco);
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,domaine_dis);
  int nb_faces = le_dom.nb_faces();
  int num_face;
  double Re;

  for (num_face=0; num_face<nb_faces  ; num_face++)
    {
      if (K_eps_Bas_Re(num_face,1)>0)
        {
          if (!is_visco_const)
            {
              int elem0 = le_dom.face_voisins(num_face,0);
              int elem1 = le_dom.face_voisins(num_face,1);
              if (elem1!=-1)
                {
                  visco = tab_visco(elem0)*le_dom.volumes(elem0)+tab_visco(elem1)*le_dom.volumes(elem1);
                  visco /= le_dom.volumes(elem0) + le_dom.volumes(elem1);
                }
              else
                visco =  tab_visco(elem0);

            }
          if (visco>0)
            {
              Re = (K_eps_Bas_Re(num_face,0)*K_eps_Bas_Re(num_face,0))/(visco*K_eps_Bas_Re(num_face,1));
              F2[num_face] = 1. - (0.3*exp(-Re*Re));
            }
        }
      else
        {
          F2[num_face] = 1.;
        }
    }
  //Cerr<<F2.mp_min_vect()<<" F2 "<<F2.mp_max_vect()<<finl;
  return F2;
}

DoubleTab&  Modele_Launder_Sharma_VEF::Calcul_Fmu( DoubleTab& Fmu,const Domaine_dis_base& domaine_dis,const Domaine_Cl_dis_base& domaine_Cl_dis,const DoubleTab& K_eps_Bas_Re,const Champ_Don& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco->valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco.valeur());
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,domaine_dis);
  int nb_faces = le_dom.nb_faces();
  int num_face;
  double Re;
  Fmu = 0;

  for (num_face=0; num_face<nb_faces  ; num_face++)
    {
      if(1) // if (K_eps_Bas_Re(num_face,1)>0)
        {
          if (!is_visco_const)
            {
              int elem0 = le_dom.face_voisins(num_face,0);
              int elem1 = le_dom.face_voisins(num_face,1);
              if (elem1!=-1)
                {
                  visco = tab_visco(elem0)*le_dom.volumes(elem0)+tab_visco(elem1)*le_dom.volumes(elem1);
                  visco /= le_dom.volumes(elem0) + le_dom.volumes(elem1);
                }
              else
                visco =  tab_visco(elem0);
            }
          Re = (K_eps_Bas_Re(num_face,0)*K_eps_Bas_Re(num_face,0))/(visco*(K_eps_Bas_Re(num_face,1)+DMINFLOAT));
          Fmu[num_face] = exp(-3.4/((1.+Re/50.)*(1.+Re/50.)));

        }
      else
        {
          Fmu[num_face] = 1.;
        }
    }
//  Cerr<<Fmu.mp_min_vect()<<" Fmu " <<Fmu.mp_max_vect()<<finl;
  return Fmu;
}


DoubleTab&  Modele_Launder_Sharma_VEF::Calcul_Fmu_BiK( DoubleTab& Fmu,const Domaine_dis_base& domaine_dis,const Domaine_Cl_dis_base& domaine_Cl_dis,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re,const Champ_Don& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco->valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco.valeur());
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,domaine_dis);
  int nb_faces = le_dom.nb_faces();
  int num_face;
  double Re;
  Fmu = 0;

  for (num_face=0; num_face<nb_faces  ; num_face++)
    {
      if(1) // if (K_eps_Bas_Re(num_face,1)>0)
        {
          if (!is_visco_const)
            {
              int elem0 = le_dom.face_voisins(num_face,0);
              int elem1 = le_dom.face_voisins(num_face,1);
              if (elem1!=-1)
                {
                  visco = tab_visco(elem0)*le_dom.volumes(elem0)+tab_visco(elem1)*le_dom.volumes(elem1);
                  visco /= le_dom.volumes(elem0) + le_dom.volumes(elem1);
                }
              else
                visco =  tab_visco(elem0);
            }
          Re = (K_Bas_Re(num_face)*K_Bas_Re(num_face))/(visco*(eps_Bas_Re(num_face)+DMINFLOAT));
          Fmu[num_face] = exp(-3.4/((1.+Re/50.)*(1.+Re/50.)));

        }
      else
        {
          Fmu[num_face] = 1.;
        }
    }
  Cerr<<Fmu.mp_min_vect()<<" Fmu " <<Fmu.mp_max_vect()<<finl;
  return Fmu;
}


DoubleTab& Modele_Launder_Sharma_VEF::Calcul_F2_BiK( DoubleTab& F2, DoubleTab& Deb, const Domaine_dis_base& domaine_dis,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re,const Champ_base& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco.valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco);
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,domaine_dis);
  int nb_faces = le_dom.nb_faces();
  int num_face;
  double Re;

  for (num_face=0; num_face<nb_faces  ; num_face++)
    {
      if (eps_Bas_Re(num_face)>0)
        {
          if (!is_visco_const)
            {
              int elem0 = le_dom.face_voisins(num_face,0);
              int elem1 = le_dom.face_voisins(num_face,1);
              if (elem1!=-1)
                {
                  visco = tab_visco(elem0)*le_dom.volumes(elem0)+tab_visco(elem1)*le_dom.volumes(elem1);
                  visco /= le_dom.volumes(elem0) + le_dom.volumes(elem1);
                }
              else
                visco =  tab_visco(elem0);

            }
          if (visco>0)
            {
              Re = (K_Bas_Re(num_face)*K_Bas_Re(num_face))/(visco*eps_Bas_Re(num_face));
              F2[num_face] = 1. - (0.3*exp(-Re*Re));
            }
        }
      else
        {
          F2[num_face] = 1.;
        }
    }
  //Cerr<<F2.mp_min_vect()<<" F2 "<<F2.mp_max_vect()<<finl;
  return F2;
}



DoubleTab& Modele_Launder_Sharma_VEF::Calcul_F1_BiK( DoubleTab& F1, const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis, const DoubleTab& P, const DoubleTab& K_Bas_Re, const DoubleTab& eps_Bas_Re,const Champ_base& ch_visco) const
{
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,domaine_dis);
  int nb_faces = le_dom.nb_faces();
  for (int num_face=0; num_face <nb_faces; num_face ++ )
    F1[num_face] = 1.;
  return F1;
}


DoubleTab& Modele_Launder_Sharma_VEF::Calcul_E_BiK(DoubleTab& E,const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis, const DoubleTab& transporte,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re,const Champ_Don& ch_visco, const DoubleTab& visco_turb ) const
{
  return Calcul_E( E,domaine_dis,domaine_Cl_dis,transporte,K_Bas_Re,ch_visco,visco_turb );
}

DoubleTab& Modele_Launder_Sharma_VEF::Calcul_D_BiK(DoubleTab& D,const Domaine_dis_base& domaine_dis, const Domaine_Cl_dis_base& domaine_Cl_dis,
                                                   const DoubleTab& vitesse,const DoubleTab& K_Bas_Re,const DoubleTab& eps_Bas_Re, const Champ_Don& ch_visco ) const
{
  double visco=-1;
  const DoubleTab& tab_visco=ch_visco->valeurs();
  int is_visco_const=sub_type(Champ_Uniforme,ch_visco.valeur());
  if (is_visco_const)
    visco=tab_visco(0,0);
  const Domaine_VEF& le_dom = ref_cast(Domaine_VEF,domaine_dis);
  //  const Domaine_Cl_VEF& le_dom_Cl = ref_cast(Domaine_Cl_VEF,domaine_Cl_dis);
  const DoubleVect& volumes = le_dom.volumes();
  //  int nb_faces = le_dom.nb_faces();
  int nb_faces_tot = le_dom.nb_faces();
  //  int nb_elem = le_dom.nb_elem();
  int nb_elem_tot = le_dom.nb_elem_tot();
  //  int nb_elem_tot = le_dom.nb_elem_tot();
  const Domaine& domaine=le_dom.domaine();
  const IntTab& elem_faces = le_dom.elem_faces();
  const int nb_faces_elem = domaine.nb_faces_elem();
  const DoubleTab& face_normales = le_dom.face_normales();
  const DoubleVect& vol_ent = le_dom.volumes_entrelaces();

  D = 0;
  //return D;
  int num_elem,i,face_loc,face_glob,num_face;
  double deriv,Vol;

  // Algo :
  //   * Boucle sur les elements
  //      * boucle locale dans l'element sur les faces -> calcul du gradient de racine de k
  //      * distribution de la valeur de l'integrale sur les faces de l'element
  //  Cela doit etre OK!!!
  // Rq : k et epsilon sont definis comme la vitesse P1NC.

  // ET le // ????? nb_elem ou nb_elem_tot
  // Pbls des C.L. en paroi???? non car on impose epsilon-D = 0!!
  for (num_elem=0; num_elem<nb_elem_tot; num_elem++)
    {
      Vol = volumes(num_elem);
      for(i=0; i<dimension; i++)
        {
          deriv = 0;
          for(face_loc=0; face_loc<nb_faces_elem; face_loc++)
            {
              face_glob = elem_faces(num_elem,face_loc);
              deriv += sqrt(K_Bas_Re(face_glob))*face_normales(face_glob,i)*le_dom.oriente_normale(face_glob,num_elem);
            }
          deriv /= Vol;
          for(face_loc=0; face_loc<nb_faces_elem; face_loc++)
            {
              face_glob = elem_faces(num_elem,face_loc);
              /*              Cerr<<"face_glob = "<<face_glob<<finl;
                              Cerr<<"K_eps_Bas_Re(face_glob,0) = "<<K_eps_Bas_Re(face_glob,0)<<finl;
                              Cerr<<"deriv = "<<deriv<<", Vol = "<<Vol<<", nb_faces_elem = "<<nb_faces_elem<<finl;*/
              if (!is_visco_const)
                visco=tab_visco(num_elem);
              // on divise par nb_faces_elem pour dire que chaque elem
              // contribue pour 1/3 de son volume en 2D.
              D(face_glob) += visco*deriv*deriv*Vol/nb_faces_elem;
            }
        }
    }
  // on divise par vol_ent car en 2d c 2/3 du volume. Donc en tout en 2D
  // chaque face recoit en tout 2 * (1/3) * vol / (vol * (2/3))=1 !
  for (num_face=0; num_face<nb_faces_tot; num_face++)
    D(num_face) *= 2./vol_ent(num_face);

  // RQ : il faut diviser par le volume entrelace car ce que nous calculons c est l integrale en V
  // Par contre dans les termes sources on multiplie par le volume entrelace donc pour rester coherent avec le reste....
  //D=0;
  //  D*=0.5;
  return D;
}

