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
#include <Scatter.h>
#include <MD_Vector_tools.h>
#include <MD_Vector_std.h>
#include <Debog.h>
#include <Domaine_ALE.h>
#include <Probleme_base.h>
#include <Zone_VDF.h>
#include <Zone_VEF.h>
#include <Motcle.h>
#include <EFichier.h>
#include <Champ_front_ALE.h>
#include <Ch_front_input_ALE.h>
#include <Champ_front_ALE_Beam.h>

Implemente_instanciable_sans_constructeur(Domaine_ALE,"Domaine_ALE",Domaine);
//XD domaine_ale domaine domaine_ale -1 Domain with nodes at the interior of the domain which are displaced in an arbitrarily prescribed way thanks to ALE (Arbitrary Lagrangian-Eulerian) description. NL2 Keyword to specify that the domain is mobile following the displacement of some of its boundaries.
Domaine_ALE::Domaine_ALE() : nb_bords_ALE(0),update_or_not_matrix_coeffs_(1)
{

}
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
  int N_som=nb_som_tot();

  update_or_not_matrix_coeffs_=0;

  //Cerr<<"Domaine_ALE::mettre_a_jour"<<finl;

  bool  check_NoZero_ALE=true; //  if ALE boundary velocity is zero, ALE mesh velocity is zero and mettre_a_jour is skipped.

  ALE_mesh_velocity = calculer_vitesse(temps,le_domaine_dis,pb, check_NoZero_ALE);

  if(check_NoZero_ALE)
    {
      for (int i=0; i<N_som; i++)
        {
          for (int k=0; k<dimension; k++)
            coord(i,k)+=ALE_mesh_velocity(i,k)*dt_;
        }

      //On recalcule les vitesses aux faces
      Zone_VF& la_zone_VF=ref_cast(Zone_VF,le_domaine_dis.zone_dis(0).valeur());
      int nb_faces=la_zone_VF.nb_faces();
      int nb_som_face=la_zone_VF.nb_som_face();
      IntTab& face_sommets=la_zone_VF.face_sommets();
      creer_mes_domaines_frontieres(la_zone_VF);//update the boundary surface domain
      calculer_vitesse_faces(ALE_mesh_velocity,nb_faces,nb_som_face,face_sommets);

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


          // PL: je trouve etonnant que le calcul des surfaces se fasse AVANT le calcul des normales
          // PL: il devrait se faire APRES

          int nb_faces_tot=face_sommets.dimension_tot(0);
          la_zone_VEF.calculer_h_carre();
          const Elem_VEF& type_elem=la_zone_VEF.type_elem();

          // Recalcul des normales
          normales=0;
          for (int num_face=0; num_face<nb_faces_tot; num_face++)
            type_elem.normale(num_face,normales, face_sommets,
                              face_voisins,elem_faces,
                              la_zone) ;
          type_elem.creer_facette_normales(la_zone, facette_normales_, rang_elem_non_standard);
          //Cerr << "carre_pas_du_maillage : " << la_zone_VEF.carre_pas_du_maillage() << finl;
          int nb_eqn=pb.nombre_d_equations();

          for(int num_eq=0; num_eq<nb_eqn; num_eq++)
            {
              Zone_Cl_dis& zcl_dis=pb.equation(num_eq).zone_Cl_dis();
              Zone_Cl_VEF& la_zcl_VEF=ref_cast(Zone_Cl_VEF, zcl_dis.valeur());
              la_zcl_VEF.remplir_volumes_entrelaces_Cl(la_zone_VEF);
              la_zcl_VEF.remplir_normales_facettes_Cl(la_zone_VEF );
            }

          // Recalcul des surfaces avec les normales:
          DoubleVect face_surfaces_(nb_faces_tot);
          for (int i=0; i<nb_faces_tot; i++)
            {
              double surf=0;
              for (int k=0; k<dimension; k++)
                surf += (la_zone_VF.face_normales(i,k)*la_zone_VF.face_normales(i,k));
              face_surfaces_(i) = sqrt(surf);
            }
          la_zone_VF.calculer_face_surfaces(face_surfaces_);


        }
      else
        {
          Cerr << "Discretisation non reconnue par ALE!" << finl;
          exit();
        }
    }
}
void Domaine_ALE::initialiser (double temps, Domaine_dis& le_domaine_dis,Probleme_base& pb)
{
  //Cerr << "Domaine_ALE::initialiser  " << finl;

  deformable_=1;
  zone(0).invalide_octree();
  bool  check_NoZero_ALE= true;
  ALE_mesh_velocity=calculer_vitesse(temps,le_domaine_dis,pb,  check_NoZero_ALE);

  //On initialise les vitesses aux faces
  Zone_VF& la_zone_VF=ref_cast(Zone_VF,le_domaine_dis.zone_dis(0).valeur());
  int nb_faces=la_zone_VF.nb_faces();
  int nb_faces_tot=la_zone_VF.nb_faces_tot();
  int nb_som_face=la_zone_VF.nb_som_face();
  IntTab& face_sommets=la_zone_VF.face_sommets();

  vf.resize(nb_faces, dimension);
  const MD_Vector& md = la_zone_VF.md_vector_faces();
  MD_Vector_tools::creer_tableau_distribue(md, vf);

  calculer_vitesse_faces(ALE_mesh_velocity,nb_faces_tot,nb_som_face,face_sommets);

  // Initializing Ch_front_input_ALE (needed if Ch_front_input_ALE is used in Imposer_vit_bords_ALE block within data-file)
  int Nr_of_vit_bords_ALE=les_champs_front.size();
  if(Nr_of_vit_bords_ALE!=0)
    {
      for(int i=0; i<Nr_of_vit_bords_ALE; i++)
        {
          Champ_front_base& firstField = les_champs_front[i].valeur();
          //Cerr <<  firstField.que_suis_je() << finl;
          if(sub_type(Ch_front_input_ALE, firstField))
            {
              Ch_front_input_ALE& aleField = ref_cast(Ch_front_input_ALE, firstField);
              const Champ_Inc_base& vit = ref_cast(Champ_Inc_base, pb.get_champ("VITESSE"));
              aleField.initialiser(temps, vit);
            }
        }
    }
  //End of initializing Ch_front_input_ALE
}

DoubleTab Domaine_ALE::calculer_vitesse(double temps, Domaine_dis& le_domaine_dis,Probleme_base& pb, bool& check_NoZero_ALE)
{
  int n; // A activer ou desactiver si on utilise le laplacien ou non
  int N_som=nb_som_tot(); //A activer ou desactiver si on utilise le laplacien ou non

  const Zone& lazone = les_zones(0);
  const MD_Vector& md = lazone.domaine().md_vector_sommets();
  DoubleTab vit_maillage(lazone.nb_som(),dimension);
  MD_Vector_tools::creer_tableau_distribue(md, vit_maillage);

  vit_maillage=0.;
  DoubleTab vit_bords(vit_maillage);
  //DoubleTab tab_champ_front(vit_maillage);
  for (n=0; n<nb_bords_ALE; n++)
    {
      const Nom& le_nom_bord_ALE=les_bords_ALE(n).le_nom();
      int rang=les_zones(0).rang_frontiere(le_nom_bord_ALE);
      const Frontiere_dis_base& la_fr_dis=le_domaine_dis.zone_dis(0).frontiere_dis(rang);
      les_champs_front[n].valeur().associer_fr_dis_base(la_fr_dis);
      const Nom& le_nom_ch_front_courant=les_champs_front[n].valeur().que_suis_je();
      if (le_nom_ch_front_courant == "Champ_front_ALE")
        {
          ref_cast(Champ_front_ALE, les_champs_front[n].valeur()).remplir_vit_som_bord_ALE(temps);
          if (nb_bords_ALE==1)
            {
              //vit_bords+=ref_cast(Champ_front_ALE, les_champs_front[n].valeur()).get_vit_som_bord_ALE();
              DoubleTab vit_bord_ale = ref_cast(Champ_front_ALE, les_champs_front[n].valeur()).get_vit_som_bord_ALE();

              for(int dim=0; dim<dimension; dim++)
                {
                  for(int i=0; i<N_som; i++)
                    {
                      vit_bords(i, dim) += vit_bord_ale(i, dim);
                    }
                }
            }
          else   //If two ALE boundaries have mutual node, then ALE boundary velocity should be taken account only once.
            {
              DoubleTab vit_bord_ale = ref_cast(Champ_front_ALE, les_champs_front[n].valeur()).get_vit_som_bord_ALE();
              for(int dim=0; dim<dimension; dim++)
                {
                  for(int i=0; i<N_som; i++)
                    {
                      if(vit_bords(i, dim) == 0.)
                        vit_bords(i, dim) += vit_bord_ale(i, dim);
                    }
                }
            }
        }
      else if (le_nom_ch_front_courant == "Ch_front_input_ALE")
        {
          ref_cast(Ch_front_input_ALE, les_champs_front[n].valeur()).remplir_vit_som_bord_ALE(temps);
          if (nb_bords_ALE==1)
            {
              DoubleTab vit_bord_ale = ref_cast(Ch_front_input_ALE, les_champs_front[n].valeur()).get_vit_som_bord_ALE();
              for(int dim=0; dim<dimension; dim++)
                {
                  for(int i=0; i<N_som; i++)
                    {
                      vit_bords(i, dim) += vit_bord_ale(i, dim);
                    }
                }
            }
          else   //If two ALE boundaries have mutual node, then ALE boundary velocity should be taken account only once.
            {
              DoubleTab vit_bord_ale = ref_cast(Ch_front_input_ALE, les_champs_front[n].valeur()).get_vit_som_bord_ALE();
              for(int dim=0; dim<dimension; dim++)
                {
                  for(int i=0; i<N_som; i++)
                    {
                      if(vit_bords(i, dim) == 0.)
                        vit_bords(i, dim) += vit_bord_ale(i, dim);
                    }
                }
            }
        }
      else if (le_nom_ch_front_courant == "Champ_front_ALE_Beam")
        {
          ref_cast(Champ_front_ALE_Beam, les_champs_front[n].valeur()).remplir_vit_som_bord_ALE(temps);
          DoubleTab vit_bord_ale = ref_cast(Champ_front_ALE_Beam, les_champs_front[n].valeur()).get_vit_som_bord_ALE();
          if (nb_bords_ALE==1)
            {
              for(int dim=0; dim<dimension; dim++)
                {
                  for(int i=0; i<N_som; i++)
                    {
                      vit_bords(i, dim) += vit_bord_ale(i, dim);
                    }
                }
            }
          else
            {
              for(int dim=0; dim<dimension; dim++)
                {
                  for(int i=0; i<N_som; i++)
                    {
                      if(vit_bords(i, dim) == 0.)
                        vit_bords(i, dim) += vit_bord_ale(i, dim);
                    }
                }
            }
          //double dt_beam = computeDtBeam(le_domaine_dis);
          //Cerr << " dt: "<< dt_beam << endl;
          //getchar();
        }
      else
        {
          Cerr << "Un champ front de type : "
               << les_champs_front[n].valeur().le_nom()
               << " ne peut etre utilise pour un probleme ALE pour le moment...."
               << finl;
          exit();
        }
    }
  vit_bords.echange_espace_virtuel();
  //If ALE boundary velocity is zero, then ALE mesh velocity is directly zero (= default initialization value). Otherwise laplacien() is used.
  if(vit_bords.mp_max_abs_vect() >=1.e-12)
    {
      laplacien(le_domaine_dis, pb, vit_bords, vit_maillage);
      check_NoZero_ALE = true;
    }
  else
    {
      check_NoZero_ALE = false;
      update_or_not_matrix_coeffs_=1;
    }

  Debog::verifier("Domaine_ALE::calculer_vitesse -vit_maillage", vit_maillage);
  return vit_maillage;

}

DoubleTab& Domaine_ALE::calculer_vitesse_faces(DoubleTab& vit_maillage,int nb_faces,int nb_som_face, IntTab& face_sommets)
{
  int i,j,s;
  vf=0.;
  //If ALE mesh velocity is zero, then face velocity is directly zero (= default initialization value). Otherwise it is calculated.
  if(vit_maillage.mp_max_abs_vect() >=1.e-12)
    {
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
    }
  vf.echange_espace_virtuel();
  Debog::verifier("Domaine_ALE::calculer_vitesse_faces -vit_maillage_aux_faces", vf);
  return vf;
}

DoubleTab& Domaine_ALE::laplacien(Domaine_dis& le_domaine_dis,Probleme_base& pb, const DoubleTab& vit_bords, DoubleTab& ch_som)
{
  //Cerr << "Domaine_ALE::laplacien" << finl;
  const Zone& lazone = les_zones(0);
  int nb_elem_tot = lazone.nb_elem_tot();
  //int nbsom = lazone.nb_som();
  int nbsom = lazone.nb_som_tot();
  int nb_som_elem = lazone.nb_som_elem();
  const Zone_VEF& zone_VEF=ref_cast(Zone_VEF,le_domaine_dis.zone_dis(0).valeur());
  const DoubleTab& normales=zone_VEF.face_normales();
  const Zone_Cl_VEF& zone_Cl_VEF = ref_cast(Zone_Cl_VEF, pb.equation(0).zone_Cl_dis().valeur());
  const IntTab& elem_som = lazone.les_elems();
  const IntTab& elem_faces=zone_VEF.elem_faces();
  double mijK;
  int nb_comp=ch_som.dimension(1);
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
  {
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
    //Cerr << "Matrice de filtrage OK" << finl;
    //Cerr << "Matrice de Laplacien P1 : " << finl;
  }
  //Debog::verifier_Mat_elems("Matrice de Laplacien", mat);
  DoubleVect secmem(lazone.nb_som());
  const MD_Vector& md = lazone.domaine().md_vector_sommets();
  MD_Vector_tools::creer_tableau_distribue(md, secmem);

  DoubleVect solution(secmem);

  for(int comp=0; comp<nb_comp; comp++)
    {
      secmem = 0.;
      solution = 0.;
      for (n_bord=0; n_bord<zone_VEF.nb_front_Cl(); n_bord++)
        {
          //for n_bord
          const Cond_lim& la_cl = zone_Cl_VEF.les_conditions_limites(n_bord);
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();

          for (int face=num1; face<num2; face++)
            {
              for(int isom=0; isom<dimension; isom++)
                {
                  int som=zone_VEF.face_sommets(face,isom);
                  som=get_renum_som_perio(som);
                  secmem(som)=1.e6*vit_bords(som,comp);
                }
            }
        }
      secmem.echange_espace_virtuel();
      Debog::verifier("Domaine_ALE::laplacien -secmem", secmem);
      // If secmem is zero, then it is zero solution. Otherwise system is solved.
      if (secmem.mp_max_abs_vect() >=1.e-15)
        {
          //Cerr << "Resolution du systeme de filtrage: ALE" << finl;
          solv.resoudre_systeme(mat, secmem, solution);
          solution.echange_espace_virtuel();
          for(int som=0; som<nbsom; som++)
            ch_som(som,comp)=solution(get_renum_som_perio(som));
        }
      else
        {
          for(int som=0; som<nbsom; som++)
            ch_som(som,comp)=0.0;
        }
    }

  ch_som.echange_espace_virtuel();
  //double x = mp_norme_vect(ch_som);
  //Cerr << "norme(c) = " << x << finl;
  Debog::verifier("Domaine_ALE::laplacien -ch_som", ch_som);
  return ch_som;
}

void Domaine_ALE::set_dt(double& dt)
{
  dt_=dt;
}

void Domaine_ALE::update_ALEjacobians(DoubleTab& NewValueOf_ALEjacobian_old,DoubleTab& NewValueOf_ALEjacobian_new, int TimeStepNr)
{

  if(TimeStepNr==0)
    {
      //Initially ALEjacobian_old= 1.0
      ALEjacobian_old=NewValueOf_ALEjacobian_new; //Give a right size
      ALEjacobian_old= 1.0; //Value to 1
      ALEjacobian_new=NewValueOf_ALEjacobian_new;
      TimeStepNr=1;
    }
  else
    {
      ALEjacobian_old=NewValueOf_ALEjacobian_old;
      ALEjacobian_new=NewValueOf_ALEjacobian_new;
    }
}
void Domaine_ALE::reading_vit_bords_ALE(Entree& is)
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
      les_bords_ALE.add(les_zones(0).faces_bord()(rang));
      is >> les_champs_front[compteur];
      compteur++;
    }
}

//  Read the solver used to solve the system giving the moving mesh velocity
void Domaine_ALE::reading_solver_moving_mesh_ALE(Entree& is)
{
  Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");
  Motcle motlu;
  Nom nomlu;
  is >> motlu;
  if (motlu != accolade_ouverte)
    {
      Cerr << "Error while reading the solveur_moving_mesh_ALE \n";
      Cerr << "We were waiting for a " << accolade_ouverte << " instead of \n"
           << motlu;
      exit();
    }
  is >>  solv;
  solv.nommer("ALE_solver");
  //Cerr << "ALE solver: " << finl;
  while(1)
    {
      // lecture d'un nom de bord ou de }
      is >> nomlu;
      motlu=nomlu;
      if (motlu == accolade_fermee)
        break;
    }
}

void Domaine_ALE::reading_beam_model(Entree& is)
{
  beam.setActivate(true);
  Cerr << "Beam activate : " <<  beam.getActivate() << finl;
  Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");
  Motcle motlu;
  Nom nomlu;
  Nom masse_and_stiffness_file_name;
  Nom phi_file_name;
  Nom absc_file_name ;
  int var_int;
  double var_double;
  is >> motlu;
  if (motlu != accolade_ouverte)
    {
      Cerr << "Erreur a la lecture des vitesses ALE aux bords\n";
      Cerr << "On attendait une " << accolade_ouverte << " a la place de \n"
           << motlu;
      exit();
    }
  while(1)
    {
      // lecture d'un nom de bord ou de }
      is >> nomlu;
      motlu=nomlu;
      if(motlu=="nb_modes")
        {
          is >> var_int;
          beam.setNbModes(var_int);
          Cerr << "Number of modes : " <<  beam.getNbModes() << finl;
        }
      if(motlu=="direction")
        {
          is >> var_int;
          beam.setDirection(var_int);
          Cerr << "Direction : " <<  beam.getDirection() << finl;
        }
      if(motlu=="Young_Module")
        {
          is >> var_double;
          beam.setYoung(var_double);
          Cerr << "Young module : " <<  beam.getYoung() << finl;
        }
      if(motlu=="Rho_beam")
        {
          is >> var_double;
          beam.setRhoBeam(var_double);
          Cerr << "Rho beam : " <<  beam.getRhoBeam() << finl;
        }
      if(motlu=="Mass_and_stiffness_file_name")
        {
          is >> nomlu;
          masse_and_stiffness_file_name=nomlu;
        }
      if(motlu=="Absc_file_name")
        {
          is >> nomlu;
          absc_file_name=nomlu;
        }

      if(motlu=="Modal_deformation_file_name")
        {
          is >> nomlu;
          phi_file_name=nomlu;
        }
      if(motlu=="NewmarkTimeScheme")
        {
          is >> nomlu;
          bool scheme=true;
          if(nomlu=="FD")
            scheme=true;
          else if(nomlu=="MA")
            scheme=false;
          else
            {
              Cerr << "NewmarkTimeScheme wrong: choose between FD and MA" <<finl;
              exit();
            }

          beam.setTimeScheme(scheme);
          Cerr << "TimeScheme: " <<  beam.getTimeScheme() << finl;
        }

      if (motlu == accolade_fermee)
        break;
    }
  beam.readInputMassStiffnessFiles(masse_and_stiffness_file_name);
  beam.readInputAbscFiles(absc_file_name);
  beam.readInputModalDeformation(phi_file_name);


}
DoubleVect Domaine_ALE::interpolationOnThe3DSurface(const double& x, const double& y, const double& z) const
{
  return beam.interpolationOnThe3DSurface(x,y,z);
}

void Domaine_ALE::initializationBeam (double velocity)
{
  beam.initialization(velocity);
}
double Domaine_ALE::computeDtBeam(Domaine_dis& le_domaine_dis)
{

  double dt = 0;
  const Zone_VEF& zone_VEF=ref_cast(Zone_VEF,le_domaine_dis.zone_dis(0).valeur());
  const DoubleVect& surfaces = zone_VEF.face_surfaces();
  double minSurf = mp_min_vect(surfaces);
  minSurf = Process::mp_min(minSurf);
  Cerr << " Surface min: "<< minSurf << endl;
  double soundSpeed=beam.soundSpeed();
  Cerr << " soundSpeed: "<< soundSpeed << endl;
  dt = 0.5*(minSurf/soundSpeed);
  /*Cerr << " dt: "<< dt << endl;
  double massBeam=beam.getMass(0);
  double stiffBeam=beam.getStiffness(0);
  dt = 0.5 * (minSurf * sqrt(massBeam)) / sqrt(stiffBeam);*/
  return dt;
}
