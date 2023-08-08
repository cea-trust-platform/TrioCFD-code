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
#include <Schema_Temps_base.h>
#include <Domaine_VDF.h>
#include <Domaine_VEF.h>
#include <Motcle.h>
#include <EFichier.h>
#include <Champ_front_ALE.h>
#include <Ch_front_input_ALE.h>
#include <Champ_front_ALE_Beam.h>
#include <Champs_front_ALE_projection.h>
#include <Noms.h>
#include <Navier_Stokes_std.h>
#include <Equation_base.h>
#include <Operateur_Diff.h>
#include <Operateur_Grad.h>
#include <communications.h>
#include <Faces.h>

Implemente_instanciable_sans_constructeur_ni_destructeur(Domaine_ALE,"Domaine_ALE",Domaine);
//XD domaine_ale domaine domaine_ale -1 Domain with nodes at the interior of the domain which are displaced in an arbitrarily prescribed way thanks to ALE (Arbitrary Lagrangian-Eulerian) description. NL2 Keyword to specify that the domain is mobile following the displacement of some of its boundaries.
Domaine_ALE::Domaine_ALE() : dt_(0.), nb_bords_ALE(0), update_or_not_matrix_coeffs_(1), resumption(0), associate_eq(false),  tempsComputeForceOnBeam(0.), meshMotionModel_(0)
{
  beam = new Beam_model();
  str_mesh_model = new Structural_dynamic_mesh_model() ;
}
Domaine_ALE::~Domaine_ALE()
{
  delete beam;
  delete str_mesh_model ;
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
  invalide_octree();
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
      Domaine_VF& le_dom_VF=ref_cast(Domaine_VF,le_domaine_dis.valeur());

      int nb_faces=le_dom_VF.nb_faces();
      int nb_som_face=le_dom_VF.nb_som_face();
      IntTab& face_sommets=le_dom_VF.face_sommets();
      //creer_mes_domaines_frontieres(le_dom_VF);//update the boundary surface domain
      calculer_vitesse_faces(ALE_mesh_velocity,nb_faces,nb_som_face,face_sommets);

      //On recalcule les metriques
      le_dom_VF.volumes()=0;
      calculer_volumes(le_dom_VF.volumes(),le_dom_VF.inverse_volumes());
      le_dom_VF.xp()=0;
      calculer_centres_gravite(le_dom_VF.xp());

      DoubleTab& xv=le_dom_VF.xv();
      xv.set_smart_resize(1);
      xv.reset();
      Type_Face type_face=type_elem().type_face();
      IntTab& elem_faces=le_dom_VF.elem_faces();
      IntTab& face_voisins=le_dom_VF.face_voisins();

      ::calculer_centres_gravite(xv, type_face,
                                 sommets_, face_sommets);
      if(sub_type(Domaine_VDF, le_dom_VF))
        {
          Domaine_VDF& le_dom_VDF=ref_cast(Domaine_VDF,le_domaine_dis.valeur());
          le_dom_VF.volumes_entrelaces()=0;
          le_dom_VDF.calculer_volumes_entrelaces();
        }
      else if(sub_type(Domaine_VEF, le_dom_VF))
        {
          Domaine_VEF& le_dom_VEF=ref_cast(Domaine_VEF,le_domaine_dis.valeur());
          DoubleTab& normales=le_dom_VEF.face_normales();
          DoubleTab& facette_normales_=le_dom_VEF.facette_normales();
          IntVect& rang_elem_non_standard=le_dom_VEF.rang_elem_non_std();
          le_dom_VF.volumes_entrelaces()=0;
          le_dom_VEF.calculer_volumes_entrelaces();

          int nb_faces_tot=face_sommets.dimension_tot(0);
          le_dom_VEF.calculer_h_carre();
          const Elem_VEF& type_elem=le_dom_VEF.type_elem();

          // Recalcul des normales
          normales=0;
          for (int num_face=0; num_face<nb_faces_tot; num_face++)
            type_elem.normale(num_face,normales, face_sommets,
                              face_voisins,elem_faces,
                              *this) ;
          type_elem.creer_facette_normales(*this, facette_normales_, rang_elem_non_standard);
          //Cerr << "carre_pas_du_maillage : " << le_dom_VEF.carre_pas_du_maillage() << finl;
          int nb_eqn=pb.nombre_d_equations();

          for(int num_eq=0; num_eq<nb_eqn; num_eq++)
            {
              Domaine_Cl_dis& zcl_dis=pb.equation(num_eq).domaine_Cl_dis();
              Domaine_Cl_VEF& la_zcl_VEF=ref_cast(Domaine_Cl_VEF, zcl_dis.valeur());
              la_zcl_VEF.remplir_volumes_entrelaces_Cl(le_dom_VEF);
              la_zcl_VEF.remplir_normales_facettes_Cl(le_dom_VEF );
            }

          // Recalcul des surfaces avec les normales:
          DoubleVect face_surfaces_(nb_faces_tot);
          for (int i=0; i<nb_faces_tot; i++)
            {
              double surf=0;
              for (int k=0; k<dimension; k++)
                surf += (le_dom_VF.face_normales(i,k)*le_dom_VF.face_normales(i,k));
              face_surfaces_(i) = sqrt(surf);
            }
          le_dom_VF.calculer_face_surfaces(face_surfaces_);


        }
      else
        {
          Cerr << "Discretisation non reconnue par ALE!" << finl;
          exit();
        }

    }
}

void Domaine_ALE::update_after_post(double temps)
{
  update_ALE_projection(temps);
}
//Compute the fluid force projected within the requested boundaries
void Domaine_ALE::update_ALE_projection(double temps,  Nom& name_ALE_boundary_projection, Champ_front_ALE_projection& field_ALE_projection, int nb_mode)
{
  Cerr<<"update_ALE_projection "<<finl;
  const Navier_Stokes_std& eqn_hydr = ref_cast(Navier_Stokes_std,getEquation());
  const Operateur_base& op_grad= eqn_hydr.operateur_gradient().l_op_base();
  const Operateur_base& op_diff= eqn_hydr.operateur_diff().l_op_base();
  const Domaine_VEF& le_dom_vef=ref_cast(Domaine_VEF,op_grad.equation().domaine_dis().valeur());
  const DoubleTab& xv=le_dom_vef.xv();
  DoubleTab& flux_bords_grad=op_grad.flux_bords();
  DoubleTab& flux_bords_diff=op_diff.flux_bords();
  double modalForce = 0.;

  for (int n=0; n<nb_bords_ALE; n++)
    {
      const Nom& le_nom_bord_ALE=les_bords_ALE(n).le_nom();
      if(name_ALE_boundary_projection==le_nom_bord_ALE)
        {
          if((flux_bords_grad.size() == flux_bords_diff.size()) && (flux_bords_grad.size() >0) )
            {
              double phi=0.;
              modalForce=0.;
              int ndeb = les_bords_ALE(n).num_premiere_face();
              int nfin = ndeb + les_bords_ALE(n).nb_faces();

              for (int face=ndeb; face<nfin; face++)
                {
                  for(int comp=0; comp<dimension; comp++)
                    {
                      phi=field_ALE_projection.evaluate(temps, xv(face,0),xv(face,1),xv(face,2), comp);
                      modalForce += (flux_bords_grad(face, comp)+ flux_bords_diff(face, comp))*phi;
                    }
                }
            }


        }

    }
  mp_sum(modalForce);
  if (je_suis_maitre()) //Write the result in the ModalForce_BoundaryName.txt file
    {
      std::string nom="ModalForce_";
      nom += name_ALE_boundary_projection;
      nom +="_";
      std::string index(std::to_string(nb_mode));
      nom +=index;
      nom+=".txt";
      std::ofstream ofs_1;
      ofs_1.open (nom, std::ofstream::out | std::ofstream::app);
      ofs_1<<temps<<" "<<modalForce;
      ofs_1<<endl;
      ofs_1.close();
    }
}
//Compute the fluid force projected within the requested boundaries
void  Domaine_ALE::update_ALE_projection(const double temps)
{

  int size_projection_boundaries=field_ALE_projection_.size();
  if(size_projection_boundaries==0)
    {
      return;
    }
  Cerr<<"update_ALE_projection "<<finl;
  const Navier_Stokes_std& eqn_hydr = ref_cast(Navier_Stokes_std,getEquation());
  const Operateur_base& op_grad= eqn_hydr.operateur_gradient().l_op_base();
  const Operateur_base& op_diff= eqn_hydr.operateur_diff().l_op_base();
  const Domaine_VEF& le_dom_vef=ref_cast(Domaine_VEF,op_grad.equation().domaine_dis().valeur());
  const DoubleTab& xv=le_dom_vef.xv();
  DoubleTab& flux_bords_grad=op_grad.flux_bords();
  DoubleTab& flux_bords_diff=op_diff.flux_bords();
  DoubleVect modalForce(size_projection_boundaries);
  modalForce=0.;
  for (int n=0; n<nb_bords_ALE; n++)
    {
      const Nom& le_nom_bord_ALE=les_bords_ALE(n).le_nom();
      for(int i=0; i<size_projection_boundaries; i++)
        {
          if(name_ALE_boundary_projection_[i]==le_nom_bord_ALE)
            {
              if((flux_bords_grad.size() == flux_bords_diff.size()) && (flux_bords_grad.size() >0) )
                {
                  double phi=0.;
                  modalForce[i]=0.;
                  int ndeb = les_bords_ALE(n).num_premiere_face();
                  int nfin = ndeb + les_bords_ALE(n).nb_faces();

                  for (int face=ndeb; face<nfin; face++)
                    {
                      for(int comp=0; comp<dimension; comp++)
                        {
                          phi=field_ALE_projection_[i].evaluate(temps, xv(face,0),xv(face,1),xv(face,2), comp);
                          modalForce[i] += (flux_bords_grad(face, comp)+ flux_bords_diff(face, comp))*phi;
                        }
                    }
                }
            }
        }

    }


  mp_sum_for_each_item(modalForce);
  if (je_suis_maitre()) //Write the result in the ModalForce_BoundaryName_[i].txt file
    {
      for(int i=0; i<size_projection_boundaries; i++)
        {
          std::string nom="ModalForce_";
          nom += name_ALE_boundary_projection_[i];
          nom +="_";
          std::string index(std::to_string(i+1));
          nom +=index;
          nom+=".txt";
          std::ofstream ofs_1;
          ofs_1.open (nom, std::ofstream::out | std::ofstream::app);
          ofs_1<<temps<<" "<<modalForce[i];
          ofs_1<<endl;
          ofs_1.close();
        }
    }

}

void Domaine_ALE::initialiser (double temps, Domaine_dis& le_domaine_dis,Probleme_base& pb)
{
  //Cerr << "Domaine_ALE::initialiser  " << finl;

  deformable_=1;
  invalide_octree();
  bool  check_NoZero_ALE= true;
  ALE_mesh_velocity=calculer_vitesse(temps,le_domaine_dis,pb,  check_NoZero_ALE);

  //On initialise les vitesses aux faces
  Domaine_VF& le_dom_VF=ref_cast(Domaine_VF,le_domaine_dis.valeur());
  int nb_faces=le_dom_VF.nb_faces();
  int nb_faces_tot=le_dom_VF.nb_faces_tot();
  int nb_som_face=le_dom_VF.nb_som_face();
  IntTab& face_sommets=le_dom_VF.face_sommets();

  if(!associate_eq)
    {
      const Equation_base& equation=pb.equation(0);
      associer_equation(equation);
      associate_eq=true;
    }

  vf.resize(nb_faces, dimension);
  const MD_Vector& md = le_dom_VF.md_vector_faces();
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

  const MD_Vector& md = md_vector_sommets();
  DoubleTab vit_maillage(nb_som(),dimension);
  MD_Vector_tools::creer_tableau_distribue(md, vit_maillage);

  vit_maillage=0.;
  DoubleTab vit_bords(vit_maillage);

  IntVect tag_nodes_bords(nb_som()) ;
  MD_Vector_tools::creer_tableau_distribue(md, tag_nodes_bords);

  //DoubleTab tab_champ_front(vit_maillage);
  for (n=0; n<nb_bords_ALE; n++)
    {
      const Nom& le_nom_bord_ALE=les_bords_ALE(n).le_nom();
      int rang=rang_frontiere(le_nom_bord_ALE);
      const Frontiere_dis_base& la_fr_dis=le_domaine_dis.frontiere_dis(rang);
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
          if(!associate_eq)
            {
              const Equation_base& equation=pb.equation(0);
              associer_equation(equation);
              associate_eq=true;
            }
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
      switch (meshMotionModel_)
        {
        case(0) : // Laplacien
          laplacien(le_domaine_dis, pb, vit_bords, vit_maillage);
          break ;

        case(1) : // Structural dynamics

          // Tag nodes with imposed velocity
          tag_nodes_bords = 0 ;
          const Domaine_VEF& domaine_VEF=ref_cast(Domaine_VEF,le_domaine_dis.valeur());
          const Domaine_Cl_VEF& domaine_Cl_VEF = ref_cast(Domaine_Cl_VEF, pb.equation(0).domaine_Cl_dis().valeur());
          for (int n_bord=0; n_bord<domaine_VEF.nb_front_Cl(); n_bord++)
            {
              //for n_bord
              const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              int num1 = le_bord.num_premiere_face();
              int num2 = num1 + le_bord.nb_faces();

              for (int face=num1; face<num2; face++)
                {
                  int elem=domaine_VEF.face_voisins(face,0);
                  for(int isom=0; isom<dimension; isom++)
                    {
                      int som=domaine_VEF.face_sommets(face,isom);
                      som=get_renum_som_perio(som);
                      tag_nodes_bords[som] = 1 ;
                    }
                }
            }

          tag_nodes_bords.echange_espace_virtuel() ;

          // Initialize x with current coordinates
          for (int i=0; i<N_som; i++)
            {
              for (int k=0; k<dimension; k++)
                str_mesh_model->x(i,k) = coord(i,k) ;
            }

          // Solve explicit dynamic problem giving mesh displacement and velocity at time "temps"
          int nbSom = nb_som() ;
          int nbElem = nb_elem() ;
          int nbSomElem= nb_som_elem() ;
          IntTab& sommets = les_elems() ;

          solveDynamicMeshProblem_(temps, vit_bords, tag_nodes_bords, vit_maillage, nbSom, nbElem, nbSomElem, sommets) ;
          break ;

        default :
          Cerr << "Unknown model for ALE grid motion" << finl ;
          exit();
        }
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
  int nb_elem_t = nb_elem_tot();
  //int nbsom = nb_som();
  int nbsom = nb_som_tot();
  int nb_som_ele = nb_som_elem();
  const Domaine_VEF& domaine_VEF=ref_cast(Domaine_VEF,le_domaine_dis.valeur());
  const DoubleTab& normales=domaine_VEF.face_normales();
  const Domaine_Cl_VEF& domaine_Cl_VEF = ref_cast(Domaine_Cl_VEF, pb.equation(0).domaine_Cl_dis().valeur());
  const IntTab& elem_som = les_elems();
  const IntTab& elem_faces=domaine_VEF.elem_faces();
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
    for (elem=0; elem<nb_elem_t; elem++)
      {
        double volume=domaine_VEF.volumes(elem);
        for (int isom=0; isom<nb_som_ele; isom++)
          {
            int facei=elem_faces(elem,isom);
            int ii=get_renum_som_perio(elem_som(elem,isom));
            for (int jsom=isom+1; jsom<nb_som_ele; jsom++)
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

                coeffij*=domaine_VEF.oriente_normale(facei,elem)*
                         domaine_VEF.oriente_normale(facej,elem);
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
    for (n_bord=0; n_bord<domaine_VEF.nb_front_Cl(); n_bord++)
      {
        //for n_bord
        const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
        const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
        int num1 = le_bord.num_premiere_face();
        int num2 = num1 + le_bord.nb_faces();

        for (int face=num1; face<num2; face++)
          {
            elem=domaine_VEF.face_voisins(face,0);
            for(int isom=0; isom<dimension; isom++)
              {
                int som=domaine_VEF.face_sommets(face,isom);
                som=get_renum_som_perio(som);
                diag[som]=1.e14;
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
  DoubleVect secmem(nb_som());
  const MD_Vector& md = md_vector_sommets();
  MD_Vector_tools::creer_tableau_distribue(md, secmem);

  DoubleVect solution(secmem);

  for(int comp=0; comp<nb_comp; comp++)
    {
      secmem = 0.;
      solution = 0.;
      for (n_bord=0; n_bord<domaine_VEF.nb_front_Cl(); n_bord++)
        {
          //for n_bord
          const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (int face=num1; face<num2; face++)
            {
              for(int isom=0; isom<dimension; isom++)
                {
                  int som=domaine_VEF.face_sommets(face,isom);
                  som=get_renum_som_perio(som);
                  secmem(som)=1.e14*vit_bords(som,comp);
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

  if(TimeStepNr==0 && resumption==0)
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

void Domaine_ALE::resumptionJacobian(DoubleTab& ValueOf_ALEjacobian_old, DoubleTab& ValueOf_ALEjacobian_new)
{

  ALEjacobian_old=ValueOf_ALEjacobian_old;
  ALEjacobian_new=ValueOf_ALEjacobian_new;
  resumption=1;
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
      int rang=rang_frontiere(nomlu);
      les_bords_ALE.add(faces_bord()(rang));
      is >> les_champs_front[compteur];
      compteur++;
    }
}
//Read the projection boundary
void Domaine_ALE::reading_projection_ALE_boundary(Entree& is)
{
  Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");
  Motcle motlu;
  Nom nomlu;
  int nb_projection;
  is >> motlu;
  if (motlu != accolade_ouverte)
    {
      Cerr << "Error when reading the 'Projection_ALE_boundary' \n";
      Cerr << "We were waiting for " << accolade_ouverte << " instead of \n"
           << motlu;
      exit();
    }
  is >> nb_projection;
  Cerr << "Number of ALE projection boundary : " <<  nb_projection << finl;
  field_ALE_projection_.dimensionner(nb_projection);
  int compteur=0;
  while(1)
    {
      // lecture d'un nom de bord ou de }
      is >> nomlu;
      motlu=nomlu;
      if (motlu == accolade_fermee)
        break;
      name_ALE_boundary_projection_.add(nomlu);
      is >> field_ALE_projection_[compteur];
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


//Read the mechanical beam model parameters. See the Beam class for details
void Domaine_ALE::reading_beam_model(Entree& is)
{
  Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");
  Motcle motlu;
  Nom nomlu;
  Nom masse_and_stiffness_file_name;
  Noms phi_file_name;
  Nom absc_file_name;
  Nom CI_file_name;
  Nom Restart_file_name="none";
  int var_int;
  int nb_modes;
  int nb_output_points_1D=0;
  DoubleVect output_position_1D(nb_output_points_1D);
  int nb_output_points_3D=0;
  DoubleTab output_position_3D(nb_output_points_3D,0);
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
          is >> nb_modes;
          beam->setNbModes(nb_modes);
          Cerr << "Number of modes : " <<  beam->getNbModes() << finl;
        }
      if(motlu=="direction")
        {
          is >> var_int;
          beam->setDirection(var_int);
          Cerr << "Direction : " <<  beam->getDirection() << finl;
        }
      if(motlu=="Young_Module")
        {
          is >> var_double;
          beam->setYoung(var_double);
          Cerr << "Young module : " <<  beam->getYoung() << finl;
        }
      if(motlu=="Rho_beam")
        {
          is >> var_double;
          beam->setRhoBeam(var_double);
          Cerr << "Rho beam : " <<  beam->getRhoBeam() << finl;
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
      if(motlu=="CI_file_name")
        {
          is >> nomlu;
          CI_file_name=nomlu;

        }

      if(motlu=="Modal_deformation_file_name")
        {
          is >> var_int;
          for (int i=0; i< var_int; i ++)
            {
              is >>nomlu;
              phi_file_name.add(nomlu);
            }
        }
      if(motlu=="NewmarkTimeScheme")
        {
          is >> nomlu;
          bool scheme=true;
          if(nomlu=="FD")
            scheme=false;
          else if(nomlu=="MA")
            scheme=true;
          else
            {
              Cerr << "NewmarkTimeScheme wrong: choose between FD and MA" <<finl;
              exit();
            }

          beam->setTimeScheme(scheme);
        }
      if(motlu=="Output_position_1D")
        {
          Cerr << "Beam Output_position_1D "<< finl;
          is>>nb_output_points_1D;
          output_position_1D.resize(nb_output_points_1D);
          double poz;
          for (int i=0; i< nb_output_points_1D; i ++)
            {
              is >>poz;
              output_position_1D[i]=poz;
            }

          beam->setOutputPosition1D(output_position_1D);
        }
      if(motlu=="Output_position_3D")
        {
          Cerr << "Beam Output_position_3D "<< finl;
          is>>nb_output_points_3D;
          output_position_3D.resize(nb_output_points_3D, 3);
          double poz=0.;
          for (int i=0; i< nb_output_points_3D; i ++)
            {
              for (int j=0; j< 3; j ++)
                {
                  is >>poz;
                  output_position_3D(i,j)=poz;
                }
            }
          beam->setOutputPosition3D(output_position_3D);
        }
      if (motlu=="Restart_file_name")
        {
          is >> nomlu;
          Restart_file_name=nomlu;
        }

      if (motlu == accolade_fermee)
        break;
    }
  beam->readInputMassStiffnessFiles(masse_and_stiffness_file_name);
  beam->readInputAbscFiles(absc_file_name);
  assert(nb_modes==phi_file_name.size());
  beam->readInputModalDeformation(phi_file_name);
  beam->readInputCIFile(CI_file_name);
  fluidForceOnBeam.resize(nb_modes);
  fluidForceOnBeam=0.;
  tempsComputeForceOnBeam=0.;
  if(Restart_file_name!="none")
    {
      beam->readRestartFile(Restart_file_name);
      tempsComputeForceOnBeam=beam->getTime();

      std::string const nomFichier(Restart_file_name);
      ifstream monFlux(nomFichier.c_str());

      if(monFlux)
        {
          double temps, displacement, speed, acceleration, force;
          for(int i=0; i<nb_modes; i++)
            {
              monFlux >> temps >> displacement >> speed >>acceleration>>force;
              fluidForceOnBeam[i]=force;
            }

          monFlux.close();
        }
      else
        {
          Cerr<< "ERROR: Unable to open the file " <<Restart_file_name<<finl;
        }
    }
  else
    {
      //delete old files if ever the simulation is released in the same folder. This will not delete the files in case of resumption of calculation
      std::remove("BeamDisplacement1D.txt");
      std::remove("BeamVelocity1D.txt");
      std::remove("BeamAcceleration1D.txt");
      std::remove("BeamDisplacement3D.txt");
      std::remove("BeamVelocity3D.txt");
      std::remove("BeamAcceleration3D.txt");
      std::remove("ModalForceFluide1D.txt");
      //end delete old files

      //prepare the headers of the output file ModalForceFluide1D
      std::ofstream ofs;
      ofs.open ("ModalForceFluide1D.txt", std::ofstream::out | std::ofstream::app);
      ofs<<"# Printing modal 1D fluid force: time mode  ";
      for(int k=0; k<nb_modes; k++) ofs<<k+1<<" ";
      ofs<<endl;
      ofs.close();
      //end prepare the headers of the output file ModalForceFluide1D

      if(nb_output_points_1D>0)
        {
          // prepare the headers of the output files of the displacement, velocity and accelerations of the beam
          std::ofstream ofs_1;
          ofs_1.open ("BeamDisplacement1D.txt", std::ofstream::out | std::ofstream::app);
          std::ofstream ofs_2;
          ofs_2.open ("BeamVelocity1D.txt", std::ofstream::out | std::ofstream::app);
          std::ofstream ofs_3;
          ofs_3.open ("BeamAcceleration1D.txt", std::ofstream::out | std::ofstream::app);

          ofs_1<<"# Printing Beam 1D displacement: time  values of x y z -component at points ";
          for(int k=0; k<nb_output_points_1D; k++)
            {
              ofs_1<<output_position_1D[k]<<" ";
            }
          ofs_1<<endl;
          ofs_1.close();

          ofs_2<<"# Printing Beam 1D velocity: time values of x y z -component at points ";
          for(int k=0; k<nb_output_points_1D; k++)
            {
              ofs_2<<output_position_1D[k]<<" ";
            }
          ofs_2<<endl;
          ofs_2.close();


          ofs_3<<"# Printing Beam 1D acceleration: time values of x y z -component at points ";
          for(int k=0; k<nb_output_points_1D; k++)
            {
              ofs_3<<output_position_1D[k]<<" ";
            }
          ofs_3<<endl;
          ofs_3.close();
          //end prepare the headers
        }
      if(nb_output_points_3D>0)
        {
          // prepare the headers of the output files of the displacement, velocity and accelerations of the beam
          std::ofstream ofs_1;
          ofs_1.open ("BeamDisplacement3D.txt", std::ofstream::out | std::ofstream::app);
          std::ofstream ofs_2;
          ofs_2.open ("BeamVelocity3D.txt", std::ofstream::out | std::ofstream::app);
          std::ofstream ofs_3;
          ofs_3.open ("BeamAcceleration3D.txt", std::ofstream::out | std::ofstream::app);

          ofs_1<<"# Printing Beam 3D displacement: time  values of x y z -component at points ";
          for(int k=0; k<nb_output_points_3D; k++)
            {
              ofs_1<<"("<< output_position_3D(k, 0)<<", "<<output_position_3D(k, 1)<<", "<<output_position_3D(k, 2)<<") ";
            }
          ofs_1<<endl;
          ofs_1.close();

          ofs_2<<"# Printing Beam 3D velocity: time values of x y z -component at points ";
          for(int k=0; k<nb_output_points_3D; k++)
            {
              ofs_2<<"("<< output_position_3D(k, 0)<<", "<<output_position_3D(k, 1)<<", "<<output_position_3D(k, 2)<<") ";
            }
          ofs_2<<endl;
          ofs_2.close();


          ofs_3<<"# Printing Beam 3D acceleration: time values of x y z -component at points ";
          for(int k=0; k<nb_output_points_3D; k++)
            {
              ofs_3<<"("<< output_position_3D(k, 0)<<", "<<output_position_3D(k, 1)<<", "<<output_position_3D(k, 2)<<") ";
            }
          ofs_3<<endl;
          ofs_3.close();
          //end prepare the headers
        }
    }

}
DoubleVect Domaine_ALE::interpolationOnThe3DSurface(const double& x, const double& y, const double& z, const DoubleTab& u, const DoubleTab& R) const
{
  return beam->interpolationOnThe3DSurface(x,y,z, u, R);
}

void Domaine_ALE::initializationBeam (double velocity)
{
  beam->initialization(velocity);
}
double Domaine_ALE::computeDtBeam(Domaine_dis& le_domaine_dis)
{

  double dt = 0.;
  const Domaine_VEF& domaine_VEF=ref_cast(Domaine_VEF,le_domaine_dis.valeur());
  const DoubleVect& surfaces = domaine_VEF.face_surfaces();
  double minSurf = mp_min_vect(surfaces);
  minSurf = Process::mp_min(minSurf);
  //Cerr << " Surface min: "<< minSurf << endl;
  double soundSpeed=beam->soundSpeed();
  //Cerr << "soundSpeed: "<< soundSpeed << endl;
  dt = 0.5*(minSurf/soundSpeed);
  return dt;
}

const DoubleTab& Domaine_ALE::getBeamDisplacement(int i) const
{
  return beam->getDisplacement(i);
}
const DoubleTab& Domaine_ALE::getBeamRotation(int i) const
{
  return beam->getRotation(i);

}
const int& Domaine_ALE::getBeamDirection() const
{
  return beam->getDirection();
}
DoubleVect& Domaine_ALE::getBeamVelocity(const double& tps, const double& dt)
{
  //if tps=tempsComputeForceOnBeam then the dynamics of the beam has already been solved. We only solve once per time step.
  if(tps!=tempsComputeForceOnBeam && dt!=0.)
    {
      computeFluidForceOnBeam();
      tempsComputeForceOnBeam = tps; // update the variable tempsComputeForceOnBeam after computing the fluid force
    }
  return beam->getVelocity(tps, dt, fluidForceOnBeam);
}
const int& Domaine_ALE::getBeamNbModes()
{
  return beam->getNbModes();
}

const DoubleVect& Domaine_ALE::getFluidForceOnBeam()
{

  return fluidForceOnBeam;
}

Equation_base& Domaine_ALE::getEquation()
{
  return eq;
}
//Compute the modal fluid force acting on the Beam: sum of a pressure (op_grad.flux_bords) and a viscous term (op_diff.flux_bords) projected on the 3D modal deformation
void  Domaine_ALE::computeFluidForceOnBeam()
{
  const Navier_Stokes_std& eqn_hydr = ref_cast(Navier_Stokes_std,getEquation());
  const Operateur_base& op_grad= eqn_hydr.operateur_gradient().l_op_base();
  const Operateur_base& op_diff= eqn_hydr.operateur_diff().l_op_base();
  const Domaine_VEF& le_dom_vef=ref_cast(Domaine_VEF,op_grad.equation().domaine_dis().valeur());
  const DoubleTab& xv=le_dom_vef.xv();
  DoubleTab& flux_bords_grad=op_grad.flux_bords();
  DoubleTab& flux_bords_diff=op_diff.flux_bords();
  const int nbModes=fluidForceOnBeam.size();

  /* if (flux_bords_grad.size()==0)
     {
       //The flux_bords are zero during the first time step following a resumption of calculation
       //and are updated only at the end of this time step. The call to the "impr" function is used to update flux_bords variables.
       SFichier os("toto_grad.txt");
       op_grad.impr(os);
     }*/

  if((flux_bords_grad.size() == flux_bords_diff.size()) && (flux_bords_grad.size() >0) )
    {
      fluidForceOnBeam=0.;
      DoubleVect phi(3);
      phi=0.;
      for (int n=0; n<nb_bords_ALE; n++)
        {
          int ndeb = les_bords_ALE(n).num_premiere_face();
          int nfin = ndeb + les_bords_ALE(n).nb_faces();

          for (int face=ndeb; face<nfin; face++)
            {
              for(int nbmodes=0; nbmodes<nbModes; nbmodes++)
                {
                  const DoubleTab& u=getBeamDisplacement(nbmodes);
                  const DoubleTab& R=getBeamRotation(nbmodes);
                  phi=interpolationOnThe3DSurface(xv(face,0),xv(face,1),xv(face,2), u, R); //compute the 3D modal deformation
                  for(int comp=0; comp<3; comp++)
                    {
                      fluidForceOnBeam[nbmodes] += (flux_bords_grad(face, comp)+ flux_bords_diff(face, comp))*phi[comp];
                    }
                }
            }

        }
    }
  mp_sum_for_each_item(fluidForceOnBeam);
  if (je_suis_maitre()) // Write the result in the ModalForceFluide1D.txt file
    {
      std::ofstream ofs_1;
      ofs_1.open ("ModalForceFluide1D.txt", std::ofstream::out | std::ofstream::app);
      ofs_1<<tempsComputeForceOnBeam<<" ";
      for(int nbmodes=0; nbmodes<nbModes; nbmodes++)
        ofs_1<<fluidForceOnBeam[nbmodes]<<" ";
      ofs_1<<endl;
      ofs_1.close();
    }
}

void Domaine_ALE::reading_structural_dynamic_mesh_model(Entree& is)
{
  Motcle accolade_ouverte("{");
  Motcle accolade_fermee("}");
  Motcle motlu;
  Nom nomlu;
  double var_double;

  std::vector<double> val_prop ;

  is >> motlu;
  if (motlu != accolade_ouverte)
    {
      Cerr << "Error reading structural dynamic mesh model\n";
      Cerr << "A " << accolade_ouverte << " was expected instead of \n"
           << motlu;
      exit();
    }
  while(1)
    {
      // lecture d'un nom variable ou de }
      is >> nomlu;
      motlu=nomlu;
      if(motlu=="Mfront_library")
        {
          is >> nomlu;
          std::string const mfront_lib_path(nomlu) ;
          str_mesh_model->setMfrontLibraryPath(mfront_lib_path) ;
        }
      if(motlu=="Mfront_model_name")
        {
          is >> nomlu;
          std::string const mfront_model_name(nomlu) ;
          str_mesh_model->setMfrontModelName(mfront_model_name) ;
        }
      if(motlu=="Mfront_material_property")
        {
          is >> motlu;
          if (motlu != accolade_ouverte)
            {
              Cerr << "Error reading Mfront material property\n";
              Cerr << "A " << accolade_ouverte << " was expected instead of \n"
                   << motlu;
              exit();
            }
          is >> motlu ;
          std::string const nom_prop(motlu) ;
          while (1)
            {
              is >> var_double ;
              val_prop.push_back(var_double) ;

              if (motlu == accolade_fermee)
                break;
            }

          str_mesh_model->addMaterialProperty(nom_prop, val_prop) ;
        }
      // Temporary: YoungModulus is needed only for stability, it should be extracted back from Mfront in some way
      if(motlu=="YoungModulus")
        {
          is >> var_double;
          str_mesh_model->setYoungModulus(var_double) ;
        }
      if(motlu=="Density")
        {
          is >> var_double;
          str_mesh_model->setDensity(var_double) ;
        }
      if(motlu=="Inertial_Damping")
        {
          is >> var_double;
          str_mesh_model->setInertialDamping(var_double) ;
        }
      if(motlu=="Time_Step_Safety_Coefficient")
        {
          is >> var_double;
          str_mesh_model->setDtSafetyCoefficient(var_double) ;
        }
      if (motlu == accolade_fermee)
        break;
    }

  if(str_mesh_model->getDensity() == 0.)
    {
      Cerr << "Error: density not provided" << finl;
      exit();
    }
  // Temporary: YoungModulus is needed only for stability, it should be extracted back from Mfront in some way
  if(str_mesh_model->getYoungModulus() == 0.)
    {
      Cerr << "Error: Young modulus not provided" << finl;
      exit();
    }

}

void Domaine_ALE::solveDynamicMeshProblem_(double temps, const DoubleTab& imposedVelocity, const IntVect& imposedVelocityTag,
                                           DoubleTab& outputMeshVelocity, int nbSom, int nbElem, int nbSomElem, IntTab& sommets)
{

  DoubleTab x0(str_mesh_model->x) ; // Copy coordinates at the beginning of the step
  double tt = str_mesh_model->gridTime ;
  double t0 = tt ;
  bool loopOnGridProblem = true ;

  while (loopOnGridProblem)
    {
      // Adjust the grid time step for smooth arrival at time = temps
      if (str_mesh_model->gridDt >= temps - tt)
        {
          str_mesh_model->gridDt = temps - tt ;
          loopOnGridProblem = false ; // final time reached after this last step
        }
      else if (str_mesh_model->gridDt >= 0.5 * (temps - tt)) // Avoid excessive time step variations
        {
          str_mesh_model->gridDt = 0.5 * (temps - tt) ;
        }

      double Dt = str_mesh_model->gridDt ;
      tt += Dt ;

      // Update mid-step velocities, displacements and coordinates
      for (int i=0; i<nbSom; i++)
        {
          for (int j=0; j<dimension; j++)
            {
              str_mesh_model->vp(i,j) = str_mesh_model->v(i,j) + 0.5 * Dt * str_mesh_model->a(i,j) ;
              if (imposedVelocityTag[j] == 1)
                {
                  str_mesh_model->vp(i,j) = imposedVelocity(j,i) ; // Apply imposed velocity from the boundary
                }
              double du = Dt * str_mesh_model->vp(i,j) ;
              str_mesh_model->u(i,j) += du ;
              str_mesh_model->x(i,j) += du ;
            }
        }

      // Compute internal forces
      str_mesh_model->ff = 0. ;
      int nbn = str_mesh_model->getNbNodesPerElem() ;
      int elnodes[nbn] ;
      double volume ;
      double xlong ;
      double E ;
      double density = str_mesh_model->getDensity() ;
      double dts ;

      str_mesh_model->gridDt = 1.E12 ; // Initialize gridDt at <huge>

      for (int elem=0; elem<nbElem; elem++)
        {
          for (int i=0; i< nbn; i++)
            {
              elnodes[i] = get_renum_som_perio(sommets(elem,i)) ;
            }

          str_mesh_model->setLocalFields(elnodes, elem) ; // x, u global to local + elem id

          str_mesh_model->computeInternalForces(volume, xlong, E) ; // local force computation on element elem

          str_mesh_model->setGlobalFields(elnodes) ; // ff back to global

          if (!str_mesh_model->isMassBuilt)
            {
              str_mesh_model->setMassElem(volume * density) ;
              for (int i=0; i< nbn; i++)
                {
                  int ii = elnodes[i] ;
                  str_mesh_model->mass[ii] += volume * density / nbn ;
                }
            }

          // Set next time step
          dts = str_mesh_model->computeCriticalDt(volume, xlong, E) ;
          str_mesh_model->gridDt = std::min(str_mesh_model->gridDt, dts) ;
        }

      MD_Vector_tools::echange_espace_virtuel(str_mesh_model->ff, MD_Vector_tools::EV_SOMME_ECHANGE) ;
      if (!str_mesh_model->isMassBuilt)
        {
          MD_Vector_tools::echange_espace_virtuel(str_mesh_model->mass, MD_Vector_tools::EV_SOMME_ECHANGE) ;
          str_mesh_model->isMassBuilt = true ;
        }

      // Compute accelerations and full step velocities
      double rhs ;
      double den ;
      double d = str_mesh_model->getDampingCoefficient() ;
      for (int i=0; i<nbSom; i++)
        {
          for (int j=0; j<dimension; j++)
            {
              rhs = -str_mesh_model->ff(i,j) - d * str_mesh_model->mass[j] * str_mesh_model->vp(i,j) ;
              den = str_mesh_model->mass[j] * (1. + 0.5 * d * Dt) ;
              str_mesh_model->a(i,j) = rhs / den ;
              str_mesh_model->v(i,j) = str_mesh_model->vp(i,j) + 0.5 * Dt * str_mesh_model->a(i,j) ;
            }
        }
    } // End time loop on grid problem

  str_mesh_model->gridTime = tt ;

  for (int i=0; i<nbSom; i++)
    {
      for (int j=0; j<dimension; j++)
        {
          outputMeshVelocity(i,j) = (str_mesh_model->x(i,j) - x0(i,j)) / (tt - t0) ;
        }
    }

}
