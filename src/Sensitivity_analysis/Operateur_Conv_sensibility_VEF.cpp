/****************************************************************************
* Copyright (c) 2020, CEA
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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Operateur_Conv_sensibility_VEF.cpp
// Directory : $$Sensitivity_analysis/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Operateur_Conv_sensibility_VEF.h>
#include <Op_Conv_VEF_base.h>
#include <Navier_Stokes_std_sensibility.h>
#include <Convection_Diffusion_Temperature_sensibility.h>
#include <Fluide_Incompressible.h>
#include <TRUSTTrav.h>
#include <Op_Conv_VEF_Face.h>
#include <Probleme_base.h>
#include <Schema_Temps_base.h>
#include <Periodique.h>
#include <Symetrie.h>
#include <Neumann_homogene.h>
#include <Neumann_paroi.h>
#include <Echange_externe_impose.h>
#include <Neumann_sortie_libre.h>
#include <Debog.h>
#include <Champ_P1NC.h>
#include <Porosites_champ.h>
#include <Convection_tools.h>
#include <CL_Types_include.h>


Implemente_instanciable( Operateur_Conv_sensibility_VEF, "Op_Conv_sensibility_VEF_P1NC",Operateur_Conv_sensibility ) ;


Sortie&  Operateur_Conv_sensibility_VEF::printOn(Sortie& os) const
{
  return os;
}


Entree&  Operateur_Conv_sensibility_VEF::readOn(Entree& is)
{
  Cerr << " Operateur_Conv_sensibility_VEF ::readOn " << finl;
  op_conv.associer_eqn(equation());
  op_conv.associer_vitesse(la_vitesse.valeur());
  op_conv.lire(is);
  Cerr << "Operateur_Conv_sensibility_VEF : " << op_conv.que_suis_je() << finl;
  return is;
}

void  Operateur_Conv_sensibility_VEF::associer (const Domaine_dis& domaine_dis ,
                                                const Domaine_Cl_dis& domaine_cl_dis,
                                                const Champ_Inc& inco )
{
  Cerr << " Operateur_Conv_sensibility_VEF::associer" << finl;
  const Domaine_VEF& zvef = ref_cast(Domaine_VEF,domaine_dis.valeur());
  const Domaine_Cl_VEF& zclvef = ref_cast(Domaine_Cl_VEF,domaine_cl_dis.valeur());

  le_dom_vef = zvef;
  la_zcl_vef = zclvef;
  le_dom_vef->creer_tableau_faces(fluent);
  Operateur_Conv_sensibility::associer(domaine_dis,domaine_cl_dis,inco);
}

DoubleTab& Operateur_Conv_sensibility_VEF::ajouter(const DoubleTab& inco, DoubleTab& resu) const
{
  Cerr << "Operateur_Conv_sensibility_VEF::ajouter" << finl;

  //Check convection scheme discretization type:
  const Op_Conv_VEF_base& opConvVEFbase = ref_cast(Op_Conv_VEF_base, op_conv.valeur());
  const Op_Conv_VEF_Face& opConvVEFFace = ref_cast(Op_Conv_VEF_Face, op_conv.valeur());
  int convectionSchemeDiscrType; // amont=0, muscl=1, centre=2
  opConvVEFFace.get_type_op(convectionSchemeDiscrType);
  const Domaine_VEF& domaine_VEF = ref_cast(Domaine_VEF, le_dom_vef.valeur());

  if(convectionSchemeDiscrType==0 || convectionSchemeDiscrType==1 || convectionSchemeDiscrType==2) // Convection scheme discr "amont".
    {
      if(equation().que_suis_je()=="Navier_Stokes_standard_sensibility")
        {
          const Navier_Stokes_std_sensibility& eq = ref_cast(Navier_Stokes_std_sensibility, equation());
          const Champ_Inc_base& state = eq.get_state();
          const DoubleTab& state_field = eq.get_state_field();
          const Motcle& uncertain_var =  eq.get_uncertain_variable_name();
          // Dimensionnement du tableau des flux convectifs au bord du domaine de calcul
          DoubleTab& flux_b = flux_bords_;
          int nb_faces_bord=domaine_VEF.nb_faces_bord();
          flux_b.resize(nb_faces_bord,inco.dimension(1));
          flux_b = 0.;
          ajouter_conv_term(state,inco ,resu, flux_b);
          ajouter_conv_term(la_vitesse,state_field ,resu, flux_b);

          //ajouter_Lstate_sensibility_Amont(state, inco, resu);
          //ajouter_Lsensibility_state_Amont(inco, state, resu);
          if(uncertain_var=="MU")
            add_diffusion_term(state.valeurs(), resu);
          opConvVEFbase.modifier_flux(*this); // Multiplication by rho

        }
      else if (equation().que_suis_je()=="Convection_Diffusion_Temperature_sensibility")
        {
          const Convection_Diffusion_Temperature_sensibility& eq = ref_cast(Convection_Diffusion_Temperature_sensibility, equation());
          const Champ_Inc_base& velocity_state = eq.get_velocity_state();
          const DoubleTab& temperature_state = eq.get_temperature_state_field();
          const Motcle& uncertain_var =  eq.get_uncertain_variable_name();

          DoubleTab& flux_b = flux_bords_;
          int nb_faces_bord=domaine_VEF.nb_faces_bord();
          flux_b.resize(nb_faces_bord,inco.dimension(0));
          flux_b = 0.;

          ajouter_conv_term(velocity_state,inco ,resu, flux_b);
          ajouter_conv_term(vitesse(), temperature_state,resu, flux_b);
          opConvVEFbase.modifier_flux(*this); // Multiplication by rhoCp
          if(uncertain_var=="LAMBDA")
            {
              add_diffusion_scalar_term(temperature_state, resu);
            }
          if(uncertain_var=="CP")
            {
              double lambda_div_Cp=-1.;
              const double Cp = eq.fluide().capacite_calorifique()->valeurs()(0, 0);
              const double lambda = eq.fluide().conductivite()->valeurs()(0, 0);
              lambda_div_Cp*=(lambda/Cp);
              add_diffusion_scalar_term(temperature_state, resu,lambda_div_Cp);
            }


          if(uncertain_var!="TEMPERATURE" && uncertain_var!="BOUSSINESQ_TEMPERATURE"
              && uncertain_var!="BETA_TH"  && uncertain_var!="LAMBDA" && uncertain_var!="CP")
            {
              Cout << "Variable "<<uncertain_var<<" Is not available yet."<<" "
                   "Try available TEMPERATURE or BOUSSINESQ_TEMPERATURE or BETA_TH or LAMBDA or CP variable." << finl;
              exit();
            }


        }
      else
        {
          Cout << "Operateur_Conv_sensibility_VEF::ajouter() =>  Sensibility cannot use currently equation " << equation().que_suis_je()
               <<". Try available Navier_Stokes_standard_sensibility or Convection_Diffusion_Temperature_sensibility." << finl;
          exit();
        }

    }
  else
    {
      Cout << "Sensibility cannot use currently convection scheme " << op_conv.type() <<". Try available Sensibility -centre -amont -muscl convection scheme." << finl;
      exit();
    }
  resu.echange_espace_virtuel();
  Debog::verifier("resu dansOperateur_Conv_sensibility_VEF::ajouter", resu);
  return resu;
}

void Operateur_Conv_sensibility_VEF::ajouter_conv_term(const Champ_Inc_base& velocity, const DoubleTab& transporte, DoubleTab& resu,  DoubleTab& flux_b) const
{
  const Domaine_Cl_VEF& domaine_Cl_VEF = la_zcl_vef.valeur();
  const Domaine_VEF& domaine_VEF = ref_cast(Domaine_VEF, le_dom_vef.valeur());
  const Op_Conv_VEF_base& opConvVEFbase = ref_cast(Op_Conv_VEF_base, op_conv.valeur());
  const Op_Conv_VEF_Face& opConvVEFFace = ref_cast(Op_Conv_VEF_Face, op_conv.valeur());

  int ordre;
  opConvVEFFace.get_ordre(ordre);

  double alpha_;
  opConvVEFFace.get_alpha(alpha_);

  Motcle type_limit;
  if(op_conv.type()=="MUSCL") //type_lumit has default value set in Op_Conv_Muscl_VEF_Face.cpp .
    {
      type_limit="vanleer";
    }
  else // In case of: "convection {  { generic muscl [limiter] [order of accuracy] [alpha] } }". op_conv.type()=="GENERIC" .
    {
      opConvVEFFace.get_type_lim(type_limit);
    }
  enum type_operateur { amont, muscl, centre };

  const DoubleTab& vitesse_face_absolue=velocity.valeurs();
  const DoubleVect& porosite_face = la_zcl_vef->equation().milieu().porosite_face();

  int marq=opConvVEFbase.phi_u_transportant(equation());
  DoubleTab transporte_face_;
  DoubleTab vitesse_face_;
  // soit on a transporte_face=phi*transporte et vitesse_face=vitesse
  // soit on a transporte_face=transporte et vitesse_face=phi*vitesse
  // cela depend si on transporte avec phi*u ou avec u.
  const DoubleTab& transporte_face = modif_par_porosite_si_flag(transporte,transporte_face_,!marq,porosite_face);
  const DoubleTab& vitesse_face    = modif_par_porosite_si_flag(vitesse_face_absolue,vitesse_face_,marq,porosite_face);

  const IntTab& elem_faces = domaine_VEF.elem_faces();
  const DoubleTab& facenormales = domaine_VEF.face_normales();
  const DoubleTab& facette_normales = domaine_VEF.facette_normales();
  const Domaine& domaine = domaine_VEF.domaine();
  const int nfa7 = domaine_VEF.type_elem().nb_facette();
  const int nb_elem_tot = domaine_VEF.nb_elem_tot();
  const IntVect& rang_elem_non_std = domaine_VEF.rang_elem_non_std();
  const IntTab& face_voisins = domaine_VEF.face_voisins();
  const DoubleTab& normales_facettes_Cl = domaine_Cl_VEF.normales_facettes_Cl();
  int premiere_face_int = domaine_VEF.premiere_face_int();
  int nfac = domaine.nb_faces_elem();
  int nsom = domaine.nb_som_elem();
  const IntTab& sommet_elem = domaine.les_elems();
  const DoubleTab& vecteur_face_facette = ref_cast_non_const(Domaine_VEF,domaine_VEF).vecteur_face_facette();
  const  DoubleTab& vecteur_face_facette_Cl = domaine_Cl_VEF.vecteur_face_facette_Cl();
  int nb_bord = domaine_VEF.nb_front_Cl();
  const IntTab& les_elems=domaine.les_elems();

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
              const Cond_lim_base& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord).valeur();
              if ( sub_type(Dirichlet,la_cl) || sub_type(Dirichlet_homogene,la_cl) )
                {
                  const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
                  int nb_faces_tot = le_bord.nb_faces_tot();
                  for (int ind_face=0; ind_face<nb_faces_tot; ind_face++)
                    {
                      int num_face = le_bord.num_face(ind_face);
                      int elem = face_voisins(num_face,0);
                      traitement_pres_bord_[elem]=1;
                    }
                }
            }
        }
      else
        {
          // Pour le muscl/centre actuels on utilise un calcul de flux a l'ordre 1
          // aux mailles de bord ou aux mailles ayant un sommet de Dirichlet
          ArrOfInt est_un_sommet_de_bord_(domaine_VEF.nb_som_tot());
          for (int n_bord=0; n_bord<nb_bord; n_bord++)
            {
              const Cond_lim_base& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord).valeur();
              if ( sub_type(Dirichlet,la_cl) || sub_type(Dirichlet_homogene,la_cl) )
                {
                  const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
                  int nb_faces_tot = le_bord.nb_faces_tot();
                  int size = domaine_VEF.face_sommets().dimension(1);
                  for (int ind_face=0; ind_face<nb_faces_tot; ind_face++)
                    for (int som=0; som<size; som++)
                      {
                        int face = le_bord.num_face(ind_face);
                        est_un_sommet_de_bord_[domaine_VEF.face_sommets(face,som)]=1;
                      }
                }
            }
          for (int elem=0; elem<nb_elem_tot; elem++)
            {
              if (rang_elem_non_std(elem)!=-1)
                traitement_pres_bord_[elem]=1;
              else
                {
                  for (int n_som=0; n_som<nsom; n_som++)
                    if (est_un_sommet_de_bord_[les_elems(elem,n_som)])
                      traitement_pres_bord_[elem]=1;
                }
            }
        }
      // Construction du tableau est_une_face_de_dirichlet_
      est_une_face_de_dirichlet_.resize_array(domaine_VEF.nb_faces_tot());
      est_une_face_de_dirichlet_=0;
      for (int n_bord=0; n_bord<nb_bord; n_bord++)
        {
          const Cond_lim_base& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord).valeur();
          if ( sub_type(Dirichlet,la_cl) || sub_type(Dirichlet_homogene,la_cl) )
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              int nb_faces_tot = le_bord.nb_faces_tot();
              for (int ind_face=0; ind_face<nb_faces_tot; ind_face++)
                {
                  int num_face = le_bord.num_face(ind_face);
                  est_une_face_de_dirichlet_[num_face] = 1;
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

  const Elem_VEF_base& type_elemvef= domaine_VEF.type_elem();
  int istetra=0;
  Nom nom_elem=type_elemvef.que_suis_je();
  if ((nom_elem=="Tetra_VEF")||(nom_elem=="Tri_VEF"))
    istetra=1;

  const DoubleVect& porosite_elem = la_zcl_vef->equation().milieu().porosite_elem();
  double psc;
  int poly,face_adj,fa7,i,j,n_bord;
  int num_face, rang;
  int num10,num20,num_som;
  int ncomp_ch_transporte=(transporte_face.nb_dim() == 1?1:transporte_face.dimension(1));

  // Traitement particulier pour les faces de periodicite
  int nb_faces_perio = 0;
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
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
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
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
  int nb_faces_ = domaine_VEF.nb_faces();
  ArrOfInt face(nfac);
  //statistiques().end_count(m1);
  //statistiques().begin_count(m2);

  // Tableau gradient base sur gradient_elem selon schema
  DoubleTab gradient_elem(nb_elem_tot,ncomp_ch_transporte,dimension);  // (du/dx du/dy dv/dx dv/dy) pour un poly
  if(op_conv.type()=="CENTRE" || op_conv.type()=="MUSCL")
    {
      Champ_P1NC::calcul_gradient(transporte_face,gradient_elem,domaine_Cl_VEF);
    }
  DoubleTab gradient;
  if (op_conv.type()=="CENTRE")
    {
      gradient.ref(gradient_elem);
    }
  else if(op_conv.type()=="MUSCL")
    {
      //  application du limiteur
      gradient.resize(0, ncomp_ch_transporte, dimension);     // (du/dx du/dy dv/dx dv/dy) pour une face
      domaine_VEF.creer_tableau_faces(gradient);
      for (n_bord=0; n_bord<nb_bord; n_bord++)
        {
          const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
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
                      //const DoubleTab& facenormales = domaine_VEF.face_normales();
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
          if (ordre==3 && (traitement_pres_bord_[elem1] || traitement_pres_bord_[elem2])) minmod_pres_du_bord = 1;
          for (comp0=0; comp0<ncomp_ch_transporte; comp0++)
            for (i=0; i<dimension; i++)
              {
                double grad1=gradient_elem(elem1, comp0, i);
                double grad2=gradient_elem(elem2, comp0, i);
                if (minmod_pres_du_bord)
                  gradient(fac, comp0, i) = minmod(grad1, grad2);
                else
                  gradient(fac, comp0, i) =application_LIMITEUR(grad1, grad2, type_limit);
              }
        } // fin du for faces
      gradient.echange_espace_virtuel();
    }// fin if(type_op==muscl)

  ArrOfDouble vs(dimension);
  ArrOfDouble vc(dimension);
  ArrOfDouble cc(dimension);
  DoubleVect xc(dimension);
  DoubleTab vsom(nsom,dimension);
  DoubleTab xsom(nsom,dimension);


  int nb_faces_bord=domaine_VEF.nb_faces_bord();
  const IntTab& KEL=type_elemvef.KEL();
  const DoubleTab& xv=domaine_VEF.xv();
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
      type_operateur type_op_boucle;
      if(op_conv.type()=="MUSCL")
        type_op_boucle = muscl;
      else if (op_conv.type()=="AMONT")
        type_op_boucle = amont;
      else
        type_op_boucle =centre;
      if (passe==2)
        {
          type_op_boucle = centre;
          gradient.ref(gradient_elem);
        }
      // Les polyedres non standard sont ranges en 2 groupes dans le Domaine_VEF:
      //  - polyedres bords et joints
      //  - polyedres bords et non joints
      // On traite les polyedres en suivant l'ordre dans lequel ils figurent
      // dans le domaine
      // boucle sur les polys
      for (poly=0; poly<nb_elem_tot; poly++)
        {
          int contrib = 0;
          // calcul des numeros des faces du polyedre
          for (face_adj=0; face_adj<nfac; face_adj++)
            {
              int face_ = elem_faces(poly,face_adj);
              face[face_adj]= face_;
              if (face_<nb_faces_) contrib=1; // Une face reelle sur l'element virtuel
            }
          //
          if (contrib)
            {
              int calcul_flux_en_un_point = (ordre != 3) && (ordre==1 || traitement_pres_bord_[poly]);
              for (j=0; j<dimension; j++)
                {
                  vs[j] = vitesse_face_absolue(face[0],j)*porosite_face[face[0]];
                  for (i=1; i<nfac; i++)
                    vs[j]+= vitesse_face_absolue(face[i],j)*porosite_face[face[i]];
                }
              // calcul de la vitesse aux sommets des polyedres
              // On va utliser les fonctions de forme implementees dans la classe Champs_P1_impl ou Champs_Q1_impl
              if (istetra==1)
                {
                  for (i=0; i<nsom; i++)
                    for (j=0; j<dimension; j++)
                      vsom(i,j) = (vs[j] - dimension*vitesse_face_absolue(face[i],j)*porosite_face[face[i]]);
                }
              else
                {
                  // pour que cela soit valide avec les hexa (c'est + lent a calculer...)
                  for (j=0; j<nsom; j++)
                    {
                      num_som = sommet_elem(poly,j);
                      for (int ncomp=0; ncomp<dimension; ncomp++)
                        vsom(j,ncomp) =velocity.valeur_a_sommet_compo(num_som,poly,ncomp);
                    }
                }

              // Determination du type de CL selon le rang
              rang = rang_elem_non_std(poly);
              int itypcl = (rang==-1 ? 0 : domaine_Cl_VEF.type_elem_Cl(rang));

              // calcul de vc (a l'intersection des 3 facettes) vc vs vsom proportionnelles a la porosite
              type_elemvef.calcul_vc(face,vc,vs,vsom,velocity,itypcl,porosite_face);

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
                  num10 = face[KEL(0,fa7)];
                  num20 = face[KEL(1,fa7)];
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
                    if (est_une_face_de_dirichlet_[num10] || est_une_face_de_dirichlet_[num20])
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
                  if (type_op_boucle==muscl && ordre==3)
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
                  if (type_op_boucle==muscl)
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
                      if (type_op_boucle==amont || appliquer_cl_dirichlet)
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
                          if (num10<nb_faces_bord) flux_b(num10,0) += flux;
                          if (num20<nb_faces_bord) flux_b(num20,0) -= flux;
                        }
                      else
                        {
                          resu(num10,comp0) -= flux;
                          resu(num20,comp0) += flux;
                          if (num10<nb_faces_bord) flux_b(num10,comp0) += flux;
                          if (num20<nb_faces_bord) flux_b(num20,comp0) -= flux;
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
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);

      if (sub_type(Neumann_sortie_libre,la_cl.valeur()))
        {
          const Neumann_sortie_libre& la_sortie_libre = ref_cast(Neumann_sortie_libre, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (num_face=num1; num_face<num2; num_face++)
            {
              psc =0;
              for (i=0; i<dimension; i++)
                psc += vitesse_face(num_face,i)*facenormales(num_face,i);
              if (psc>0)
                if (ncomp_ch_transporte == 1)
                  {
                    resu(num_face) -= psc*transporte_face(num_face);
                    flux_b(num_face,0) = -psc*transporte_face(num_face);
                  }
                else
                  for (i=0; i<ncomp_ch_transporte; i++)
                    {
                      resu(num_face,i) -= psc*transporte_face(num_face,i);
                      flux_b(num_face,i) = -psc*transporte_face(num_face,i);
                    }
              else
                {
                  if (ncomp_ch_transporte == 1)
                    {
                      resu(num_face) -= psc*la_sortie_libre.val_ext(num_face-num1);
                      flux_b(num_face,0) = -psc*la_sortie_libre.val_ext(num_face-num1);
                    }
                  else
                    for (i=0; i<ncomp_ch_transporte; i++)
                      {
                        resu(num_face,i) -= psc*la_sortie_libre.val_ext(num_face-num1,i);
                        flux_b(num_face,i) = -psc*la_sortie_libre.val_ext(num_face-num1,i);
                      }
                }
            }
        }
      else if (sub_type(Periodique,la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
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
                      /* On ne doit pas ajouter a flux_b, c'est deja calcule au dessus
                         flux_b(voisine,0) += diff1;
                         flux_b(num_face,0) += diff2;*/
                    }
                  else
                    for (int comp=0; comp<ncomp_ch_transporte; comp++)
                      {
                        diff1 = resu(num_face,comp)-tab(nb_faces_perio,comp);
                        diff2 = resu(voisine,comp)-tab(nb_faces_perio+voisine-num_face,comp);
                        resu(voisine,comp)  += diff1;
                        resu(num_face,comp) += diff2;
                        /* On ne doit pas ajouter a flux_b, c'est deja calcule au dessus
                           flux_b(voisine,comp) += diff1;
                           flux_b(num_face,comp) += diff2; */
                      }

                  fait[num_face-num1]= 1;
                  fait[voisine-num1] = 1;
                }
              nb_faces_perio++;
            }
        }
    }
}


//state_field grad_inco
void Operateur_Conv_sensibility_VEF::ajouter_Lstate_sensibility_Amont(const DoubleTab& state_field, const DoubleTab& inco, DoubleTab& resu ) const
{
  Cout<<"ajouter_Lstate_sensibility "<<finl;
  const Domaine_Cl_VEF& domaine_Cl_VEF = la_zcl_vef.valeur();
  const Domaine_VEF& domaine_VEF = ref_cast(Domaine_VEF, le_dom_vef.valeur());
  const IntTab& elem_faces = domaine_VEF.elem_faces();
  const DoubleTab& facenormales = domaine_VEF.face_normales();
  const DoubleTab& facette_normales = domaine_VEF.facette_normales();
  const Domaine& domaine = domaine_VEF.domaine();
  const int nfa7 = domaine_VEF.type_elem().nb_facette();
  const int nb_elem_tot = domaine_VEF.nb_elem_tot();
  const IntVect& rang_elem_non_std = domaine_VEF.rang_elem_non_std();
  const DoubleTab& normales_facettes_Cl = domaine_Cl_VEF.normales_facettes_Cl();
  int nfac = domaine.nb_faces_elem();
  int nsom = domaine.nb_som_elem();
  int nb_bord = domaine_VEF.nb_front_Cl();
  const IntTab& les_elems=domaine.les_elems();
  int option_appliquer_cl_dirichlet = 0 ;

  // Definition d'un tableau pour un traitement special des schemas pres des bords
  if (traitement_pres_bord_.size_array()!=nb_elem_tot)
    {
      traitement_pres_bord_.resize_array(nb_elem_tot);
      traitement_pres_bord_=0;
      ArrOfInt est_un_sommet_de_bord_(domaine_VEF.nb_som_tot());
      for (int n_bord=0; n_bord<nb_bord; n_bord++)
        {
          const Cond_lim_base& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord).valeur();
          if ( sub_type(Dirichlet,la_cl) || sub_type(Dirichlet_homogene,la_cl) )
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              int nb_faces_tot = le_bord.nb_faces_tot();
              int size = domaine_VEF.face_sommets().dimension(1);
              for (int ind_face=0; ind_face<nb_faces_tot; ind_face++)
                for (int som=0; som<size; som++)
                  {
                    int face = le_bord.num_face(ind_face);
                    est_un_sommet_de_bord_[domaine_VEF.face_sommets(face,som)]=1;
                  }
            }
        }
      for (int elem=0; elem<nb_elem_tot; elem++)
        {
          if (rang_elem_non_std(elem)!=-1)
            traitement_pres_bord_[elem]=1;
          else
            {
              for (int n_som=0; n_som<nsom; n_som++)
                if (est_un_sommet_de_bord_[les_elems(elem,n_som)])
                  traitement_pres_bord_[elem]=1;
            }
        }
      // Construction du tableau est_une_face_de_dirichlet_
      est_une_face_de_dirichlet_.resize_array(domaine_VEF.nb_faces_tot());
      est_une_face_de_dirichlet_=0;
      for (int n_bord=0; n_bord<nb_bord; n_bord++)
        {
          const Cond_lim_base& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord).valeur();
          if ( sub_type(Dirichlet,la_cl) || sub_type(Dirichlet_homogene,la_cl) )
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              int nb_faces_tot = le_bord.nb_faces_tot();
              for (int ind_face=0; ind_face<nb_faces_tot; ind_face++)
                {
                  int num_face = le_bord.num_face(ind_face);
                  est_une_face_de_dirichlet_[num_face] = 1;
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

  const Elem_VEF_base& type_elemvef= domaine_VEF.type_elem();
  int istetra=0;
  Nom nom_elem=type_elemvef.que_suis_je();
  if ((nom_elem=="Tetra_VEF")||(nom_elem=="Tri_VEF"))
    istetra=1;

  double psc, psc_inco;
  int poly,face_adj,fa7,i,j,n_bord;
  int num_face, rang;
  int num10,num20;
  int ncomp_ch_transporte=(inco.nb_dim() == 1?1:inco.dimension(1));

  // Traitement particulier pour les faces de periodicite
  int nb_faces_perio = 0;
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
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
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
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
  int nb_faces_ = domaine_VEF.nb_faces();
  ArrOfInt face(nfac);
  ArrOfDouble vs(dimension);
  ArrOfDouble vs_inco(dimension);
  ArrOfDouble vc(dimension);
  ArrOfDouble vc_inco(dimension);
  ArrOfDouble cc(dimension);
  //DoubleVect xc(dimension);
  DoubleTab vsom(nsom,dimension);
  DoubleTab vsom_inco(nsom,dimension);
  //DoubleTab xsom(nsom,dimension);


  // Dimensionnement du tableau des flux convectifs au bord du domaine de calcul
  DoubleTab& flux_b = flux_bords_;
  int nb_faces_bord=domaine_VEF.nb_faces_bord();
  flux_b.resize(nb_faces_bord,ncomp_ch_transporte);
  flux_b = 0.;

  const IntTab& KEL=type_elemvef.KEL();



  // Les polyedres non standard sont ranges en 2 groupes dans le Domaine_VEF:
  //  - polyedres bords et joints
  //  - polyedres bords et non joints
  // On traite les polyedres en suivant l'ordre dans lequel ils figurent
  // dans le domaine
  // boucle sur les polys

  for (poly=0; poly<nb_elem_tot; poly++)
    {
      int contrib = 0;
      // calcul des numeros des faces du polyedre
      for (face_adj=0; face_adj<nfac; face_adj++)
        {
          int face_ = elem_faces(poly,face_adj);
          face[face_adj]= face_;
          if (face_<nb_faces_) contrib=1; // Une face reelle sur l'element virtuel
        }
      //
      if (contrib)
        {
          for (j=0; j<dimension; j++)
            {
              vs[j] = state_field(face[0],j);
              vs_inco[j] = inco(face[0],j);
              for (i=1; i<nfac; i++)
                {
                  vs[j]+= state_field(face[i],j);
                  vs_inco[j] = inco(face[i],j);
                }
            }
          // calcul de la vitesse aux sommets des polyedres
          if (istetra==1)
            {
              for (i=0; i<nsom; i++)
                for (j=0; j<dimension; j++)
                  {
                    vsom(i,j) = (vs[j] - dimension*state_field(face[i],j));
                    vsom_inco(i,j) = (vs_inco[j] - dimension*inco(face[i],j));
                  }
            }
          else
            {
              Cout << "Sensibility is currently working only with Tetra_VEF (3D) or Tri_VEF (2D)." << finl;
              exit();
            }

          // Determination du type de CL selon le rang
          rang = rang_elem_non_std(poly);
          int itypcl = (rang==-1 ? 0 : domaine_Cl_VEF.type_elem_Cl(rang));

          // calcul de vc (a l'intersection des 3 facettes) vc vs vsom
          calcul_vc(face,vc,vs,vsom,state_field,itypcl);
          calcul_vc(face,vc_inco,vs_inco,vsom_inco,inco,itypcl);


          // Boucle sur les facettes du polyedre non standard:
          for (fa7=0; fa7<nfa7; fa7++)
            {
              num10 = face[KEL(0,fa7)];
              num20 = face[KEL(1,fa7)];

              // normales aux facettes
              if (rang==-1)
                for (i=0; i<dimension; i++)
                  cc[i] = facette_normales(poly, fa7, i);
              else
                for (i=0; i<dimension; i++)
                  cc[i] = normales_facettes_Cl(rang,fa7,i);

              // Calcul des vitesses en C,S,S2 les 3 extremites de la fa7 et M le centre de la fa7
              double psc_c=0,psc_s=0,psc_m,psc_s2=0;
              double psc_c_inco=0,psc_s_inco=0,psc_m_inco=0;
              if (dimension==2)
                {
                  for (i=0; i<2; i++)
                    {
                      psc_c+=vc[i]*cc[i];
                      psc_s+=vsom(KEL(2,fa7),i)*cc[i];
                      psc_c_inco+=vc_inco[i]*cc[i];
                      psc_s_inco+=vsom_inco(KEL(2,fa7),i)*cc[i];

                    }
                  psc_m=(psc_c+psc_s)/2.;
                  psc_m_inco=(psc_c_inco+psc_s_inco)/2.;
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
              if (option_appliquer_cl_dirichlet)
                if (est_une_face_de_dirichlet_[num10] || est_une_face_de_dirichlet_[num20])
                  {
                    psc_m = psc_c;
                    psc_m_inco = psc_c_inco;
                  }
              // Determination de la face amont pour M
              int face_amont_m;
              if(false) Cout<<psc_m_inco<<finl;
              if (psc_m >= 0)
                //if (psc_m_inco >= 0)
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
                  double inco_m = (ncomp_ch_transporte==1?inco(face_amont_m):inco(face_amont_m,comp0));
                  flux = inco_m*psc_m;

                  if (ncomp_ch_transporte == 1)
                    {
                      resu(num10) -= flux;
                      resu(num20) += flux;
                      if (num10<nb_faces_bord) flux_b(num10,0) += flux;
                      if (num20<nb_faces_bord) flux_b(num20,0) -= flux;
                    }
                  else
                    {
                      resu(num10,comp0) -= flux;
                      resu(num20,comp0) += flux;
                      if (num10<nb_faces_bord) flux_b(num10,comp0) += flux;
                      if (num20<nb_faces_bord) flux_b(num20,comp0) -= flux;
                    }

                }// boucle sur comp
            } // fin de la boucle sur les facettes
        }
    } // fin de la boucle
  resu.echange_espace_virtuel();
  int voisine;
  nb_faces_perio = 0;
  double diff1,diff2;

  // Boucle sur les bords pour traiter les conditions aux limites
  // il y a prise en compte d'un terme de convection pour les
  // conditions aux limites de Neumann_sortie_libre seulement
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);

      if (sub_type(Neumann_sortie_libre,la_cl.valeur()))
        {
          const Neumann_sortie_libre& la_sortie_libre = ref_cast(Neumann_sortie_libre, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (num_face=num1; num_face<num2; num_face++)
            {
              psc =0.;
              psc_inco = 0.;
              for (i=0; i<dimension; i++)
                {
                  psc += state_field(num_face,i)*facenormales(num_face,i);
                  psc_inco += inco(num_face,i)*facenormales(num_face,i);
                }
              if (psc>0)
                if (ncomp_ch_transporte == 1)
                  {
                    resu(num_face) -= psc*inco(num_face);
                    flux_b(num_face,0) -= psc*inco(num_face);
                  }
                else
                  for (i=0; i<ncomp_ch_transporte; i++)
                    {
                      resu(num_face,i) -= psc*inco(num_face,i);
                      flux_b(num_face,i) -= psc*inco(num_face,i);
                    }
              else
                {
                  if (ncomp_ch_transporte == 1)
                    {
                      resu(num_face) -= psc_inco*la_sortie_libre.val_ext(num_face-num1);
                      flux_b(num_face,0) -= psc_inco*la_sortie_libre.val_ext(num_face-num1);
                    }
                  else
                    for (i=0; i<ncomp_ch_transporte; i++)
                      {
                        resu(num_face,i) -= psc_inco*la_sortie_libre.val_ext(num_face-num1,i);
                        flux_b(num_face,i) -= psc_inco*la_sortie_libre.val_ext(num_face-num1,i);
                      }
                }
            }
        }
      else if (sub_type(Periodique,la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
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
//inco grad_state_field
void Operateur_Conv_sensibility_VEF::ajouter_Lsensibility_state_Amont(const DoubleTab& inco, const DoubleTab& state_field, DoubleTab& resu ) const
{
  Cout<<"ajouter_Lsensibility_state "<<finl;
  const Domaine_Cl_VEF& domaine_Cl_VEF = la_zcl_vef.valeur();
  const Domaine_VEF& domaine_VEF = ref_cast(Domaine_VEF, le_dom_vef.valeur());
  const IntTab& elem_faces = domaine_VEF.elem_faces();
  const DoubleTab& facenormales = domaine_VEF.face_normales();
  const DoubleTab& facette_normales = domaine_VEF.facette_normales();
  const Domaine& domaine = domaine_VEF.domaine();
  const int nfa7 = domaine_VEF.type_elem().nb_facette();
  const int nb_elem_tot = domaine_VEF.nb_elem_tot();
  const IntVect& rang_elem_non_std = domaine_VEF.rang_elem_non_std();
  const DoubleTab& normales_facettes_Cl = domaine_Cl_VEF.normales_facettes_Cl();
  int nfac = domaine.nb_faces_elem();
  int nsom = domaine.nb_som_elem();
  int nb_bord = domaine_VEF.nb_front_Cl();
  const IntTab& les_elems=domaine.les_elems();
  int option_appliquer_cl_dirichlet = 0 ;


  // Definition d'un tableau pour un traitement special des schemas pres des bords
  if (traitement_pres_bord_.size_array()!=nb_elem_tot)
    {
      traitement_pres_bord_.resize_array(nb_elem_tot);
      traitement_pres_bord_=0;
      ArrOfInt est_un_sommet_de_bord_(domaine_VEF.nb_som_tot());
      for (int n_bord=0; n_bord<nb_bord; n_bord++)
        {
          const Cond_lim_base& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord).valeur();
          if ( sub_type(Dirichlet,la_cl) || sub_type(Dirichlet_homogene,la_cl) )
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              int nb_faces_tot = le_bord.nb_faces_tot();
              int size = domaine_VEF.face_sommets().dimension(1);
              for (int ind_face=0; ind_face<nb_faces_tot; ind_face++)
                for (int som=0; som<size; som++)
                  {
                    int face = le_bord.num_face(ind_face);
                    est_un_sommet_de_bord_[domaine_VEF.face_sommets(face,som)]=1;
                  }
            }
        }
      for (int elem=0; elem<nb_elem_tot; elem++)
        {
          if (rang_elem_non_std(elem)!=-1)
            traitement_pres_bord_[elem]=1;
          else
            {
              for (int n_som=0; n_som<nsom; n_som++)
                if (est_un_sommet_de_bord_[les_elems(elem,n_som)])
                  traitement_pres_bord_[elem]=1;
            }
        }

      // Construction du tableau est_une_face_de_dirichlet_
      est_une_face_de_dirichlet_.resize_array(domaine_VEF.nb_faces_tot());
      est_une_face_de_dirichlet_=0;
      for (int n_bord=0; n_bord<nb_bord; n_bord++)
        {
          const Cond_lim_base& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord).valeur();
          if ( sub_type(Dirichlet,la_cl) || sub_type(Dirichlet_homogene,la_cl) )
            {
              const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
              int nb_faces_tot = le_bord.nb_faces_tot();
              for (int ind_face=0; ind_face<nb_faces_tot; ind_face++)
                {
                  int num_face = le_bord.num_face(ind_face);
                  est_une_face_de_dirichlet_[num_face] = 1;
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

  const Elem_VEF_base& type_elemvef= domaine_VEF.type_elem();
  int istetra=0;
  Nom nom_elem=type_elemvef.que_suis_je();
  if ((nom_elem=="Tetra_VEF")||(nom_elem=="Tri_VEF"))
    istetra=1;

  double psc;
  int poly,face_adj,fa7,i,j,n_bord;
  int num_face, rang;
  int num10,num20;
  int ncomp_ch_transporte=(state_field.nb_dim() == 1?1:state_field.dimension(1));

  // Traitement particulier pour les faces de periodicite
  int nb_faces_perio = 0;
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
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
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
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
  int nb_faces_ = domaine_VEF.nb_faces();
  ArrOfInt face(nfac);
  ArrOfDouble vs(dimension);
  ArrOfDouble vs_state(dimension);
  ArrOfDouble vc(dimension);
  ArrOfDouble vc_state(dimension);
  ArrOfDouble cc(dimension);
  //DoubleVect xc(dimension);
  DoubleTab vsom(nsom,dimension);
  DoubleTab vsom_state(nsom,dimension);
  //DoubleTab xsom(nsom,dimension);


  // Dimensionnement du tableau des flux convectifs au bord du domaine de calcul
  DoubleTab& flux_b = flux_bords_;
  int nb_faces_bord=domaine_VEF.nb_faces_bord();
  flux_b.resize(nb_faces_bord,ncomp_ch_transporte);
  flux_b = 0.;

  const IntTab& KEL=type_elemvef.KEL();

  // Les polyedres non standard sont ranges en 2 groupes dans le Domaine_VEF:
  //  - polyedres bords et joints
  //  - polyedres bords et non joints
  // On traite les polyedres en suivant l'ordre dans lequel ils figurent
  // dans le domaine
  // boucle sur les polys

  for (poly=0; poly<nb_elem_tot; poly++)
    {
      int contrib = 0;
      // calcul des numeros des faces du polyedre
      for (face_adj=0; face_adj<nfac; face_adj++)
        {
          int face_ = elem_faces(poly,face_adj);
          face[face_adj]= face_;
          if (face_<nb_faces_) contrib=1; // Une face reelle sur l'element virtuel
        }
      //
      if (contrib)
        {
          for (j=0; j<dimension; j++)
            {
              vs[j] = inco(face[0],j);
              vs_state[j] = state_field(face[0],j);
              for (i=1; i<nfac; i++)
                {
                  vs[j]+= inco(face[i],j);
                  vs_state[j] = state_field(face[i],j);
                }
            }
          // calcul de la vitesse aux sommets des polyedres
          if (istetra==1)
            {
              for (i=0; i<nsom; i++)
                for (j=0; j<dimension; j++)
                  {
                    vsom(i,j) = (vs[j] - dimension*inco(face[i],j));
                    vsom_state(i,j) = (vs_state[j] - dimension*state_field(face[i],j));
                  }
            }
          else
            {
              Cout << "Sensibility is currently working only with Tetra_VEF (3D) or Tri_VEF (2D)." << finl;
              exit();
            }

          // Determination du type de CL selon le rang
          rang = rang_elem_non_std(poly);
          int itypcl = (rang==-1 ? 0 : domaine_Cl_VEF.type_elem_Cl(rang));

          // calcul de vc (a l'intersection des 3 facettes) vc vs vsom
          calcul_vc(face,vc,vs,vsom,inco,itypcl);
          calcul_vc(face,vc_state,vs_state,vsom_state,state_field,itypcl);

          // Boucle sur les facettes du polyedre non standard:
          for (fa7=0; fa7<nfa7; fa7++)
            {
              num10 = face[KEL(0,fa7)];
              num20 = face[KEL(1,fa7)];

              // normales aux facettes
              if (rang==-1)
                for (i=0; i<dimension; i++)
                  cc[i] = facette_normales(poly, fa7, i);
              else
                for (i=0; i<dimension; i++)
                  cc[i] = normales_facettes_Cl(rang,fa7,i);

              // Calcul des vitesses en C,S,S2 les 3 extremites de la fa7 et M le centre de la fa7
              double psc_c=0,psc_s=0,psc_m,psc_s2=0;
              double psc_m_state=0.,psc_c_state=0.,psc_s_state=0.;
              if (dimension==2)
                {
                  for (i=0; i<2; i++)
                    {
                      psc_c+=vc[i]*cc[i];
                      psc_s+=vsom(KEL(2,fa7),i)*cc[i];
                      psc_c_state+=vc_state[i]*cc[i];
                      psc_s_state+=vsom_state(KEL(2,fa7),i)*cc[i];

                    }
                  psc_m=(psc_c+psc_s)/2.;
                  psc_m_state=(psc_c_state+psc_s_state)/2.;
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
              //int appliquer_cl_dirichlet=0;
              if (option_appliquer_cl_dirichlet)
                if (est_une_face_de_dirichlet_[num10] || est_une_face_de_dirichlet_[num20])
                  {
                    //appliquer_cl_dirichlet = 1;
                    psc_m = psc_c;
                    psc_m_state = psc_c_state;
                  }

              // Determination de la face amont pour M
              int face_amont_m;
              if(false) Cout<<psc_m_state<<finl;
              if (psc_m >= 0)
                //if (psc_m_state >= 0)
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
                  double inco_m = (ncomp_ch_transporte==1?state_field(face_amont_m):state_field(face_amont_m,comp0));

                  flux = inco_m*psc_m;

                  if (ncomp_ch_transporte == 1)
                    {
                      resu(num10) -= flux;
                      resu(num20) += flux;
                      if (num10<nb_faces_bord) flux_b(num10,0) += flux;
                      if (num20<nb_faces_bord) flux_b(num20,0) -= flux;
                    }
                  else
                    {
                      resu(num10,comp0) -= flux;
                      resu(num20,comp0) += flux;
                      if (num10<nb_faces_bord) flux_b(num10,0) += flux;
                      if (num20<nb_faces_bord) flux_b(num20,0) -= flux;
                    }

                }// boucle sur comp
            } // fin de la boucle sur les facettes
        }
    } // fin de la boucle
  resu.echange_espace_virtuel();
  int voisine;
  nb_faces_perio = 0;
  double diff1,diff2;

  // Boucle sur les bords pour traiter les conditions aux limites
  // il y a prise en compte d'un terme de convection pour les
  // conditions aux limites de Neumann_sortie_libre seulement
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);

      if (sub_type(Neumann_sortie_libre,la_cl.valeur()))
        {
          const Neumann_sortie_libre& la_sortie_libre = ref_cast(Neumann_sortie_libre, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int num1 = le_bord.num_premiere_face();
          int num2 = num1 + le_bord.nb_faces();
          for (num_face=num1; num_face<num2; num_face++)
            {
              psc =0;
              for (i=0; i<dimension; i++)
                psc += inco(num_face,i)*facenormales(num_face,i);
              if (psc>0)
                if (ncomp_ch_transporte == 1)
                  {
                    resu(num_face) -= psc*state_field(num_face);
                    flux_b(num_face,0) -= psc*state_field(num_face);
                  }
                else
                  for (i=0; i<ncomp_ch_transporte; i++)
                    {
                      resu(num_face,i) -= psc*state_field(num_face,i);
                      flux_b(num_face,i) -= psc*state_field(num_face,i);
                    }
              else
                {
                  if (ncomp_ch_transporte == 1)
                    {
                      resu(num_face) -= psc*la_sortie_libre.val_ext(num_face-num1);
                      flux_b(num_face) -= psc*la_sortie_libre.val_ext(num_face-num1);
                    }
                  else
                    for (i=0; i<ncomp_ch_transporte; i++)
                      {
                        resu(num_face,i) -= psc*la_sortie_libre.val_ext(num_face-num1,i);
                        flux_b(num_face,i) -= psc*la_sortie_libre.val_ext(num_face-num1,i);
                      }
                }
            }
        }
      else if (sub_type(Periodique,la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
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


void Operateur_Conv_sensibility_VEF::calcul_vc(const ArrOfInt& Face, ArrOfDouble& vc, const ArrOfDouble& vs, const DoubleTab& vsom,
                                               const DoubleTab& vitesse_face,int type_cl) const
{
// Operateur_Conv_sensibility_VEF::calcul_vc(...) is based on Tetra_VEF::calcul_vc(...) and on Tri_VEF::calcul_vc(...).
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
            vc[0]= vitesse_face(Face[2],0);
            vc[1]= vitesse_face(Face[2],1);
            break;
          }

        case 2: // le triangle a une Face de Dirichlet :la Face 1
          {
            vc[0]= vitesse_face(Face[1],0);
            vc[1]= vitesse_face(Face[1],1);
            break;
          }

        case 4: // le triangle a une Face de Dirichlet :la Face 0
          {
            vc[0]= vitesse_face(Face[0],0);
            vc[1]= vitesse_face(Face[0],1);
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
              vc[comp] = vitesse_face(Face[3],comp);
            break;
          }

        case 2: // le tetraedre a une Face de Dirichlet : KEL2
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = vitesse_face(Face[2],comp);
            break;
          }

        case 4: // le tetraedre a une Face de Dirichlet : KEL1
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = vitesse_face(Face[1],comp);
            break;
          }

        case 8: // le tetraedre a une Face de Dirichlet : KEL0
          {
            for (comp=0; comp<3; comp++)
              vc[comp] = vitesse_face(Face[0],comp);
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
void Operateur_Conv_sensibility_VEF::remplir_fluent(DoubleVect& tab_fluent) const
{

  const Domaine_Cl_VEF& domaine_Cl_VEF = la_zcl_vef.valeur();
  const Domaine_VEF& domaine_VEF = ref_cast(Domaine_VEF, le_dom_vef.valeur());
  const Champ_Inc_base& velocity=vitesse();
  const DoubleTab& vitesse_face=velocity.valeurs();
  const IntTab& elem_faces = domaine_VEF.elem_faces();
  const DoubleTab& face_normales = domaine_VEF.face_normales();
  const DoubleTab& facette_normales = domaine_VEF.facette_normales();
  const Domaine& domaine = domaine_VEF.domaine();
  const int nfa7 = domaine_VEF.type_elem().nb_facette();
  const int nb_elem_tot = domaine_VEF.nb_elem_tot();
  const IntVect& rang_elem_non_std = domaine_VEF.rang_elem_non_std();
  const DoubleTab& normales_facettes_Cl = domaine_Cl_VEF.normales_facettes_Cl();
  int nfac = domaine.nb_faces_elem();
  int nsom = domaine.nb_som_elem();
  const IntTab& sommet_elem = domaine.les_elems();
  DoubleVect& fluent_ = tab_fluent;


  const Elem_VEF_base& type_elemvef= domaine_VEF.type_elem();
  int istetra=0;
  Nom nom_elem=type_elemvef.que_suis_je();
  if ((nom_elem=="Tetra_VEF")||(nom_elem=="Tri_VEF"))
    istetra=1;


  double psc;
  int poly,face_adj,fa7,i,j,n_bord;
  int num_face, rang ,itypcl;
  int num1,num2,num_som;

  ArrOfInt face(nfac);
  ArrOfDouble vs(dimension);
  ArrOfDouble vc(dimension);
  DoubleTab vsom(nsom,dimension);
  ArrOfDouble cc(dimension);

  // Dimensionnement du tableau des flux convectifs au bord du domaine de calcul
  const IntTab& KEL=type_elemvef.KEL();

  // On remet a zero le tableau qui sert pour
  // le calcul du pas de temps de stabilite
  fluent_ = 0;

  // Les polyedres non standard sont ranges en 2 groupes dans le Domaine_VEF:
  //  - polyedres bords et joints
  //  - polyedres bords et non joints
  // On traite les polyedres en suivant l'ordre dans lequel ils figurent
  // dans le domaine

  // boucle sur les polys
  for (poly=0; poly<nb_elem_tot; poly++)
    {
      rang = rang_elem_non_std(poly);
      // On cherche, pour un elem qui n'est pas de bord (rang==-1),
      // si un des sommets est sur un bord (tableau des sommets) (C MALOD 17/07/2007)

      if (rang==-1)
        itypcl=0;
      else
        itypcl=domaine_Cl_VEF.type_elem_Cl(rang);

      // calcul des numeros des faces du polyedre
      for (face_adj=0; face_adj<nfac; face_adj++)
        face[face_adj]= elem_faces(poly,face_adj);

      for (j=0; j<dimension; j++)
        {
          vs[j] = vitesse_face(face[0],j);
          for (i=1; i<nfac; i++)
            vs[j]+= vitesse_face(face[i],j);
        }
      // calcul de la vitesse aux sommets des polyedres
      if (istetra==1)
        {
          for (i=0; i<nsom; i++)
            for (j=0; j<dimension; j++)
              vsom(i,j) = (vs[j] - dimension*vitesse_face(face[i],j));
        }
      else
        {
          // pour que cela soit valide avec les hexa (c'est + lent a calculer...)
          int ncomp;
          for (j=0; j<nsom; j++)
            {
              num_som = sommet_elem(poly,j);
              for (ncomp=0; ncomp<dimension; ncomp++)
                vsom(j,ncomp) =velocity.valeur_a_sommet_compo(num_som,poly,ncomp);
            }
        }

      // calcul de vc (a l'intersection des 3 facettes) vc vs vsom proportionnelles a la prosite
      calcul_vc(face,vc,vs,vsom,vitesse().valeurs(),itypcl);
      // Boucle sur les facettes du polyedre non standard:
      for (fa7=0; fa7<nfa7; fa7++)
        {
          num1 = face[KEL(0,fa7)];
          num2 = face[KEL(1,fa7)];
          // normales aux facettes
          if (rang==-1)
            {
              for (i=0; i<dimension; i++)
                cc[i] = facette_normales(poly, fa7, i);
            }
          else
            {
              for (i=0; i<dimension; i++)
                cc[i] = normales_facettes_Cl(rang,fa7,i);
            }
          // On applique le schema de convection a chaque sommet de la facette

          double psc_c=0,psc_s=0,psc_m,psc_s2=0;
          if (dimension==2)
            {
              for (i=0; i<dimension; i++)
                {
                  psc_c+=vc[i]*cc[i];
                  psc_s+=vsom(KEL(2,fa7),i)*cc[i];
                }
              psc_m=(psc_c+psc_s)/2.;
            }
          else
            {
              for (i=0; i<dimension; i++)
                {
                  psc_c+=vc[i]*cc[i];
                  psc_s+=vsom(KEL(2,fa7),i)*cc[i];
                  psc_s2+=vsom(KEL(3,fa7),i)*cc[i];
                }
              psc_m=(psc_c+psc_s+psc_s2)/3.;
            }

          // int amont,dir;
          if (psc_m >= 0)
            {
              // amont = num1;
              fluent_(num2)  += psc_m;
              //dir=0;
            }
          else
            {
              //amont = num2;
              fluent_(num1)  -= psc_m;
              //dir=1;
            }

        } // fin de la boucle sur les facettes
    } // fin de la boucle

  // Boucle sur les bords pour traiter les conditions aux limites
  // il y a prise en compte d'un terme de convection pour les
  // conditions aux limites de Neumann_sortie_libre seulement
  int nb_bord = domaine_VEF.nb_front_Cl();
  for (n_bord=0; n_bord<nb_bord; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Neumann_sortie_libre,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int num1b = le_bord.num_premiere_face();
          int num2b = num1b + le_bord.nb_faces();
          for (num_face=num1b; num_face<num2b; num_face++)
            {
              psc = 0;
              for (i=0; i<dimension; i++)
                psc += vitesse_face(num_face,i)*face_normales(num_face,i);
              if (psc>0)
                ;
              else
                fluent_(num_face) -= psc;
            }
        }
    }
}
double Operateur_Conv_sensibility_VEF::calculer_dt_stab() const
{

  return DMAXFLOAT; //on resout la sensibilite avec le meme dt que l'etat

  const Domaine_Cl_VEF& domaine_Cl_VEF = la_zcl_vef.valeur();
  const Domaine_VEF& domaine_VEF = le_dom_vef.valeur();
  const DoubleVect& volumes_entrelaces = domaine_VEF.volumes_entrelaces();
  const DoubleVect& volumes_entrelaces_Cl = domaine_Cl_VEF.volumes_entrelaces_Cl();
  remplir_fluent(fluent);
  double dt_face,dt_stab =1.e30;

  // On traite les conditions aux limites
  // Si une face porte une condition de Dirichlet on n'en tient pas compte
  // dans le calcul de dt_stab
  for (int n_bord=0; n_bord<domaine_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      if ( (sub_type(Dirichlet,la_cl.valeur()))
           ||
           (sub_type(Dirichlet_homogene,la_cl.valeur()))
         )
        ;
      else
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
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
  int ndeb = domaine_VEF.premiere_face_int();
  int nfin = domaine_VEF.premiere_face_std();

  for (int num_face=ndeb; num_face<nfin; num_face++)
    {
      dt_face = volumes_entrelaces_Cl(num_face)/(fluent[num_face]+DMINFLOAT);
      dt_stab =(dt_face < dt_stab) ? dt_face : dt_stab;
    }

  // On traite les faces internes standard
  ndeb = nfin;
  nfin = domaine_VEF.nb_faces();
  for (int num_face=ndeb; num_face<nfin; num_face++)
    {
      dt_face = volumes_entrelaces(num_face)/(fluent[num_face]+DMINFLOAT);
      dt_stab =(dt_face < dt_stab) ? dt_face : dt_stab;
    }

  dt_stab = Process::mp_min(dt_stab);
  // astuce pour contourner le type const de la methode
  Operateur_Conv_sensibility_VEF& op = ref_cast_non_const(Operateur_Conv_sensibility_VEF,*this);
  op.fixer_dt_stab_conv(dt_stab);
  return dt_stab;
}



void  Operateur_Conv_sensibility_VEF::add_diffusion_term(const DoubleTab& state, DoubleTab& resu) const
{
  Cerr << "add_diffusion_term" << finl;
  const Domaine_Cl_VEF& domaine_Cl_VEF = la_zcl_vef.valeur();
  const Domaine_VEF& domaine_VEF = le_dom_vef.valeur();

  int nb_comp = 1;
  int nb_dim = resu.nb_dim();
  if(nb_dim==2)
    nb_comp=resu.dimension(1);

  DoubleTab& tab_flux_bords = flux_bords_;

  const IntTab& elemfaces = domaine_VEF.elem_faces();
  const IntTab& face_voisins = domaine_VEF.face_voisins();
  int i0,j,num_face;
  int nb_faces = domaine_VEF.nb_faces();
  int nb_faces_elem = domaine_VEF.domaine().nb_faces_elem();
  int n_bord0;
  double valA;//,flux;
  //DoubleVect n(Objet_U::dimension);
  //DoubleTrav Tgrad(Objet_U::dimension,Objet_U::dimension);

  // On dimensionne et initialise le tableau des bilans de flux:
  tab_flux_bords.resize(domaine_VEF.nb_faces_bord(),nb_comp);
  tab_flux_bords=0.;

  assert(nb_comp>1);
  int nb_bords=domaine_VEF.nb_front_Cl();
  int ind_face;

  for (n_bord0=0; n_bord0<nb_bords; n_bord0++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord0);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
      // const IntTab& elemfaces = domaine_VEF.elem_faces();
      int num1 = 0;
      int num2 = le_bord.nb_faces_tot();
      int nb_faces_bord_reel = le_bord.nb_faces();

      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique,la_cl.valeur());
          int fac_asso;
          for (ind_face=num1; ind_face<nb_faces_bord_reel; ind_face++)
            {
              fac_asso = la_cl_perio.face_associee(ind_face);
              fac_asso = le_bord.num_face(fac_asso);
              num_face = le_bord.num_face(ind_face);
              for (int kk=0; kk<2; kk++)
                {
                  int elem = face_voisins(num_face, kk);
                  for (i0=0; i0<nb_faces_elem; i0++)
                    {
                      if ( ( (j= elemfaces(elem,i0)) > num_face ) && (j != fac_asso ) )
                        {
                          valA = viscA(num_face,j,elem);
                          for (int nc=0; nc<nb_comp; nc++)
                            {
                              resu(num_face,nc)+=valA*state(j,nc);
                              resu(num_face,nc)-=valA*state(num_face,nc);
                              if(j<nb_faces) // face reelle
                                {
                                  resu(j,nc)+=0.5*valA*state(num_face,nc);
                                  resu(j,nc)-=0.5*valA*state(j,nc);
                                }
                            }
                        }
                    }
                }

            }
        }// fin if periodique
      else
        {
          for (ind_face=num1; ind_face<num2; ind_face++)
            {
              num_face = le_bord.num_face(ind_face);
              int elem=face_voisins(num_face,0);

              // Boucle sur les faces :
              for (int i=0; i<nb_faces_elem; i++)
                if (( (j= elemfaces(elem,i)) > num_face ) || (ind_face>=nb_faces_bord_reel))
                  {
                    valA = viscA(num_face,j,elem);
                    for (int nc=0; nc<nb_comp; nc++)
                      {
                        double flux=valA*(state(j,nc)-state(num_face,nc));
                        if (ind_face<nb_faces_bord_reel)
                          {
                            resu(num_face,nc)+=flux;
                            tab_flux_bords(num_face,nc)+=flux;
                          }

                        if(j<nb_faces) // face reelle
                          {
                            resu(j,nc)-=flux;
                          }
                      }
                  }
            }
        }
    }//Fin for n_bord

  // On traite les faces internes

  for (num_face=domaine_VEF.premiere_face_int(); num_face<nb_faces; num_face++)
    {
      for (int k=0; k<2; k++)
        {
          int elem = face_voisins(num_face,k);
          for (i0=0; i0<nb_faces_elem; i0++)
            {
              if ( (j= elemfaces(elem,i0)) > num_face )
                {
                  int el1,el2;
                  int contrib=1;
                  if(j>=nb_faces) // C'est une face virtuelle
                    {
                      el1 = face_voisins(j,0);
                      el2 = face_voisins(j,1);
                      if((el1==-1)||(el2==-1))
                        contrib=0;
                    }
                  if(contrib)
                    {
                      valA = viscA(num_face,j,elem);
                      for (int nc=0; nc<nb_comp; nc++)
                        {
                          resu(num_face,nc)+=valA*state(j,nc);
                          resu(num_face,nc)-=valA*state(num_face,nc);
                          if(j<nb_faces) // On traite les faces reelles
                            {
                              resu(j,nc)+=valA*state(num_face,nc);
                              resu(j,nc)-=valA*state(j,nc);
                            }
                          else
                            {
                              // La face j est virtuelle
                            }
                        }
                    }
                }
            }
        }
    }// Fin faces internes

  for (int n_bord=0; n_bord<nb_bords; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Symetrie,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();
          for (int face=ndeb; face<nfin; face++)
            tab_flux_bords(face,0) = 0.;
        }
    }


}


inline double Operateur_Conv_sensibility_VEF::viscA(int i, int j, int num_elem, double diffu) const
{
  const Domaine_VEF& domaine=le_dom_vef.valeur();
  const IntTab& face_voisins=domaine.face_voisins();
  const DoubleTab& face_normales=domaine.face_normales();
  const DoubleVect& inverse_volumes=domaine.inverse_volumes();
  double pscal = face_normales(i,0)*face_normales(j,0)
                 + face_normales(i,1)*face_normales(j,1);
  if (Objet_U::dimension == 3)
    pscal += face_normales(i,2)*face_normales(j,2);


  if ( (face_voisins(i,0) == face_voisins(j,0)) ||
       (face_voisins(i,1) == face_voisins(j,1)) )
    return -(pscal*diffu)*inverse_volumes(num_elem);
  else
    return (pscal*diffu)*inverse_volumes(num_elem);
}

void Operateur_Conv_sensibility_VEF::add_diffusion_scalar_term(const DoubleTab& inconnue, DoubleTab& resu, double diffu) const
{

  Cerr << "add_diffusion_scalar_term" << finl;
  const Domaine_Cl_VEF& domaine_Cl_VEF = la_zcl_vef.valeur();
  const Domaine_VEF& domaine_VEF = le_dom_vef.valeur();
  DoubleTab& tab_flux_bords = flux_bords_;

  const IntTab& elemfaces = domaine_VEF.elem_faces();
  const IntTab& face_voisins = domaine_VEF.face_voisins();
  int i,j,num_face;
  int nb_faces = domaine_VEF.nb_faces();
  int nb_faces_elem = domaine_VEF.domaine().nb_faces_elem();
  double valA,flux;
  int n_bord, ind_face;
  int nb_bords=domaine_VEF.nb_front_Cl();
  // On dimensionne et initialise le tableau des bilans de flux:
  tab_flux_bords.resize(domaine_VEF.nb_faces_bord(),1);
  tab_flux_bords=0.;
  const int premiere_face_int=domaine_VEF.premiere_face_int();

  // On traite les faces bord
  for (n_bord=0; n_bord<nb_bords; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
      //const IntTab& elemfaces = domaine_VEF.elem_faces();
      int num1=0;
      int num2=le_bord.nb_faces_tot();
      int nb_faces_bord_reel = le_bord.nb_faces();

      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique,la_cl.valeur());
          int fac_asso;
          for (ind_face=num1; ind_face<nb_faces_bord_reel; ind_face++)
            {
              num_face = le_bord.num_face(ind_face);
              fac_asso = la_cl_perio.face_associee(ind_face);
              fac_asso = le_bord.num_face(fac_asso);
              for (int kk=0; kk<2; kk++)
                {
                  int elem = face_voisins(num_face,kk);
                  for (i=0; i<nb_faces_elem; i++)
                    {
                      if ( ( (j= elemfaces(elem,i)) > num_face ) && (j != fac_asso) )
                        {
                          valA = viscA(num_face,j,elem, diffu);
                          resu(num_face)+=valA*inconnue(j);
                          resu(num_face)-=valA*inconnue(num_face);
                          if(j<nb_faces) // face reelle
                            {
                              resu(j)+=0.5*valA*inconnue(num_face);
                              resu(j)-=0.5*valA*inconnue(j);
                            }
                        }
                    }
                }
            }
        }
      else   // Il n'y a qu'une seule composante, donc on traite
        // une equation scalaire (pas la vitesse) on a pas a utiliser
        // le tau tangentiel (les lois de paroi thermiques ne calculent pas
        // d'echange turbulent a la paroi pour l'instant
        {
          for (ind_face=num1; ind_face<num2; ind_face++)
            {
              num_face = le_bord.num_face(ind_face);
              int elem = face_voisins(num_face,0);

              for (i=0; i<nb_faces_elem; i++)
                {
                  if (( (j= elemfaces(elem,i)) > num_face ) || (ind_face>=nb_faces_bord_reel))
                    {
                      valA = viscA(num_face,j,elem,diffu);

                      if (ind_face<nb_faces_bord_reel)
                        {
                          flux=valA*(inconnue(j)-inconnue(num_face));
                          // PL : c'est bien un - ici pour flux_bords. Cette valeur est ensuite
                          // ecrasee pour les bords avec Neumann et ne sert donc que pour les bords
                          // avec Dirichlet ou le volume de controle est nul
                          tab_flux_bords(num_face,0)-=flux;
                          resu(num_face)+=flux;
                        }

                      if(j<nb_faces) // face reelle
                        {
                          flux=valA*(inconnue(num_face)-inconnue(j));
                          if (j<premiere_face_int)
                            tab_flux_bords(j,0)-=flux;
                          resu(j)+=flux;
                        }
                    }
                }
            }
        }
    }

  // Faces internes :
  for (num_face=premiere_face_int; num_face<nb_faces; num_face++)
    {
      for (int k=0; k<2; k++)
        {
          int elem = face_voisins(num_face,k);
          {
            for (i=0; i<nb_faces_elem; i++)
              {
                j=elemfaces(elem,i);
                if ( j  > num_face )
                  {
                    int el1,el2;
                    int contrib=1;
                    if(j>=nb_faces) // C'est une face virtuelle
                      {
                        el1 = face_voisins(j,0);
                        el2 = face_voisins(j,1);
                        if((el1==-1)||(el2==-1))
                          contrib=0;
                      }
                    if(contrib)
                      {
                        valA = viscA(num_face,j,elem,diffu);
                        resu(num_face)+=valA*inconnue(j);
                        resu(num_face)-=valA*inconnue(num_face);
                        if(j<nb_faces) // On traite les faces reelles
                          {
                            resu(j)+=valA*inconnue(num_face);
                            resu(j)-=valA*inconnue(j);
                          }
                      }
                  }
              }
          }
        }
    }

  // Neumann :
  for (n_bord=0; n_bord<nb_bords; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);

      if (sub_type(Neumann_paroi,la_cl.valeur()))
        {
          const Neumann_paroi& la_cl_paroi = ref_cast(Neumann_paroi, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();
          for (int face=ndeb; face<nfin; face++)
            {
              flux=la_cl_paroi.flux_impose(face-ndeb)*domaine_VEF.surface(face);
              resu[face] += flux;
              tab_flux_bords(face,0) = flux;
            }
        }
      else if (sub_type(Echange_externe_impose,la_cl.valeur()))
        {
          const Echange_externe_impose& la_cl_paroi = ref_cast(Echange_externe_impose, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();
          for (int face=ndeb; face<nfin; face++)
            {
              flux=la_cl_paroi.h_imp(face-ndeb)*(la_cl_paroi.T_ext(face-ndeb)-inconnue(face))*domaine_VEF.surface(face);
              resu[face] += flux;
              tab_flux_bords(face,0) = flux;
            }
        }
      else if (sub_type(Neumann_homogene,la_cl.valeur())
               || sub_type(Symetrie,la_cl.valeur())
               || sub_type(Neumann_sortie_libre,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl->frontiere_dis());
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();
          for (int face=ndeb; face<nfin; face++)
            tab_flux_bords(face,0) = 0.;
        }
    }
}

double Operateur_Conv_sensibility_VEF::application_LIMITEUR(double grad1, double grad2, Motcle& type_limit) const
{
  double gradlim=0.;
  if (type_limit=="minmod")
    gradlim=minmod(grad1,grad2);
  else  if (type_limit=="vanleer")
    gradlim=vanleer(grad1,grad2);
  else  if (type_limit=="vanalbada")
    gradlim=vanalbada(grad1,grad2);
  else  if (type_limit=="chakravarthy")
    gradlim=chakravarthy(grad1,grad2);
  else if (type_limit=="superbee")
    gradlim=superbee(grad1,grad2);
  else
    {
      Cerr << type_limit << " is not implemented. " << finl;
      Cerr << " Choose from: minmod - vanleer - vanalbada - chakravarthy - superbee " << finl;
      exit();
    }

  return gradlim;
}
