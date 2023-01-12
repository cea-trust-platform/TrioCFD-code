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
// File:        Navier_Stokes_QC_impl.cpp
// Directory:   $TRIO_U_ROOT/src/ThHyd/Quasi_Compressible
// Version:     /main/47
//
//////////////////////////////////////////////////////////////////////////////

#include <Navier_Stokes_QC_impl.h>
#include <Navier_Stokes_std.h>
#include <Fluide_Quasi_Compressible.h>
#include <Debog.h>
#include <Dirichlet.h>
#include <Zone_VF.h>
#include <Schema_Temps_base.h>
#include <Schema_Temps.h>
#include <TRUSTTrav.h>
#include <Champs_compris.h>
#include <Transport_Interfaces_base.h>
#include <Probleme_base.h>
#include <Loi_Etat_GP.h>
#include <Champ_Fonc_Face_VDF.h>
#include <Assembleur_P_VDF.h>
#include <Matrice_Morse.h>
#include <Multigrille_Adrien.h>
#include <Statistiques.h>

#ifdef CHECKOPERATEUR_CENTRE4_IJK
#include <OpCentre4IJK.h>
#include <IJK_discretization.h>
#include <Interprete_bloc.h>
#include <VDF_to_IJK.h>

#endif

extern Stat_Counter_Id assemblage_sys_counter_;


DoubleTab& rho_vitesse_impl(const DoubleTab& tab_rho,const DoubleTab& vitesse,DoubleTab& rhovitesse)
{
  int i,j, n = vitesse.dimension(0);
  if (vitesse.nb_dim()==1)
    {
      for (i=0 ; i<n ; i++)
        {
          rhovitesse(i) = tab_rho(i)*vitesse(i);
        }
    }
  else
    {
      for (i=0 ; i<n ; i++)
        {
          for (j=0 ; j<Objet_U::dimension ; j++)
            {
              rhovitesse(i,j) = tab_rho(i)*vitesse(i,j);
            }
        }
    }
  rhovitesse.echange_espace_virtuel();
  Debog::verifier("Navier_Stokes_QC::rho_vitesse : ", rhovitesse);

  return rhovitesse;
}


/*! @brief Calcule la derivee en temps de l'inconnue vitesse, i.
 *
 * e. l'acceleration dU/dt et la renvoie.
 *     Appelle Equation_base::derivee_en_temps_inco(DoubleTab& )
 *     Calcule egalement la pression.
 *
 * @param (DoubleTab& vpoint) le tableau des valeurs de l'acceleration dU/dt
 * @return (DoubleTab&) le tableau des valeurs de l'acceleration (derivee de la vitesse)
 */
DoubleTab& Navier_Stokes_QC_impl::derivee_en_temps_inco_impl(Navier_Stokes_std& eqn,DoubleTab& vpoint, Fluide_Incompressible& le_fluide, Matrice& matrice_pression_,Assembleur& assembleur_pression_)
{
  return derivee_en_temps_inco(eqn,vpoint,le_fluide,matrice_pression_,assembleur_pression_,0);
}

DoubleTab& Navier_Stokes_QC_impl::derivee_en_temps_inco_diffusion_impl(Navier_Stokes_std& eqn,DoubleTab& vpoint, Fluide_Incompressible& le_fluide, Matrice& matrice_pression_,Assembleur& assembleur_pression_)
{
  return derivee_en_temps_inco(eqn,vpoint,le_fluide,matrice_pression_,assembleur_pression_,1);
}
#ifdef CHECKOPERATEUR_CENTRE4_IJK
void compare(const IJK_Field_double& ref, const IJK_Field_double& x, const char *msg)
{
  int ni = x.ni();
  int nj = x.nj();
  int nk = x.nk();
  int i0 = x.get_splitting().get_offset_local(DIRECTION_I);
  int j0 = x.get_splitting().get_offset_local(DIRECTION_J);
  int k0 = x.get_splitting().get_offset_local(DIRECTION_K);
  double maxd = 0.;
  for(int k=0; k<nk; k++)
    for(int j=0; j<nj; j++)
      for(int i=0; i<ni; i++)
        {
          double d = fabs(x(i,j,k)-ref(i,j,k));
          maxd = max(d, maxd);
          if (d > 1e-15)
            {
              char s[1000];
              snprintf(s, 1000, " local_ijk=( %3d %3d %3d )  global_ijk=( %3d %3d %3d ) ref= %10g  current= %10g  delta=%5g",
                      i, j, k, i+i0, j+j0, k+k0, ref(i,j,k), x(i,j,k), d);
              Process::Journal() << msg << s << finl;
            }
        }
  maxd = Process::mp_max(maxd);
  Process::Journal() << msg << " max_error= " << maxd << finl;
  Cerr << msg << " max_error= " << maxd << finl;
}
#endif

DoubleTab&  Navier_Stokes_QC_impl::derivee_en_temps_impl_p1(Navier_Stokes_std& eqn,DoubleTab& vpoint, Fluide_Incompressible& le_fluide, Matrice& matrice_pression_,Assembleur& assembleur_pression_ ,int diffusion_implicite)
{
  int i, n = vpoint.dimension(0);
  DoubleTab& pression=eqn.pression().valeurs();
  DoubleTab& vitesse=eqn.vitesse().valeurs();

  Fluide_Quasi_Compressible& fluide_QC=ref_cast(Fluide_Quasi_Compressible,le_fluide);

  //F.A Passage a Rho_n+1 pour essaye d'augmenter l'odre du shema en temps
  //const DoubleTab& tab_rho_face_n=fluide_QC.rho_face_n();
  const DoubleTab& tab_rho_face_np1=fluide_QC.rho_face_np1(); // on fera propre si ca marche

  vpoint=0;

  //Resolution vitesse
  //ajout diffusion (avec la viscosite dynamique)
  //Debog::verifier("Navier_Stokes_QC::derivee_en_temps_inco, vpoint DEBUT : ", vpoint);
  if (!diffusion_implicite)
    eqn.operateur(0).ajouter(vpoint);
  //Debog::verifier("Navier_Stokes_QC::derivee_en_temps_inco, vpoint OP0 : ", vpoint);

  vpoint.echange_espace_virtuel();
  Debog::verifier("Navier_Stokes_QC::derivee_en_temps_inco, OP0 2 : ", vpoint);

  DoubleTab& rhovitesse = ref_cast_non_const(DoubleTab,eqn.rho_la_vitesse().valeurs());
  rho_vitesse_impl(tab_rho_face_np1,vitesse,rhovitesse);

  //ajout convection utilise rhovitesse
  if (!diffusion_implicite)
    eqn.operateur(1).ajouter(rhovitesse,vpoint);
  else
    {

      DoubleTrav trav(vpoint);
      eqn.derivee_en_temps_conv(trav,rhovitesse);
      vpoint=trav;
    }
  vpoint.echange_espace_virtuel();
  Debog::verifier("Navier_Stokes_QC::derivee_en_temps_inco, CONV vpoint : ", vpoint);

  //ajout source
  eqn.sources().ajouter(vpoint);
  vpoint.echange_espace_virtuel(); //les sources ne semblent rien changer ds le cas source de chaleur volumique
  Debog::verifier("Navier_Stokes_QC::derivee_en_temps_inco,  SOURCES vpoint : ", vpoint);

  // ajout de u div(rho u)

  DoubleTrav divrhou(pression);
  DoubleTrav divrhou_face(vitesse);

  eqn.operateur_divergence().ajouter(rhovitesse,divrhou);
  divrhou.echange_espace_virtuel();
  ref_cast(Fluide_Quasi_Compressible,le_fluide).divu_discvit(divrhou,divrhou_face);

  for (i=0 ; i<n ; i++)
    vpoint(i)+= vitesse(i)*divrhou_face(i);

#ifdef CHECKOPERATEUR_CENTRE4_IJK
  {
    DoubleTrav resu(vpoint);
    eqn.operateur(1).ajouter(rhovitesse,resu);
    resu.echange_espace_virtuel();
    for (i=0 ; i<n ; i++)
      resu(i)+= vitesse(i) * divrhou_face(i);

    // Calcul avec l'operateur centre4:
    Cerr << "Check operateur centre4:" << finl;
    // fetch the vdf_to_ijk translator (assume there is one unique object, with conventional name)
    const char * ijkdis_name = IJK_discretization::get_conventional_name();
    const IJK_discretization& ijkdis = ref_cast(IJK_discretization, Interprete_bloc::objet_global(ijkdis_name));
    const IJK_Splitting& split = ijkdis.get_IJK_splitting();
    OpConvCentre4IJK_double op_ijk_;
    op_ijk_.initialize(split);
    IJK_Field_double inputx_, inputy_, inputz_, vx_, vy_, vz_;
    IJK_Field_double dvx_, dvy_, dvz_, divrhou_, divrhou_ref_;


    const VDF_to_IJK& vdf_to_ijk_i = ijkdis.get_vdf_to_ijk(IJK_Splitting::FACES_I);
    const VDF_to_IJK& vdf_to_ijk_j = ijkdis.get_vdf_to_ijk(IJK_Splitting::FACES_J);
    const VDF_to_IJK& vdf_to_ijk_k = ijkdis.get_vdf_to_ijk(IJK_Splitting::FACES_K);
    const VDF_to_IJK& vdf_to_ijk = ijkdis.get_vdf_to_ijk(IJK_Splitting::ELEM);
    vx_.allocate(split, IJK_Splitting::FACES_I, 1);
    vy_.allocate(split, IJK_Splitting::FACES_J, 1);
    vz_.allocate(split, IJK_Splitting::FACES_K, 1);
    inputx_.allocate(split, IJK_Splitting::FACES_I, 2);
    inputy_.allocate(split, IJK_Splitting::FACES_J, 2);
    inputz_.allocate(split, IJK_Splitting::FACES_K, 2);
    dvx_.allocate(split, IJK_Splitting::FACES_I, 0);
    dvy_.allocate(split, IJK_Splitting::FACES_J, 0);
    dvz_.allocate(split, IJK_Splitting::FACES_K, 0);
    divrhou_.allocate(split, IJK_Splitting::ELEM, 2);
    divrhou_ref_.allocate(split, IJK_Splitting::ELEM, 2);
    vdf_to_ijk_i.convert_to_ijk(rhovitesse, inputx_);
    vdf_to_ijk_j.convert_to_ijk(rhovitesse, inputy_);
    vdf_to_ijk_k.convert_to_ijk(rhovitesse, inputz_);
    vdf_to_ijk_i.convert_to_ijk(vitesse, vx_);
    vdf_to_ijk_j.convert_to_ijk(vitesse, vy_);
    vdf_to_ijk_k.convert_to_ijk(vitesse, vz_);
    vx_.echange_espace_virtuel(vx_.ghost());
    vy_.echange_espace_virtuel(vy_.ghost());
    vz_.echange_espace_virtuel(vz_.ghost());
    inputx_.echange_espace_virtuel(inputx_.ghost());
    inputy_.echange_espace_virtuel(inputy_.ghost());
    inputz_.echange_espace_virtuel(inputz_.ghost());
    dvx_.data() = 0;
    dvy_.data() = 0;
    dvz_.data() = 0;
    divrhou_.data() = 0;
    op_ijk_.ajouter_avec_u_div_rhou(inputx_, inputy_, inputz_,
                                    vx_, vy_, vz_,
                                    dvx_, dvy_, dvz_ ,divrhou_);
    vdf_to_ijk_i.convert_to_ijk(resu, vx_);
    vdf_to_ijk_j.convert_to_ijk(resu, vy_);
    vdf_to_ijk_k.convert_to_ijk(resu, vz_);
    vdf_to_ijk.convert_to_ijk(divrhou, divrhou_ref_);
    Nom msg_("comparaison centre4 ");
    compare(vx_, dvx_, msg_ + "vx");
    compare(vy_, dvy_, msg_ + "vy");
    compare(vz_, dvz_, msg_ + "vz");
    compare(divrhou_ref_, divrhou_, msg_ + "divrhou");
  }
#endif

// F.A 5/04/11 division par v* rho pour revenir a dudt homogene a du/dt
  const Champ_base& rho_vit=eqn.get_champ("rho_comme_v");
  ref_cast_non_const(DoubleTab,rho_vit.valeurs())=tab_rho_face_np1;
  if (diffusion_implicite)
    {
      abort();
    }
  else
    {
      // on va diviser par M* rho
      eqn.solv_masse()->set_name_of_coefficient_temporel("rho_comme_v");

      eqn.solv_masse().appliquer(vpoint);

      eqn.solv_masse()->set_name_of_coefficient_temporel("no_coeff");
    }

  const Conds_lim& lescl=eqn.zone_Cl_dis().les_conditions_limites();
  const IntTab& face_voisins = eqn.zone_dis().valeur().face_voisins();
  int nbcondlim=lescl.size();
  int taille=vpoint.nb_dim();
  if (taille==1)
    {
      if (orientation_VDF_.size()==0)
        orientation_VDF_.ref(ref_cast(Zone_VF,eqn.zone_dis().valeur()).orientation());
    }

  double dt_reel= eqn.inconnue().valeur().recuperer_temps_futur(1)-eqn.schema_temps().temps_courant();
  double dt_ = eqn.schema_temps().pas_de_temps();

  // on revient au pas de temps du schema, c.a.d au pas de temps intermediaire

  dt_reel=dt_;
  nbcondlim=0;
  for (int icl=0; icl<nbcondlim; icl++)
    {
      const Cond_lim_base& la_cl_base = lescl[icl].valeur();
      if (sub_type(Dirichlet,la_cl_base))
        {
          const Front_VF& la_front_dis = ref_cast(Front_VF,la_cl_base.frontiere_dis());
          const Dirichlet& diri=ref_cast(Dirichlet,la_cl_base);
          int ndeb = la_front_dis.num_premiere_face();
          int nfin = ndeb + la_front_dis.nb_faces();

          if (taille==1)
            {
              for (int num_face=ndeb; num_face<nfin; num_face++)
                {
                  int n0 = face_voisins(num_face, 0);
                  if (n0 == -1)
                    n0 = face_voisins(num_face, 1);

                  vpoint(num_face)=(diri.val_imp(num_face-ndeb,orientation_VDF_(num_face))- vitesse(num_face))/dt_reel;
                }
            }
          else
            {
              for (int num_face=ndeb; num_face<nfin; num_face++)
                for (int jj=0; jj<Objet_U::dimension; jj++)
                  {
                    vpoint(num_face,jj)=(diri.val_imp(num_face-ndeb,jj) - vitesse(num_face,jj))/dt_reel;
                  }
            } // fin else (taille==1)
        } // fin if (sub_type(Dirichlet,la_cl_base))
    } // fin  for (int icl=0;icl<nbcondlim;icl++)





// fin ajout

  return vpoint;
}

DoubleTab&  Navier_Stokes_QC_impl::derivee_en_temps_impl_projection(Navier_Stokes_std& eqn,DoubleTab& vpoint, Fluide_Incompressible& le_fluide, Matrice& matrice_pression_,Assembleur& assembleur_pression_ ,int diffusion_implicite)
{
  DoubleTab& pression=eqn.pression().valeurs();
  DoubleTab& vitesse=eqn.vitesse().valeurs();
  DoubleTab secmem(pression);
  REF(Champ_base) gradient_pression;

  try
    {
      gradient_pression = eqn.get_champ("gradient_pression");
    }
  catch (Champs_compris_erreur)
    {
      Cerr<<" l'equation ne comprend pas gradient_pression "<<finl;
      Process::exit();
    }

  DoubleTab& gradP=gradient_pression.valeur().valeurs();
  eqn.operateur_gradient().calculer(pression, gradP); // F.A 8/04/11

  DoubleTrav inc_pre(pression);
  DoubleTrav Mmoins1grad(gradP);

  if (tab_W.size()==0) tab_W=secmem;

  Fluide_Quasi_Compressible& fluide_QC=ref_cast(Fluide_Quasi_Compressible,le_fluide);
  // F.A Passage a Rho_n+1 pour essayer d'augmenter l'odre du schema en temps
  // const DoubleTab& tab_rho_face_n=fluide_QC.rho_face_n();
  const DoubleTab& tab_rho_face_np1=fluide_QC.rho_face_np1();
  double dt_ = eqn.schema_temps().pas_de_temps();
  const Champ_base& rho_vit=eqn.get_champ("rho_comme_v");
  ref_cast_non_const(DoubleTab,rho_vit.valeurs())=tab_rho_face_np1;

  assert(diffusion_implicite==eqn.schema_temps().diffusion_implicite());

  // ajout de gradP
  // F.A 5/04/11
// eqn.corriger_derivee_expl(vpoint);

  eqn.solv_masse()->set_name_of_coefficient_temporel("rho_comme_v");
  eqn.solv_masse().appliquer(gradP);
  vpoint -= gradP;
  eqn.solv_masse()->set_name_of_coefficient_temporel("no_coeff");

  //Resolution pression
  vpoint.echange_espace_virtuel();
  Debog::verifier("Navier_Stokes_QC::derivee_en_temps_inco, MASSE vpoint : ", vpoint);

  tab_W=0;

  eqn.probleme().equation(1).operateur(0).ajouter(fluide_QC.temperature(),tab_W);  //devrait calculer div.(lambda grad(T))
  eqn.probleme().equation(1).sources().ajouter(tab_W);//devrait ajouter la source volumique pour projection div(u)



  const Loi_Etat_GP& la_loi = ref_cast(Loi_Etat_GP,fluide_QC.loi_etat().valeur());
  const double R = la_loi.R();
  const double Cp = la_loi.Cp();
  const double Pth = fluide_QC.pression_th();
  //const  double Pthn = fluide_QC.pression_thn();
  double gamma_moins_1=R/(Cp-R);
  // Corrections pour equation
  //double dpth=(Pth-Pthn)/dt_;
  //dpth=(Pth-Pthn)/dt_;
  //tab_W-=dpth/gamma_moins_1;
  const DoubleTab& flux_bord =  eqn.probleme().equation(1).operateur(0).l_op_base().flux_bords();
  int size=flux_bord.dimension(0);
  double bilan=0;
  for (int f=0; f<size; f++)
    {
      bilan+=flux_bord(f,0);
    }
  bilan=Process::mp_sum(bilan);

  const DoubleVect& volumes=ref_cast(Zone_VF,eqn.zone_dis().valeur()).volumes();
  double vtot=0;
  int i;

  for ( i=0 ; i<secmem.dimension(0) ; i++)
    {
      vtot+=volumes(i);
    }
  vtot=Process::mp_sum(vtot);

  double dpth2=(gamma_moins_1/vtot)*bilan;

  for (i=0 ; i<secmem.dimension_tot(0) ; i++)
    {
      tab_W(i)= R/(Cp*Pth)*(tab_W(i)-dpth2/gamma_moins_1*volumes(i)); // volumes ici !!
      secmem(i) = -tab_W(i);
    }

  Debog::verifier("C::derivee_en_temps_inco, tab_W : ", tab_W);

  secmem.echange_espace_virtuel();
  Debog::verifier("Navier_Stokes_QC::derivee_en_temps_inco, secmem0 : ", secmem);

  //ici vpoint=d/dt(u)_1
  eqn.operateur_divergence().ajouter(vitesse,secmem);
  secmem /= dt_;
  eqn.operateur_divergence().ajouter(vpoint,secmem);


  secmem*=-1;
  secmem.echange_espace_virtuel();
  Debog::verifier("Navier_Stokes_QC::derivee_en_temps_inco, DIV secmem : ", secmem);

  //On ne fait appel qu une seule fois a assembler dans preparer calcul (au lieu de assembler_QC)
  // Correction du second membre d'apres les conditions aux limites :
  assembleur_pression_.modifier_secmem(secmem);
  Debog::verifier("Navier_Stokes_QC::derivee_en_temps_inco, modifier secmem : ", secmem);
  Debog::verifier("Navier_Stokes_QC::derivee_en_temps_inco, pression avant : ", pression);

  //on modifie la matrice de pression avec div(gradP/rho)
  statistiques().begin_count(assemblage_sys_counter_);

  Assembleur_P_VDF& l_assembleur = ref_cast(Assembleur_P_VDF,assembleur_pression_.valeur());

  l_assembleur.assembler_rho_variable(matrice_pression_, ref_cast(Champ_Don_base,rho_vit));
  l_assembleur.set_resoudre_increment_pression(1);
  l_assembleur.set_resoudre_en_u(0);

  statistiques().end_count(assemblage_sys_counter_);
// On a modifie la matrice, il faut reinitialiser le solveur
// (recalcul des preconditionnement, factorisation pour Cholesky, etc...)
  eqn.solveur_pression().valeur().reinit();

  if (sub_type(Multigrille_Adrien, eqn.solveur_pression().valeur()))
    {
      Multigrille_Adrien& mg = ref_cast(Multigrille_Adrien, eqn.solveur_pression().valeur());
      mg.set_rho(fluide_QC.rho_np1()); // on passe a np1
    }

  eqn.solveur_pression().resoudre_systeme(matrice_pression_.valeur(), secmem,inc_pre);
  pression+=inc_pre;
  l_assembleur.modifier_solution(pression);

// Correction de la vitesse en pression
// M-1 Bt P

  vpoint += gradP;

  // On conserve Bt P pour la prochaine fois.
  eqn.operateur_gradient().calculer(pression, gradP);
  Mmoins1grad=gradP;

  eqn.solv_masse()->set_name_of_coefficient_temporel("rho_comme_v");
  eqn.solv_masse().appliquer(Mmoins1grad);
  eqn.solv_masse()->set_name_of_coefficient_temporel("no_coeff");

  // Correction en pression
  vpoint -= Mmoins1grad;
  vpoint.echange_espace_virtuel();

  return vpoint;
}
DoubleTab& Navier_Stokes_QC_impl::derivee_en_temps_inco(Navier_Stokes_std& eqn,DoubleTab& vpoint, Fluide_Incompressible& le_fluide, Matrice& matrice_pression_,Assembleur& assembleur_pression_ ,int diffusion_implicite)
{
  derivee_en_temps_impl_p1( eqn, vpoint, le_fluide,  matrice_pression_, assembleur_pression_ , diffusion_implicite);
  derivee_en_temps_impl_projection( eqn, vpoint, le_fluide,  matrice_pression_, assembleur_pression_ , diffusion_implicite);
  return vpoint;

}

//void calculer_rho_np1_vitesse(const Equation_base& eqn,const DoubleTab& vitesse,DoubleTab& resu);
void Navier_Stokes_QC_impl::assembler_impl(const Navier_Stokes_std& eqn,Matrice_Morse& matrice, const DoubleTab& present, DoubleTab& resu)
{
  Cerr<<"  Navier_Stokes_QC_impl::assembler_impl ne doit pas etre utilise. COntatctez la maintenance."<<finl;
  Process::exit();
}
void Navier_Stokes_QC_impl::assembler_avec_inertie_impl(const Navier_Stokes_std& eqn,Matrice_Morse& matrice, const DoubleTab& present, DoubleTab& resu)
{

  // assemblage special NS_QC
  // avant inertie

  // diffusion en div(mu grad u ) or on veut impliciter en rho u
  // on divise les contributions par le rho_face associe
  // GF on ajoute apres avoir contribuer pour avoir les bons flux bords
  DoubleTrav rhovitesse(present);
  eqn.operateur(0).l_op_base().contribuer_a_avec(present,matrice);

  eqn.operateur(0).ajouter(resu);
  const Fluide_Quasi_Compressible& fluide_QC=ref_cast(Fluide_Quasi_Compressible,eqn.milieu());
  const DoubleTab& tab_rho_face_np1=fluide_QC.rho_face_np1();
  const DoubleTab& tab_rho_face_n=fluide_QC.rho_face_n();
  int nb_compo = 1;
  if (present.nb_dim() == 2) nb_compo = present.dimension(1);
  const IntVect& tab1= matrice.get_tab1();
  const IntVect& tab2= matrice.get_tab2();
  DoubleVect& coeff=matrice.get_set_coeff();
  for (int i=0; i<matrice.nb_lignes(); i++)
    {
      for (int k=tab1(i)-1; k<tab1(i+1)-1; k++)
        {
          int j=tab2(k)-1;
          double rapport=tab_rho_face_np1(j/nb_compo);
          coeff(k)/=rapport;
        }
    }
  // calculer_rho_np1_vitesse(eqn,present,rhovitesse);
  rho_vitesse_impl(tab_rho_face_np1,present,rhovitesse);

  eqn.operateur(1).l_op_base().contribuer_a_avec(rhovitesse,matrice);
  eqn.operateur(1).ajouter(rhovitesse,resu);
  eqn.sources().ajouter(resu);
  eqn.sources().contribuer_a_avec(present,matrice);
  // on resout en rho u on stocke donc rho u dans present
  // calculer_rho_np1_vitesse(eqn,present,present);
  rho_vitesse_impl(tab_rho_face_np1,present,ref_cast_non_const(DoubleTab,present));
  matrice.ajouter_multvect(present,resu);

  // contribution a la matrice de l'inertie
  // en attenddant de faire mieux
  // on divisie la diagonale par rhon+1 face
  // on ajoute l'inertiede facon standard
  // on remultiplie la diagonale par rhon+1
  // ajout de l'inertie
  const double& dt=eqn.schema_temps().pas_de_temps();
  eqn.solv_masse().ajouter_masse(dt,matrice,0);

  // calculer_rho_n_vitesse(eqn,passe,rhovitesse);
  rho_vitesse_impl(tab_rho_face_n,eqn.inconnue().passe(),rhovitesse);
  eqn.solv_masse().ajouter_masse(dt,resu,rhovitesse,0);


  // blocage_cl faux si dirichlet u!=0 !!!!!! manque multiplication par rho
  for (int op=0; op< eqn.nombre_d_operateurs(); op++)
    eqn.operateur(op).l_op_base().modifier_pour_Cl(matrice,resu);
  // correction finale pour les dirichlets
  // on ne doit pas imposer un+1 mais rho_un+1
  // on multiplie dons le resu par rho_face_np1
  const Conds_lim& lescl=eqn.zone_Cl_dis().les_conditions_limites();
  int nbcondlim=lescl.size();

  for (int icl=0; icl<nbcondlim; icl++)
    {
      const Cond_lim_base& la_cl_base = lescl[icl].valeur();
      if (sub_type(Dirichlet,la_cl_base))
        {

          const Front_VF& la_front_dis = ref_cast(Front_VF,la_cl_base.frontiere_dis());
          int ndeb = la_front_dis.num_premiere_face();
          int nfin = ndeb + la_front_dis.nb_faces();
          for (int num_face=ndeb; num_face<nfin; num_face++)
            {
              if (present.nb_dim()==1)
                resu(num_face)*=tab_rho_face_np1(num_face);
              else
                for (int dir=0; dir<Objet_U::dimension; dir++)
                  resu(num_face,dir)*=tab_rho_face_np1(num_face);
            }
        }
    }

}
int Navier_Stokes_QC_impl::impr_impl(const Navier_Stokes_std& eqn ,Sortie& os) const
{
  {
    const DoubleTab& vitesse=eqn.vitesse().valeurs();

    const Fluide_Quasi_Compressible& fluide_QC=ref_cast(Fluide_Quasi_Compressible,eqn.fluide());
    DoubleTab prov(vitesse);
    const DoubleTab& rho=fluide_QC.rho_face_np1();
    rho_vitesse_impl(rho,vitesse,prov);
    DoubleTab res;
    res.copy(eqn.div().valeurs(), Array_base::NOCOPY_NOINIT); // init structure uniquement
    if (tab_W.get_md_vector().non_nul())
      {
        // tab_W initialise si on passe dans
        operator_egal(res, tab_W ); //, VECT_REAL_ITEMS);
        res*=-1;
      }
    else
      {
        // remarque B.M.: certaines implementations de cette methode ne font pas echange espace virtuel:
        fluide_QC.secmembre_divU_Z(res);
        res*=-1;
      }

    res.echange_espace_virtuel();
    eqn.operateur_divergence().ajouter(prov, res);
    double bilan_massique = mp_max_abs_vect(res);
    os << "-------------------------------------------------------------------"<< finl;
    os << "Cell balance mass flow rate control for the problem " << eqn.probleme().le_nom() << " : " << finl;
    os << "Absolute value : " << bilan_massique << " kg/s" << finl; ;
    eqn.operateur_divergence().volumique(res); // On divise res par vol(i)
    // On divise par un rho moyen
    double rho_moyen = mp_moyenne_vect(rho);
    double dt = eqn.probleme().schema_temps().pas_de_temps();
    double bilan_massique_relatif = mp_max_abs_vect(res) * dt / rho_moyen;
    os << "Relative value : " << bilan_massique_relatif << finl; // =max|bilan_massique/(rho_moyen*Vol(i)/dt)|
    /*
    // Nouveau 1.6.1, arret si bilans de masse mauvais
    // Commente car plusieurs cas tests plantent en implicite, peut etre un probleme
    // du calcul du bilan massique en implicite
    if (bilan_massique_relatif>0.01)
    {
    Cerr << "The mass balance is too bad (relative value>1%)." << finl;
    Cerr << "Please check and lower the convergence value of the pressure solver." << finl;
    Process::exit();
    }
    */
  }
  return 1;
}



