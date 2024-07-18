/****************************************************************************
* Copyright (c) 2019, CEA
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
// File:        Op_Diff_K_Omega_VEF_Face.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VEF/Operateurs
//
//////////////////////////////////////////////////////////////////////////////

#include <Op_Diff_K_Omega_VEF_Face.h>
#include <Modele_turbulence_hyd_K_Omega.h>
#include <Champ_P1NC.h>
#include <Periodique.h>
#include <Neumann_paroi.h>
#include <Paroi_hyd_base_VEF.h>
#include <Champ_Don.h>
#include <Champ_Uniforme.h>
#include <Fluide_Incompressible.h>

Implemente_instanciable(Op_Diff_K_Omega_VEF_Face,
                        "Op_Diff_K_Omega_VEF_P1NC",
                        Op_Diff_K_Omega_VEF_base);

Sortie& Op_Diff_K_Omega_VEF_Face::printOn(Sortie& s ) const
{
  return s << que_suis_je() ;
}
Entree& Op_Diff_K_Omega_VEF_Face::readOn(Entree& s )
{
  return s ;
}
/*! @brief complete l'evaluateur
 *
 */
void Op_Diff_K_Omega_VEF_Face::associer(const Domaine_dis& domaine_dis,
                                        const Domaine_Cl_dis& domaine_cl_dis,
                                        const Champ_Inc& ch_diffuse)
{
  inconnue_  = ref_cast(Champ_P1NC, ch_diffuse.valeur());
  le_dom_vef = ref_cast(Domaine_VEF, domaine_dis.valeur());
  la_zcl_vef = ref_cast(Domaine_Cl_VEF, domaine_cl_dis.valeur());
}

void Op_Diff_K_Omega_VEF_Face::associer_diffusivite_turbulente()
{
  assert(mon_equation.non_nul());

  const Transport_K_Omega& eqn_transport = ref_cast(Transport_K_Omega,
                                                    mon_equation.valeur());
  const Modele_turbulence_hyd_K_Omega& mod_turb = ref_cast(Modele_turbulence_hyd_K_Omega,
                                                           eqn_transport.modele_turbulence());
  turbulence_model = mod_turb;
  const Champ_Fonc& diff_turb = mod_turb.viscosite_turbulente();
  Op_Diff_K_Omega_VEF_base::associer_diffusivite_turbulente(diff_turb);
}

double Op_Diff_K_Omega_VEF_Face::blender(double const val1, double const val2,
                                         int const face) const
{
  const DoubleTab& F1 = turbulence_model->get_blenderF1();
  return F1(face)*val1 + (1 - F1(face))*val2;
}


void  Op_Diff_K_Omega_VEF_Face::calc_visc(ArrOfDouble& diffu_tot, const Domaine_VEF& le_dom,
                                          int num_face, int num2, int dimension_inut,
                                          int num_elem, double diffu_turb,
                                          const DoubleTab& diffu, int is_mu_unif,
                                          const ArrOfDouble& inv_Prdt) const
{
  double valA = viscA(num_face, num2, num_elem, diffu_turb);
  double nu = 0;
  if (is_mu_unif)
    nu = diffu(0, 0);
  else
    nu = diffu(num_elem);
  diffu_tot[0] = nu + valA*inv_Prdt[0];
  diffu_tot[1] = nu + valA*inv_Prdt[1];
}


DoubleTab& Op_Diff_K_Omega_VEF_Face::ajouter(const DoubleTab& inconnue, DoubleTab& resu) const
{
  const Domaine_Cl_VEF& domaine_Cl_VEF = la_zcl_vef.valeur();
  const Domaine_VEF& domaine_VEF = le_dom_vef.valeur();

  const IntTab& elem_faces = domaine_VEF.elem_faces();
  const IntTab& face_voisins = domaine_VEF.face_voisins();

  ArrOfDouble diffu_tot(2);
  ArrOfDouble inv_Prdt(2);
  inv_Prdt[0] = 1./Prdt_K;
  inv_Prdt[1] = 1./Prdt_Omega;

  double invPrdtOmega = 1./Prdt_Omega;
  const bool is_SST = turbulence_model->is_SST();

  int is_mu_unif = sub_type(Champ_Uniforme, diffusivite_.valeur());
  const DoubleTab& mu = diffusivite_->valeurs();

  const int nb_faces_elem = domaine_VEF.domaine().nb_faces_elem();
  const int nb_front = domaine_VEF.nb_front_Cl();
  const DoubleTab& mu_turb = diffusivite_turbulente_->valeurs();
  // double inverse_Prdt_K = 1.0/Prdt_K;
  //double inverse_Prdt_Omega = 1.0/Prdt_Omega;

  // On traite les faces bord :
  for (int n_bord = 0; n_bord < nb_front; n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF, la_cl.frontiere_dis());;
      int num1 = 0;
      int num2 = le_bord.nb_faces_tot();

      if (sub_type(Periodique, la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique, la_cl.valeur());

          for (int ind_face = num1; ind_face < num2; ind_face++)
            {
              int fac_asso = la_cl_perio.face_associee(ind_face);
              fac_asso = le_bord.num_face(fac_asso);
              int num_face = le_bord.num_face(ind_face);

              inv_Prdt[1] = is_SST
                            ? blender(SIGMA_OMEGA1, SIGMA_OMEGA2, num_face)
                            : invPrdtOmega;

              for (int k = 0; k < 2; k++)
                {
                  int elem = face_voisins(num_face, k);
                  double d_mu = mu_turb(elem);

                  for (int i = 0; i < nb_faces_elem; i++)
                    {
                      int j = elem_faces(elem, i);
                      if ((j > num_face) && (j != fac_asso))
                        {
                          calc_visc(diffu_tot, domaine_VEF, num_face, j,
                                    dimension, elem, d_mu, mu, is_mu_unif, inv_Prdt);

                          resu(num_face, 0) += diffu_tot[0]*inconnue(j, 0);
                          resu(num_face, 0) -= diffu_tot[0]*inconnue(num_face, 0);
                          resu(num_face, 1) += diffu_tot[1]*inconnue(j, 1);
                          resu(num_face, 1) -= diffu_tot[1]*inconnue(num_face, 1);

                          resu(j, 0) += 0.5*diffu_tot[0]*inconnue(num_face, 0);
                          resu(j, 0) -= 0.5*diffu_tot[0]*inconnue(j, 0);
                          resu(j, 1) += 0.5*diffu_tot[1]*inconnue(num_face, 1);
                          resu(j, 1) -= 0.5*diffu_tot[1]*inconnue(j, 1);
                        }
                    }
                }
            }
        }
      else   //  conditions aux limites autres que periodicite, paroi_fixe, symetrie
        {
          for (int ind_face = num1; ind_face < num2; ind_face++)
            {
              int num_face = le_bord.num_face(ind_face);
              int elem = face_voisins(num_face,0);
              double d_mu = mu_turb(elem);
              inv_Prdt[1] = is_SST
                            ? blender(SIGMA_OMEGA1, SIGMA_OMEGA2, num_face)
                            : invPrdtOmega;

              for (int i = 0; i < nb_faces_elem; i++)
                {
                  int j = elem_faces(elem, i);
                  if (j > num_face)
                    {
                      // valA = viscA(domaine_VEF,num_face,j,dimension,elem,d_mu);
                      calc_visc(diffu_tot, domaine_VEF, num_face, j,
                                dimension, elem, d_mu, mu, is_mu_unif, inv_Prdt);

                      resu(num_face, 0) += diffu_tot[0]*inconnue(j, 0);
                      resu(num_face, 0) -= diffu_tot[0]*inconnue(num_face, 0);
                      resu(num_face, 1) += diffu_tot[1]*inconnue(j, 1);
                      resu(num_face, 1) -= diffu_tot[1]*inconnue(num_face, 1);

                      resu(j, 0) += diffu_tot[0]*inconnue(num_face, 0);
                      resu(j, 0) -= diffu_tot[0]*inconnue(j, 0);
                      resu(j, 1) += diffu_tot[1]*inconnue(num_face, 1);
                      resu(j, 1) -= diffu_tot[1]*inconnue(j, 1);
                    }
                }
            }
        }
    }

  // On traite les faces internes
  int n0 = domaine_VEF.premiere_face_int();
  int n1 = domaine_VEF.nb_faces();

  for (int num_face = n0; num_face < n1; num_face++)
    {
      inv_Prdt[1] = is_SST
                    ? blender(SIGMA_OMEGA1, SIGMA_OMEGA2, num_face)
                    : invPrdtOmega;

      for (int k = 0; k < 2; k++)
        {
          int elem = face_voisins(num_face, k);
          double d_mu = mu_turb(elem);

          for (int i = 0; i < nb_faces_elem; i++)
            {
              int j = elem_faces(elem, i);
              if (j > num_face)
                {
                  //   double valA = viscA(num_face,j,elem,d_mu);
                  calc_visc(diffu_tot, domaine_VEF, num_face, j,
                            dimension, elem, d_mu, mu, is_mu_unif, inv_Prdt);

                  resu(num_face, 0) += diffu_tot[0]*inconnue(j, 0);
                  resu(num_face, 0) -= diffu_tot[0]*inconnue(num_face, 0);
                  resu(num_face, 1) += diffu_tot[1]*inconnue(j, 1);
                  resu(num_face, 1) -= diffu_tot[1]*inconnue(num_face, 1);

                  resu(j, 0) += diffu_tot[0]*inconnue(num_face, 0);
                  resu(j, 0) -= diffu_tot[0]*inconnue(j, 0);
                  resu(j, 1) += diffu_tot[1]*inconnue(num_face, 1);
                  resu(j, 1) -= diffu_tot[1]*inconnue(j, 1);
                }
            }
        }
    }

  // Prise en compte des conditions aux limites de Neumnan a la paroi
  for (int n_bord = 0; n_bord < domaine_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);

      if (sub_type(Neumann_paroi, la_cl.valeur()))
        {
          const Neumann_paroi& la_cl_paroi = ref_cast(Neumann_paroi, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF, la_cl.frontiere_dis());
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();

          for (int face = ndeb; face < nfin; face++)
            {
              resu(face, 0) += la_cl_paroi.flux_impose(face-ndeb, 0)*domaine_VEF.surface(face);
              resu(face, 1) += la_cl_paroi.flux_impose(face-ndeb, 1)*domaine_VEF.surface(face);
            }
        }
    }

  modifier_flux(*this);
  return resu;
}


DoubleTab& Op_Diff_K_Omega_VEF_Face::calculer(const DoubleTab& inconnue, DoubleTab& resu) const
{
  resu = 0;
  return ajouter(inconnue, resu);
}


/////////////////////////////////////////
// Methode pour l'implicite
/////////////////////////////////////////

void Op_Diff_K_Omega_VEF_Face::ajouter_contribution(const DoubleTab& transporte,
                                                    Matrice_Morse& matrice) const

{
  modifier_matrice_pour_periodique_avant_contribuer(matrice,equation());
  const Domaine_Cl_VEF& domaine_Cl_VEF = la_zcl_vef.valeur();
  const Domaine_VEF& domaine_VEF = le_dom_vef.valeur();

  const IntTab& elem_faces = domaine_VEF.elem_faces();
  const IntTab& face_voisins = domaine_VEF.face_voisins();

  int i,j,nc,l,num_face;
  int nb_faces = domaine_VEF.nb_faces();
  int elem;
  int nb_faces_elem = domaine_VEF.domaine().nb_faces_elem();
  double  d_mu;
  ArrOfDouble inv_Prdt(2),diffu_tot(2);
  inv_Prdt[0]=1./Prdt_K;
  inv_Prdt[1]=1./Prdt_Omega;
  const DoubleTab& mu_turb=diffusivite_turbulente_->valeurs();

  const int is_mu_unif=sub_type(Champ_Uniforme,diffusivite_.valeur());
  const DoubleTab& mu=diffusivite_->valeurs();
  int nb_comp = 1;

  nb_comp=transporte.line_size();
  assert(nb_comp==2);

  // On traite les faces bord
  for (int n_bord=0; n_bord<domaine_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
      int num1 = 0;
      int num2 = le_bord.nb_faces_tot();
      int nb_faces_bord_reel = le_bord.nb_faces();
      if (sub_type(Periodique,la_cl.valeur()))
        {
          const Periodique& la_cl_perio = ref_cast(Periodique,la_cl.valeur());
          // on ne parcourt que la moitie des faces volontairement...
          // GF il ne faut pas s'occuper des faces virtuelles
          num2=nb_faces_bord_reel/2;
          for (int ind_face=num1; ind_face<num2; ind_face++)
            {
              int fac_asso = la_cl_perio.face_associee(ind_face);
              fac_asso = le_bord.num_face(fac_asso);
              int num_face_b = le_bord.num_face(ind_face);
              for (int ll=0; ll<2; ll++)
                {
                  int elem2 = face_voisins(num_face_b,ll);
                  double d_mu2 = mu_turb(elem2);
                  for (int ii=0; ii<nb_faces_elem; ii++)
                    {
                      if ( (j= elem_faces(elem2,ii)) > num_face_b )
                        {
                          //double valAp = viscA(num_face_b,j,elem2,d_mu2);
                          calc_visc(diffu_tot,domaine_VEF,num_face_b,j,dimension,elem2,d_mu2,mu,is_mu_unif,inv_Prdt);
                          int fac_loc=0;
                          int ok=1;
                          while ((fac_loc<nb_faces_elem) && (elem_faces(elem2,fac_loc)!=num_face_b)) fac_loc++;
                          if (fac_loc==nb_faces_elem)
                            ok=0;
                          int contrib=1;
                          if(j>=nb_faces) // C'est une face virtuelle
                            {
                              int el1 = face_voisins(j,0);
                              int el2 = face_voisins(j,1);
                              if((el1==-1)||(el2==-1))
                                contrib=0;
                            }
                          if (contrib)
                            {
                              for (int nc2=0; nc2<nb_comp; nc2++)
                                {
                                  int n0=num_face_b*nb_comp+nc2;
                                  int n0perio=fac_asso*nb_comp+nc2;
                                  int j0=j*nb_comp+nc2;
                                  double valA=diffu_tot[nc2];
                                  matrice(n0,n0)+=valA;
                                  matrice(n0,j0)-=valA;
                                  if(j<nb_faces) // On traite les faces reelles
                                    {
                                      //if (l==0)
                                      if (ok==1)
                                        matrice(j0,n0)-=valA;
                                      else
                                        matrice(j0,n0perio)-=valA;
                                      matrice(j0,j0)+=valA;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
      else   //  conditions aux limites autres que periodicite
        {
          for (int ind_face=num1; ind_face<num2; ind_face++)
            {
              int num_face_b = le_bord.num_face(ind_face);
              elem = face_voisins(num_face_b,0);
              d_mu = mu_turb(elem);
              for (i=0; i<nb_faces_elem; i++)
                if (( (j= elem_faces(elem,i)) > num_face_b ) || (ind_face>=nb_faces_bord_reel))
                  {
                    //double valAp = viscA(num_face_b,j,elem,d_mu);
                    calc_visc(diffu_tot,domaine_VEF,num_face_b,j,dimension,elem,d_mu,mu,is_mu_unif,inv_Prdt);

                    for (int nc2=0; nc2<nb_comp; nc2++)
                      {
                        int n0=num_face_b*nb_comp+nc2;
                        int j0=j*nb_comp+nc2;
                        double valA=diffu_tot[nc2];

                        if (ind_face<nb_faces_bord_reel)
                          {
                            // double valA=diffu_tot[nc];
                            matrice(n0,n0)+=valA;
                            matrice(n0,j0)-=valA;
                          }
                        if(j<nb_faces)
                          {
                            matrice(j0,n0)-=valA;
                            matrice(j0,j0)+=valA;
                          }
                      }
                  }
            }
        }
    }

  // On traite les faces internes
  int n0 = domaine_VEF.premiere_face_int();
  for (num_face=n0; num_face<nb_faces; num_face++)
    {
      for ( l=0; l<2; l++)
        {
          int elem2 = face_voisins(num_face,l);
          double d_mu2 = mu_turb(elem2);
          for (i=0; i<nb_faces_elem; i++)
            {
              int jj = elem_faces(elem2,i);
              if ( jj > num_face )
                {
                  int contrib=1;
                  if(jj>=nb_faces) // C'est une face virtuelle
                    {
                      int el1 = face_voisins(jj,0);
                      int el2 = face_voisins(jj,1);
                      if((el1==-1)||(el2==-1))
                        contrib=0;
                    }
                  if (contrib)
                    {
                      //double valAp = viscA(num_face,jj,elem2,d_mu2);
                      calc_visc(diffu_tot,domaine_VEF,num_face,jj,
                                dimension,elem2,d_mu2,mu,
                                is_mu_unif,inv_Prdt);

                      for (nc=0; nc<nb_comp; nc++)
                        {
                          double valA=diffu_tot[nc];
                          int num0=num_face*nb_comp+nc;
                          int j0=jj*nb_comp+nc;
                          matrice(num0,num0)+=valA;
                          matrice(num0,j0)-=valA;
                          if (jj < nb_faces) // On traite les faces reelles
                            {
                              matrice(j0, num0) -= valA;
                              matrice(j0, j0) += valA;
                            }
                        }
                    }
                }
            }
        }
    }

  modifier_matrice_pour_periodique_apres_contribuer(matrice,equation());
}

void Op_Diff_K_Omega_VEF_Face::contribue_au_second_membre(DoubleTab& resu) const
{
  const Domaine_Cl_VEF& domaine_Cl_VEF = la_zcl_vef.valeur();
  const Domaine_VEF& domaine_VEF = le_dom_vef.valeur();

  // Prise en compte des conditions aux limites de Neumnan a la paroi
  for (int n_bord = 0; n_bord < domaine_VEF.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = domaine_Cl_VEF.les_conditions_limites(n_bord);
      if (sub_type(Neumann_paroi, la_cl.valeur()))
        {
          const Neumann_paroi& la_cl_paroi = ref_cast(Neumann_paroi, la_cl.valeur());
          const Front_VF& le_bord = ref_cast(Front_VF, la_cl.frontiere_dis());
          int ndeb = le_bord.num_premiere_face();
          int nfin = ndeb + le_bord.nb_faces();
          for (int face = ndeb; face < nfin; face++)
            {
              resu(face, 0) += la_cl_paroi.flux_impose(face-ndeb, 0)*domaine_VEF.surface(face);
              resu(face, 1) += la_cl_paroi.flux_impose(face-ndeb, 1)*domaine_VEF.surface(face);
            }
        }
    }
}

/*! @brief On modifie le second membre et la matrice dans le cas des conditions de dirichlet.
 *
 */
void Op_Diff_K_Omega_VEF_Face::modifier_pour_Cl(Matrice_Morse& matrice, DoubleTab& secmem) const
{
  // on recupere le tableau
  const Transport_K_Omega& eqn_k_omega = ref_cast(Transport_K_Omega, equation());
  const DoubleTab& val = equation().inconnue().valeurs();
  const Turbulence_paroi& mod = eqn_k_omega.modele_turbulence().loi_paroi();
  const Paroi_hyd_base_VEF& paroi = ref_cast(Paroi_hyd_base_VEF, mod.valeur());
  const ArrOfInt& face_komega_imposee = paroi.face_keps_imposee();
  int size = secmem.dimension(0);
  const IntVect& tab1 = matrice.get_tab1();
  DoubleVect& coeff = matrice.get_set_coeff();
  int nb_comp = 2;

  Op_VEF_Face::modifier_pour_Cl(le_dom_vef.valeur(), la_zcl_vef.valeur(), matrice, secmem);

  if (face_komega_imposee.size_array()>0)
    // en plus des dirichlets ????
    // on change la matrice et le resu sur toutes les lignes ou k_omega_ est imposee....
    for (int face = 0; face < size; face++)
      {
        if (face_komega_imposee[face] != -2)
          {
            for (int comp = 0; comp < nb_comp; comp++)
              {
                // on doit remettre la ligne a l'identite et le secmem a l'inconnue
                int idiag = tab1[face*nb_comp + comp]-1;
                coeff[idiag] = 1;
                // pour les voisins
                int nbvois = tab1[face*nb_comp + 1 + comp] - tab1[face*nb_comp + comp];
                for (int k = 1; k < nbvois; k++)
                  coeff[idiag + k] = 0;
                secmem(face, comp) = val(face, comp);
              }
          }
      }
}
