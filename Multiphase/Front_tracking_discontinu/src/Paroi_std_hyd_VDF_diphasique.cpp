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
/////////////////////////////////////////////////////////////////////////////
//
// File      : Paroi_std_hyd_VDF_diphasique.cpp
// Directory : $TrioCFD_ROOT/Front_tracking_discontinu/src
//
/////////////////////////////////////////////////////////////////////////////

#include <Paroi_std_hyd_VDF_diphasique.h>
#include <Zone_Cl_VDF.h>
#include <Dirichlet_paroi_fixe.h>
#include <Dirichlet_paroi_defilante.h>
#include <Fluide_Quasi_Compressible.h>
#include <Equation_base.h>
#include <Fluide_Diphasique.h>
#include <Transport_Interfaces_FT_Disc.h>
#include <Probleme_base.h>
#include <Champ_Uniforme.h>

Implemente_instanciable( Paroi_std_hyd_VDF_diphasique, "loi_standard_hydr_diphasique_VDF", Paroi_std_hyd_VDF ) ;

Sortie& Paroi_std_hyd_VDF_diphasique::printOn( Sortie& os ) const
{
  return os << que_suis_je() << " " << le_nom();
}

Entree& Paroi_std_hyd_VDF_diphasique::readOn( Entree& is )
{
  Paroi_std_hyd_VDF::readOn(is);
  return is;
}

// tab1 = tab_nu_t and tab2 = tab_k
int Paroi_std_hyd_VDF_diphasique::calculer_hyd(DoubleTab& tab1, DoubleTab& tab2)
{
  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const IntVect& orientation = zone_VDF.orientation();
  const IntTab& face_voisins = zone_VDF.face_voisins();
  const Equation_base& eqn_hydr = mon_modele_turb_hyd->equation();
  const DoubleVect& vit = eqn_hydr.inconnue().valeurs();

  // Physical properties of both phases
  const Fluide_Diphasique& le_fluide = ref_cast(Fluide_Diphasique, eqn_hydr.milieu());
  const Fluide_Incompressible& phase_1 = le_fluide.fluide_phase(1);
  const Fluide_Incompressible& phase_0 = le_fluide.fluide_phase(0);
  const Champ_Don& ch_visco_cin_ph1 = phase_1.viscosite_cinematique();
  const Champ_Don& ch_visco_cin_ph0 = phase_0.viscosite_cinematique();
  const DoubleTab& tab_visco_ph1 = phase_1.viscosite_cinematique().valeur().valeurs();
  const DoubleTab& tab_visco_ph0 = phase_0.viscosite_cinematique().valeur().valeurs();
  const double delta_nu = tab_visco_ph1(0,0) - tab_visco_ph0(0,0);

  // One way to get the Transport equation to pass the indicator DoubleTab
  const Zone_Cl_dis_base& zone_Cl_dis_base = eqn_hydr.zone_Cl_dis().valeur();
  const Equation_base& eqn_trans = zone_Cl_dis_base.equation().probleme().equation("Transport_Interfaces_FT_Disc");
  const Transport_Interfaces_FT_Disc& eqn_interf = ref_cast(Transport_Interfaces_FT_Disc, eqn_trans);
  const DoubleTab& indic = eqn_interf.inconnue().valeurs();

  double visco_ph0=-1;
  int l_unif;

  if (sub_type(Champ_Uniforme,ch_visco_cin_ph1.valeur()) && sub_type(Champ_Uniforme,ch_visco_cin_ph0.valeur()) )
    {
      visco_ph0 = max(tab_visco_ph0(0,0),DMINFLOAT);
      l_unif = 1;
    }
  else
    l_unif = 0;
  if ((!l_unif) && ((tab_visco_ph1.local_min_vect()<DMINFLOAT) || (tab_visco_ph0.local_min_vect()<DMINFLOAT) ))
    {
      Cerr << "Negative viscosity !!!" << finl;
      Process::exit(-1);
    }

  int ndeb,nfin;
  int elem,ori;
  double norm_v;
  double dist;
  double u_plus_d_plus,d_visco;
  //double val,val1,val2;

  //****************************************************************
  // Modifs du 19/10 pour etre coherent avec la diffusion turbulente
  double signe;
  //****************************************************************

  // Boucle sur les bords

  for (int n_bord=0; n_bord<zone_VDF.nb_front_Cl(); n_bord++)
    {

      // pour chaque condition limite on regarde son type
      // On applique les lois de paroi uniquement
      // aux voisinages des parois

      const Cond_lim& la_cl = la_zone_Cl_VDF->les_conditions_limites(n_bord);
      const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
      ndeb = le_bord.num_premiere_face();
      nfin = ndeb + le_bord.nb_faces();

      if (sub_type(Dirichlet_paroi_fixe,la_cl.valeur()) || sub_type(Dirichlet_paroi_defilante,la_cl.valeur() ))
        {
          int isdiri=0;
          DoubleTab vitesse_imposee_face_bord(le_bord.nb_faces(),dimension);
          if (sub_type(Dirichlet_paroi_defilante,la_cl.valeur()) )
            {
              isdiri=1;
              const Dirichlet_paroi_defilante& cl_diri = ref_cast(Dirichlet_paroi_defilante,la_cl.valeur());
              for (int face=ndeb; face<nfin; face++)
                for (int k=0; k<dimension; k++)
                  vitesse_imposee_face_bord(face-ndeb,k) = cl_diri.val_imp(face-ndeb,k);
            }
          ArrOfDouble vit_paroi(dimension);
          ArrOfDouble val(dimension-1);


          for (int num_face=ndeb; num_face<nfin; num_face++)
            {
              if (isdiri)
                {
                  int rang = num_face-ndeb;
                  for (int k=0; k<dimension; k++)
                    vit_paroi[k]=vitesse_imposee_face_bord(rang,k);
                }
              ori = orientation(num_face);
              if ( (elem =face_voisins(num_face,0)) != -1)
                {
                  //norm_v=norm_2D_vit(vit,elem,ori,zone_VDF,val);
                  norm_v=norm_vit(vit,elem,ori,zone_VDF,vit_paroi,val);
                  signe = -1.;
                }
              else
                {
                  elem = face_voisins(num_face,1);
                  //norm_v=norm_2D_vit(vit,elem,ori,zone_VDF,val);
                  norm_v=norm_vit(vit,elem,ori,zone_VDF,vit_paroi,val);
                  signe = 1.;
                }
              if (axi)
                dist=zone_VDF.dist_norm_bord_axi(num_face);
              else
                dist=zone_VDF.dist_norm_bord(num_face);
              if (l_unif)
                d_visco = visco_ph0 + indic(elem) * delta_nu;
              else
                d_visco = (tab_visco_ph0.nb_dim()==1 ? (tab_visco_ph0(elem) + indic(elem) * delta_nu) :
                           (tab_visco_ph0(elem,0) + indic(elem) * delta_nu));

              u_plus_d_plus = norm_v*dist/d_visco;

              // Calcul de u* et des grandeurs turbulentes:


              // tab1=tab_nu_t
              // tab2=tab_k
              //calculer_local(u_plus_d_plus,d_visco,tab_nu_t,tab_k,norm_v,dist,elem,num_face);
              calculer_local(u_plus_d_plus,d_visco,tab1,tab2,norm_v,dist,elem,num_face);

              // Calcul de la contrainte tangentielle

              double vit_frot = tab_u_star(num_face)*tab_u_star(num_face);

              if (ori == 0)
                {
                  Cisaillement_paroi_(num_face,1) = vit_frot*val[0]*signe;
                  if (dimension==3)
                    Cisaillement_paroi_(num_face,2) = vit_frot*val[1]*signe;
                  Cisaillement_paroi_(num_face,0) = 0.;
                }
              else if (ori == 1)
                {
                  Cisaillement_paroi_(num_face,0) = vit_frot*val[0]*signe;
                  if (dimension==3)
                    Cisaillement_paroi_(num_face,2) = vit_frot*val[1]*signe;
                  Cisaillement_paroi_(num_face,1) = 0.;
                }
              else
                {
                  Cisaillement_paroi_(num_face,0) = vit_frot*val[0]*signe;
                  Cisaillement_paroi_(num_face,1) = vit_frot*val[1]*signe;
                  Cisaillement_paroi_(num_face,2) = 0.;
                }

              // Calcul de u+ d+
              calculer_uplus_dplus(uplus_, tab_d_plus_, tab_u_star_, num_face, dist, d_visco, norm_v) ;
            }


        }
    }
  Cisaillement_paroi_.echange_espace_virtuel();
  tab1.echange_espace_virtuel();
  tab2.echange_espace_virtuel();

  return 1;

}

int Paroi_std_hyd_VDF_diphasique::calculer_hyd(DoubleTab& tab_k_eps)
{
  Cerr << " Paroi_std_hyd_VDF_diphasique::calculer_hyd(DoubleTab& tab_k_eps) " << finl;
  Cerr <<  "on ne doit pas entrer dans cette methode" << finl;
  Cerr << " car elle est definie uniquement pour la LES " << finl ;
  Process::exit(-1);
  return 1 ;
} // fin de calcul_hyd (K-eps)
