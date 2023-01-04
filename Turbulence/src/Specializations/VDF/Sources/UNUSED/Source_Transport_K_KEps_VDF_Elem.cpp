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
// File:        Source_Transport_K_KEps_VDF_Elem.cpp
// Directory:   $TURBULENCE_ROOT/src/Specializations/VDF/Sources/UNUSED
//
//////////////////////////////////////////////////////////////////////////////

#include <Source_Transport_K_KEps_VDF_Elem.h>
#include <Convection_Diffusion_Temperature.h>
#include <Modele_turbulence_scal_base.h>
#include <Dirichlet_paroi_defilante.h>
#include <Dirichlet_paroi_fixe.h>
#include <Paroi_std_hyd_VDF.h>
#include <Transport_K_KEps.h>
#include <Champ_Uniforme.h>
#include <Probleme_base.h>
#include <Fluide_base.h>
#include <Zone_Cl_VDF.h>
#include <Champ_Face_VDF.h>
#include <TRUSTTrav.h>

Implemente_instanciable_sans_constructeur(Source_Transport_K_KEps_VDF_Elem,"Source_Transport_K_KEps_VDF_P0_VDF",Source_Transport_VDF_Elem_base);

Sortie& Source_Transport_K_KEps_VDF_Elem::printOn(Sortie& s) const { return s << que_suis_je() ; }
Entree& Source_Transport_K_KEps_VDF_Elem::readOn(Entree& is) { return Source_Transport_VDF_Elem_base::readOn(is); }

void Source_Transport_K_KEps_VDF_Elem::associer_pb(const Probleme_base& pb)
{
  Source_Transport_VDF_Elem_base::associer_pb(pb);
  mon_eq_transport_K_Eps = ref_cast(Transport_K_KEps,equation());
}

void Source_Transport_K_KEps_VDF_Elem::ajouter_blocs(matrices_t matrices, DoubleTab& secmem, const tabs_t& semi_impl) const
{
  const Zone_VDF& zone_VDF = la_zone_VDF.valeur();
  const Zone_VDF& zone_VDF_NS = ref_cast(Zone_VDF,eq_hydraulique->zone_dis().valeur());
  const Zone_Cl_VDF& zone_Cl_VDF_NS = ref_cast(Zone_Cl_VDF,eq_hydraulique->zone_Cl_dis().valeur());
  const DoubleTab& K_eps = mon_eq_transport_K_Eps->inconnue().valeurs();
  const DoubleTab& visco_turb = mon_eq_transport_K_Eps->modele_turbulence().viscosite_turbulente().valeurs();
  const DoubleTab& vit = eq_hydraulique->inconnue().valeurs();
  const DoubleVect& volumes = zone_VDF.volumes();
  const DoubleVect& porosite_vol = la_zone_Cl_VDF->equation().milieu().porosite_elem();
  const IntTab& face_voisins = zone_VDF.face_voisins();
  const IntTab& elem_faces = zone_VDF.elem_faces();
  const int nbcouches = mon_eq_transport_K_Eps->get_nbcouches();
  int ndeb,nfin,elem,face_courante,elem_courant;
  double dist, d_visco,y_etoile, critere_switch=0.;
  Loi_2couches_base& loi2couches =ref_cast_non_const(Transport_K_KEps,mon_eq_transport_K_Eps.valeur()).loi2couches();
  int valswitch, icouche;
  const int impr2 =  mon_eq_transport_K_Eps->get_impr();

  const int typeswitch = mon_eq_transport_K_Eps->get_switch();
  if ( typeswitch == 0 ) valswitch = mon_eq_transport_K_Eps->get_yswitch();
  else valswitch = mon_eq_transport_K_Eps->get_nutswitch();

  const int nb_elem = zone_VDF.nb_elem(), nb_elem_tot = zone_VDF.nb_elem_tot();
  DoubleTrav P(nb_elem_tot), Eps(nb_elem_tot), tab_couches(nb_elem_tot);

  //Extraction de la viscosite moleculaire
  const Fluide_base& le_fluide = ref_cast(Fluide_base,eq_hydraulique->milieu());
  const Champ_Don& ch_visco_cin = le_fluide.viscosite_cinematique();
  const DoubleTab& tab_visco = ch_visco_cin->valeurs();
  double visco=-1;
  int l_unif;
  if (sub_type(Champ_Uniforme,ch_visco_cin.valeur()))
    {
      visco = std::max(tab_visco(0,0),DMINFLOAT);
      l_unif = 1;
    }
  else
    l_unif = 0;

  //essayer de calculer u+d+ ici pour chaque elt du bord?!
  for (elem=0; elem<nb_elem; elem++)
    {
      Eps[elem] = K_eps(elem,1);
      tab_couches[elem]=0;
    }

  // Boucle sur les bords
  for (int n_bord=0; n_bord<zone_VDF_NS.nb_front_Cl(); n_bord++)
    {
      const Cond_lim& la_cl = zone_Cl_VDF_NS.les_conditions_limites(n_bord);
      if (sub_type(Dirichlet_paroi_fixe,la_cl.valeur()) || sub_type(Dirichlet_paroi_defilante,la_cl.valeur()))
        {
          const Front_VF& le_bord = ref_cast(Front_VF,la_cl.frontiere_dis());
          ndeb = le_bord.num_premiere_face();
          nfin = ndeb + le_bord.nb_faces();
          if (dimension == 2 )
            for (int num_face=ndeb; num_face<nfin; num_face++)
              {
                if ( (elem =face_voisins(num_face,0)) != -1) ;
                else
                  elem = face_voisins(num_face,1);
                if (l_unif)
                  d_visco = visco;
                else
                  d_visco = tab_visco[elem];
                dist=zone_VDF_NS.distance_normale(num_face);
                //Cerr << "dist = " << dist << finl;

                face_courante = num_face;
                elem_courant = elem;


                for(icouche=0; (icouche<nbcouches) && (elem_courant != -1); icouche++)
                  {
                    double kSqRt=1.e-3;
                    if (K_eps(elem_courant,0)>0)
                      kSqRt=sqrt(K_eps(elem_courant,0));
                    y_etoile = kSqRt*dist/d_visco;
                    if (typeswitch == 0)        //selon le critere de switch choisit:y* ou nu_t
                      critere_switch = y_etoile;
                    else
                      critere_switch = visco_turb[elem_courant]/d_visco;
                    if (critere_switch < valswitch)        // si on est en proche paroi
                      {
                        tab_couches[elem_courant]=1;
                        double Leps,Lmu,vvSqRt;
                        if (K_eps(elem_courant,0)>0)
                          loi2couches.LepsLmu(K_eps(elem_courant,0),d_visco,dist,y_etoile,Leps,Lmu,vvSqRt);
                        else
                          loi2couches.LepsLmu(1.e-3,d_visco,dist,y_etoile,Leps,Lmu,vvSqRt);

                        if ( Leps!=0.)
                          Eps[elem_courant] = K_eps(elem_courant,0)*kSqRt/Leps;
                        else Eps[elem_courant] = 1.e-3;
                      }
                    else
                      {
                        break;
                      }
                    if ( elem_faces(elem_courant,0) == face_courante)
                      face_courante = elem_faces(elem_courant,2);
                    else if ( elem_faces(elem_courant,1) == face_courante)
                      face_courante = elem_faces(elem_courant,3);
                    else if ( elem_faces(elem_courant,2) == face_courante)
                      face_courante = elem_faces(elem_courant,0);
                    else if ( elem_faces(elem_courant,3) == face_courante)
                      face_courante = elem_faces(elem_courant,1);
                    if ( face_voisins(face_courante,0) != elem_courant)
                      elem_courant = face_voisins(face_courante,0);
                    else
                      elem_courant = face_voisins(face_courante,1);
                    dist+=zone_VDF_NS.distance_normale(face_courante);

                  }
                if ((eq_hydraulique->schema_temps().limpr()) && (impr2 == 1) )
                  {
                    if ( (typeswitch == 0) )        //selon le critere de switch choisit:y* ou nu_t
                      Cout << "Changement de couche a la maille " << icouche <<  " (" << zone_VDF_NS.xp(elem_courant,0) << ";" << zone_VDF_NS.xp(elem_courant,1) << ") y* = " << critere_switch << finl;
                    else
                      Cout << "Changement de couche a la maille " << icouche << " (" << zone_VDF_NS.xp(elem_courant,0) << ";" << zone_VDF_NS.xp(elem_courant,1)  << ") nu_t/nu = " << critere_switch << finl;
                  }

              }

          else if (dimension == 3)
            for (int num_face=ndeb; num_face<nfin; num_face++)
              {
                if ( (elem =face_voisins(num_face,0)) != -1) ;
                else
                  elem = face_voisins(num_face,1);
                if (l_unif)
                  d_visco = visco;
                else
                  d_visco = tab_visco[elem];
                dist=zone_VDF_NS.distance_normale(num_face);

                face_courante = num_face;
                elem_courant = elem;

                for(icouche=0; (icouche<nbcouches) && (elem_courant != -1); icouche++)
                  {
                    double kSqRt=1.e-3;
                    if (K_eps(elem_courant,0)>0)
                      kSqRt=sqrt(K_eps(elem_courant,0));
                    y_etoile = kSqRt*dist/d_visco;
                    if (typeswitch == 0)        //selon le critere de switch choisit:y* ou nu_t
                      critere_switch = y_etoile;
                    else
                      critere_switch = visco_turb[elem_courant]/d_visco;
                    if (critere_switch < valswitch)        // si on est en proche paroi
                      {
                        tab_couches[elem_courant]=1;
                        double Leps,Lmu,vvSqRt;
                        if (K_eps(elem_courant,0)>0)
                          loi2couches.LepsLmu(K_eps(elem_courant,0),d_visco,dist,y_etoile,Leps,Lmu,vvSqRt);
                        else
                          loi2couches.LepsLmu(1.e-3,d_visco,dist,y_etoile,Leps,Lmu,vvSqRt);
                        if ( Leps!=0.)
                          Eps[elem_courant] = K_eps(elem_courant,0)*kSqRt/Leps;
                        else Eps[elem_courant] = 1.e-3;
                      }
                    else
                      {
                        break;
                      }


                    if ( elem_faces(elem_courant,0) == face_courante)
                      face_courante = elem_faces(elem_courant,3);
                    else if ( elem_faces(elem_courant,1) == face_courante)
                      face_courante = elem_faces(elem_courant,4);
                    else if ( elem_faces(elem_courant,2) == face_courante)
                      face_courante = elem_faces(elem_courant,5);
                    else if ( elem_faces(elem_courant,3) == face_courante)
                      face_courante = elem_faces(elem_courant,0);
                    else if ( elem_faces(elem_courant,4) == face_courante)
                      face_courante = elem_faces(elem_courant,1);
                    else if ( elem_faces(elem_courant,5) == face_courante)
                      face_courante = elem_faces(elem_courant,2);

                    if ( face_voisins(face_courante,0) != elem_courant)
                      elem_courant = face_voisins(face_courante,0);
                    else
                      elem_courant = face_voisins(face_courante,1);
                    dist+=zone_VDF_NS.distance_normale(face_courante);
                  }
                if ((eq_hydraulique->schema_temps().limpr()) && (impr2 == 1) && (elem_courant != -1))
                  {
                    if ( (typeswitch == 0) )        //selon le critere de switch choisit:y* ou nu_t
                      Cout << "Changement de couche a la maille " << icouche << " (" << zone_VDF_NS.xp(elem_courant,0) << ";" << zone_VDF_NS.xp(elem_courant,1) << ";" << zone_VDF_NS.xp(elem_courant,2) << ") y* = " << critere_switch << finl;
                    else
                      Cout << "Changement de couche a la maille " << icouche << " (" << zone_VDF_NS.xp(elem_courant,0) << ";" << zone_VDF_NS.xp(elem_courant,1) << ";" << zone_VDF_NS.xp(elem_courant,2) << ") nu_t/nu = " << critere_switch << finl;
                  }
              }
        }
    }

  if (axi)
    {
      Champ_Face_VDF& vitesse = ref_cast_non_const(Champ_Face_VDF,eq_hydraulique->inconnue().valeur());
      calculer_terme_production_K_Axi(zone_VDF,vitesse,P,K_eps,visco_turb);
    }
  else
    {
      Champ_Face_VDF& vitesse = ref_cast_non_const(Champ_Face_VDF,eq_hydraulique->inconnue().valeur());
      calculer_terme_production_K(zone_VDF,zone_Cl_VDF_NS,P,K_eps,vit,vitesse,visco_turb);
    }

  for (elem=0; elem<nb_elem; elem++)
    {
      secmem(elem,0) += (P(elem)-Eps(elem))*volumes(elem)*porosite_vol(elem);
      if (K_eps(elem,0) >= 10.e-10)
        secmem(elem,1) += (C1*P(elem)- C2*Eps(elem))*volumes(elem)*porosite_vol(elem)*Eps(elem)/K_eps(elem,0);
    }

}
