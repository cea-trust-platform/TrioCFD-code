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
// File:        Flux_radiatif_VDF.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src/VDF
// Version:     /main/18
//
//////////////////////////////////////////////////////////////////////////////

#include <Flux_radiatif_VDF.h>
#include <Neumann_paroi_rayo_semi_transp_VDF.h>
#include <Echange_contact_rayo_semi_transp_VDF.h>
#include <Echange_externe_impose_rayo_semi_transp.h>
#include <Echange_global_impose_rayo_semi_transp.h>
#include <Frontiere_ouverte_temperature_imposee_rayo_semi_transp.h>
#include <Frontiere_ouverte_rayo_semi_transp.h>
#include <Eq_rayo_semi_transp_VDF.h>
#include <Modele_rayo_semi_transp.h>
#include <Fluide_base.h>
#include <Champ_Uniforme.h>
#include <Champ_front_uniforme.h>
#include <Domaine_VDF.h>
#include <Schema_Temps_base.h>

Implemente_instanciable(Flux_radiatif_VDF,"Flux_radiatif_VDF",Flux_radiatif_base);

/*! @brief Imprime le type de l'equation sur un flot de sortie.
 *
 * @param (Sortie& s) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Flux_radiatif_VDF::printOn(Sortie& s ) const
{
  return s << que_suis_je() << "\n";
}



/*! @brief
 *
 * @param (Entree& s) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Flux_radiatif_VDF::readOn(Entree& s )
{
  return Flux_radiatif_base::readOn(s);
}


/*! @brief
 *
 */
void Flux_radiatif_VDF::evaluer_cl_rayonnement(Champ_front_base& Tb, const Champ_Don_base& coeff_abs, const Champ_Don_base& longueur_rayo,
                                               const Champ_Don_base& indice,const Domaine_VF& zvf, const double sigma, double temps)
{
  //  Cerr<<"Flux_radiatif_VDF::evaluer_cl_rayonnement() : Debut"<<finl;
  const DoubleTab& l_rayo = longueur_rayo.valeurs();
  const DoubleTab& n = indice.valeurs();
  const DoubleTab& epsilon = emissivite().valeurs();
  const DoubleTab& kappa = coeff_abs.valeurs();

  const Domaine_VDF& zvdf = ref_cast(Domaine_VDF,zvf);
  const Front_VF& le_bord = ref_cast(Front_VF,frontiere_dis());
  const IntTab& face_voisins = zvdf.face_voisins();

  //  const DoubleTab& xv = zvdf.xv();

  // On dimensionne le_champ_front
  assert(le_champ_front->nb_comp() == 1);
  DoubleTab& tab = le_champ_front->valeurs_au_temps(temps);
  tab.resize(le_bord.nb_faces(),le_champ_front->nb_comp());

  // Boucle sur les faces de le_bord
  int ndeb = le_bord.num_premiere_face();
  int nfin = ndeb + le_bord.nb_faces();
  int face = ndeb;

  //  Cerr<<"nom du bord : "<<le_bord.le_nom()<<finl;
  for(face = ndeb; face<nfin; face++)
    {
      int elem = face_voisins(face,0);

      // Mise en commentaire de la determination de signe car on a
      // corrige le bug sur la condition a la limite Neumann_paroi
      double signe = 1;
      //      if(elem < 0)
      //        {
      //          elem = face_voisins(face,1);
      //          signe *= -1;
      //        }

      double eF = zvdf.dist_norm_bord(face);

      double nn;
      assert(indice.nb_comp() == 1);
      if(sub_type(Champ_Uniforme,indice))
        nn = n(0,0);
      else
        nn = n(elem,0);

      double l_r;
      assert(longueur_rayo.nb_comp() == 1);
      if(sub_type(Champ_Uniforme,longueur_rayo))
        l_r = l_rayo(0,0);
      else
        l_r = l_rayo(elem,0);

      double k;
      assert(coeff_abs.nb_comp() == 1);
      if(sub_type(Champ_Uniforme,coeff_abs))
        k = kappa(0,0);
      else
        k = kappa(elem,0);

      // determination de la temperature de paroi en fonction
      // de la face consideree
      double T;
      assert(Tb.nb_comp() == 1);
      if (sub_type(Champ_front_uniforme,Tb))
        T = Tb.valeurs()(0,0);
      else
        T = Tb.valeurs_au_temps(temps)(face-ndeb,0);

      // Determination de l'emissivite de paroi en fonction
      // de la face consideree
      double epsi;
      assert(emissivite().nb_comp() == 1);
      if (sub_type(Champ_front_uniforme,emissivite()))
        epsi = epsilon(0,0);
      else
        epsi = epsilon(face-ndeb,0);

      // Remplissage de la condition a la limite pour l'equation de
      // rayonnement
      double numer_coeff = l_r;
      numer_coeff *= 4*nn*nn*sigma*pow(T,4);

      double denum_coeff = 3*k*epsi;
      denum_coeff = 1/denum_coeff;
      denum_coeff *= A_*(2-epsi);
      denum_coeff = denum_coeff + eF;

      if(epsi<DMINFLOAT)
        tab(face-ndeb,0) = 0;
      else
        tab(face-ndeb,0) = signe*numer_coeff/denum_coeff;
    }
  //  Cerr<<"Flux_radiatif_VDF::evaluer_cl_rayonnement() : Fin"<<finl;
  tab.echange_espace_virtuel();
}


void Flux_radiatif_VDF::calculer_flux_radiatif(const Equation_base& eq_temp)
{
  //  Cerr<<"Flux_radiatif_VDF::calculer_flux_radiatif : Debut"<<finl;
  // On doit recuperer la temperature de bord
  const Front_VF& le_bord = ref_cast(Front_VF,frontiere_dis());
  int nb_faces = le_bord.nb_faces();
  OBS_PTR(Champ_front_base) Tb;
  const Conds_lim& les_cl_temp = eq_temp.domaine_Cl_dis().les_conditions_limites();
  int num_cl_temp = 0;

  int test_nom = 0;
  for(num_cl_temp = 0; num_cl_temp<les_cl_temp.size(); num_cl_temp++)
    {
      const Cond_lim& la_cl_temp = eq_temp.domaine_Cl_dis().les_conditions_limites(num_cl_temp);
      Nom nom_cl_temp = la_cl_temp->frontiere_dis().le_nom();
      if(nom_cl_temp == frontiere_dis().le_nom())
        {
          test_nom = 1;
          if (sub_type(Neumann_paroi_rayo_semi_transp_VDF,la_cl_temp.valeur()))
            {
              const Neumann_paroi_rayo_semi_transp_VDF& la_cl_temper
                = ref_cast(Neumann_paroi_rayo_semi_transp_VDF,la_cl_temp.valeur());
              Tb = la_cl_temper.temperature_bord();
            }
          else if (sub_type(Echange_contact_rayo_semi_transp_VDF,la_cl_temp.valeur()))
            {
              Echange_contact_rayo_semi_transp_VDF& la_cl_temper
                = ref_cast_non_const(Echange_contact_rayo_semi_transp_VDF,la_cl_temp.valeur());
              Tb = la_cl_temper.temperature_bord();
            }
          else if (sub_type(Echange_externe_impose_rayo_semi_transp,la_cl_temp.valeur()))
            {
              Echange_externe_impose_rayo_semi_transp& la_cl_temper
                = ref_cast_non_const(Echange_externe_impose_rayo_semi_transp,la_cl_temp.valeur());
              Tb = la_cl_temper.temperature_bord();
            }
          else if (sub_type(Frontiere_ouverte_temperature_imposee_rayo_semi_transp,la_cl_temp.valeur()))
            {
              const Frontiere_ouverte_temperature_imposee_rayo_semi_transp& la_cl_temper
                = ref_cast(Frontiere_ouverte_temperature_imposee_rayo_semi_transp,la_cl_temp.valeur());
              Tb = la_cl_temper.temperature_bord();
            }
          else if (sub_type(Frontiere_ouverte_rayo_semi_transp,la_cl_temp.valeur()))
            {
              const Frontiere_ouverte_rayo_semi_transp& la_cl_temper
                = ref_cast(Frontiere_ouverte_rayo_semi_transp,la_cl_temp.valeur());
              Tb = la_cl_temper.temperature_bord();
            }
          else if (sub_type(Echange_global_impose_rayo_semi_transp,la_cl_temp.valeur()))
            {
              Echange_global_impose_rayo_semi_transp& la_cl_temper
                = ref_cast_non_const(Echange_global_impose_rayo_semi_transp,la_cl_temp.valeur());
              Tb = la_cl_temper.temperature_bord();
            }
          else
            {
              Cerr <<"Coder pour les autres condition limites de l'equation de temperature 1 "<<finl;
              exit();
            }
        }
    }
  if(test_nom == 0)
    {
      Cerr<<"Erreur : il n'y a pas de condition limite sur une frontiere portant le nom : "<<le_nom()<<finl;
      exit();
    }
  //  Tb contient les temperatures de bord
  //  Cerr<<"Tb = "<<Tb.valeurs()<<finl;
  // Calcul du flux radiatif
  DoubleTab& Flux = flux_radiatif_->valeurs();
  Flux.resize(le_bord.nb_faces(),1);
  Eq_rayo_semi_transp_VDF& eq_rayo = ref_cast( Eq_rayo_semi_transp_VDF,domaine_Cl_dis().equation());
  Fluide_base& fluide = eq_rayo.fluide();
  DoubleTab& kappa = fluide.kappa().valeurs();
  DoubleTab& indice = fluide.indice().valeurs();

  DoubleTab& irradiance = eq_rayo.inconnue().valeurs();

  const Domaine_VDF& zvdf = ref_cast(Domaine_VDF,domaine_Cl_dis().domaine_dis());
  const IntTab& face_voisins = zvdf.face_voisins();
  const DoubleVect& face_surfaces = zvdf.face_surfaces();

  double bilan_flux=0;
  // On fait une boucle sur les faces
  int face=0;
  int ndeb = le_bord.num_premiere_face();
  //  Cerr<<"ndeb = "<<ndeb<<finl;

  for(face=0; face<nb_faces; face++)
    {
      int elem = face_voisins(face+ndeb,0);
      if (elem < 0)
        elem = face_voisins(face+ndeb,1);

      double G_F = irradiance(elem);
      double kappa_F;
      assert(fluide.kappa().nb_comp()==1);
      if(sub_type(Champ_Uniforme,fluide.kappa()))
        kappa_F = kappa(0,0);
      else
        kappa_F = kappa(elem,0);
      double epsi;
      assert(emissivite().nb_comp() == 1);
      if (sub_type(Champ_front_uniforme,emissivite()))
        epsi = emissivite().valeurs()(0,0);
      else
        epsi = emissivite().valeurs()(face,0);
      double eF = zvdf.dist_norm_bord(face+ndeb);
      double n;
      assert(fluide.indice().nb_comp() == 1);
      if(sub_type(Champ_Uniforme,fluide.indice()))
        n = indice(0,0);
      else
        n = indice(elem,0);

      double sigma = eq_rayo.Modele().valeur_sigma();
      double Tbord;

      assert(Tb.valeur().nb_comp() == 1);
      if(sub_type(Champ_front_uniforme,Tb.valeur()))
        Tbord = Tb.valeur().valeurs()(0,0);
      else
        Tbord = Tb.valeur().valeurs()(face,0);

      double denum = A()*(2-epsi);
      denum /= 3*kappa_F*epsi;
      denum += eF;

      double numer = G_F - 4*n*n*sigma*pow(Tbord,4);
      double grad_G = numer/denum;

      if(epsi<DMINFLOAT)
        Flux(face,0) = 0;
      else
        Flux(face,0) = -grad_G/(3*kappa_F);

      bilan_flux += face_surfaces(face+ndeb)*Flux(face,0);
      //      Cerr<<"G_F = "<<G_F<<", Tbord = "<<Tbord<<", Flux("<<face<<",0) = "<<Flux(face,0)<<", eF = "<<eF<<"kappa_F = "<<kappa_F<<finl;
    }
  if (eq_rayo.schema_temps().limpr())
    Cout << "Flux radiatif sur le bord "<<le_bord.le_nom()<<" : "<<bilan_flux<<finl;
}
