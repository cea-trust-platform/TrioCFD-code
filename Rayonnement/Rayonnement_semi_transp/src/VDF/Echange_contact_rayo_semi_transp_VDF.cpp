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
// File:        Echange_contact_rayo_semi_transp_VDF.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement_semi_transp/src/VDF
// Version:     /main/19
//
//////////////////////////////////////////////////////////////////////////////

#include <Echange_contact_rayo_semi_transp_VDF.h>
#include <Pb_Conduction.h>
#include <Modele_rayo_semi_transp.h>
#include <Champ_front_calc.h>
#include <Milieu_base.h>
#include <Zone_VDF.h>

Implemente_instanciable_sans_constructeur(Echange_contact_rayo_semi_transp_VDF,"Paroi_echange_contact_rayo_semi_transp_VDF",Echange_contact_VDF);


Echange_contact_rayo_semi_transp_VDF::Echange_contact_rayo_semi_transp_VDF():num_premiere_face_dans_pb_fluide(-1) {}

/*! @brief
 *
 * @param (Sortie& os) un flot de sortie
 * @return (Sortie&) le flot de sortie modifie
 */
Sortie& Echange_contact_rayo_semi_transp_VDF::printOn(Sortie& os) const
{
  return os;
}


/*! @brief Lecture et construction du champ_front_calc T_autre_Pb qui servira pour calculer les conditions d'echange a la paroi
 *
 *     Lecture de l'emissivite de la paroi
 *     Lecture du coefficient A
 *        typage et dimenssionnement du champ_front T_ext() qui contient
 *        les temperatures de paroi et qui est utilise lors du calcul des
 *        flux diffusifs.
 *
 *
 * @param (Entree& is) un flot d'entree
 * @return (Entree&) le flot d'entree modifie
 */
Entree& Echange_contact_rayo_semi_transp_VDF::readOn(Entree& is)
{
  return Echange_contact_VDF::readOn(is);
}


const Cond_lim_base& Echange_contact_rayo_semi_transp_VDF::la_cl() const
{
  return (*this);
}



/*! @brief Renvoie le Champ_front des temperatures de bord.
 *
 */
Champ_front& Echange_contact_rayo_semi_transp_VDF::temperature_bord()
{
  return T_ext();
}


/*! @brief la methode mettre_a_jour(temps) a pour role de remplire le tableau T_ext avec les temperatures de paroi.
 *
 * Ces temperatures de paroi sont utilisees pour imposer la condition a la limite du
 *       probleme de rayonnement.
 *       Les temperatures de paroi sont calculees en tenant compte du caractere rayonnant
 *       de la paroi.
 *
 *
 */
void Echange_contact_rayo_semi_transp_VDF::calculer_temperature_bord(double temps)
{
  if (T_autre_pb()->valeurs_au_temps(temps).size()==0)
    {
      T_autre_pb().associer_fr_dis_base(T_ext().frontiere_dis());
      Zone_dis_base& zone_dis1 = zone_Cl_dis().zone_dis().valeur();
      Nom nom_racc1=frontiere_dis().frontiere().le_nom();
      if (zone_dis1.zone().raccord(nom_racc1).valeur().que_suis_je() !="Raccord_distant_homogene")
        verifier_correspondance();
    }


  int is_pb_fluide=0;
  DoubleTab& mon_h= h_imp_->valeurs();

  int opt=0;
  assert(h_paroi!=0.);
  double invhparoi=1./h_paroi;
  calculer_h_autre_pb( autre_h, invhparoi, opt);
  calculer_h_mon_pb(mon_h,0.,opt);

  calculer_Teta_paroi(T_paroi->valeurs_au_temps(temps),mon_h,autre_h,is_pb_fluide,temps);
}


/*! @brief la methode mettre_a_jour(temps) a pour role de remplire le tableau T_ext avec les temperatures de paroi.
 *
 * T_ext sera utilise pour calculer le flux impose a la limite de chaque  domaine du
 *       probleme couple.
 *       Les temperatures de paroi sont calculees en tenant compte du caractere rayonnant
 *       de la paroi.
 *       Appel de la methode Echange_externe_impose::mettre_a_jour(temps);
 *
 *
 * @param (double temps) temps courant
 */
void Echange_contact_rayo_semi_transp_VDF::mettre_a_jour(double temps)
{
  //  Cerr<<"Echange_contact_rayo_semi_transp_VDF::mettre_a_jour : Debut"<<finl;

  Champ_front_calc& ch=ref_cast(Champ_front_calc, T_autre_pb().valeur());
  const Milieu_base& le_milieu=ch.milieu();
  int nb_comp = le_milieu.conductivite()->nb_comp();

  if (num_premiere_face_dans_pb_fluide==-1)
    {
      Cerr<<"fin de construction dans "<<que_suis_je()<<finl;
      T_autre_pb().associer_fr_dis_base(T_ext().frontiere_dis());
      // on regarde qui est le pb fluide
      const Equation_base* eqn=NULL;
      const Equation_base& mon_eqn = zone_Cl_dis().equation();

      const Champ_front_calc& chcal=ref_cast(Champ_front_calc,T_autre_pb().valeur());
      const  Equation_base& autre_eqn=chcal.inconnue().equation();
      int m=0;
      ////if (mon_eqn.probleme().equation(0).comprend_mot("vitesse"))
      /*
        Noms liste_noms_mon_pb;
        mon_eqn.probleme().get_noms_champs_postraitables(liste_noms_mon_pb);
        //if (liste_noms_mon_pb.contient("vitesse"))
        for (int i=0; i<liste_noms_mon_pb.size();i++)
        if (liste_noms_mon_pb[i]=="vitesse")
        {
        // on est du cote fluide
        eqn=&mon_eqn;
        m++;
        }
      */
      if (mon_eqn.probleme().milieu().is_rayo_semi_transp())
        {
          // on est du cote fluide
          eqn=&mon_eqn;
          m++;
        }

      ////if (autre_eqn.probleme().equation(0).comprend_mot("vitesse"))
      /*
        Noms liste_noms_autre_pb;
        autre_eqn.probleme().get_noms_champs_postraitables(liste_noms_autre_pb);
        //if (liste_noms_autre_pb.contient("vitesse"))
        for (int i=0; i<liste_noms_autre_pb.size();i++)
        if (liste_noms_autre_pb[i]=="vitesse")
        {
        eqn=&autre_eqn;
        m++;

        }
      */
      if (autre_eqn.probleme().milieu().is_rayo_semi_transp())
        {
          eqn=&autre_eqn;
          m++;

        }
      if (m!=1)
        {
          Cerr<<"gros pb "<<que_suis_je()<<finl;
          assert(0);
          exit();
        }

      const Front_VF& frontvf=ref_cast(Front_VF,eqn->zone_dis().valeur().frontiere_dis(frontiere_dis().le_nom()));
      num_premiere_face_dans_pb_fluide=frontvf.num_premiere_face();

      Zone_dis_base& zone_dis1 = zone_Cl_dis().zone_dis().valeur();
      Nom nom_racc1=frontiere_dis().frontiere().le_nom();
      if (zone_dis1.zone().raccord(nom_racc1).valeur().que_suis_je() !="Raccord_distant_homogene")
        verifier_correspondance();
    }

  T_autre_pb().mettre_a_jour(temps);
  assert(nb_comp==1);

  DoubleTab& mon_h= h_imp_->valeurs();

  int opt=0;
  assert(h_paroi!=0.);
  double invhparoi=1./h_paroi;
  calculer_h_autre_pb( autre_h, invhparoi, opt);
  calculer_h_mon_pb(mon_h,0.,opt);

  int is_pb_fluide=0;
  calculer_Teta_equiv(T_ext()->valeurs_au_temps(temps),mon_h,autre_h,is_pb_fluide,temps);

  // on a calcule T_ext, on peut calculer heq dans himp (= mon_h)
  int taille=mon_h.dimension(0);
  for (int ii=0; ii<taille; ii++)
    for (int jj=0; jj<nb_comp; jj++)
      {
        mon_h(ii,jj)=1./(1./autre_h(ii,jj)+1./mon_h(ii,jj));
      }

  Echange_global_impose::mettre_a_jour(temps);
}


int Echange_contact_rayo_semi_transp_VDF::compatible_avec_eqn(const Equation_base& eqn) const
{
  Motcle dom_app=eqn.domaine_application();
  Motcle Thermique="Thermique";
  Motcle indetermine="indetermine";
  if ( (dom_app==Thermique) || (dom_app==indetermine) )
    return 1;
  else
    {
      err_pas_compatible(eqn);
      return 0;
    }
}


void  Echange_contact_rayo_semi_transp_VDF::completer()
{
  Echange_contact_VDF::completer();

  // On dimenssionne le champ_front T_ext qui va contenir les temperatures de paroi
  int nb_comp = 1;

  // T_ext et Tp ne sont pas les memes, on type T_paroi
  T_paroi.typer("Champ_front_fonc");
  T_paroi->fixer_nb_comp(nb_comp);
}


void Echange_contact_rayo_semi_transp_VDF::calculer_Teta_paroi(DoubleTab& Teta_p,const DoubleTab& mon_h,const DoubleTab& lautre_h,int i,double temps)
{
  // i servira en radiatif
  //
  //mon_h(monT-Tw)=autre_h*(Tw-Tautre)
  // Tautre=Text;
  //  Tw=(autr_h*Tautre+mon_h*monT)/(autr_h+mon_h)
  const Equation_base& mon_eqn = zone_Cl_dis().equation();
  const DoubleTab& mon_inco=mon_eqn.inconnue().valeurs();
  const Zone_VDF& ma_zvdf = ref_cast(Zone_VDF,zone_Cl_dis().zone_dis().valeur());
  const Front_VF& ma_front_vf = ref_cast(Front_VF,frontiere_dis());

  //DoubleTab& mon_h=h_imp_->valeurs();
  int ndeb = ma_front_vf.num_premiere_face();
  //int ndeb2 = l_autre_front_vf.num_premiere_face();
  int nb_faces_bord = ma_front_vf.nb_faces();
  int ind_fac,elem;
  //DoubleTab& Teta_i=T_ext().valeurs();
  Teta_p.resize(nb_faces_bord,1);
  DoubleTab& t_autre=T_autre_pb()->valeurs_au_temps(temps);

  const Modele_rayo_semi_transp& le_modele = modele();
  const DoubleTab& flux_radiatif = le_modele.flux_radiatif(frontiere_dis().le_nom()).valeurs();
  for (int numfa=0; numfa<nb_faces_bord; numfa++)
    {
      ind_fac = numfa+ndeb;
      if (ma_zvdf.face_voisins(ind_fac,0)!= -1)
        elem = ma_zvdf.face_voisins(ind_fac,0);
      else
        elem = ma_zvdf.face_voisins(ind_fac,1);

      double flux_radia = flux_radiatif(numfa,0);
      Teta_p(numfa,0) = (mon_h(numfa,0)*mon_inco(elem) + lautre_h(numfa,0)*t_autre(numfa,0)-flux_radia)/(mon_h(numfa,0)+lautre_h(numfa,0));

    }

  const Probleme_base& pb =  ma_zone_cl_dis->equation().probleme();
  if (sub_type(Pb_Conduction,pb))
    {
      Cerr<<"On ne devrait pas a avoir a calculer la temperature de bord pour le probleme"<<finl;
      Cerr<<"solide "<<finl;
      exit();
    }
}


void Echange_contact_rayo_semi_transp_VDF::calculer_Teta_equiv(DoubleTab& Teta_eq,const DoubleTab& mon_h,const DoubleTab& lautre_h,int i,double temps)
{

  // on ne peut plus simple!!
  // Teta_equiv=T_autre_pb

  // i servira en radiatif
  //
  // mon_h(monT-Tw)=autre_h*(Tw-Tautre)
  // Tautre=Text;
  // Tw=(autr_h*Tautre+mon_h*monT)/(autr_h+mon_h)
  const Front_VF& ma_front_vf = ref_cast(Front_VF,frontiere_dis());
  int nb_faces_bord = ma_front_vf.nb_faces();

  assert(Teta_eq.dimension(0)==nb_faces_bord);
  assert(Teta_eq.dimension(1)==1);
  const Modele_rayo_semi_transp& le_modele = modele();
  const DoubleTab& flux_radiatif = le_modele.flux_radiatif(frontiere_dis().le_nom()).valeurs();
  DoubleTab& t_autre=T_autre_pb()->valeurs_au_temps(temps);
  for (int numfa=0; numfa<nb_faces_bord; numfa++)
    Teta_eq(numfa,0) = t_autre(numfa,0) - (1/lautre_h(numfa,0))*flux_radiatif(numfa,0);

  Teta_eq.echange_espace_virtuel();

}

/*! @brief Renvoie la CL portee par le probleme de l'autre cote de la frontiere.
 *
 */
Echange_contact_rayo_semi_transp_VDF& Echange_contact_rayo_semi_transp_VDF::la_Cl_opposee()
{
  Champ_front_calc& ch = ref_cast(Champ_front_calc, T_autre_pb().valeur());
  const Zone_Cl_dis_base& zcld = ch.zone_Cl_dis();
  for (int i=0; i<zcld.nb_cond_lim(); i++)
    {
      const Cond_lim& lacl=zcld.les_conditions_limites(i);
      if (lacl.frontiere_dis().le_nom()==frontiere_dis().le_nom())
        return ref_cast_non_const(Echange_contact_rayo_semi_transp_VDF,lacl.valeur());
    }

  Cerr << "Erreur lors de la recherche de la CL opposee." << finl;
  exit();

  // Pour le compilo
  return *this;
}
