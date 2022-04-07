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
// File:        Echange_contact_Rayo_transp_VDF.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement/src/VDF
// Version:     /main/25
//
//////////////////////////////////////////////////////////////////////////////

#include <Echange_contact_Rayo_transp_VDF.h>
#include <Probleme_base.h>
#include <Zone_VDF.h>
#include <Champ_front_uniforme.h>
#include <TRUSTTabs.h>
#include <Champ_Uniforme.h>
#include <Neumann_paroi.h>
#include <Modele_Rayonnement_Milieu_Transparent.h>
#include <Debog.h>

Implemente_instanciable_sans_constructeur(Echange_contact_Rayo_transp_VDF,"Echange_contact_Rayo_transp_VDF",Echange_contact_VDF);

Echange_contact_Rayo_transp_VDF::Echange_contact_Rayo_transp_VDF():num_premiere_face_dans_pb_fluide(-1),alpha_(0.5) {}
// printOn et readOn

Sortie& Echange_contact_Rayo_transp_VDF::printOn(Sortie& is ) const
{
  return is;
}


// Echange_contact_Rayo_transp_sans_relax_VDF comme  Echange_contact_Rayo_transp_VDF avec alpha_=0

Implemente_instanciable_sans_constructeur(Echange_contact_Rayo_transp_sans_relax_VDF,"Echange_contact_Rayo_transp_sans_relax_VDF",Echange_contact_Rayo_transp_VDF);
Sortie& Echange_contact_Rayo_transp_sans_relax_VDF::printOn(Sortie& is ) const
{
  return Echange_contact_Rayo_transp_VDF::printOn(is);
}
Entree& Echange_contact_Rayo_transp_sans_relax_VDF::readOn(Entree& s )
{
  return  Echange_contact_Rayo_transp_VDF::readOn(s);
}
Echange_contact_Rayo_transp_sans_relax_VDF::Echange_contact_Rayo_transp_sans_relax_VDF()
{
  num_premiere_face_dans_pb_fluide=-1;
  alpha_=0;
}



Entree& Echange_contact_Rayo_transp_VDF::readOn(Entree& s )
{
  return  Echange_contact_VDF::readOn(s);
}

void Echange_contact_Rayo_transp_VDF::completer()
{
  Echange_contact_VDF::completer();
  preparer_surface(frontiere_dis(),zone_Cl_dis());

  // calcul de Teta_i 0

  const Equation_base& mon_eqn = zone_Cl_dis().equation();
  const DoubleTab& mon_inco=mon_eqn.inconnue().valeurs();
  const Zone_VDF& ma_zvdf = ref_cast(Zone_VDF,zone_Cl_dis().zone_dis().valeur());
  const Front_VF& ma_front_vf = ref_cast(Front_VF,frontiere_dis());

  int ndeb = ma_front_vf.num_premiere_face();

  int nb_faces_bord = ma_front_vf.nb_faces();


  for (int numfa=0; numfa<nb_faces_bord; numfa++)
    {

      int ind_fac = numfa+ndeb;
      int elem;
      if (ma_zvdf.face_voisins(ind_fac,0)!= -1)
        elem = ma_zvdf.face_voisins(ind_fac,0);
      else
        elem = ma_zvdf.face_voisins(ind_fac,1);
      Teta_i(numfa)=mon_inco(elem);
    }
  return;

}



void Echange_contact_Rayo_transp_VDF::mettre_a_jour(double temps)
{
  //Cerr << " Echange_contact_Rayo_transp_VDF::mettre_a_jour au temps " << temps << finl;
  //  Echange_contact_VDF::mettre_a_jour(temps);
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
      if (mon_eqn.probleme().milieu().is_rayo_transp())
        {
          // on est du cote fluide
          eqn=&mon_eqn;
          m++;
        }

      if (autre_eqn.probleme().milieu().is_rayo_transp())
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

      if (frontiere_dis().le_nom()!=chcal.front_dis().le_nom())
        {
          Cerr<<"Le nom de bord doit etre le meme pour les deux domaines au niveau du raccord"<<finl;
          exit();
        }

      const Front_VF& frontvf=ref_cast(Front_VF,eqn->zone_dis().valeur().frontiere_dis(frontiere_dis().le_nom()));
      num_premiere_face_dans_pb_fluide=frontvf.num_premiere_face();
      Zone_dis_base& zone_dis1 = zone_Cl_dis().zone_dis().valeur();
      Nom nom_racc1=frontiere_dis().frontiere().le_nom();
      if (zone_dis1.zone().raccord(nom_racc1).valeur().que_suis_je() !="Raccord_distant_homogene")
        verifier_correspondance();
    }
  else
    {
      //Cerr<<T_autre_pb().valeurs()<<finl;
      // exit();
    }
  T_autre_pb().mettre_a_jour(temps);

  assert(nb_comp==1);
  int is_pb_fluide=0;
  // ATTENTION
  //regarder equation pour determiner is_pb_fluide
  // necessaire pour le rayonnement
  DoubleTab& mon_h= h_imp_->valeurs();
  int opt=0;
  assert(h_paroi!=0.);
  double invhparoi=1./h_paroi;

  calculer_h_autre_pb( autre_h, invhparoi, opt);

  calculer_h_mon_pb(mon_h,0.,opt);

  calculer_Teta_paroi(T_ext()->valeurs_au_temps(temps),mon_h,autre_h,is_pb_fluide,temps);

  // on a calculer Teta paroi, on peut calculer htot dans himp (= mon_h)

  int taille=mon_h.dimension(0);
  for (int ii=0; ii<taille; ii++)
    for (int jj=0; jj<nb_comp; jj++)
      {
        mon_h(ii,jj)=1./(1./autre_h(ii,jj)+1./mon_h(ii,jj));
      }

  Echange_global_impose::mettre_a_jour(temps);
}

void Echange_contact_Rayo_transp_VDF::calculer_Teta_equiv(DoubleTab& Teta_eq,const DoubleTab& mon_h,const DoubleTab& lautre_h,int i,double temps)
{

  Cerr<<"On ne doit pas passer par la Echange_contact_Rayo_transp_VDF::calculer_Teta_equiv "<<__FILE__<<finl;
  exit();

}


void Echange_contact_Rayo_transp_VDF::calculer_Teta_paroi(DoubleTab& Teta_equiv,const DoubleTab& mon_h,const DoubleTab& lautre_h,int i,double temps)
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
  //Teta_equiv.resize(nb_faces_bord,1);
  Modele_Rayonnement_Milieu_Transparent& modrayo=ref_cast(Modele_Rayonnement_Milieu_Transparent,le_modele_rayo.valeur());
  if (modrayo.relaxation()==0)
    alpha_=0;

  // on est oblige de tracer le flux_radiatif en // pour le pb solide
  Champ_front_calc& ch=ref_cast(Champ_front_calc, T_autre_pb().valeur());
  DoubleTab& t_autre=ch.valeurs_au_temps(temps);
  const Front_VF& autre_front_vf=ref_cast(Front_VF, ch.front_dis());

  int ispbfluide=mon_eqn.probleme().milieu().is_rayo_transp();

  const Front_VF& frontfluide=(ispbfluide?ma_front_vf:autre_front_vf);
  int nb_faces_bord_fluide = frontfluide.nb_faces();
  DoubleVect flux_radia_fluide(nb_faces_bord_fluide);
  int ndebfluide=frontfluide.num_premiere_face();
  assert(ndebfluide==num_premiere_face_dans_pb_fluide);
  for (int nfluide=0; nfluide<nb_faces_bord_fluide; nfluide++)
    {
      flux_radia_fluide(nfluide)=modrayo.flux_radiatif(nfluide+ndebfluide);
    }
  DoubleVect fluxradia;
  Zone_dis_base& zone_dis1 = zone_Cl_dis().zone_dis().valeur();
  Nom nom_racc1=frontiere_dis().frontiere().le_nom();
  if ((ispbfluide)||(zone_dis1.zone().raccord(nom_racc1).valeur().que_suis_je() !="Raccord_distant_homogene"))
    {
      fluxradia=flux_radia_fluide;
    }
  else
    {
      autre_front_vf.frontiere().trace_face_distant(flux_radia_fluide,fluxradia);
    }
  for (int numfa=0; numfa<nb_faces_bord; numfa++)
    {
      ind_fac = numfa+ndeb;
      if (ma_zvdf.face_voisins(ind_fac,0)!= -1)
        elem = ma_zvdf.face_voisins(ind_fac,0);
      else
        elem = ma_zvdf.face_voisins(ind_fac,1);

      //double flux_radia =modrayo.flux_radiatif(numfa+num_premiere_face_dans_pb_fluide);
      //assert(flux_radia==modrayo.flux_radiatif(frontiere_dis().le_nom(),numfa));
      double new_Teta_i = (mon_h(numfa,0)*mon_inco(elem) + lautre_h(numfa,0)*t_autre(numfa,0)-fluxradia(numfa))/(mon_h(numfa,0)+lautre_h(numfa,0));
      //double alpha=0.5;

      double tetisa=Teta_i(numfa);
      //alpha_=0.;

      Teta_i(numfa)=(1-alpha_)*new_Teta_i+alpha_*Teta_i(numfa);
      Teta_equiv(numfa,0)= (-mon_h(numfa,0)*mon_inco(elem) + (lautre_h(numfa,0)+mon_h(numfa,0))*Teta_i(numfa))/lautre_h(numfa,0);
      //      Teta_equiv(numfa,0)= Teta_i(numfa);
      if (Teta_i(numfa)<=0)
        {
          Cerr<<"Teta_equiv fin "<<Teta_equiv(numfa,0)<< " Teta_i_old "<<tetisa<<" Teta_i "<< Teta_i(numfa);
          Cerr<<" Teta_equiv avt modif "<<(mon_h(numfa,0)*mon_inco(elem) + lautre_h(numfa,0)*t_autre(numfa)-fluxradia(numfa))/(mon_h(numfa,0)+lautre_h(numfa,0));
          Cerr<<" h "<<mon_h(numfa,0)<<" autreh "<<lautre_h(numfa,0)<<" flux radia "<<fluxradia(numfa);
          Cerr<<" mon T "<<mon_inco(elem)<<" autre T "<< t_autre(numfa)<<finl;
          exit();
        }
      //        assert(Teta_equiv(numfa,0)>0);
    }
  Debog::verifier("Teta_equiv dans Echange_contact_Rayo_transp_VDF::calculer_Teta_paroi",Teta_equiv);
}
