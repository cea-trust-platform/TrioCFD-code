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
// File:        Face_Rayonnante.cpp
// Directory:   $TRUST_ROOT/src/Rayonnement
// Version:     /main/22
//
//////////////////////////////////////////////////////////////////////////////

#include <Face_Rayonnante.h>
#include <Zone_VF.h>
#include <Zone_Cl_dis_base.h>
#include <EcrFicPartage.h>

Implemente_instanciable_sans_destructeur(Face_Rayonnante,"Face_Rayonnante",Objet_U);


Face_Rayonnante::~Face_Rayonnante()
{
  //ecrire_temperature_bord();
}
Entree& Face_Rayonnante::readOn(Entree& is)
{

  is >> nom_bord_rayo_lu_;
  nom_bord_rayo_=nom_bord_rayo_lu_;

  /* A debugger: sur gcc 3.2, les ecritures sur la sortie erreur ne marchent plus
     apres un appel a Face_Rayonnante::readOn !
  */
  // on recupere le nom du bord (c.a.d la partie a gauche du %)
  const char* marq=strchr(nom_bord_rayo_lu_,'%');
  if (marq)
    {
      nom_bord_rayo_.prefix(marq);
      Cerr<< nom_bord_rayo_lu_ << " associe au bord " << nom_bord_rayo_ <<finl;
    }
  is >> surf_;
  is >> emissivite_;
  // dimensionnement des ensemble_face
  nb_ensembles_faces_ =1;
  Les_ensembles_faces_bord.dimensionner(nb_ensembles_faces());

  flux_radiatif_ = 0;
  // on ne  peut pas le faire avant de typer!!!
  //  if (nom_bord_rayo_lu_!=nom_bord_rayo_)
  //   {
  //       Ensemble_Faces_base& faces_j =  ensembles_faces_bord(0).valeur();
  //       faces_j.lire(nom_bord_rayo_lu_);
  //     }
  return is;
}


Sortie& Face_Rayonnante::printOn(Sortie& os) const
{
  return os;
}

double  Face_Rayonnante::calculer_temperature()
{
  //Cerr << "Face_Rayonnante::calculer_temperature sur le bord de nom " << nom_bord_rayo_ << finl;
  // cas particulier ou Ensemble_Faces coincide avec un bord: nb_ensembles_faces() = 1
  for (int j=0; j<nb_ensembles_faces(); j++)
    {
      Ensemble_Faces_base& faces_j =  ensembles_faces_bord(j);
      //Cond_Lim_Rayo& la_cl_rayon= faces_j.cond_lim_rayo();
      double sum_surf=0;
      double sum_T=0;
      if (faces_j.nb_faces_bord()!=0)
        {
          const Frontiere& le_bord =  faces_j.la_cl_base().frontiere_dis().frontiere();
          double surf;
          int nbfacesbord=le_bord.nb_faces();
          const IntVect& tables_faces =faces_j.Table_faces();
          int nbfaces_ens=tables_faces.size();
          if ((nbfacesbord==nbfaces_ens)||(nbfaces_ens==0))
            {
              for (int face=0; face<nbfacesbord; face++)
                {
                  surf = faces_j.surface(face);
                  sum_surf += surf;
                  sum_T += surf*pow(faces_j.teta_i(face),4);
                }
            }
          else
            {
              for (int face=0; face<nbfaces_ens; face++)
                {
                  int face2=tables_faces(face);
                  surf = faces_j.surface(face2);
                  sum_surf += surf;
                  sum_T += surf*pow(faces_j.teta_i(face2),4);
                }
            }
        }
      double sum_T_tot=sum_T;
      double sum_surf_tot=sum_surf;
      sum_surf_tot=mp_sum(sum_surf_tot);
      // test important mais gene validation Rayo
      //if (!est_egal(sum_surf_tot,surf_,1e-4)) { Cerr<< " Pb surface calculee "<<sum_surf_tot<<" surfaces du fichier facesrayo "<<surf_<<finl;abort();}
      //surf_=sum_surf_tot;
      sum_T_tot=mp_sum(sum_T_tot);
      double coef = sum_T_tot/sum_surf_tot;
      T_face_rayo_ =  pow(coef,0.25);
    }
  return T_face_rayo_;
}

double Face_Rayonnante::imprimer_flux_radiatif(Sortie& os, Sortie& os1, Sortie& os2 ) const
{

  double flux_=flux_radiatif_*surface_rayo();
  Nom espace="\t";
  os << nom_bord_rayo_lu() << "\t: " <<  flux_ << " W (temperature bord:"<< T_face_rayo_ <<" K)" << finl;
  os1 << espace << flux_;
  os2 << espace << T_face_rayo_;
  return flux_;
}

int Face_Rayonnante::chercher_ensemble_faces(const Nom& nom_bord) const
{

  // recherche si le nom_bord passe est une face rayonnante.

  if (nom_bord==nom_bord_rayo_)
    {
      if (emissivite_ != -1)
        {
          return 1;
        }
      else
        {
          return -1;
        }
    }
  // else
  //     {
  //       return -1;
  //     }
  return -1;
}
void  Face_Rayonnante::ecrire_temperature_bord() const
{
  if (emissivite_==-1) return;
  Nom fic(nom_bord_rayo_lu());
  fic+="_Temp";
  EcrFicPartage sortie(fic);
  if (je_suis_maitre())
    {
      Cerr<<"ecriture des temperatures de bord dans "<<fic<<finl;
      sortie <<"# temperatures bord du bord "<<nom_bord_rayo_lu()<<" emissivites "<<emissivite_<<finl;
      sortie <<"# x y (z) temperature surface"<<finl;
    }
  sortie.lockfile();
  for (int j=0; j<nb_ensembles_faces(); j++)
    {
      // cast en dur a cause de teta_i
      // a nettoyer quand teta_i sera une methode const...
      Ensemble_Faces_base& faces_j =  ref_cast_non_const(Ensemble_Faces_base,ensembles_faces_bord(j));
      double T=0;
      if (faces_j.nb_faces_bord()!=0)
        {
          const Frontiere& le_bord =  faces_j.la_cl_base().frontiere_dis().frontiere();
          const Zone_VF& zone_dis=ref_cast(Zone_VF, faces_j.la_cl_base().zone_Cl_dis().zone_dis().valeur());
          const DoubleTab& xv=zone_dis.xv();
          double surf;
          int nbfacesbord=le_bord.nb_faces();
          const IntVect& tables_faces =faces_j.Table_faces();
          int nbfaces_ens=tables_faces.size();
          if ((nbfacesbord==nbfaces_ens)||(nbfaces_ens==0))
            {
              for (int face=0; face<nbfacesbord; face++)
                {
                  surf = faces_j.surface(face);
                  T=faces_j.teta_i(face);
                  sortie<<xv(face,0)<<" "<<xv(face,1)<<" ";
                  if (dimension==3)
                    sortie<<xv(face,1)<<" ";
                  sortie<<T<<" "<<surf<<finl;
                }
            }
          else
            {
              for (int face=0; face<nbfaces_ens; face++)
                {
                  int face2=tables_faces(face);
                  surf = faces_j.surface(face2);
                  T=faces_j.teta_i(face2);
                  sortie<<xv(face,0)<<" "<<xv(face,1)<<" ";
                  if (dimension==3)
                    sortie<<xv(face,1)<<" ";
                  sortie<<T<<" "<<surf<<finl;
                }
            }
        }

    }
  sortie.unlockfile();
  sortie.syncfile();
}
