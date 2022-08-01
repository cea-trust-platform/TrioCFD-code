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
// File:        ConcatAnsys.cpp
// Directory:   $TRUST_ROOT/../Composants/TrioCFD/Rayonnement/src
// Version:     /main/19
//
//////////////////////////////////////////////////////////////////////////////

#include <ConcatAnsys.h>
#include <EFichier.h>
#include <SFichier.h>

Implemente_instanciable(ConcatAnsys,"ConcatAnsys",Interprete_geometrique_base);

Sortie& ConcatAnsys::printOn(Sortie& os) const
{
  return Interprete::printOn(os);
}

Entree& ConcatAnsys::readOn(Entree& is)
{
  return Interprete::readOn(is);
}

Entree& ConcatAnsys::interpreter_(Entree& is)
{
  Nom nom,nom2,nomS,nom_dom;
  is >> nom_dom >> nom;
  associer_domaine(nom_dom);
  Domaine& dom=domaine();
  nomS=nom_du_cas();
  nom2=nom_du_cas();
  nom2+=Nom(".factforme");
  Nom nomS2=nom_du_cas();
  nomS2+=Nom(".facesrayo");
  Zone& zone=dom.zone(0);
  // Definition du nombre de faces rayonnantes
  int nombre_faces_rayonnantes=zone.nb_bords()+zone.nb_raccords();
  // Tableau des facteurs de forme regroupes
  DoubleTab FIJ(nombre_faces_rayonnantes,nombre_faces_rayonnantes);
  DoubleVect SI(nombre_faces_rayonnantes);
  int j,J,i,I;

  Cerr << "Ouverture du fichier des bords " << nomS2 << finl;
  SFichier fr(nomS2);
  fr << nombre_faces_rayonnantes << " " << nombre_faces_rayonnantes << finl;
  Cerr << "Ouverture du fichier Ansys " << nom << finl;
  EFichier premiere_lecture_ansys(nom);
  Nom motlu;
  int nombre_faces_ansys;
  //  int meth=1;
  premiere_lecture_ansys>>motlu;
  if (motlu=="Number")
    {
      // meth=2;
      // Correction car si le nombre de faces depasse 9999 il est illisible
      // dans le fichier Ansys .vf !
      /* while(motlu!="Surfaces")
         premiere_lecture_ansys>>motlu;
         premiere_lecture_ansys>>motlu; */
      // On recupere donc le nombre de faces en comptant
      // le nombre de valeurs lues de la premiere ligne de la matrice
      // des facteurs de formes:
      nombre_faces_ansys=-1;
      while(motlu!="TOTAL=")
        premiere_lecture_ansys>>motlu;
      premiere_lecture_ansys>>motlu;
      while(motlu!="Element" && motlu!="TOTAL=" && motlu!="*****")
        {
          premiere_lecture_ansys>>motlu;
          if (motlu=="WARNING")
            {
              Cerr << "Warning trouve dans le fichier " << nom << ". Calcul des facteurs de " << finl;
              Cerr << "forme interrompu. Contacter le support TRUST." << finl;
              exit();
            }
          nombre_faces_ansys++;
        }
    }
  else
    {
      while(motlu!="ELEMENTS=")
        premiere_lecture_ansys>>motlu;
      premiere_lecture_ansys >> nombre_faces_ansys;
    }
  Cerr << "Nombre de faces qui vont etre lus:" << nombre_faces_ansys << finl;
  EFichier ansys(nom);
  //
  // Cas ou Ansys a regroupe les faces sur les bords (non_hidden==3)
  //
  int deja_regroupe_dans_ansys=0;
  if (nombre_faces_rayonnantes==nombre_faces_ansys)
    deja_regroupe_dans_ansys=1;

  double min_total=1, max_total=0;
  //
  // Calcul du facteurs FIJ en fonction des Fij
  //
  for (I=0; I<nombre_faces_rayonnantes; I++)
    {
      Faces& faces=(I<zone.nb_bords()?zone.bord(I).faces():zone.raccord(I-zone.nb_bords()).valeur().faces());
      int nb_faces=faces.nb_faces();
      // Calcul des surfaces de chaque face
      DoubleVect Si;
      faces.associer_zone(zone);
      faces.calculer_surfaces(Si);
      for (i=0; i<nb_faces; i++)
        {
          SI(I)+=Si(i);
        }
      // Ecriture dans le fichier des bords pour TRUST (avec par defaut une emissivite de 1)
      fr << (I<zone.nb_bords()?zone.bord(I).le_nom():zone.raccord(I-zone.nb_bords()).le_nom()) << " " << SI(I) << " 1." << finl;
      if (deja_regroupe_dans_ansys)
        nb_faces=1;
      for (i=0; i<nb_faces; i++)
        {
          double fij,total_lu;
          // Lecture du fichier Ansys
          ansys>>motlu;
          while(motlu!="TOTAL=")
            ansys>>motlu;
          ansys>>total_lu;
          Cerr << "\rTotal des facteurs de forme sur la ligne " <<I*nb_faces+i+1<< " : " << total_lu << finl;
          for (J=0; J<nombre_faces_rayonnantes; J++)
            {
              double FiJ=0;
              int nombre_faces=(J<zone.nb_bords()?zone.bord(J).faces().nb_faces():zone.raccord(J-zone.nb_bords()).valeur().faces().nb_faces());
              if (deja_regroupe_dans_ansys)
                nombre_faces=1;
              for (j=0; j<nombre_faces; j++)
                {
                  ansys>>fij;
                  FiJ+=fij;
                }
              if (deja_regroupe_dans_ansys)
                FIJ(I,J)=FiJ;
              else
                FIJ(I,J)+=Si(i)*FiJ/SI(I);
            }
          if (sup_strict(total_lu,1.001))
            {
              Cerr << "La somme des facteurs de forme (" << total_lu << ") de la ligne " << I+1 << " depasse 1 !" << finl;
              Cerr << "Il se peut que vous ayez fait un calcul Ansys avec" << finl;
              Cerr << "option faces non cachees qui ne soit pas justifie." << finl;
              exit();
            }
          else if (est_egal(total_lu,0))
            {
              Cerr << "La somme des facteurs de forme de la ligne " << I+1 << " vaut 0 !" << finl;
              Cerr << "Cela concerne le bord " << (I<zone.nb_bords()?zone.bord(I).le_nom():zone.raccord(I-zone.nb_bords()).le_nom()) << finl;
              if (est_egal(SI(I),0))
                {
                  Cerr << "La surface de ce bord est nulle..." << finl;
                }
              if (Objet_U::bidim_axi)
                {
                  Cerr << "Ce bord est peut etre sur l'axe de revolution de votre calcul 2D axisymetrique."  << finl;
                  Cerr << "Retirez ce bord de votre maillage .geom et relancer le calcul des facteurs de forme." << finl;
                }
              else
                {
                  Cerr << "Une face du maillage .geom est mal orientee." << finl;
                  Cerr << "Contactez le support TRUST." << finl;
                }
              exit();
            }
          else
            {
              min_total=std::min(min_total,total_lu);
              max_total=std::max(max_total,total_lu);
            }
        }
    }
  //Cerr << FIJ << finl;
  Cerr << "Somme des facteurs de forme sur la matrice Ansys lue: Min=" << min_total << " Max=" << max_total << finl;
  if (min_total<0.9) Cerr << "Les facteurs de forme ne sont pas bien calcules. Essayer d'ameliorer la precision du calcul Ansys." << finl;
  if (min_total<0.5) Cerr << "Il se peut qu'il y'ait un probleme dans le convertisseur TRUST->Ansys. Contacter le support TRUST." << finl;

  min_total=1;
  max_total=0;
  for (I=0; I<nombre_faces_rayonnantes; I++)
    {
      double total_calcule=0;
      for (J=0; J<nombre_faces_rayonnantes; J++)
        total_calcule+=FIJ(I,J);
      min_total=std::min(min_total,total_calcule);
      max_total=std::max(max_total,total_calcule);
    }
  Cerr << "Somme des facteurs de forme apres regroupement: Min=" << min_total << " Max=" << max_total << finl;
  //Cerr << FIJ << finl;
  //
  // On normalise
  //
  Cerr << "Debut symetrisation et normalisation des facteurs de forme." << finl;
  DoubleTab Fij(FIJ);
  DoubleTab Err(FIJ);
  Fij=FIJ;
  int inter=0;
  double x=0,err=1;
  while (inter++<2000 && err>Objet_U::precision_geom)
    {
      for (i=0; i<nombre_faces_rayonnantes; i++)
        for (j=i+1; j<nombre_faces_rayonnantes; j++)
          {
            err=SI(i)*Fij(i,j)-SI(j)*Fij(j,i);
            Fij(i,j)-=(err/2./SI(i));
            Fij(j,i)+=(err/2./SI(j));
          }
      Err=Fij;
      Err-=FIJ;
      err=local_max_abs_vect(Err);
      // on arrange pour que somme=1
      //FIJ(i,j)=FIJ(i,j)/Sum(FIJ(i,j),j=1..N);
      min_total=1.;
      for (i=0; i<nombre_faces_rayonnantes; i++)
        {
          x=0;
          for (j=0; j<nombre_faces_rayonnantes; j++)
            x+=Fij(i,j);
          for (j=0; j<nombre_faces_rayonnantes; j++)
            Fij(i,j)/=x;
          min_total=std::min(x,min_total);
        }
      // On affiche apres la symetrisation...
      Cerr << "\rIteration " << inter << ": Residu="<<err<< " Min(somme normalisee des facteurs de forme)=" << min_total << "         ";
      if (inter==1)
        Cerr << finl;

      Err=Fij;
      Err-=FIJ;
      // FiJ=FIJ de l'etape precedente - FIJ de l'etape courante apres normalisattion
      // on affecte la matrice modifie a FIJ
      FIJ=Fij;
      err=local_max_abs_vect(Err);
    }
  //Cerr << FIJ << finl;
  if (min_total<0.999)
    {
      Cerr << "Erreur dans le calcul des facteurs de forme apres le regroupement." << finl;
      Cerr << "Contacter le support TRUST." << finl;
      exit();
    }
  Cerr << finl << "Ecriture du fichier " << nom2 << " contenant les facteurs de forme pour TRUST." << finl;
  SFichier ff(nom2);
  ff << nombre_faces_rayonnantes << finl;
  for (I=0; I<nombre_faces_rayonnantes; I++)
    {
      for (J=0; J<nombre_faces_rayonnantes; J++)
        ff << FIJ(I,J) << " ";
      ff << finl;
    }
  Cerr << "Ecriture du fichier " << nomS2 << " contenant les faces rayonnantes." << finl;
  return is;
}

