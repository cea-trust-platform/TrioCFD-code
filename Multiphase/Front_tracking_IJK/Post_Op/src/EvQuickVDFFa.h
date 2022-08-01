//////////////////////////////////////////////////////////////////////////////
//
// File:        EvQuickVDFFa.h
// Directory:   $TRIO_U_ROOT/VDF/Operateurs/Evaluateurs
// Version:     /main/9
//
//////////////////////////////////////////////////////////////////////////////


// .NOM Eval_Quick_VDF_Face : classe instanciable des evaluateurs de convection VDF
// .ENTETE Trio_U VDF
// .LIBRAIRIE libevalvdf
// .FILE EvQuickVDFFa.h
// .FILE EvQuickVDFFa.cpp

#ifndef Eval_Quick_VDF_Face_H
#define Eval_Quick_VDF_Face_H

#include <EvAmontVDFFa.h>
#include <EvalVDFFace.h>
#include <Zone_VDF.h>
//
// .DESCRIPTION class Eval_Quick_VDF_Face
//
// Evaluateur VDF pour la convection
// Le champ convecte est un Champ_Face
// Schema de convection Quick
// Rq:Les evaluateurs de flux convectifs calculent en fait le terme 
// convectif qui figure au second membre de l'equation d'evolution 
// c.a.d l'oppose du flux convectif

//
// .SECTION voir aussi Eval_Conv_VDF


class Eval_Quick_VDF_Face : public Eval_Amont_VDF_Face
{


public:
  inline Eval_Quick_VDF_Face();
   
  inline entier calculer_arete_fluide() const ;
  inline entier calculer_arete_paroi() const ;
  inline entier calculer_arete_paroi_fluide() const ;
  inline entier calculer_arete_symetrie() const ;
  inline entier calculer_arete_interne() const ;
  inline entier calculer_arete_mixte() const ;
  inline entier calculer_fa7_sortie_libre() const ;
  inline entier calculer_arete_periodicite() const;
  inline entier calculer_arete_symetrie_paroi() const;
  inline entier calculer_arete_symetrie_fluide() const;
 
  // Fonctions qui servent a calculer le flux de grandeurs scalaires
  // Elles sont de type double et renvoient le flux     
   
  inline double flux_fa7_sortie_libre(const DoubleTab&, entier , const Neumann_sortie_libre&, entier ) const;
  inline double flux_fa7_elem(const DoubleTab&, entier, entier, entier) const ; 
  inline double flux_arete_interne(const DoubleTab&, entier, entier, entier, entier) const ;
  inline double flux_arete_mixte(const DoubleTab&, entier, entier, entier, entier) const ;
  inline double flux_arete_paroi(const DoubleTab&, entier, entier, entier, entier ) const ;   
  inline double flux_arete_symetrie(const DoubleTab&, entier, entier, entier, entier) const ;
  inline void flux_arete_fluide(const DoubleTab&, entier, entier, entier, entier,
                                double& ,double&) const ;
  inline void flux_arete_paroi_fluide(const DoubleTab&, entier, entier, entier, entier,
                                      double& ,double&) const ; 
  inline void flux_arete_periodicite(const DoubleTab&, entier, entier, entier, entier,
                                     double&, double&) const ; 
  inline void flux_arete_symetrie_fluide(const DoubleTab&, entier, entier, entier, entier, 
                                         double&, double&) const ; 
  inline double flux_arete_symetrie_paroi(const DoubleTab&, entier, entier, entier, entier) const ;

  // Fonctions qui servent a calculer le flux de grandeurs vectorielles
  // Elles sont de type void et remplissent le tableau flux 
   
  inline void flux_fa7_elem(const DoubleTab&, entier, entier, entier, DoubleVect& flux) const;
  inline void flux_fa7_sortie_libre(const DoubleTab&, entier , const Neumann_sortie_libre&,
                                    entier, DoubleVect& flux) const;  
  inline void flux_arete_interne(const DoubleTab&, entier, entier, entier,
                                 entier, DoubleVect& flux) const ;
  inline void flux_arete_mixte(const DoubleTab&, entier, entier, entier, 
                               entier, DoubleVect& flux) const ;   
  inline void flux_arete_symetrie(const DoubleTab&, entier, entier, entier, 
                                  entier, DoubleVect& flux) const ;
  inline void flux_arete_paroi(const DoubleTab&, entier, entier, entier, 
                               entier, DoubleVect& flux) const ;
  inline void flux_arete_fluide(const DoubleTab&, entier, entier, entier, 
                                entier, DoubleVect&, DoubleVect&) const;   
  inline void flux_arete_paroi_fluide(const DoubleTab&, entier, entier,
                                      entier, entier, DoubleVect& , DoubleVect&) const ;      
  inline void flux_arete_periodicite(const DoubleTab&, entier, entier,
                                     entier, entier, DoubleVect&, DoubleVect& ) const ; 
  inline void flux_arete_symetrie_fluide(const DoubleTab&, entier, entier, entier, 
                                         entier, DoubleVect&, DoubleVect& ) const;
  inline void flux_arete_symetrie_paroi(const DoubleTab&, entier, entier, entier, 
                                        entier, DoubleVect& flux) const ;     
protected:

private:
   
  inline double dist_face(entier n1,entier n2,entier k) const                { return la_zone->dist_face(n1,n2,k); };
  inline double dist_face_period(entier n1,entier n2,entier k) const         { return la_zone->dist_face_period(n1,n2,k); };
  inline double dim_face(entier ,entier ) const;
  inline entier face_amont_conj(entier ,entier ,entier ) const;
  inline entier face_amont_princ(entier ,entier ) const;
  inline double dist_elem(entier n1,entier n2,entier k) const                { return la_zone->dist_elem(n1,n2,k); };
  inline double dist_elem_period(entier n1, entier n2, entier k) const        { return la_zone->dist_elem_period(n1,n2,k); };
  inline double dim_elem(entier ,entier ) const;          

};

//
// Fonctions inline de la classe Eval_Quick_VDF_Face
//

inline Eval_Quick_VDF_Face::Eval_Quick_VDF_Face()
{
}
   
 
//// calculer_arete_fluide
//

inline entier Eval_Quick_VDF_Face::calculer_arete_fluide() const
{
  return 1;
}


//// calculer_arete_paroi
//

inline entier Eval_Quick_VDF_Face::calculer_arete_paroi() const
{
  return 0;
}


//// calculer_arete_paroi_fluide
//

inline entier Eval_Quick_VDF_Face::calculer_arete_paroi_fluide() const
{
  return 1;
}


//// calculer_arete_symetrie
//

inline entier Eval_Quick_VDF_Face::calculer_arete_symetrie() const
{
  return 0;
}

//// calculer_arete_interne
//

inline entier Eval_Quick_VDF_Face::calculer_arete_interne() const 
{
  return 1;
}

//// calculer_arete_mixte
//

inline entier Eval_Quick_VDF_Face::calculer_arete_mixte() const 
{
  return 1;
}

//// calculer_fa7_sortie_libre
//

inline entier Eval_Quick_VDF_Face::calculer_fa7_sortie_libre() const 
{
  return 1;
}

//// calculer_arete_periodicite
//

inline entier Eval_Quick_VDF_Face::calculer_arete_periodicite() const 
{
  return 1;
}

//// calculer_arete_symetrie_paroi
//

inline entier Eval_Quick_VDF_Face::calculer_arete_symetrie_paroi() const
{
  return 0;
}

//// calculer_arete_symetrie_fluide
//

inline entier Eval_Quick_VDF_Face::calculer_arete_symetrie_fluide() const
{
  return 1;
}  

// Fonction de calcul de cf(limiteur de pente) dans le schema Quick-sharp

inline double sharp(const double& utc) 
{
  double cf;
  if ( (utc <= -1) || (utc >= 1.5) )
    cf = 0.125;  
  else if ((utc > -1) && (utc <= 0) )
    cf = 0.5 + 0.375*utc;
  else if ((utc <= 0.25) && (utc > 0) )
    cf = 0.5 - 0.625*sqrt(utc);
  else if ((utc > 0.25) && (utc <= 1.) )
    cf = 0.25*(1.-utc);
  else
    cf = 0.25*(utc-1.);
  return cf;
}

inline double conv_quick_sharp_plus(const double& psc,const double& vit_0, const double& vit_1,
                                    const double& vit_0_0, const double& dx,
                                    const double& dm, const double& dxam)
{ 
  double cf;
  double curv;
  double delta_0 = vit_0 - vit_0_0;
  double delta = vit_1 - vit_0;
  double dd1,utc;
  double delta_delta;
   
  curv = (delta/dx - delta_0/dxam)/dm ;
   
  // Calcul de cf:
   
  delta_delta = delta_0+delta;
  dd1 = dabs(delta_delta);
  if (dd1 < 1.e-5)
    cf = 0.125;
  else {
    utc = delta_0/delta_delta;
    cf = sharp(utc);
  }
   
  return (0.5*(vit_0 + vit_1) - cf*(dx*dx)*curv)*psc;
   
}

inline double conv_quick_sharp_moins(const double& psc,const double& vit_0,const double& vit_1,
                                     const double& vit_1_1,const double& dx,
                                     const double& dm,const double& dxam)
{
  double cf;
  double curv;
  double delta_1 = vit_1_1 - vit_1;
  double delta = vit_1 - vit_0;
  double dd1,utc;
  double delta_delta;
   
  curv = ( delta_1/dxam - delta/dx )/dm ;
   
  // Calcul de cf:
   
  delta_delta = delta_1+delta;
  dd1 = dabs(delta_delta);
  if (dd1 < 1.e-5)
    cf = 0.125;
  else {
    utc = delta_1/delta_delta;
    cf = sharp(utc);
  }
   
  return (0.5*(vit_0 + vit_1) - cf*(dx*dx)*curv)*psc;

}

inline double Eval_Quick_VDF_Face::dim_elem(entier n1, entier k) const
{
  return la_zone->dim_elem(n1,k);
}

inline double Eval_Quick_VDF_Face::dim_face(entier n1, entier k) const
{
  return la_zone->dim_face(n1,k);
}

inline entier Eval_Quick_VDF_Face::face_amont_conj(entier num_face, entier i, entier k) const
{
  return la_zone->face_amont_conj(num_face, i, k);
}

inline entier Eval_Quick_VDF_Face::face_amont_princ(entier num_face, entier i) const
{
  return la_zone->face_amont_princ(num_face, i);
}

// Fonctions de calcul des flux pour une inconnue scalaire

//// flux_arete_fluide
//

inline void Eval_Quick_VDF_Face::flux_arete_fluide(const DoubleTab& inco, 
                                                   entier fac1, entier fac2, 
                                                   entier fac3, entier signe,
                                                   double& flux3,double& flux1_2) const
{
  double flux,psc;
  entier ori;
  psc = 0.25 * ((dt_vitesse[fac1]*porosite[fac1]    + 
                 dt_vitesse[fac2]*porosite[fac2] )  * 
                (surface(fac1) + surface(fac2))  );     
   
  if ((psc*signe)>0) 
    flux=inco[fac3]*psc ;
  else {
    ori =orientation(fac3);
    flux = 0.5*(inconnue->val_imp_face_bord(fac1,ori) 
                +  inconnue->val_imp_face_bord(fac2,ori)) * psc ;
  }        
  flux3 = -flux;
   
  // Calcul de flux1_2:
   
  psc = 0.5*dt_vitesse(fac3)*surface(fac3)*porosite(fac3);
  if (psc>0)
    flux = psc*inco[fac1];
  else
    flux = psc*inco[fac2];
   
  flux1_2 = -flux;
}


//// flux_arete_interne
//

inline double Eval_Quick_VDF_Face::flux_arete_interne(const DoubleTab& inco, entier fac1, 
                                                      entier fac2, entier fac3, entier fac4) const
{
  double flux;
  double psc = 0.25* (dt_vitesse[fac1]*porosite(fac1) + dt_vitesse[fac2]*porosite(fac2))
    * (surface(fac1) + surface(fac2));   
  entier ori = orientation(fac1);
  entier num0_0 = face_amont_conj(fac3,ori,0);
  entier num1_1 = face_amont_conj(fac4,ori,1);
   
  if (psc > 0) 
    if (num0_0 == -1) // Schema amont
      flux = psc*inco[fac3];
    else // Schema quick
      {
        double dx = dist_face_period(fac3,fac4,ori);
        double dm = dim_face(fac3,ori);
        double dxam = dist_face_period(num0_0,fac3,ori);
        double vit_0 = inco[fac3];
        double vit_0_0 = inco[num0_0];
        double vit_1 = inco[fac4];            
	//		Cerr<<dx <<" tailles1 "<<dm<<" "<<dxam<<finl;
        flux = conv_quick_sharp_plus(psc,vit_0,vit_1,vit_0_0,dx,dm,dxam);
      } 
      
  else // (psc <= 0)
    if (num1_1 == -1) // Schema amont
      flux = psc*inco[fac4];
    else // Schema quick
      {
        double dx = dist_face_period(fac3,fac4,ori);

        double dm = dim_face(fac4,ori);
        double dxam = dist_face_period(fac4,num1_1,ori);   
	//		Cerr<<dx <<" tailles "<<dm<<" "<<dxam<<finl;
        double vit_0 = inco[fac3];
        double vit_1 = inco[fac4];
        double vit_1_1 = inco[num1_1];            
        flux = conv_quick_sharp_moins(psc,vit_0,vit_1,vit_1_1,dx,dm,dxam);
      } 

  return -flux;
}


//// flux_arete_mixte
//

inline double Eval_Quick_VDF_Face::flux_arete_mixte(const DoubleTab& inco , 
                                                    entier fac1, entier fac2, 
                                                    entier fac3, entier fac4) const
{
  double flux;
  double psc = 0.25 * ( dt_vitesse[fac1]*porosite[fac1] + 
                        dt_vitesse[fac2]*porosite[fac2] )  * 
    ( surface[fac1] + surface[fac2] );
  if (psc>0)
    flux = psc*inco[fac3];
  else
    flux = psc*inco[fac4];
 
  return -flux;
}


//// flux_arete_paroi
//

inline double Eval_Quick_VDF_Face::flux_arete_paroi(const DoubleTab&, entier, entier, entier, entier) const
{
  return 0;
}


//// flux_arete_paroi_fluide
//

inline void Eval_Quick_VDF_Face::flux_arete_paroi_fluide(const DoubleTab& inco, 
                                                         entier fac1, entier fac2, 
                                                         entier fac3, entier signe,
                                                         double& flux3, double& flux1_2) const
{
  double flux,psc;
  entier ori;
  psc = 0.25 * ((dt_vitesse[fac1]*porosite[fac1]    + 
                 dt_vitesse[fac2]*porosite[fac2] )  * 
                (surface(fac1) + surface(fac2))  );     
   
  if ((psc*signe)>0) 
    flux=inco[fac3]*psc ;
  else {
    ori =orientation(fac3);
    flux = 0.5*(inconnue->val_imp_face_bord(fac1-premiere_face_bord, ori) 
                +  inconnue->val_imp_face_bord(fac2-premiere_face_bord, ori))* psc ;
  }
   
  flux3 = -flux;        
   
  // Calcul de flux1_2:
   
  psc = 0.5*dt_vitesse(fac3)*surface(fac3)*porosite(fac3);
  if (psc>0)
    flux = psc*inco[fac1];
  else
    flux = psc*inco[fac2];
   
  flux1_2 = -flux;
}

//// flux_arete_periodicite
//

void Eval_Quick_VDF_Face::flux_arete_periodicite(const DoubleTab& inco, 
                                                 entier fac1, entier fac2, 
                                                 entier fac3, entier fac4,
                                                 double& flux3_4 , double& flux1_2) const
{
  flux3_4=flux_arete_interne(inco,fac1, fac2,fac3, fac4) ;
  flux1_2=flux_arete_interne(inco,fac3, fac4,fac1, fac2) ;
  return;
  double flux;
  double psc;
   
  // On calcule le flux convectif entre les volumes de controle associes a fac3 et fac4:
   
  psc = 0.25 * ( dt_vitesse[fac1]*porosite[fac1] + 
                 dt_vitesse[fac2]*porosite[fac2] ) *
    ( surface[fac1] + surface[fac2] );
   
  if (psc>0)
    flux = psc*inco[fac3];
  else
    flux = psc*inco[fac4];
   
  flux3_4 = -flux;
   
  // On calcule le flux convectif entre les volumes de controle associes a fac1 et fac2:
   
  psc = 0.25 * ( dt_vitesse[fac3]*porosite[fac3] + 
                 dt_vitesse[fac4]*porosite[fac4] ) *
    ( surface[fac3] + surface[fac4] );
   
  if (psc>0)
    flux = psc*inco[fac1];
  else
    flux = psc*inco[fac2];
   
  flux1_2 = -flux;
   
}

//// flux_arete_symetrie
//

inline double Eval_Quick_VDF_Face::flux_arete_symetrie(const DoubleTab&, entier, entier, entier, entier) const
{
  return 0;
}

//// flux_arete_symetrie_fluide
//

inline void Eval_Quick_VDF_Face::flux_arete_symetrie_fluide(const DoubleTab& inco, 
                                                            entier fac1, entier fac2, 
                                                            entier fac3, entier signe,
                                                            double& flux3,double& flux1_2) const
{
  double flux,psc;
  entier ori;
  psc = 0.25 * ((dt_vitesse[fac1]*porosite[fac1]    + 
                 dt_vitesse[fac2]*porosite[fac2] )  * 
                (surface(fac1) + surface(fac2))  );     
   
  if ((psc*signe)>0) 
    flux=inco[fac3]*psc ;
  else {
    ori =orientation(fac3);
    flux = 0.5*(inconnue->val_imp_face_bord(fac1,ori) 
                +  inconnue->val_imp_face_bord(fac2,ori)) * psc ;
  }        
  flux3 = -flux;
   
  // Calcul de flux1_2:
   
  psc = 0.5*dt_vitesse(fac3)*surface(fac3)*porosite(fac3);
  if (psc>0)
    flux = psc*inco[fac1];
  else
    flux = psc*inco[fac2];
   
  flux1_2 = -flux;
}

//// flux_arete_symetrie_paroi
//

inline double Eval_Quick_VDF_Face::flux_arete_symetrie_paroi(const DoubleTab&, entier, entier, entier, entier) const
{
  return 0;
}


//// flux_fa7_elem
//

inline double Eval_Quick_VDF_Face::flux_fa7_elem(const DoubleTab& inco, entier num_elem, 
                                                 entier fac1, entier fac2) const
{
  double flux;
  entier num0_0 = face_amont_princ(fac1,0);
  entier num1_1 = face_amont_princ(fac2,1);
  double psc = 0.25 * (dt_vitesse[fac1]+dt_vitesse[fac2])*(surface(fac1)+surface(fac2));
  if  (psc > 0) 
    {
      if (num0_0 == -1) // Schema amont
        flux = psc*inco[fac1]*porosite(fac1);
      else {
        entier ori = orientation(fac1);
        double dx = dim_elem(num_elem,ori);
        entier elem_amont = elem(fac1,0);     
        double dm = dist_elem_period(elem_amont,num_elem,ori);
        double dxam = dim_elem(elem_amont,ori);                         
        double vit_0 = inco[fac1]*porosite(fac1);
        double vit_0_0 = inco[num0_0]*porosite(num0_0);
        double vit_1 = inco[fac2]*porosite(fac2);
         
        flux = conv_quick_sharp_plus(psc,vit_0,vit_1,vit_0_0,dx,dm,dxam);
      }
    }
  else // (psc < 0)
    {
      if (num1_1 == -1) // Schema amont 
        flux = psc*inco[fac2]*porosite(fac2);
      else {
        entier ori = orientation(fac2);
        double dx = dim_elem(num_elem,ori);
        entier elem_amont = elem(fac2,1);
        double dm = dist_elem_period(num_elem,elem_amont,ori);
        double dxam = dim_elem(elem_amont,ori);
        double vit_0 = inco[fac1]*porosite(fac1);
        double vit_1 = inco[fac2]*porosite(fac2);
        double vit_1_1 = inco[num1_1]*porosite(num1_1);
         
        flux = conv_quick_sharp_moins(psc,vit_0,vit_1,vit_1_1,dx,dm,dxam);
      }
    }
   
  return -flux;
}


inline double Eval_Quick_VDF_Face::flux_fa7_sortie_libre(const DoubleTab& inco, entier face, 
                                                         const Neumann_sortie_libre& la_cl, 
                                                         entier num1) const
{
  double flux;
  double psc = dt_vitesse[face]*surface(face);
  entier elem1 = elem(face,0);
  if (elem1 != -1)
    if (psc > 0) 
      flux = psc*inco[face];
    else 
      flux = psc*la_cl.val_ext(face-num1,orientation(face));
  else // (elem2 != -1)
    if (psc < 0) 
      flux = psc*inco[face];
    else 
      flux = psc*la_cl.val_ext(face-num1,orientation(face));
   
  return -flux;
}


// Fonctions de calcul des flux pour une inconnue vectorielle

//// flux_arete_fluide
//

inline void Eval_Quick_VDF_Face::flux_arete_fluide(const DoubleTab& inco,entier fac1,
                                                   entier fac2,entier fac3,entier signe,
                                                   DoubleVect& flux3, DoubleVect& flux1_2) const
{
  entier k;
  double psc = 0.25 * ((dt_vitesse[fac1]*porosite[fac1]    + 
                        dt_vitesse[fac2]*porosite[fac2] )  * 
                       (surface(fac1) + surface(fac2))  );     
   
  if ((psc*signe)>0) 
    for (k=0; k<flux3.size(); k++)
      flux3(k) = -inco(fac3,k)*psc ;
  else {
    entier ori = orientation(fac3);
    for (k=0; k<flux3.size(); k++)
      flux3(k) = -0.5*(inconnue->val_imp_face_bord(fac1,k,ori) 
                       +  inconnue->val_imp_face_bord(fac2,k,ori)) * psc ;
  }
 
  // Calcul de flux1_2:
   
  psc = 0.5*dt_vitesse(fac3)*surface(fac3)*porosite(fac3);
  if (psc>0)
    for (k=0; k<flux1_2.size(); k++)
      flux1_2(k) = -psc*inco(fac1,k);
  else
    for (k=0; k<flux1_2.size(); k++)
      flux1_2(k) = -psc*inco(fac2,k);
   
}


//// flux_arete_interne
//

inline void Eval_Quick_VDF_Face::flux_arete_interne(const DoubleTab& inco, entier fac1, 
                                                    entier fac2, entier fac3, entier fac4,
                                                    DoubleVect& flux) const
{
  entier k;
  entier ncomp = flux.size();
  double psc = 0.25* (dt_vitesse[fac1]*porosite(fac1) + dt_vitesse[fac2]*porosite(fac2))
    * (surface(fac1) + surface(fac2));   
  entier ori = orientation(fac1);
  entier num0_0 = face_amont_conj(fac3,ori,0);
  entier num1_1 = face_amont_conj(fac4,ori,1);
   
  if (psc > 0) 
    if (num0_0 == -1) // Schema amont
      for (k=0; k<ncomp; k++)
        flux(k) = -psc*inco(fac3,k);
    else // Schema quick
      {
        double dx = dist_face(fac3,fac4,ori);
        double dm = dim_face(fac3,ori);
        double dxam = dist_face_period(num0_0,fac3,ori);
        ArrOfDouble vit_0(ncomp);
        ArrOfDouble vit_0_0(ncomp);
        ArrOfDouble vit_1(ncomp);
        for (k=0; k<ncomp; k++) {
          vit_0(k) = inco(fac3,k);
          vit_0_0(k) = inco(num0_0,k);
          vit_1(k) = inco(fac4,k);
          flux(k) = -conv_quick_sharp_plus(psc,vit_0(k),vit_1(k),vit_0_0(k),dx,dm,dxam);
        }
      } 
   
  else // (psc <= 0)
    if (num1_1 == -1) // Schema amont
      for (k=0; k<ncomp; k++)
        flux(k) = -psc*inco(fac4,k);
    else // Schema quick
      {
        double dx = dist_face(fac3,fac4,ori);
        double dm = dim_face(fac4,ori);
        double dxam = dist_face_period(fac4,num1_1,ori);
        ArrOfDouble vit_0(ncomp);
        ArrOfDouble vit_1(ncomp);
        ArrOfDouble vit_1_1(ncomp);
        for (k=0; k<ncomp; k++) {
          vit_0(k) = inco(fac3,k);
          vit_1(k) = inco(fac4,k);
          vit_1_1(k) = inco(num1_1,k);            
          flux(k) = -conv_quick_sharp_moins(psc,vit_0(k),vit_1(k),vit_1_1(k),dx,dm,dxam);
        }
      } 
}


//// flux_arete_mixte
//

inline void Eval_Quick_VDF_Face::flux_arete_mixte(const DoubleTab& inco,entier fac1,
                                                  entier fac2,entier fac3, entier fac4,
                                                  DoubleVect& flux) const
{
  entier k;
  double psc = 0.25 * ( dt_vitesse[fac1]*porosite[fac1] + 
                        dt_vitesse[fac2]*porosite[fac2] )  * 
    ( surface[fac1] + surface[fac2] );
  if (psc>0)
    for (k=0; k<flux.size(); k++)
      flux(k) = -psc*inco(fac3,k);
  else
    for (k=0; k<flux.size(); k++)
      flux(k) = -psc*inco(fac4,k);
}


//// flux_arete_paroi
//

inline void Eval_Quick_VDF_Face::flux_arete_paroi(const DoubleTab&, entier, entier,
                                                  entier, entier,  DoubleVect& ) const
{
  ;
}


//// flux_arete_paroi_fluide
//

inline void Eval_Quick_VDF_Face::flux_arete_paroi_fluide(const DoubleTab& inco,entier fac1,
                                                         entier fac2,entier fac3,entier signe,
                                                         DoubleVect& flux3, DoubleVect& flux1_2) const
{
  entier k;
  double psc = 0.25 * ((dt_vitesse[fac1]*porosite[fac1]    + 
                        dt_vitesse[fac2]*porosite[fac2] )  * 
                       (surface(fac1) + surface(fac2))  );     
   
  if ((psc*signe)>0) 
    for (k=0; k<flux3.size(); k++)
      flux3(k) = -inco(fac3,k)*psc;
  else {
    entier ori = orientation(fac3);
    for (k=0; k<flux3.size(); k++)
      flux3(k) = -0.5*(inconnue->val_imp_face_bord(fac1-premiere_face_bord,k,ori) 
                       +  inconnue->val_imp_face_bord(fac2-premiere_face_bord,k,ori))* psc ;        
  }
   
  // Calcul de flux1_2:
   
  psc = 0.5*dt_vitesse(fac3)*surface(fac3)*porosite(fac3);
  if (psc>0)
    for (k=0; k<flux1_2.size(); k++)
      flux1_2(k) = -psc*inco(fac1,k);
  else
    for (k=0; k<flux1_2.size(); k++)
      flux1_2(k) = -psc*inco(fac2,k);
   
}

//// flux_arete_periodicite
//

inline void Eval_Quick_VDF_Face::flux_arete_periodicite(const DoubleTab& inco,
                                                        entier fac1, entier fac2 , entier fac3, entier fac4, 
                                                        DoubleVect& flux3_4, DoubleVect& flux1_2) const
{
  double psc;
  entier k;
   
  // On calcule le flux convectif entre les volumes de controle associes a fac3 et fac4:
   
  psc = 0.25 * ( dt_vitesse[fac1]*porosite[fac1] + 
                 dt_vitesse[fac2]*porosite[fac2] ) *
    ( surface[fac1] + surface[fac2] );
   
  if (psc>0)
    for (k=0; k<flux3_4.size(); k++)
      flux3_4(k) = -psc*inco(fac3,k);
  else
    for (k=0; k<flux3_4.size(); k++)
      flux3_4(k) = -psc*inco(fac4,k);
   
   
  // On calcule le flux convectifs entre les volumes de controle associes a fac1 et fac2:
   
  psc = 0.25 * ( dt_vitesse[fac3]*porosite[fac3] + 
                 dt_vitesse[fac4]*porosite[fac4] ) *
    ( surface[fac3] + surface[fac4] );
   
  if (psc>0)
    for (k=0; k<flux1_2.size(); k++)
      flux1_2(k) = -psc*inco(fac1,k);
  else
    for (k=0; k<flux3_4.size(); k++)
      flux1_2(k) = -psc*inco(fac2,k);
}   

//// flux_arete_symetrie
//

inline void Eval_Quick_VDF_Face::flux_arete_symetrie(const DoubleTab&, entier, entier,
                                                     entier, entier, DoubleVect& ) const
{
  ;
}

//// flux_arete_symetrie_fluide
//

inline void Eval_Quick_VDF_Face::flux_arete_symetrie_fluide(const DoubleTab& inco,entier fac1,
                                                            entier fac2,entier fac3,entier signe,
                                                            DoubleVect& flux3, DoubleVect& flux1_2) const
{
  entier k;
  double psc = 0.25 * ((dt_vitesse[fac1]*porosite[fac1]    + 
                        dt_vitesse[fac2]*porosite[fac2] )  * 
                       (surface(fac1) + surface(fac2))  );     
   
  if ((psc*signe)>0) 
    for (k=0; k<flux3.size(); k++)
      flux3(k) = -inco(fac3,k)*psc ;
  else {
    entier ori = orientation(fac3);
    for (k=0; k<flux3.size(); k++)
      flux3(k) = -0.5*(inconnue->val_imp_face_bord(fac1,k,ori) 
                       +  inconnue->val_imp_face_bord(fac2,k,ori)) * psc ;
  }
 
  // Calcul de flux1_2:
   
  psc = 0.5*dt_vitesse(fac3)*surface(fac3)*porosite(fac3);
  if (psc>0)
    for (k=0; k<flux1_2.size(); k++)
      flux1_2(k) = -psc*inco(fac1,k);
  else
    for (k=0; k<flux1_2.size(); k++)
      flux1_2(k) = -psc*inco(fac2,k);
   
}

//// flux_arete_symetrie_paroi
//

inline void Eval_Quick_VDF_Face::flux_arete_symetrie_paroi(const DoubleTab&, entier, entier, entier, 
                                                           entier, DoubleVect& ) const 
{
  ;
}

//// flux_fa7_elem
//

inline void Eval_Quick_VDF_Face::flux_fa7_elem(const DoubleTab& inco,entier num_elem, 
                                               entier fac1,entier fac2,DoubleVect& flux) const
{
  entier k;
  entier num0_0 = face_amont_princ(fac1,0);
  entier num1_1 = face_amont_princ(fac2,1);
  entier ncomp = flux.size();
  double psc = 0.25 * (dt_vitesse[fac1]+dt_vitesse[fac2])*(surface(fac1)+surface(fac2));
  if  (psc > 0) 
    {
      if (num0_0 == -1) // Schema amont
        for (k=0; k<ncomp; k++)
          flux(k) = -psc*inco(fac1,k)*porosite(fac1);
      else {
        entier ori = orientation(fac1);
        double dx = dim_elem(num_elem,ori);
        entier elem_amont = elem(fac1,0);     
        double dm = dist_elem_period(elem_amont,num_elem,ori);
        double dxam = dim_elem(elem_amont,ori);
        ArrOfDouble vit_0(ncomp);
        ArrOfDouble vit_0_0(ncomp);
        ArrOfDouble vit_1(ncomp);
        for (k=0; k<ncomp; k++) {
          vit_0(k) = inco(fac1,k)*porosite(fac1);
          vit_0_0(k) = inco(num0_0,k)*porosite(num0_0);
          vit_1(k) = inco(fac2,k)*porosite(fac2);
          flux(k) = -conv_quick_sharp_plus(psc,vit_0(k),vit_1(k),vit_0_0(k),dx,dm,dxam);            
        }
      }
    }
  else // (psc < 0)
    {
      if (num1_1 == -1) // Schema amont 
        for (k=0; k<ncomp; k++)
          flux(k) = -psc*inco(fac2,k)*porosite(fac2);
      else {
        entier ori = orientation(fac2);
        double dx = dim_elem(num_elem,ori);
        entier elem_amont = elem(fac2,1);
        double dm = dist_elem_period(num_elem,elem_amont,ori);
        double dxam = dim_elem(elem_amont,ori);
        ArrOfDouble vit_0(ncomp);
        ArrOfDouble vit_1(ncomp);
        ArrOfDouble vit_1_1(ncomp);         
        for (k=0; k<ncomp; k++) {
          vit_0(k) = inco(fac1,k)*porosite(fac1);
          vit_1(k) = inco(fac2,k)*porosite(fac2);
          vit_1_1(k) = inco(num1_1,k)*porosite(num1_1);
          flux(k) = -conv_quick_sharp_moins(psc,vit_0(k),vit_1(k),vit_1_1(k),dx,dm,dxam);
        }
      }
    }
}


inline void Eval_Quick_VDF_Face::flux_fa7_sortie_libre(const DoubleTab& inco, entier face, 
                                                       const Neumann_sortie_libre& la_cl, 
                                                       entier num1,DoubleVect& flux) const
{
  entier k;
  double psc = dt_vitesse[face]*surface(face);
  entier elem1 = elem(face,0);
  if (elem1 != -1)
    if (psc > 0) 
      for (k=0; k<flux.size(); k++)
        flux(k) = -psc*inco(face,k);
    else 
      for (k=0; k<flux.size(); k++)
        flux(k) = -psc*la_cl.val_ext(face-num1,k);
  else // (elem2 != -1)
    if (psc < 0) 
      for (k=0; k<flux.size(); k++)
        flux(k) = -psc*inco(face,k);
    else 
      for (k=0; k<flux.size(); k++)
        flux(k) = -psc*la_cl.val_ext(face-num1,k);
}

#endif
