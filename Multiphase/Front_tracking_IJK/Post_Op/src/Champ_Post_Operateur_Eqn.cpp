//////////////////////////////////////////////////////////////////////////////
//
// File:        Champ_Post_Operateur_Eqn.cpp
// Directory:   $TRIO_U_ROOT/Kernel/Champs
// Version:     /main/2
//
//////////////////////////////////////////////////////////////////////////////

#include <Champ_Post_Operateur_Eqn.h>
#include <Pb_base.h>
#include <Post.h>
#include <Domaine_VF.h>
#include <Nom.h>
#include <Champ_Generique_refChamp.h>
#include <Discr_base.h>
#include <Eqn_base.h>
#include <Op_base.h>
#include <Param_lecture.h>

Implemente_instanciable_sans_constructeur(Champ_Post_Operateur_Eqn,"Champ_Post_Operateur_Eqn",Champ_Generique_Operateur_base);
   
Sortie& Champ_Post_Operateur_Eqn::printOn(Sortie& s ) const {
  return s << que_suis_je() << " " << le_nom();
}
Champ_Post_Operateur_Eqn::Champ_Post_Operateur_Eqn() { numero_op_=-1;numero_source_=-1; sans_solveur_masse_=0; }

void Champ_Post_Operateur_Eqn::set_param(Param& param)
{
 Champ_Generique_Operateur_base::set_param(param);
 param.ajouter("numero_source",&numero_source_);
 param.ajouter("numero_op",&numero_op_);
 param.ajouter_flag("sans_solveur_masse",&sans_solveur_masse_);
}
/*
Entree &   Champ_Post_Operateur_Eqn::lire(const Motcle & motcle, Entree & is)
{

  Motcles motcles(3);
  //Motcle motlu;
  motcles[0]="numero_source";
  motcles[1]="numero_op";
  motcles[2]="sans_solveur_masse";
   if (motcle == description()) {
     Cerr << motcles << finl;
     Champ_Generique_Operateur_base::lire(motcle, is);
  }
  if (motcle==motcles[0])
    {
      is>> numero_source_;
      return is;
    }
  else if (motcle==motcles[1])
    {
      is>> numero_op_;
      return is;
    } 
  else if (motcle==motcles[2])
    {
      sans_solveur_masse_=1;
      return is;
    } 
  
    return Champ_Generique_Operateur_base::lire(motcle, is);
}
*/
Entree& Champ_Post_Operateur_Eqn::readOn(Entree& s ) {
  Champ_Generique_Operateur_base::readOn(s);
  
  return s ;
}

void Champ_Post_Operateur_Eqn::completer(const Postraitement_base& post)
{
  
  Champ_Gen_de_Champs_Gen::completer(post);
  
  const Probleme_base& Pb = ref_cast(Postraitement,post).probleme();
  int numero_eq_=-1;
  if (sub_type(Champ_Generique_refChamp,get_source(0))) {
  
   Champ espace_stockage;
   const Champ_base& mon_champ = get_source(0).get_champ(espace_stockage);
    if (sub_type(Champ_Inc_base,mon_champ)) {
      const Champ_Inc_base& mon_champ_inc = ref_cast(Champ_Inc_base,mon_champ);
     
    {						   
    
     //On recupere l equation alors qu elle n est pas encore associee au Champ_Inc
     //On parcours les equatiosn du probleme et on identifie celle qui correspond au champ inconnu
     
     entier nb_eq = Pb.nombre_d_equations();
     entier i=0;
     
     while (i<nb_eq)
      {
       const Equation_base& eq_test = Pb.equation(i);
       if ((eq_test.inconnue().le_nom() == mon_champ_inc.le_nom())) 
	 { numero_eq_=i;
	   break;
	 }
	i++;
      }
    }
    }
  } 
  if (numero_eq_==-1)
  {
   Cerr<<"Can only apply a Champ_Post_Operateur_Eqn to an unknown field of the problem"<<finl;
   exit();
  }
 
   ref_eq_=Pb.equation(numero_eq_);

   int ok=0;
   const MD_Vector & md = ref_eq_.valeur().inconnue().valeurs().get_md_vector();
   const Domaine_VF& zvf= ref_cast( Domaine_VF,ref_eq_.valeur().domaine_dis().valeur());
   if (md== zvf.face_sommets().get_md_vector())
     {
       localisation_inco_=FACE;
       ok=1;
     }
   if (md== zvf.domaine().les_elems().get_md_vector())
     {
       localisation_inco_=ELEMENT;
       ok=1;
     }
   if (md == zvf.domaine().domaine().les_sommets().get_md_vector())
     {
       ok=1;
       localisation_inco_=NODE;
     }
   if (ok==0)
     {
       Cerr<<"Error in "<<que_suis_je()<<" unkonown localisation"<<finl;
       exit();
     }
     
}

const Champ_base& Champ_Post_Operateur_Eqn::get_champ(Champ & espace_stockage) const
{
  
  Champ_Fonc espace_stockage_fonc;
  Champ source_espace_stockage;
  //const Champ_base& source = get_source(0).get_champ(source_espace_stockage);
 
  //  const Champ_Inc_base& ch_inc=ref_cast(Champ_Inc_base,source);
  double temps=0.;
  Nom directive; directive=ref_eq_->inconnue().le_nom();
  // bidouille EF
  if (directive=="enthalpie") directive="temperature";
  ref_eq_.valeur().discretisation().discretiser_champ(directive,ref_eq_->domaine_dis().valeur(),"oooo","unit", -1,temps,espace_stockage_fonc);
  espace_stockage=espace_stockage_fonc;
  DoubleTab& es =(espace_stockage.valeurs());
  //if (ref_eq_->schema_temps().temps_courant()!=0)
  {
  if (numero_op_!=-1)
    {
      // certains calculer  sont faux !!!! il faudrait tous les recoder en res =0 ajouter();
      es=0;
      Operateur().calculer(ref_eq_->inconnue().valeurs(),es);
    }
  else
    ref_eq_->sources()(numero_source_).calculer(es);
  if (!sans_solveur_masse_)
   ref_eq_->solv_masse().appliquer(es);
  // Hack: car Masse_PolyMAC_Face::appliquer_impl ne divise par le volume (matrice de masse)....
  if (ref_eq_->solv_masse().valeur().que_suis_je()=="Masse_PolyMAC_Face")
    {
      //Cerr << "Volumic source terms on faces with PolyMAC can't be post-processed yet." << finl;
      Cerr << "Warning, source terms on faces with PolyMAC are post-processed as S*dV not as volumic source terms S." << finl;
      Cerr << "Cause Masse_PolyMAC_Face::appliquer_impl do not divide per volume." << finl;
      //Process::exit();
    }
  }
  // espace_stockage.valeurs()=es;
  return espace_stockage.valeur();
   

}

const Noms Champ_Post_Operateur_Eqn::get_property(const Motcle & query) const
{
 //Creation des composantes serait a faire de maniere dynamique (Eqn_...)
 
  Motcles motcles(2);  
  motcles[0] = "composantes";
  motcles[1] = "unites";
  
  entier rang = motcles.search(query);
  switch(rang) {
  case 0:
  {
  
    entier nb_comp = dimension;
    Noms compo(nb_comp);
    for (entier i=0; i<nb_comp; i++) {
     Nom nume(i);
     compo[i] = nom_post_+nume; 
    }
     
    return compo;
    break;
  }
   
  case 1:
  {
    Noms unites(dimension);
    //Noms source_unites = get_source(0).get_property("unites");
    
    for (entier i=0; i<dimension; i++) {
     unites[i] = "unit";
    }
 
    return unites;  
    break;
  }

 }
  return Champ_Gen_de_Champs_Gen::get_property(query);
}

Entity Champ_Post_Operateur_Eqn::get_localisation(const entier index) const
{

 return localisation_inco_;
 
}
//Nomme le champ en tant que source par defaut
//"Eqn_" + nom_champ_source
void Champ_Post_Operateur_Eqn::nommer_source()
{
 if (nom_post_=="??")
  {
   Nom nom_post_source, nom_champ_source;
   //const Noms nom = get_source(0).get_property("nom");
   nom_champ_source = "coucou";
   nom_post_source =  "oooooooooooEqn_";
   nom_post_source +=  nom_champ_source;
   nommer(nom_post_source);
  }
  
}
const Operateur_base& Champ_Post_Operateur_Eqn::Operateur() const
{
  return ref_eq_->operateur(numero_op_).l_op_base();
}

Operateur_base& Champ_Post_Operateur_Eqn::Operateur()
{
  return ref_eq_->operateur(numero_op_).l_op_base();
} 

