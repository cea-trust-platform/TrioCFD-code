//
// Warning : DO NOT EDIT !
// To update this file, run: make depend
//
#include <verifie_pere.h>
#include <ConcatAnsys.h>
#include <Ecr_fic_Ansys.h>
#include <Ensemble_Faces_base.h>
#include <Face_Rayonnante.h>
#include <Frontiere_Ouverte_Rayo_transp.h>
#include <Frontiere_Ouverte_temperature_imposee_Rayo_transp.h>
#include <Modele_Rayonnement.h>
#include <Modele_Rayonnement_Milieu_Transparent.h>
#include <Pb_Couple_Rayonnement.h>
#include <VECT_Ensemble_Faces.h>
void instancie_src_Rayonnement() {
Cerr << "src_Rayonnement" << finl;
ConcatAnsys inst1;verifie_pere(inst1);
Ecr_fic_Ansys inst2;verifie_pere(inst2);
Ensemble_Faces_base inst3;verifie_pere(inst3);
Face_Rayonnante inst4;verifie_pere(inst4);
Frontiere_Ouverte_Rayo_transp inst5;verifie_pere(inst5);
Frontiere_Ouverte_temperature_imposee_Rayo_transp inst6;verifie_pere(inst6);
Modele_Rayonnement inst7;verifie_pere(inst7);
Modele_Rayonnement_Milieu_Transparent inst8;verifie_pere(inst8);
Pb_Couple_Rayonnement inst9;verifie_pere(inst9);
VECT_Ensemble_Faces inst10;verifie_pere(inst10);
}
