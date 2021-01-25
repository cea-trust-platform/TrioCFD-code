//
// Warning : DO NOT EDIT !
// To update this file, run: make depend
//
#include <verifie_pere.h>
#include <Champ_front_rayonnement.h>
#include <Equation_rayonnement.h>
#include <Flux_radiatif.h>
#include <Frontiere_ouverte_rayo_semi_transp.h>
#include <Frontiere_ouverte_temperature_imposee_rayo_semi_transp.h>
#include <Modele_rayo_semi_transp.h>
#include <Pb_Couple_rayo_semi_transp.h>
#include <Source_rayo_semi_transp.h>
void instancie_src_Rayonnement_semi_transp() {
Cerr << "src_Rayonnement_semi_transp" << finl;
Champ_front_rayonnement inst1;verifie_pere(inst1);
Equation_rayonnement inst2;verifie_pere(inst2);
Flux_radiatif inst3;verifie_pere(inst3);
Frontiere_ouverte_rayo_semi_transp inst4;verifie_pere(inst4);
Frontiere_ouverte_temperature_imposee_rayo_semi_transp inst5;verifie_pere(inst5);
Modele_rayo_semi_transp inst6;verifie_pere(inst6);
Pb_Couple_rayo_semi_transp inst7;verifie_pere(inst7);
Source_rayo_semi_transp inst8;verifie_pere(inst8);
}
