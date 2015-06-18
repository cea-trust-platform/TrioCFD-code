//
// Warning : DO NOT EDIT !
// To update this file, run: make depend
//
#include <verifie_pere.h>
#include <Champ_front_contact_rayo_semi_transp_VEF.h>
#include <Eq_rayo_semi_transp_VEF.h>
#include <Flux_radiatif_VEF.h>
#include <Neumann_paroi_rayo_semi_transp_VEF.h>
#include <Source_rayo_semi_transp_QC_VEF_P1NC.h>
#include <Source_rayo_semi_transp_VEF_P1NC.h>
#include <Temperature_imposee_paroi_rayo_semi_transp.h>
void instancie_src_Rayonnement_semi_transp_VEF() {
Cerr << "src_Rayonnement_semi_transp_VEF" << finl;
Champ_front_contact_rayo_semi_transp_VEF inst1;verifie_pere(inst1);
Eq_rayo_semi_transp_VEF inst2;verifie_pere(inst2);
Flux_radiatif_VEF inst3;verifie_pere(inst3);
Neumann_paroi_rayo_semi_transp_VEF inst4;verifie_pere(inst4);
Source_rayo_semi_transp_QC_VEF_P1NC inst5;verifie_pere(inst5);
Source_rayo_semi_transp_VEF_P1NC inst6;verifie_pere(inst6);
Temperature_imposee_paroi_rayo_semi_transp inst7;verifie_pere(inst7);
}
