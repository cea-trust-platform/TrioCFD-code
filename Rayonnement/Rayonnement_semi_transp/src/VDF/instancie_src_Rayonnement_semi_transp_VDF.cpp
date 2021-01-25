//
// Warning : DO NOT EDIT !
// To update this file, run: make depend
//
#include <verifie_pere.h>
#include <Echange_contact_rayo_semi_transp_VDF.h>
#include <Echange_externe_impose_rayo_semi_transp.h>
#include <Echange_global_impose_rayo_semi_transp.h>
#include <Eq_rayo_semi_transp_VDF.h>
#include <Flux_radiatif_VDF.h>
#include <Neumann_paroi_rayo_semi_transp_VDF.h>
#include <Source_rayo_semi_transp_QC_VDF_P0_VDF.h>
#include <Source_rayo_semi_transp_VDF_P0_VDF.h>
void instancie_src_Rayonnement_semi_transp_VDF() {
Cerr << "src_Rayonnement_semi_transp_VDF" << finl;
Echange_contact_rayo_semi_transp_VDF inst1;verifie_pere(inst1);
Echange_externe_impose_rayo_semi_transp inst2;verifie_pere(inst2);
Echange_global_impose_rayo_semi_transp inst3;verifie_pere(inst3);
Eq_rayo_semi_transp_VDF inst4;verifie_pere(inst4);
Flux_radiatif_VDF inst5;verifie_pere(inst5);
Neumann_paroi_rayo_semi_transp_VDF inst6;verifie_pere(inst6);
Source_rayo_semi_transp_QC_VDF_P0_VDF inst7;verifie_pere(inst7);
Source_rayo_semi_transp_VDF_P0_VDF inst8;verifie_pere(inst8);
}
