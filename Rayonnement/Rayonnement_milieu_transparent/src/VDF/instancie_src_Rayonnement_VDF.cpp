//
// Warning : DO NOT EDIT !
// To update this file, run: make depend
//
#include <verifie_pere.h>
#include <Echange_contact_Rayo_transp_VDF.h>
#include <Echange_externe_impose_rayo_transp.h>
#include <Paroi_flux_impose_Rayo_transp.h>
void instancie_src_Rayonnement_VDF() {
Cerr << "src_Rayonnement_VDF" << finl;
Echange_contact_Rayo_transp_VDF inst1;verifie_pere(inst1);
Echange_contact_Rayo_transp_sans_relax_VDF inst2;verifie_pere(inst2);
Echange_externe_impose_rayo_transp inst3;verifie_pere(inst3);
Paroi_flux_impose_Rayo_transp inst4;verifie_pere(inst4);
}
