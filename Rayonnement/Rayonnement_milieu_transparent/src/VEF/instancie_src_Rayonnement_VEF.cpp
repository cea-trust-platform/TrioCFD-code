//
// Warning : DO NOT EDIT !
// To update this file, run: make depend
//
#include <verifie_pere.h>
#include <Champ_front_contact_rayo_transp_VEF.h>
#include <Temperature_imposee_paroi_rayo_transp.h>
void instancie_src_Rayonnement_VEF() {
Cerr << "src_Rayonnement_VEF" << finl;
Champ_front_contact_rayo_transp_VEF inst1;verifie_pere(inst1);
Temperature_imposee_paroi_rayo_transp inst2;verifie_pere(inst2);
}
