//
// Warning : DO NOT EDIT !
// To update this file, run: make depend
//
#include <verifie_pere.h>
#include <Echange_contact_VDF_VEF_Zoom.h>
#include <Echange_contact_VDF_Zoom_fin.h>
#include <Echange_contact_VDF_Zoom_grossier.h>
#include <Echange_contact_VEF_VDF_Zoom.h>
void instancie_src_Zoom_VDF_Cond_Lim() {
Cerr << "src_Zoom_VDF_Cond_Lim" << finl;
Echange_contact_VDF_VEF_Zoom inst1;verifie_pere(inst1);
Echange_contact_VDF_Zoom_fin inst2;verifie_pere(inst2);
Echange_contact_VDF_Zoom_grossier inst3;verifie_pere(inst3);
Echange_contact_VEF_VDF_Zoom inst4;verifie_pere(inst4);
}
