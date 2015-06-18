//
// Warning : DO NOT EDIT !
// To update this file, run: make depend
//
#include <verifie_pere.h>
#include <AssembleurPVDF_PF.h>
#include <Source_Con_Phase_field.h>
#include <Source_Gravite_PF_VDF.h>
#include <Source_Qdm_VDF_Phase_field.h>
void instancie_src_Phase_field_VDF() {
Cerr << "src_Phase_field_VDF" << finl;
AssembleurPVDF_PF inst1;verifie_pere(inst1);
Source_Con_Phase_field inst2;verifie_pere(inst2);
Source_Gravite_PF_VDF inst3;verifie_pere(inst3);
Source_Qdm_VDF_Phase_field inst4;verifie_pere(inst4);
}
