//
// Warning : DO NOT EDIT !
// To update this file, run: make depend
//
#include <verifie_pere.h>
#include <Prolongement_elem_elem_DWF.h>
#include <Prolongement_elem_elem_FMG.h>
#include <Prolongement_elem_elem_complet.h>
#include <Prolongement_elem_face.h>
#include <Prolongement_face_face.h>
#include <Prolongement_face_face_DWF.h>
#include <Prolongement_face_face_FMG.h>
#include <Prolongement_face_face_complet.h>
#include <Prolongement_face_face_lin.h>
#include <Prolongement_identite.h>
#include <Restriction_elem_elem.h>
#include <Restriction_elem_elem_2.h>
#include <Restriction_face_elem.h>
#include <Restriction_face_face.h>
#include <Restriction_face_face_vect_VDF.h>
void instancie_src_Zoom_Operateurs() {
Cerr << "src_Zoom_Operateurs" << finl;
Prolongement_elem_elem_DWF inst1;verifie_pere(inst1);
Prolongement_elem_elem_FMG inst2;verifie_pere(inst2);
Prolongement_elem_elem_complet inst3;verifie_pere(inst3);
Prolongement_elem_face inst4;verifie_pere(inst4);
Prolongement_face_face inst5;verifie_pere(inst5);
Prolongement_face_face_DWF inst6;verifie_pere(inst6);
Prolongement_face_face_FMG inst7;verifie_pere(inst7);
Prolongement_face_face_complet inst8;verifie_pere(inst8);
Prolongement_face_face_lin inst9;verifie_pere(inst9);
Prolongement_identite inst10;verifie_pere(inst10);
Restriction_elem_elem inst11;verifie_pere(inst11);
Restriction_elem_elem_2 inst12;verifie_pere(inst12);
Restriction_face_elem inst13;verifie_pere(inst13);
Restriction_face_face inst14;verifie_pere(inst14);
Restriction_face_face_vect_VDF inst15;verifie_pere(inst15);
}
