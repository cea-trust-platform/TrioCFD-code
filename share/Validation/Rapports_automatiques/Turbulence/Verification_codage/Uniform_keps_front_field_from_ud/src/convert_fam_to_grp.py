""" In a MED file, add groups based on existing families.
This allows avoiding the former option 'no_family_names_from_group_names' in the LireMED keyword
of TRUST
"""

import medcoupling as mc

fInName = "downcomer.med"
fOutName = "downcomer_NEW.med"
 
msh = mc.MEDFileUMesh(fInName)

famNams = set(msh.getFamiliesNames())
grpNams = set(msh.getGroupsNames())

for f in famNams:
    if f in ["FAMILLE_ZERO", "F_3D_8"]: continue
    print(f"Handling family '{f}'")
    if f in grpNams:
        msh.removeGroup(f)
        print(f"MED file {fInName} has a family named {f} but there is already a group with the same name - group has been replaced!")
    # Get cell Ids for this family
    fam = msh.getFamilyArr(-1, f, True)
    # Create corresponding group
    fam.setName(f)
    msh.addGroup(-1, fam)
    
msh.write(fOutName, 2)

