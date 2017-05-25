#! /bin/sh
Loc=$PWD
## run mask
/hdata/fou/NorgeIsModelling/Code/pre_3.01c/mask $Loc/Control/control_mask.txt
## run pre
/hdata/fou/NorgeIsModelling/Code/pre_3.01c/predew $Loc/Control/control_pre.txt
mv $Loc/dew_landscape.txt $Loc/During/
mv $Loc/dew_grid_index.txt $Loc/NoVet/
mv $Loc/dew_waterland.txt $Loc/During/
