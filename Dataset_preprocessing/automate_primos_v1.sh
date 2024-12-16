#!/bin/bash
#Shell script to run PRIMOS calculations for the binding sites

while read fname
do
	#mv ${fname}.ent ${fname}.pdb
	bash primos.sh <<EOF
S
./Pocket_wormbase/${fname}
./PDB_complete_wormbase
EOF
	echo ${fname}
done <./Probes_list_v1.txt
