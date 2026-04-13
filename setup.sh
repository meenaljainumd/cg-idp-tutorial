#!/bin/bash
#set -e

peptides=("SS14" "CSS14" "PACAP27")
replicas=("with_hp" "no_hp")
posre_files=("ss14" "css14" "pacap27")
#gro_files=("8pep_1hp.gro" "8pep.gro")

for r in ${!peptides[@]}; do
    pep=${peptides[$r]}
    posre=${posre_files[$r]}

    cd $pep/
    for i in "${!replicas[@]}"; do
        rep=${replicas[$i]}
        #gro=${gro_files[$i]}
	
	#mkdir $rep/
        cd $rep/
	#cp ../../md.mdp .
	#cp ../$gro .
	#cp ../topol.top .	

        # Solvate
	#sed -i "s:${posre}_posre.itp:posre_itps/${posre}_posre.itp:g" topol.top
	#sed -i "s:heparin_posre.itp:posre_itps/heparin_posre.itp:g" topol.top
	#sed -i "s:water.md.itp:water.em.itp:g" topol.top
        #gmx solvate -cp $gro -cs ../../water_8.gro -p topol.top -o solv.gro

        # Add ions
        #gmx grompp -f ../em.mdp -c solv.gro -r solv.gro -p topol.top -o ions.tpr -maxwarn 1
        #gmx genion -s ions.tpr -o solv_neutral.gro -p topol.top -neutral

        # Energy minimization
	#gmx grompp -f ../em.mdp -c solv_neutral.gro -r solv_neutral.gro -p topol.top -o em.tpr -maxwarn 1
        #gmx mdrun -v -s em.tpr -deffnm em -tableb ../../ITP/Tables/table*

        # Equilibration
	#sed -i "s:water.em.itp:water.md.itp:g" topol.top
	#gmx make_ndx -f em.gro -o index.ndx 
        #grompp -f eq.mdp -c em.gro -n index.ndx -p topol.top -o eq.tpr -maxwarn 1

	# Production run 
        #sed -i "s:water.em.itp:water.md.itp:g" topol.top
        #grompp -f md.mdp -c eq.gro -n index.ndx -p topol.top -o run.tpr -maxwarn 1
	#mdrun -s run.tpr -g run_trial.log -x run_trial.xtc -tableb ../../forcefield/Tables/table*
        rm \#*
	rm run_trial*
	cd ../   # back to $pep/
    done
    #scp -rp with_hp/ no_hp/ mjain123@login.zaratan.umd.edu:/home/mjain123/scratch/Sims_SM/$pep 
    cd ../
done
