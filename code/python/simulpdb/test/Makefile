PDB=1VII.pdb
SIMULPDB=python ../simulpdb.py

v5.0::
	$(SIMULPDB) $(PDB) \
	  --gmxexe=~/work/gmx/gromacs5.0/buildgcc32 \
	  --gmxsrc=~/work/gmx/gromacs5.0

v4.6::
	$(SIMULPDB) $(PDB) \
	  --gmxexe=~/work/gmx/gromacs4.6/buildgcc32 \
	  --gmxsrc=~/work/gmx/gromacs4.6

v4.5::
	$(SIMULPDB) $(PDB) \
	  --gmxexe=~/work/gmx/gromacs4.5/buildgcc32 \
	  --gmxsrc=~/work/gmx/gromacs4.5

v4.0::
	../../../gmx4.0/amberize $(PDB)
	$(SIMULPDB) am_$(PDB) \
	  --gmxexe=~/work/gmx/gromacs4.0/buildgcc32 \
	  --gmxsrc=~/work/gmx/gromacs4.0 \
	  --ff=amber99sb

clean::
	rm -rf \#* *.log *.gro *.edr *.top *.tpr *.mdp TMP* *.trr \
	  am_*.pdb *.itp *.cpt *.xtc *~ prep*
