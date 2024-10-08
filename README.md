Welcome to mcif_OUTCAR program for VASP. It converts OUTCAR and POSCAR files of a VASP calculation into a visualizable magnetic structure file i.e. mcif file. mcif file can be visualised in VESTA, Fullprof,expo etc. 
1. I need POSCAR and OUTCAR.
2. Make sure the POSCAR is in fractional coordinates.
3. Use zmcif.sh for spin-polarized only calculations
4. mcif_outcar.sh  for non-collinear calculations to view the spin moments and also enlist the spin moments at the output.
5. orbital_mcif_outcar.sh for non-collinear calculations to view the orbital moments and also enlist the total moments (J) and orbital moments (L) at the output. This should be run only after you have ran mcif_outcar.sh (this is because it needs to sum the spin and orbital moment to get the total moment)

NOTE: For non-orthogonal lattice vectors, the magnetic moments may not point exactly along the angle you would expect if it points other than only x,y,z. They need to be transformed according to the angles between the non-orthogonal lattice vectors. This update is in process.

![Alt text](example/CrI3_mag_struc_image.png)
