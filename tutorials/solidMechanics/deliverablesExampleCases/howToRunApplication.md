<hr>

<h1><p align="center"> How to run application  
</p></h1>	



#### 1. Spherical particles in a 2D compaction:
   This is a setup to demonstrate elastic contact between multiple spherical 
particles. This has been illustrated with the case 
"1_sphericalParticlesIn2DCompaction".

i. Meshing
 - Before running this application, we need to first import the mesh from the case "23ParticlesMesh". 
This mesh case runs with the OpenFOAM-5.0.
 - The input for this mesh case can be handled from the "system/blocKMeshDict" file.
 - Clean and run the case with "`./Allclean`" and "`./Allrun`" commmand correspondingly.
 - Now the mesh can be visualized in paraview by loading the dummy file m.foam 
(from the case directory).

 

ii. Simulation
 - Now, go back to the case ("1_sphericalParticlesIn2DCompaction") directory 
and make sure the path to the mesh case directory is correct in the copyMesh.sh 
 - Check if the "solidGeneralContact" boundary condition is specified at the 
correct boundary patches in the "0/U" (here it is done with contactFieldSub) file.
 - The boundary condition representing the motion of the piston (i.e for the 
patch "topPunch" is mentioned in the "0/timeVaryingFixedSub") can be  specified 
with the "constant/timeVsTopDisp" file.   
 - The application will be running with the "elasticSolidFoam" solver and the 
mechanical properties can be specified in the "constant/rheologyProperties" 
file.
 - From the case directory clean the case with "`./Allclean`" commmand and then 
run the case with "`./Allrun`" commmand

<br/>
<hr>  

#### 2. Thermal contact between 2 particles with different temperatures: 
   This is a setup to demonstrate thermal contact between two hexahedral 
particles. This has been illustrated with the case 
"2_thermalContactOf2HexaParticles".  
 - Check if the "solid4Contact" boundary condition is specified in the 0/U file 
and "thermalContact" boundary condition in 0/T file for the correct boundary 
patches.    
 - The application will be running with the "elasticThermalSolidFoam" solver 
and the mechanical properties can be specified in the 
"constant/rheologyProperties" file while the thermal properties can be 
specified in the "constant/thermalProperties" file.
 - From the case directory clean the case with "./Allclean" commmand and then 
run the case with "`./Allrun`" commmand

<br/>
<hr>  


Author
------------------

- Ranjan Dhakal 
- IPPE, Graz University of Technology 
- June 2022
