<hr>

<h1><p align="center"> Solid Contact Analysis Tool With Foam-Extend 
</p></h1>

<br/>
		 

## 1. Introduction

This code repository is based on the Finite Volume Method and it is forked from " https://github.com/Unofficial-Extend-Project-Mirror/foam-extend-foam-extend-4.0 ". 
The improved code library " [solidModel](src/solidModels/) " is an extension of the tools to predict deformation due to contact in solids. This 
includes the new tool " [solidGeneralContact](src/solidModels/fvPatchFields/solidContact/solidGeneralContactFvPatchVectorField.C) " for contact analysis to allow the simulation of 
multiple solid bodies (e.g. particles) in contact. It also provides the coupling of thermal contact boundary condition. Please note that this is limited to one-to-one contact 
and no multiple solid bodies are allowed currently. Furthermore, the " [elasticThermalSolidFoam](applications/solvers/solidMechanics/elasticThermalSolidFoam) " solver was improved to make it functional with the  
" [solidContact](src/solidModels/fvPatchFields/solidContact/) " library. 

<br/>
<hr>


## 2. Documentation
The documention below provides guidelines on how to install this foam-extend code repository. Please follow the link to the instruction for installation given below:

* [Installation](howToInstall.md)

Also, there is a link to the tutorial cases that have been developed for the purpose of demonstrating mechanical and thermal contact behavior of solids. Please follow the link to the example cases given below: 

* [Example Cases](tutorials/solidMechanics/deliverablesExampleCases/howToRunExampleCases.md)

<br/>
<hr>

## 4. How to cite this work?
If you are using this repository for your research, especially the library 
" [solidModel](src/solidModels/) " and the solver " [elasticThermalSolidFoam](applications/solvers/solidMechanics/elasticThermalSolidFoam) ", then please cite it as follows.

      author      = {Dhakal, Ranjan},
      title       = {Solid contact analysis tool with Foam-extend},
      year        = {2022},
      institution = {Institute of Process and Particle Engineering, Graz University of Technology},
      url         = {https://github.com/dllrun/foam-extend-general-contact.git}

<br/>
<hr>

## 5. License
This toolkit is released under the GNU General Public License (version 3). More 
details can be found in the "[LICENSE](LICENSE)" file.

<br/>
<hr>

## 6. Disclaimer
This offering is not approved or endorsed by OpenCFD Limited, producer and 
distributor of the OpenFOAM software via www.openfoam.com, and owner of the 
OPENFOAM®  and OpenCFD®  trade marks.

<br/>
<hr>

## 7. Acknowledgement
This contribution is part of the Innovative Training Network (MSCA-ITN) “MATHEGRAM” 
Project. The project has received funding from the European Union’s Horizon 
2020 research and innovation programme under the Marie Skłodowska-Curie grant 
agreement No 812638


<h1><p align="center"> <img src="mathegram_logo.png" alt="drawing" width="600"/>


<br/>
<hr>


Copyright Notice
------------------

- Copyright 2022 - Graz University of Technology (R. Dhakal, S. Radl).

- Version: foam-extend-4.0

- Licence: GPL-3.0 licence

- Last commit: June 2022
