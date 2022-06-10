<hr>

<h1><p align="center"> How to install this Code Repository
</p></h1>	           



#### 1. Before installation: 
- Get the source code
  - Go to home directory in your terminal using,
  - `cd ~`
- Make a directory and download the Git repository copy into this directory,
  - `mkdir foam`
  - `cd foam`
  - `git clone https://github.com/dllrun/foam-extend-general-contact.git `

<br/>
<hr>  


#### 2. (Install missing packages if any)
#### 3. Adjust the `~/.bashrc` file: 
- In the home directory, type 
  - `vim .bashrc`
- (Set the environment variables if any)
- Set the path correctly
  - `alias switchToExtend='source $HOME/foam/foam-extend-4.0/etc/bashrc'`

#### 4. Source the .bashrc file
- Type,
  - `source .bashrc`
  - `switchToExtend`
  
<br/>
<hr>

#### 5. Start Compilation
- Go to installation directory,
  - `cd foam/foam-extend-4.0/`
  - `./Allwmake`
XXX detail what happens here

<br/>
<hr>

OpenFOAM-5.0:
---------------- 
Since, foam-extend-4.0 does not support the mesh generation for spherical 
geometry, we need to import the mesh from OpenFOAM-5.0. Therefore, it is 
recommended to have a OpenFOAM-5.0 installed.
- For the official release, details can be found in this page: 
https://openfoam.org/release/5-0/ 
- If you are using an OpenSUSE Leap with Windows 10, a further detail on 
downloading the precomplied version can be found here: 
https://wiki.openfoam.com/Win10OpenSUSEShell_by_Stefan_Radl 

XXX mention that it might also work with other versions of OpenFOAM, but was tested with OpenFOAM 5.0

<br/>
<hr>


Note: foam-extend in windows 
--------------------------------
i. If we want to use the foam-extend installed in windows directory,
- Always move the foam-extend in windows
- Never copy it to windows directory

ii. Changing the bashrc: 
- Make a change of relative path in the etc/bashrc file to 
  - `foamInstall=/mnt/c/Users/<yourUserName>/Documents/foam`
- Source the bashrc file once again.

<br/>
<hr>


Author
------------------

- Ranjan Dhakal 
- IPPE, Graz University of Technology 
- June 2022
 
