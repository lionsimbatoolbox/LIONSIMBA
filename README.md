# LIONSIMBA - Lithium-ION SIMulation BAttery Toolbox
 A Matlab framework based on a finite volume model suitable for Li-ion battery design, simulation, and control

Consumer electronics, wearable and personal health devices, power networks, microgrids, and hybrid electric vehicles (HEVs) are some of the many applications of lithium-ion batteries. Their optimal design and management are important for safe and profitable operations. The use of accurate mathematical models can help in achieving the best performance. This page provides a parametrizable Matlab framework for battery design, simulation, and control of Li-ion cells or battery packs. Based on the well known, theory-based pseudo-two-dimensional (P2D) model, the framework has been coded by reformulating the set of PDEs describing the cell behavior into a set of DAEs. The time domain is left continuous, while the spatial domain is discretized according to the Finite Volume Method (FVM). The time-adaptive DAE solver IDA is used to solve the resulting set of DAEs.

##### **Read the Journal paper** [here](http://jes.ecsdl.org/content/163/7/A1192.abstract)
-----------------------------------------------------------------
## Official Web Page

Connect to the official web page to get the latest news
 
[http://sisdin.unipv.it/labsisdin/lionsimba.php](http://sisdin.unipv.it/labsisdin/lionsimba.php)
-----------------------------------------------------------------
## Authors

+ [Marcello Torchio](https://www.linkedin.com/in/marcello-torchio-4176368a)
+ [Lalo Magni](http://sisdin.unipv.it/labsisdin/people/maglal/maglal.php)
+ [Bhushan R. Gopaluni](http://www.chbe.ubc.ca/profile/bhushan-gopaluni/)
+ [Richard D.Braatz](http://web.mit.edu/cheme/people/profile.html?id=48)
+ [Davide M. Raimondo](http://sisdin.unipv.it/labsisdin/raimondo/raimondo.php)

-----------------------------------------------------------------
## Citations

If LIONSIMBA Toolbox is used for research purposes, the authors would like to have it mentioned. Here below the necessary information can be found
+ **Title:** LIONSIMBA: A Matlab framework based on a finite volume model suitable for Li-ion battery design, simulation, and control

+ **Journal:** The Electrochemical Society

+ **Volume:** 163

+ **Number:** 7

+ **Pages:** A1192-A1205

+ **Year:** 2016

##### **Download here the [BibTeX](http://sisdin.unipv.it/labsisdin/mtorchio/lionsimba.bib) file**



-----------------------------------------------------------------

## Usage
 
 Download the latest zip package of [LIONSIMBA](https://github.com/lionsimbatoolbox/LIONSIMBA/archive/master.zip) or clone the Git Repository using the following command
 
 ```sh
$ git clone https://github.com/lionsimbatoolbox/LIONSIMBA.git
```
## Bugs report

Please feel free to use the 'issue' section on GitHub or write me at

marcello (**dot**) torchio01 (**at**) ateneopv (**dot**) it

## Changelog

### Last Update 09/23/2016 - V 1.022 Released

+ Added support for analytical Jacobian. LIONSIMBA is now able to derive automatically the analytical form of the Jacobian describing the P2D dynamics. This knowledge is the exploited from the integration process to speed up the resolution of the DAEs. (Thanks to Dr. Sergio Lucia and Prof. Rolf Findeisen for pointing us out the automatic differentiation provided by CasADi toolbox)
+ Minor fixes in the examples.

### 08/27/2016

+ Fixed bug in multicell simulation (Thanks to Chintan Pathak for pointing out the bug)

### V 1.021b
+ Fixed SOC calculation bug for Fick's diffusion
+ Minor fixes
