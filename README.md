# LIONSIMBA - Lithium-ION SIMulation BAttery Toolbox
 A Matlab framework based on a finite volume model suitable for Li-ion battery design, simulation, and control

Consumer electronics, wearable and personal health devices, power networks, microgrids, and hybrid electric vehicles (HEVs) are some of the many applications of lithium-ion batteries. Their optimal design and management are important for safe and profitable operations. The use of accurate mathematical models can help in achieving the best performance. This page provides a parametrizable Matlab framework for battery design, simulation, and control of Li-ion cells or battery packs. Based on the well known, theory-based pseudo-two-dimensional (P2D) model, the framework has been coded by reformulating the set of PDEs describing the cell behavior into a set of DAEs. The time domain is left continuous, while the spatial domain is discretized according to the Finite Volume Method (FVM). The time-adaptive DAE solver IDA is used to solve the resulting set of DAEs.
 
## Usage
 
 Download the latest zip package of [LIONSIMBA](https://github.com/lionsimbatoolbox/LIONSIMBA/archive/master.zip) or clone the Git Repository using the following command
 
 ```sh
$ git clone https://github.com/lionsimbatoolbox/LIONSIMBA.git
```
 
## Changelog

### Last Update 08/27/2016

+ Fixed bug in multicell simulation (Thanks to Chintan Pathak for pointing out the bug)

### V 1.021b
+ Fixed SOC calculation bug for Fick's diffusion
+ Minor fixes
 
 
