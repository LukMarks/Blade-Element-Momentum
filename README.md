# Blade-Element-Momentum
A propeller design using BEM

## Introduction

 These scripts were part of my gradruation's project, it's consist in a implementation of the algorithm propoused by **Hensen**. Since the original method was used to design a **turbine** it was necessery to modify some parameters, like the signal of the air speed, to convert it into a **propeller** design. A disclamer here, if you want something more robust i suggest to you to check the [JBLADE](https://sites.google.com/site/joaomorgado23/Home) application. 

#### WARNING: Before start it's recommended that you have xfoil, for a greater precision in your design. If you are using a **Windows** operating system you'll have to keep the executable version of **Xfoil** in the same folder of BEM.py, you can donwload it in the official [Xfoil](https://web.mit.edu/drela/Public/web/xfoil/) website, otherwise if you are using a **Linux distribution** you can use the next command to donwload and install it.

```
sudo apt-get install xfoil
```
## How to Start

First things first, the main file here is the [BEM.py](https://github.com/LukMarks/Blade-Element-Momentum/blob/master/BEM.py). Inside of it a couple of options will be seted, like:

- Diameter;

- Rotational Speed;

- Airfoil of each section;

- etc.

 The method allows to determine the type of geometre the blade will have, for this the script [corda.py](https://github.com/LukMarks/Blade-Element-Momentum/blob/master/corda.py) have some built-in generator, the options are the following: 
 
 - Constant Chord; 
 
 - Trapezoidal;
 
 - [Koch's Suggestion](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19980017535.pdf).


The gif right below shows the program up and running using **Xfoil** functionatilies.

![runner](https://github.com/LukMarks/Blade-Element-Momentum/blob/master/src%20images/bem.gif)

## Results

At the end of the calculation your console will show the final values for your design, like the next few lines shows:

```
Time elapsed (m:s):  0.0 : 7.0

Results: 

Thurst:  59.99 [N]

Momentum:  4.05 [Nm]

Total Mass:  0.672 [kg]

Required Power:  762.46 [W]

Flight Power:  720.0 [W]

Efficiency:  94.42 [%]

Section(s) failed to converge 0 

Maximum external diameter:  10.0 mm

```
With it severel plots will be displayed, The following images illutrates some examples.
![angle](https://github.com/LukMarks/Blade-Element-Momentum/blob/master/src%20images/twist_angle.png)

![shape](https://github.com/LukMarks/Blade-Element-Momentum/blob/master/src%20images/shape.png)

![Cl_Cd](https://github.com/LukMarks/Blade-Element-Momentum/blob/master/src%20images/cl_cd_ratio.png)

You can also use the [solid.py](https://github.com/LukMarks/Blade-Element-Momentum/blob/master/solid.py) to export each section as a dat extension. This files are built rotating and positioning every section, like the next image shows. 

![section](https://github.com/LukMarks/Blade-Element-Momentum/blob/master/src%20images/blade_section.png)

The file follows the nomenclature: Radius_distance_from_hub_mm.tx eg.R1175.txt, which means that section had been placed in a plane with 1175 mm of distance from the propeller's hub. With all those .dat files in hand it's possible to import these curves in a CAD software, like SolidWorks, a create a loft operation. A brief example of one those dat files is showed in the next few lines.

```
0.13499931674694216 	-0.1497776659082142 	1.175267116361462 
0.1345821137585913 	-0.14889393212862786 	1.175267116361462 
0.13349693961943118 	-0.14637500647734134 	1.175267116361462 
0.1318837501613349 	-0.1425686220947054 	1.175267116361462 
0.12959037350145597 	-0.13773955886428949 	1.175267116361462 
0.1264655854301031 	-0.1318968054196098 	1.175267116361462 
0.12255326892483118 	-0.12500873084659958 	1.175267116361462 
0.1179460944998143 	-0.11715756752421341 	1.175267116361462 
0.11271830530751524 	-0.10845448795952943 	1.175267116361462 
0.10694124716504534 	-0.09902631556727995 	1.175267116361462 
0.10068621678290388 	-0.08900824695232043 	1.175267116361462 
0.09402136014250122 	-0.078556065903427 	1.175267116361462 
0.08700961481751225 	-0.0678157709930124 	1.175267116361462 

```

## References
- [Hansen](https://www.amazon.com/Aerodynamics-Wind-Turbines-Martin-Hansen/dp/1844074382/ref=sr_1_4?dchild=1&keywords=Aerodynamics+of+Wind+Turbines&qid=1584484238&sr=8-4)

- [Xfoil](https://web.mit.edu/drela/Public/web/xfoil/)

- [Koch](https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/19980017535.pdf)

- [Bézier Curves](https://ocw.mit.edu/courses/electrical-engineering-and-computer-science/6-837-computer-graphics-fall-2012/lecture-notes/MIT6_837F12_Lec01.pdf)
