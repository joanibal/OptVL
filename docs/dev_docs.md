# Development Overview

OptVL is an ongoing project done entirely in my free time. 
Pull requests are welcome if you have something that you would like to contribute.
Suggestions are also welcome, but I can't promise that I will ever have time to address them. 
I'm using Github to track issues and bugs so open a new issue there with any concerns.

## Code philosophy
OptVL should be easy to install and use. 
To make it easy to install, OptVL should always be pip installable and only rely on core scientific Python packages like NumPy, and Matplotlib.
The user could have additional programs like Tecplot and ParaView, but they should never be required.
To keep OptVL easy to use, it is import to keep the required level of Python knowledge low as most users will be aerospace engineers first and may not have much Python experience.
 

## Development roadmap 
Below are some features I'm considering implementing.

- Make it easier to use Camber and thickness as design variables
- Flush out OpenMDAO geometry component that makes it easy to set geometric variables via splines
- Add body components to geometry plots and Cp surface output
- Add actuator disk model of Conway similar to VSPAERO
- Add interface to AeroSandbox

