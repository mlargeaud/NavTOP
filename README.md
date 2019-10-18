# NavTOP
Orbit propagation code

NavTOP is an Earth-bound orbit determination code developed for research purposes within the Navigation Laboratory at the Illinois Institute of Technology. This program is developed using MATLAB, and relies heavily on object-oriented programming. 

The present features of the code are the following:
  - Generate a spacecraft object using intrinsic properties such as mass, reflectivity and reflective area.
  - Generate the initial orbit using osculating orbital elements.
  - Compute the spacecraft's trajectory through time using an integrating method (Runge-Kutta 4 for now) on equations of           motion. The perturbing forces taken into account are the gravitational perturbations due to Earth's mass distribution         around its center, the solar radiation pressure and the third body gravitational pull exerted by the Sun (the Moon's           influence is to be added in future implementations).
  - Generate a projection of the successive positions taken by the spacecraft around the Earth on a map after propagation.
  - Generate plots of the orbital elements and perturbing forces evolution during propagation.
  - Export plots and data file in a results folder if prompted by user.
  
 
 INPUT FILES:
 
The parameters needed to run the code must be indicated in an input file. Examples of such files can be found in the            'input_files' folder. 
 
The 3 first lines give a description of the case studied, the file's readiness status ("parameters need checking" or "ready    for use", for example) and the source of the parameters used. This part of the file can be modified as the user wants, provided it remains above the "case_name" section.
 
HOW TO USE THE CODE:
 
  - Create or edit an input file. 
  - Open the script called "main" in MATLAB, and run the entire code.
