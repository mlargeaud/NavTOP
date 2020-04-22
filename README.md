# NavTOP ver. 1.0
## Orbit propagation code

NavTOP is an Earth-bound orbit determination code developed for research purposes within the Navigation Laboratory at the Illinois Institute of Technology. This program is developed using MATLAB, and relies heavily on Object-Oriented Programming. 

The present features of the code are the following:
  - Generate a spacecraft object using intrinsic properties such as mass, reflectivity and reflective area.
  - Generate the initial orbit using osculating orbital elements. 
  - Compute the spacecraft's trajectory through time using an integrating method (Runge-Kutta 4 for now) on equations of           motion. The perturbing forces taken into account are the gravitational perturbations due to Earth's mass distribution         around its center, the solar radiation pressure and the third body gravitational pull exerted by the Sun (the Moon's           influence is to be added in future implementations).
  - Generate a ground track plot
  - Generate plots of the orbital elements and perturbing forces evolution during propagation.
  - Export plots and data file in a results folder if prompted by user.
  
 
## Input files
 
The parameters needed to run the code must be indicated in an input file. Examples of such files can be found in the            'input_files' folder. 
 
The lines above "Case name" give a description of the case studied, the file's readiness status ('parameters need checking' or 'ready for use', for example) and the source of the parameters used. This part of the file can be modified as the user wants, provided it remains above the "Case name" entry. The entry parameters are then defined as follows:

  - Case name: string parameter referring to the scenario's name the user wants to define.
  
  Spacecraft properties:
  
  - Main gravitating body: string parameter indicating around which celestial body the spacecraft is orbiting (only 'earth'       can be chosen for now, but additional bodies may be added through future work),
  - Spacecraft mass (kg),
  - Spacecraft reflective area (m2): necessary to compute the solar radiation pressure,
  - Spacecraft reflective coefficient (unitless): also necessary to compute the solar radiation pressure. 
  
  Main body's attitude angles with respect to ECI frame (3-1-3 rotation defined by right ascension, declination and sidereal     time angles) and initial orbit around the Sun:
  
  - Right ascension (deg),
  - Declination (deg),
  - Sidereal time (deg),
  
  - Main body's orbital inclination (deg),           							
  - Main body's longitude of ascending node (deg),  								
  - Main body's argument of periapsis (deg),      								
  - Main body's true anomaly (deg),                                
  - Main body's semi-major axis (m),                							
  - Main body's eccentricity (unitless),                 							       
  - Main body's initial true anomaly (deg).
  
  Spacecraft's initial conditions:
  
  - Initial parameters: "Orbital elements" (initial radius and velocity can not be chosen as entry parameters yet)
  
  - Inclination (deg),
  - Longitude of ascending node (deg),                            
  - Argument of periapsis (deg),                                   
  - True anomaly (deg),                                          
  - Semi-major axis (m),                                         
  - Eccentricity (unitless).
  
  Propagation parameters:
  
  - Potential maximum degree (integer, greater than or equal to 2): gravitational potential maximum degree in the series           computation,
  - Time step (sec): time step used for numerical integration,
  - Final time increment (sec): propagation timespan.

The values associated with each parameter are given on the same corresponding line. The parameters' names (left hand side of the file) should not be changed, only the corresponding values (right hand side of the file) can be modified according to the case studied.

 
## How to run the code
 
  - Create or edit an input file.  
  - Open the script called 'main' in MATLAB, and run the entire script either by hitting "run" or entering "main" in console.
  - When prompted, provide the address to the input file through the shell, such as: 
    'input_files/test_cases/molniya/molniya.txt'
  - Propagation begins after providing an input file. When propagation is complete, 3 results figures appear (successive           positions plotted on main body's 2D map, orbital elements, and perturbing forces values through time)
  - When asked by the shell, indicate whether the results need to be saved or not ('y' for yes, 'n' for 'no').
  - If the results need to be saved, a folder having the same name as the case name is created in the 'results' directory,         where the input file, the 3 results figures and a file containing the entire propagation data are saved. In case the           folder already exists, the user will be asked if the results can be overwritten. If not, the results will be saved in an       'autosave' folder, and the user will be asked to modify the case name.
  
  Should you have any questions or concerns, please contact me at maxime.largeaud@gmail.com.
