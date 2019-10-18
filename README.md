# NavTOP
Orbit propagation code

NavTOP is an Earth-bound orbit determination code developed for research purposes within the Navigation Laboratory at the Illinois Institute of Technology. This program is developed using MATLAB, and relies heavily on object-oriented programming. 

The present functionalities of the code are as follows :
  - Generate a spacecraft object using intrinsic properties such as mass, reflectivity or reflective area.
  - Generate the initial orbit using the osculating orbital elements.
  - Compute the spacecraft's trajectory through time using an integrating method (Runge-Kutta 4 for now) on equations of           motion. The perturbing forces taken into account are the gravitational perturbations due to Earth's mass distribution         around its center, the solar radiation pressure and the third body gravitational pull exerted by the Sun (the Moon's           influence is to be added in future implementations).
  - Generate a projection of the successive positions taken by the spacecraft around the Earth on a map after propagation.
  - Generate plots of the orbital elements and perturbing forces evolution during propagation.
  - Export plots and data file in a results folder if prompted by user. 
