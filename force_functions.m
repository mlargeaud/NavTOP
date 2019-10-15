classdef force_functions
    
    methods(Static)
        
        % ACCELERATION PRODUCED BY SOLAR RADIATION PRESSURE
        % 
        % Inputs :
        %   r_sc (scalar, meters) - gives the mean distance between the Sun 
        %                   and the spacecraft
        %
        %   lon_asc_node_sc (scalar, radians) - gives the longitude of 
        %                   ascending node of the spacecraft's orbit
        %
        %   inc_sc (scalar, radians) - gives the inclination of the
        %                   spacecraft's orbit
        %
        %   arg_per_sc (scalar, radians) - gives the argument of periapsis
        %                   of the spacecraft's orbit
        %
        %   tr_anom_sc (scalar, radians) - gives the true anomaly of the
        %                   spacecraft
        %
        %   A (scalar, meters squared) - reflective area of the spacecraft
        %
        %   C (scalar, unitless) - coefficient associated with reflective 
        %                   property of the spacecraft
        %                   0 < C < 2
        %                   0 all radiation is absorbed 
        %                   2 all radiation is reflected
        %
        %   m (scalar, kilograms) - mass of the spacecraft
        %
        %   e_s (vector, unitless) - unit vector from spacecraft to the Sun 
        %                   in inertial frame
        %
        % Outputs :
        %   f_SRP (vector, meters/seconds squared) - acceleration vector due 
        %   to solar radiation pressure, expressed in local orbital frame
        function f_SRP = SRP(r_sc, lon_asc_node_sc, inc_sc, arg_per_sc,  ...
                tr_anom_sc, A, C, m, e_s)
           
            % Sun-Earth mean distance (1 U.A)
            AU = 149597870700;
            
            % Solar radiation pressure at Earth
            P = 1361/299792458;
            
            % rotation matrix from inertial to local frame
            R_OI = (rotation_functions.compute_rot_matrix(lon_asc_node_sc, ...
                inc_sc, arg_per_sc + tr_anom_sc, '313'))';
            
            e_s_loc = R_OI*e_s;
            
            % acceleration perturbation
            f_SRP = -((P*A*C/m)*(AU/r_sc)^2)*e_s_loc;
            
        end
        
        
        
        % THIRD BODY GRAVITATIONAL ACCELERATION
        %
        % Inputs :
        %   mu_s (scalar, meters cubed/seconds squared) - gravitational 
        %                   parameter of third body 
        %   
        %   r_b (scalar, meters) - distance between the closest body and the 
        %                   spacecraft
        %   
        %   r_s (scalar, meters) - distance between spacecraft and third body
        %   
        %   rho (vector, unitless) - unit vector from first to third body 
        %                   expressed in inertial frame
        %   
        %   lon_asc_node_sc (scalar, radians) - gives the longitude 
        %                   of ascending node of the spacecraft's orbit
        %
        %   inc_sc (scalar, radians) - gives the inclination of the
        %                   spacecraft's orbit
        %
        %   arg_per_sc (scalar, radians) - gives the argument of periapsis
        %                   of the spacecraft's orbit
        %
        %   tr_anom_sc (scalar, radians) - gives the true anomaly of the
        %                   spacecraft
        %
        % Outputs :
        %   f_th_grav (vector, meters/seconds squared) - acceleration due 
        %                   to third body gravitational pull, expressed in 
        %                   local orbital frame
        function f_th_grav = third_body_gravity(r_b, mu_s, r_s, rho, ...
                lon_asc_node_sc, inc_sc, arg_per_sc, tr_anom_sc) 
            
            % rotation matrix from inertial to local frame
            R_OI = (rotation_functions.compute_rot_matrix(lon_asc_node_sc, ...
                inc_sc, arg_per_sc + tr_anom_sc, '313'))';
            
            % third body gravity
            f_th_grav = (mu_s/(r_s^3))*R_OI*(3*(rho*rho') - eye(3))*R_OI' ...
                *[r_b; 0; 0];
            
        end
        
        
        
        % COMPUTE GRAVITATIONAL PERTURBATION IN LOCAL ORBITAL FRAME
        % 
        % Inputs : 
        %   orb_elems (vector) contains the orbital elements describing 
        %                   the orbit of the spacecraft in this order : 
        %                   (1) inclination (radians)
        %                   (2) longitude of ascending node (radians)
        %                   (3) argument of periapsis (radians)
        %                   (4) true anomaly (radians)
        %                   (5) semi-major axis (meters)
        %                   (6) eccentricity (unitless)
        %
        %   mu (scalar, meters cubed/seconds squared) - gravitational 
        %                   parameter of orbited body
        %
        %   r0 (scalar, meters) - mean radius of orbited body
        %
        %   rasc (scalar, radians) - right ascension of orbited body
        %
        %   dec (scalar, radians) - declination of orbit body
        %
        %   stime (scalar, radians) - sidereal time of orbited body
        %
        %   harms (matrix) - contains the spherical harmonic coefficients
        %                   associated with the orbited body. Use the
        %                   function force_functions.read_harm_coeff_from_file 
        %                   with the raw data from space agencies' servers 
        %                   to generate this matrix
        %
        %   max_deg (integer) - gravitational potential maximum degree 
        %
        %   par (string) - parameter used to consider J2 harmonic (or not) 
        %                   in the computations
        %                   "oblateness" -> with J2
        %                   "llumpiness" -> without J2
        function g = GravAcc(orb_elems, mu, r0, rasc, dec, stime, harms, ...
                max_deg, par)
            
            i = orb_elems(1);
            O = orb_elems(2);
            o = orb_elems(3);
            n = orb_elems(4);
            a = orb_elems(5);
            e = orb_elems(6);

            p = a*(1 - e^2);
            r = p/(1 + e*cos(n));
            
            % rotation matrix from planetary frame to inertial frame
            R_PI = rotation_functions.compute_rot_matrix(rasc, dec, stime, '313');

            % rotation matrix from local orbit frame to inertial frame
            R_IO = rotation_functions.compute_rot_matrix ...
                (O, i, o + n, '313');
            
            % position in planetary frame
            pos = R_PI'*R_IO*r*[1;0;0];
            
            lon = atan2(pos(2), pos(1));
            lat = acos(pos(3)/r);
            
            dU = force_functions.GravPotentialGradO4 ...
                    (mu, r0, r, lat, lon, 1000, 1000, 1000, ...
                    harms, max_deg, par);
                
            g = R_IO'*R_PI*dU;
                   
        end
        
        
        
        % COMPUTE POTENTIAL GRADIENT
        % Gradient is order 1 on the bounds, order 2 one step out of the 
        % bounds, order 4 inside the domain
        %
        % Inputs : 
        %   mu (scalar, meters cubed/seconds squared) - gravitational 
        %                   parameter of orbited body
        %   
        %   r_0 (scalar, meters) - mean radius of orbited body
        %   
        %   r (scalar, meters) - planet-centric radius
        %   
        %   lon, lat (scalars, radians) - longitude and latitude of the 
        %                   spacecraft in body-fixed frame 
        %
        %   dx, dy, dz (scalars, meters) - space steps used in finite 
        %                   difference method, along body-fixed frame's axes
        % 
        %   HarmCoeff (matrix) - contains the spherical harmonic coefficients
        %                   associated with the orbited body. Use the
        %                   function force_functions.read_harm_coeff_from_file 
        %                   with the raw data from space agencies' servers 
        %                   to generate this matrix
        %
        %   max_deg (integer) - gravitational potential maximum degree 
        %
        %   HarmParameter (string) - parameter used to consider J2 harmonic 
        %                   (or not) in the computations
        %                   "oblateness" -> with J2
        %                   "llumpiness" -> without J2
        %   
        % Outputs : 
        %   dU (scalar, meters/second squared) - gravitational potential 
        %   perturbation gradient, expressed body-fixed frame
        function dU = GravPotentialGradO4(mu, r0, r, lat, lon, dx, dy, dz, ...
                HarmCoeffs, max_deg, HarmParameter)
            
            U = zeros(5, 5, 5);
            
            % Body-fixed frame coordinates
            x = r*sin(lat)*cos(lon);
            y = r*sin(lat)*sin(lon);
            z = r*cos(lat);
            
            for l = 1:5

                for m = 1:5 

                    for n = 1:5
                        
                        U(l, m, n) = force_functions.PotentialFromBF...
                            (mu, r0, x + dx*(l-3), ...
                            y + dy*(m-3), ...
                            z + dz*(n-3), ...
                            HarmCoeffs, max_deg, HarmParameter);

                    end

                end 

            end
            
            % computation of the gradient
            [gradX, gradY, gradZ] = finite_difference_schemes.gradient_o4 ...
                (U, dx, dy, dz);
            
            dU = [gradX(3,3,3); gradY(3,3,3); gradZ(3,3,3)];
            
        end
        
        
        
        % COMPUTE GRAVITATIONAL POTENTIAL PERTURBATION
        % From body-fixed frame coordinates (cartesian)
        % 
        % Inputs :
        %   mu (scalar, meters cubed/seconds squared) - gravitational 
        %                   parameter of orbited body
        %   
        %   r_0 (scalar, meters) - mean radius of orbited body
        %
        %   x, y, z (scalars, meters) - location coordinates in body-fixed frame
        %   
        %   harmonic_coeff (matrix) - containing the spherical harmonic coefficients
        %                   associated with the orbited body. Use the
        %                   function force_functions.read_harm_coeff_from_file 
        %                   with the raw data from space agencies' servers 
        %                   to generate this matrix
        %
        %   max_deg (integer) - gravitational potential maximum degree 
        %
        %   harm_parameter (string) - parameter used to consider J2 harmonic 
        %                   (or not) in the computations
        %                   "oblateness" -> with J2
        %                   "llumpiness" -> without J2
        %   
        % Outputs : 
        %   U (Joules/kilograms.meters) - gravitational potential perturbation       
        function U = PotentialFromBF(mu, r_0, x, y, z, harmonic_coeff, ...
                max_deg, harm_parameter) 
            
            r = sqrt(x^2 + y^2 + z^2);
            lambda = atan2(y, x);
            phi = acos(z/r);
            
            temp_TS = 0.;

            start_j = 2; %default value
            if (harm_parameter == 'oblateness') %#ok<BDSCA>
                start_j = 1;
            end 

            for i = 2:max_deg
                
                % associated Legendre polynomials of degree i and order j = 0, ..., i  
                P = legendre(i, cos(phi));

                N = i+1;
                for j = start_j:N
                    
                    temp_TS = temp_TS + ((r_0/r)^i)*P(j)* ...
                        (harmonic_coeff(i,j,2)*cos((j-1)*(-pi+lambda)) ...
                        + harmonic_coeff(i,j,3)*sin((j-1)*(-pi+lambda)));             

                end 

            end  

            % add 1 to temp_TS to have whole potential
            % default value : (mu/r)*(temp_TS), PERTURBATION potential  
            U = (mu/r)*(temp_TS);

        end
        
        
        
        % READING SPHERICAL HARMONIC COEFFICIENTS FROM DATA FILE AND DENORMALIZING
        %
        % Inputs :
        % Can be used as :
        %       read_harm_coeff_from_file(raw_coeff)
        %   OR  read_harm_coeff_from_file(raw_coeff, n_max)
        %
        %   raw_coeff (matrix) - contains the normalized spherical harmonic
        %                   coefficients describing a body's gravitational
        %                   potential. To generate this matrix :
        %                   1) Download the raw data file from a space 
        %                   agency's website. The spherical harmonics should
        %                   be normalized. 
        %                   2) Import the raw data directly to MATLAB 
        %                   (should be imported as a matrix).
        %                   3) The resulting matrix from 2) is the one that 
        %                   is used in this function.
        %
        %   n_max (scalar) - maximum spherical harmonic degree to consider,
        %                   default value is 20
        %
        % Outputs : 
        %   harm_coeff (matrix) - contains the denormalized spherical 
        %                   harmonic coefficients associated with the 
        %                   orbited body. This matrix can then be used in 
        %                   potential computations.
        function harm_coeff = read_harm_coeff_from_file(varargin)
            
            if nargin == 1
                raw_coeff = varargin{1};
                n_max = 20; % default value
            else
                raw_coeff = varargin{1};
                n_max = varargin{2};
            end 

            % check if degree required is between 2 and maximum degree
            % provided by the raw data file
            if n_max < 2 || n_max > raw_coeff(size(raw_coeff, 1), 1)  
                'Warning : spherical harmonics maximum degree required is out of range' %#ok<NOPRT>           
            end 
             
            harm_coeff = zeros(n_max, n_max+1, 3);

            i = 1;
            while raw_coeff(i, 1) <= n_max
                
                n = raw_coeff(i, 1);
                j = raw_coeff(i, 2);

                C = sqrt(factorial(n+j)/(factorial(n-j)*(2-eq(j,0))*(2*n+1)));
                harm_coeff(n, j+1, 2) = raw_coeff(i, 3)/C;
                harm_coeff(n, j+1, 3) = raw_coeff(i, 4)/C;
                
                i = i+1;
                
            end 

        end
        
    end 
    
end 