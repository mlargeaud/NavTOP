% SPACECRAFT CLASS
classdef spacecraft
    
    properties

        mb;                 % main_body object, main orbited body
        
        mass;               % (scalar, kilograms), spacecraft's mass
        radius;             % (3 by 1 vector, meters), inertial planet-centric radius
        velocity;           % (3 by 1 vector, meters/second), inertial velocity
        radBF;              % (3 by 1 vector, meters), radius vector in body-fixed frame
        velBF;              % (3 by 1 vector, meters/second), velocity vector in body-fixed frame
        r_2sun;             % (3 by 1 vector, meters), vector from s/c to Sun in inertial frame
        reflarea;           % (scalar, meters squared), reflective area
        reflcoeff;          % (scalar, unitless), reflective coefficient
        
        orb_elems;          % (6 by 1 vector), contains orbital elements : 
                            % (1) inclination (scalar, radians)
                            % (2) longitude of ascending node (scalar, radians)
                            % (3) argument of periapsis (scalar, radians)
                            % (4) true anomaly (scalar, radians)
                            % (5) semi-major axis (scalar, meters)
                            % (6) eccentricity (scalar, unitless)
                            
        h;                  % (scalar, meters squared/second), specific angular momentum
        gamma;              % (scalar, radians), flight path angle
        
        IsInUmbra;          % (boolean), indicates whether the spacecraft is in the body's umbra
            
    end
    
    methods(Static)
        
        % CONSTRUCTOR METHOD
        % Creates a satellite object when called
        %
        % Can be used as :
        %       spacecraft(mb, mass, orb_elems)
        %   OR  spacecraft(mb, mass, rad, vel)
        %
        % Inputs :
        %   mb (main_body object) - main orbited body generated
        %                   using main_body.m 
        %   
        %   mass (scalar, kilograms) - spacecraft's mass
        %
        %   orb_elems - vector containing the orbital elements describing 
        %                   the orbit of the spacecraft in this order : 
        %                   (1) inclination (scalar, radians)
        %                   (2) longitude of ascending node (scalar, radians)
        %                   (3) argument of periapsis (scalar, radians)
        %                   (4) true anomaly (scalar, radians)
        %                   (5) semi-major axis (scalar, meters)
        %                   (6) eccentricity (scalar, unitless)
        %
        %   rad (scalar, meters) - planet-centric radius in inertial frame
        %
        %   vel (scalar, meters/second) - spacecraft's inertial velocity
        function sc = spacecraft(varargin)
            
            if nargin == 5
                
                sc.mb = varargin{1};
                sc.mass = varargin{2};       
                sc.reflarea = varargin{3};
                sc.reflcoeff = varargin{4};
                sc.orb_elems = varargin{5};
                
                sc.gamma = spacecraft.setGammaFromOrbElems ...
                    (sc.orb_elems(4), sc.orb_elems(6));
                
                sc.radius = spacecraft.setR(sc.orb_elems);
                
                sc.velocity = spacecraft.setV ... 
                    (sc.mb.mu, sc.radius, sc.gamma, sc.orb_elems);
                
                sc.radBF = sc.mb.BF2I'*sc.radius;
                sc.velBF = sc.mb.BF2I'*sc.velocity;
                
                sc.h = cross(sc.radius, sc.velocity);
                                
            else 

                sc.mb = varargin{1};
                sc.mass = varargin{2};
                sc.reflarea = varargin{3};
                sc.reflcoeff = varargin{4};
                sc.radius = varargin{5};
                sc.velocity = varargin{6};
                sc.orb_elems = spacecraft.compute_orb_elems ...
                    (sc.mb.mu, sc.radius, sc.velocity);
                sc.radBF = sc.mb.BF2I'*sc.radius;
                sc.velBF = sc.mb.BF2I'*sc.velocity;
                
            end 
            
            sc.r_2sun = spacecraft.vec_sc2sun(sc.orb_elems, sc.mb.OrbElemsMB);
            sc.IsInUmbra = spacecraft.InUmbra(sc.radius, sc.mb.OrbElemsMB, ...
                sc.mb.SCONE2I, sc.mb.hCone,sc.mb.ConeApert);
            
        end
        
        
        
        % COMPUTE ORBITAL ELEMENTS FROM INERTIAL RADIUS AND VELOCITY
        %
        % Inputs :
        %   mu (scalar, meters cubed/second squared) - gravitating body's 
        %                   gravitational parameter
        %
        %   rad (3 by 1 vector, meters) - radius vector from body's center 
        %                   to s/c expressed in inertial frame
        %
        %   vel (3 by 1 vector, meters/second) - s/c velocity vector expressed 
        %                   in inertial frame 
        %
        % Outputs :
        %   orb_elems - vector containing the orbital elements describing 
        %                   the orbit of the spacecraft in this order : 
        %                   (1) inclination (scalar, radians)
        %                   (2) longitude of ascending node (scalar, radians)
        %                   (3) argument of periapsis (scalar, radians)
        %                   (4) true anomaly (scalar, radians)
        %                   (5) semi-major axis (scalar, meters)
        %                   (6) eccentricity (scalar, unitless)
        function orb_elems = compute_orb_elems(mu, rad, vel)
            
            orb_elems = zeros(6, 1);
     
            r = norm(rad);
            v = norm(vel);
            
            I = [1;0;0];
            J = [0;1;0];
            K = [0;0;1];
            
            % angular momentum
            h = cross(rad, vel);
            
            % specific energy
            ksi = (v^2/2) - (mu/r);
            
            % node vector 
            n = cross(K, h);
            
            % semi-major axis
            orb_elems(5) = -mu/(2*ksi);
            
            % eccentricity
            e = (1/mu)*((v^2 - (mu/r))*rad - dot(rad, vel)*vel);
            orb_elems(6) = norm(e);
            
            % inclination
            orb_elems(1) = acos(dot(K, h)/norm(h));
            
            % longitude of ascending node 
            CosOmega = dot(I, n)/norm(n);
            SinOmega = dot(J, n)/norm(n);
            orb_elems(2) = atan2(SinOmega, CosOmega);
            
            % argument of periapsis
            ez = dot(e, K);
            if (ez >= 0)
                orb_elems(3) = acos(dot(n, e)/(norm(n)*norm(e)));
            else
                orb_elems(3) = 2*pi - acos(dot(n, e)/(norm(n)*norm(e)));
            end 
            
            % true anomaly
            if (dot(rad, vel) >= 0)
                orb_elems(4) = acos(dot(e, rad)/(norm(e)*r));
            else 
                orb_elems(4) = 2*pi - acos(dot(e, rad)/(norm(e)*r));
            end 
          
        end
        
        
        
        % COMPUTE INERTIAL PLANET-CENTRIC RADIUS FROM ORBITAL ELEMENTS
        %
        % Inputs :
        %   OrbElemsSC - vector containing the orbital elements describing 
        %                   the orbit of the spacecraft in this order : 
        %                   (1) inclination (scalar, radians)
        %                   (2) longitude of ascending node (scalar, radians)
        %                   (3) argument of periapsis (scalar, radians)
        %                   (4) true anomaly (scalar, radians)
        %                   (5) semi-major axis (scalar, meters)
        %                   (6) eccentricity (scalar, unitless)
        %
        % Outputs :
        %   r (scalar, meters) - spacecraft's inertial planet-centric radius
        function r = setR(OrbElemsSC)
            
            i = OrbElemsSC(1);
            O = OrbElemsSC(2);
            o = OrbElemsSC(3);
            n = OrbElemsSC(4);
            a = OrbElemsSC(5);
            e = OrbElemsSC(6);
            
            p = a*(1 - e^2);
            
            r = p/(1 + e*cos(n)); 
            
            % compute rotation matrix from local orbit frame to inertial frame
            R_IO = rotation_functions.compute_rot_matrix(O, i, o + n, '313');

            r = r*R_IO*[1; 0; 0];
            
        end
        
        
        
        % COMPUTE INERTIAL VELOCITY FROM ORBITAL ELEMENTS
        %
        % Inputs :
        %   mu (scalar, meters cubed/second squared) - gravitating body's 
        %                   gravitational parameter
        %
        %   rad (3 by 1 vector, meters) - radius vector from body's center 
        %                   to s/c expressed in inertial frame
        %
        %   gamma (scalar, radians) - spacecraft's flight path angle
        %
        %   OrbElemsSC - vector containing the orbital elements describing 
        %                   the orbit of the spacecraft in this order : 
        %                   (1) inclination (scalar, radians)
        %                   (2) longitude of ascending node (scalar, radians)
        %                   (3) argument of periapsis (scalar, radians)
        %                   (4) true anomaly (scalar, radians)
        %                   (5) semi-major axis (scalar, meters)
        %                   (6) eccentricity (scalar, unitless)
        %
        % Outputs :
        %   vel (scalar, meters/second) - spacecraft's inertial velocity 
        function vel = setV(mu, rad, gamma, OrbElemsSC)
            
            r = norm(rad);
            
            i = OrbElemsSC(1);
            O = OrbElemsSC(2);
            o = OrbElemsSC(3);
            n = OrbElemsSC(4);
            a = OrbElemsSC(5); 
            
            vel = sqrt(mu*((2/r) - (1/a)));
            vec = vel*[sin(gamma); cos(gamma); 0];
            
            % compute rotation matrix from local orbit frame to inertial frame
            R_IO = rotation_functions.compute_rot_matrix(O, i, o + n, '313');
            
            vel = R_IO*vec;
            
        end
        
        
        
        % COMPUTE FLIGHT PATH ANGLE FROM ORBITAL ELEMENTS
        %
        % Inputs :
        %   nu (scalar, radians) - spacecraft's true anomaly
        %
        %   ecc (scalar, unitless) - spacecraft orbit's eccentricity
        %
        % Outputs :
        %   gam (scalar, radians) - spacecraft's flight path angle
        function gam = setGammaFromOrbElems(nu, ecc)
            
            gam = atan2((ecc*sin(nu)), (1 + ecc*cos(nu)));
            
        end
        
        
        % COMPUTE ECCENTRIC ANOMALY FROM MEAN ANOMALY
        %
        % Inputs :
        %   M (scalar, radians) - spacecraft's mean anomaly
        %
        %   ecc (scalar, unitless) - spacecraft orbit's eccentricity
        %
        % Outputs :
        %   E (scalar, radians) - spacecraft's eccentric anomaly
        %
        % Comments :
        %   The eccentric anomaly is computed using the Newton-Raphson 
        %   root-finding method. This works as long as the root is far 
        %   from a potential stationary point on the function. 
        %   In this specific case, the function can be written as 
        %   
        %       f(E) = E - ecc*sin(E) - M,
        %
        %   which implies
        %   
        %       f'(E) = 1 - ecc*cos(E).
        %
        %   As long as ecc is inferior to 1, f'(E) can not reach 0, and the
        %   Newton-Raphson method can be used flawlessly. However, if
        %   ecc >= 1, f'(E) can have potential roots, and the method may 
        %   fail. This is not likely to lead to a division by 0 (because of
        %   the round-off errors), but the method may converge to an incorrect
        %   value. Moreover, the concept of mean anomaly can not be used
        %   for parabolic or hyperbolic orbits, since it refers to an angle
        %   computed from a fictitious circular orbit having the same period
        %   as the actual one. Therefore, the use of an eccentricity higher
        %   or equal to 1 is prevented. 
        function E = Mean2EccAnom(M, ecc)
            
            % Controls if the eccentricity is superior or equal to 1
            if ecc >= 1
                error('Eccentricity must remain below 1')
            end 
            
            % Convergence and number of iterations tolerance for Newton-
            % Raphson method 
            threshold = 10^(-6);
            it_max = 1000000;
            
            % Initial guess 
            E_i = pi;
            
            % First iteration of Newton-Raphson method
            count = 1;
            E = E_i - (E_i - ecc*sin(E_i) - M)/(1 - ecc*cos(E_i));
            
            % Main loop
            while (abs(E - E_i) >= threshold)
                
                % No convergence 
                if count == it_max
                    warning('on', 'backtrace');
                    warning ...
                    ('Newton-Raphson did not converge to eccentric anomaly value');
                    return;
                end
                
                count = count + 1;
                E_i = E;
                
                % Compute next value
                E = E_i - (E_i - ecc*sin(E_i) - M)/(1 - ecc*cos(E_i));
                
            end  
            
        end 
        
        
        % COMPUTE TRUE ANOMALY FROM ECCENTRIC ANOMALY
        %
        % Inputs :
        %   E (scalar, radians) - spacecraft's eccentric anomaly
        %
        %   ecc (scalar, unitless) - spacecraft orbit's eccentricity
        %
        % Outputs :
        %   nu (scalar, radians) - spacecraft's true anomaly
        function nu = Ecc2TrueAnom(E, ecc)
            
            % Test if E < pi (spacecraft is before apogee)
            if E < pi
                
                % Then tan(E/2) > 0 and nu is given in the correct quadrant 
                % by the expression used 
                nu = 2*atan(sqrt((1+ecc)/(1-ecc))*tan(E/2));
            
            else 
                
                % Here, tan(E/2) < 0 and the expression above gives nu < 0. 
                % Shifting by 2*pi
                nu = 2*pi + 2*atan(sqrt((1+ecc)/(1-ecc))*tan(E/2));
                
            end
            
        end 
        
        
        % COMPUTE RADIUS VECTOR FROM SPACECRAFT TO SUN IN INERTIAL FRAME
        % Considering first (orbited) body orbits the Sun. 
        %
        % Input :
        %   OrbElemsSC - vector containing the orbital elements describing 
        %                   the orbit of the spacecraft in this order : 
        %                   (1) inclination (scalar, radians)
        %                   (2) longitude of ascending node (scalar, radians)
        %                   (3) argument of periapsis (scalar, radians)
        %                   (4) true anomaly (scalar, radians)
        %                   (5) semi-major axis (scalar, meters)
        %                   (6) eccentricity (scalar, unitless)
        %
        %   OrbElemsMB - vector containing the orbital elements describing 
        %                   the orbit of the main body around the sun in this order : 
        %                   (1) inclination (scalar, radians)
        %                   (2) longitude of ascending node (scalar, radians)
        %                   (3) argument of periapsis (scalar, radians)
        %                   (4) true anomaly (scalar, radians)
        %                   (5) semi-major axis (scalar, meters)
        %                   (6) eccentricity (scalar, unitless)
        %
        % Output : 
        %   r_b (vector, meters) - radius vector from spacecraft to the Sun
        %                   expressed in inertial frame
        function r_b = vec_sc2sun(OrbElemsSC, OrbElemsMB)
        
            i_bd = OrbElemsMB(1);
            O_bd = OrbElemsMB(2);
            o_bd = OrbElemsMB(3);
            n_bd = OrbElemsMB(4);
            a_bd = OrbElemsMB(5);
            e_bd = OrbElemsMB(6);
            
            i_sc = OrbElemsSC(1);
            O_sc = OrbElemsSC(2);
            o_sc = OrbElemsSC(3);
            n_sc = OrbElemsSC(4);
            a_sc = OrbElemsSC(5);
            e_sc = OrbElemsSC(6);
            
            % compute rotation matrix from local orbit frame to inertial
            % frame for main body
            R_IO_bd = rotation_functions.compute_rot_matrix(O_bd, i_bd, ...
                o_bd + n_bd, '313');
            
            % compute rotation matrix from local orbit frame to inertial
            % frame for spacecraft
            R_IO_sc = rotation_functions.compute_rot_matrix(O_sc, i_sc, ...
                o_sc + n_sc, '313');
            
            % radius vector from s/c to main body in local frame
            p_sc = a_sc*(1 - e_sc^2);
            r_sc = -p_sc/(1 + e_sc*cos(n_sc))*[1; 0; 0];
            
            % radius vector from main to third body in third body's local frame
            p_bd = a_bd*(1 - e_bd^2);
            r_bd = p_bd/(1 + e_bd*cos(n_bd))*[1; 0; 0];
            
            % expressing r_sc and r_bd in inertial frame (fixed to main body)
            r_sc = R_IO_sc*r_sc;
            r_bd = -R_IO_bd*r_bd;
            
            % final vector from spacecraft to third body
            r_b = r_sc + r_bd;
            
        end
        
        
        
        % INDICATE WHETHER THE SPACECRAFT IS IN GRAVITATING BODY'S UMBRA
        %
        % Input : 
        %   rad (3 by 1 vector, meters) - radius vector from body's center 
        %                   to s/c expressed in inertial frame
        %
        %   OrbElemsMB - vector containing the orbital elements describing 
        %                   the orbit of the main body around the Sun in this order : 
        %                   (1) inclination (scalar, radians)
        %                   (2) longitude of ascending node (scalar, radians)
        %                   (3) argument of periapsis (scalar, radians)
        %                   (4) true anomaly (scalar, radians)
        %                   (5) semi-major axis (scalar, meters)
        %                   (6) eccentricity (scalar, unitless)
        %
        %   SCONE2I (3 by 3 matrix) - rotation matrix from umbra cone to inertial frame
        %
        %   hCone (scalar, meters) - shadow cone's height
        %
        %   ConeApert (scalar, radians) - cone's aperture angle
        %
        % Output :
        %   IsInUmbra (logical) - indicates whether the spacecraft is in
        %                   the umbra or not
        %                   True -> s/c is inside the umbra
        %                   False -> s/c is outside the umbra
        function IsInUmbra = InUmbra(rad, OrbElemsMB, SCONE2I, hCone, ConeApert)
            
            i_bd = OrbElemsMB(1);
            O_bd = OrbElemsMB(2);
            o_bd = OrbElemsMB(3);
            n_bd = OrbElemsMB(4);
            a_bd = OrbElemsMB(5);
            e_bd = OrbElemsMB(6);
            
            % compute rotation matrix from local orbit frame to inertial
            % frame for main body
            R_IO_bd = rotation_functions.compute_rot_matrix(O_bd, i_bd, ...
                o_bd + n_bd, '313');
            
            % body-centric radius vector from body to Sun in local frame
            p_bd = a_bd*(1 - e_bd^2);
            r_bd = -p_bd/(1 + e_bd*cos(n_bd))*[1; 0; 0];
            
            % body-centric radius vector from body to Sun in inertial frame
            r_bd = R_IO_bd*r_bd;
            
            % unit vector from body to Sun
            e_bd = r_bd/norm(r_bd); 
            
            % compute spacecraft position in umbra frame
            Pos = SCONE2I'*(hCone*(e_bd) + rad);
            
            F = (Pos(1)^2 + Pos(2)^2)*(cos(ConeApert)^2) ...
                - (Pos(3)^2)*(sin(ConeApert)^2);
            
            % test whether the position coordinates verify the cone's equation
            if F <= 0 && Pos(3) >= 0 && Pos(3) <= hCone 
                
                IsInUmbra = true;
                return
                
            else 
                
                IsInUmbra = false;
                return
                
            end      
            
        end 
        
        
        
        % PROPAGATE ORBIT USING RUNGE-KUTTA 4 METHOD
        %
        % Inputs :
        %   sc (spacecraft) - spacecraft object
        %
        %   ti (scalar, seconds) - initial time of propagation
        %
        %   tf (scalar, seconds) - time at which the propagation stops
        %
        %   dt (scalar, seconds) - time step
        %
        %   raw_coeff - matrix containg the normalized spherical harmonic
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
        %   max_deg (integer) - gravitational potential maximum degree 
        %
        % Outputs :
        %   sc (spacecraft) - spacecraft object with properties corresponding
        %                   to final instant tf
        %
        %   tab - matrix containing the values taken by the spacrecraft's 
        %                   properties at each instant during propagation. 
        %                   Each row corresponds to a particular instant :
        %                   tab(:, 1) -> time increment
        %                   tab(:, 2:4) -> cartesian coordinates (ECI)
        %                   tab(:, 5:7) -> cartesian velocity (ECI)
        %                   tab(:, 8) -> planet-centric radius
        %                   tab(:, 9) -> latitude
        %                   tab(:, 10) -> longitude
        %                   tab(:, 11) -> main body's right ascension
        %                   tab(:, 12) -> main body's declination
        %                   tab(:, 13) -> main body's sidereal time
        %                   tab(:, 14:16) -> gravitational perturbation
        %                   tab(:, 17:19) -> sol. rad. pressure 
        %                   tab(:, 20:22) -> third body gravity
        %                   tab(:, 23) -> inside/outside umbra index
        %                   tab(:, 24:29) -> orbital elements
        %                   tab(:, 30) -> RK4 computation error on orb. elems
        function [sc, tab] = propagate_orbit_orb_elems_RK4 ...
                (sc, ti, tf, dt, raw_coeff, max_deg) 
            
            harmonic_coeff = force_functions.read_harm_coeff_from_file ...
                (raw_coeff, max_deg);
                
            OrbElemsSC = sc.orb_elems;

            i = OrbElemsSC(1);
            O = OrbElemsSC(2);
            o = OrbElemsSC(3);
            n = OrbElemsSC(4);
            a = OrbElemsSC(5);
            e = OrbElemsSC(6);  

            OrbElemsMB = sc.mb.OrbElemsMB;
            SCONE2I = sc.mb.SCONE2I;
            hCone = sc.mb.hCone;
            ConeApert = sc.mb.ConeApert;

            rasc = sc.mb.rasc;
            dec = sc.mb.dec;
            stime = sc.mb.stime;
            
            % compute spacecraft-to-Sun initial radius and unit vectors
            e_sc2sun = sc.r_2sun/norm(sc.r_2sun);
            
            % compute gravitational acceleration at initial location
            Fg = force_functions.GravAcc(OrbElemsSC, sc.mb.mu, sc.mb.r0, ... 
                rasc, dec, stime, harmonic_coeff, max_deg, 'oblateness');
            
            % compute solar radiation pressure at initial location 
            Fsrp = [0;0;0]; % default value
              
            if sc.IsInUmbra == false
                Fsrp = force_functions.SRP ...
                    (norm(sc.r_2sun), O, i, o, n, sc.reflarea, ...
                    sc.reflcoeff, sc.mass, e_sc2sun);
            end 
            
            % compute third body gravity from Sun at initial location
            Ftb = force_functions.third_body_gravity(norm(sc.radius), ... 
                sc.mb.muSun, norm(sc.r_2sun), e_sc2sun, O, i, o, n);
           
            % initial values for time increment and error on orbital elements
            t = ti;
            err = NaN;
            
            % initializing result table
            N = ((tf - ti)/dt) + 1;
            tab = zeros(N, 36);
               
            % propagating orbit for each time increment until tf is reached 
            for j = 1:N

                % computing the location of the s/c in spherical coordinates
                sph_coord = rotation_functions.LocalToPlanetaryCoord ...
                    (rasc, dec, stime, e, a, O, i, o, n);
                
                % computing spacecraft location in body-fixed frame
                sc.radBF = sc.mb.BF2I'*sc.radius;
                
                % spacecraft velocity w.r.t body-fixed frame, expressed in
                % body fixed frame as well
                ECI_vel_BF = sc.velocity - cross(sc.mb.RotVecIn, sc.radius);
                sc.velBF = sc.mb.BF2I'*ECI_vel_BF; 
                
                % saving data in array
                tab(j, :) = [t, sc.radius', sc.velocity', sc.radBF', ...
                    sc.velBF', sph_coord', rasc, dec, stime, Fg', Fsrp', ...
                    Ftb', sc.IsInUmbra, OrbElemsSC', err];
                
                % solve variational equations using RK4 method                                 
                [sc.orb_elems, Fg, Fsrp, Ftb] = ...
                    variational_eq_solvers.compute_next_orb_elems_RK4(...
                    dt, sc.mb.mu, sc.mb.r0, rasc, dec, stime, ...
                    OrbElemsMB, harmonic_coeff, max_deg, SCONE2I, ...
                    hCone, ConeApert, sc.mb.muSun, ...
                    sc.reflarea, sc.reflcoeff, sc.mass, OrbElemsSC);
                
                % compute error
                err = log(norm(OrbElemsSC - tab(j, 10:15)'))/log(10);
                
                % update remaining orbital parameters
                sc.gamma = spacecraft.setGammaFromOrbElems ...
                    (sc.orb_elems(4), sc.orb_elems(6));
                
                sc.radius = spacecraft.setR(sc.orb_elems);
                
                sc.velocity = spacecraft.setV ...
                    (sc.mb.mu, sc.radius, sc.gamma, sc.orb_elems);
                
                sc.h = cross(sc.radius, sc.velocity);
                
                sc.r_2sun = spacecraft.vec_sc2sun(OrbElemsSC, OrbElemsMB);
                
                % update spacecraft position w.r.t. umbra 
                sc.IsInUmbra = spacecraft.InUmbra ...
                    (sc.radius, OrbElemsMB, SCONE2I, hCone, ConeApert);
                
                % update main body's orientation
                sc.mb = main_body.update_angles(sc.mb, dt);
                
                % update time increment
                t = t + dt;
                
                % collecting new parameters                
                OrbElemsSC = sc.orb_elems;
            
                i = OrbElemsSC(1);
                O = OrbElemsSC(2);
                o = OrbElemsSC(3);
                n = OrbElemsSC(4);
                a = OrbElemsSC(5);
                e = OrbElemsSC(6);  
                
                OrbElemsMB = sc.mb.OrbElemsMB;
                SCONE2I = sc.mb.SCONE2I;
                hCone = sc.mb.hCone;
                ConeApert = sc.mb.ConeApert;
                
                rasc = sc.mb.rasc;
                dec = sc.mb.dec;
                stime = sc.mb.stime;

                %j/N
                
            end 
                        
        end
        
    end 
    
end 