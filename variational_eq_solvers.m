classdef variational_eq_solvers
    
    methods(Static)

        % PROPAGATE ORBITAL ELEMENTS ONE STEP FORWARD USING RUNGE-KUTTA 4 METHOD
        %
        % Inputs :
        %   dt (scalar, seconds) - time step
        %
        %   mu (scalar, meters cubed/second squared) - gravitational parameter 
        %                   of orbited body
        %
        %   r0 (scalar, meters) - mean radius of orbited body
        %
        %   rasc (scalar, radians) - right ascension of body-fixed frame
        %                   with respect to inertial frame 
        %
        %   dec (scalar, radians) - declination of body-fixed frame with
        %                   respect to inertial frame 
        %
        %   stime (scalar, radians) - sidereal time of body-fixed frame
        %                   with respect to inertial frame
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
        %   harms - matrix containing the spherical harmonic coefficients
        %                   associated with the orbited body. Use the
        %                   function force_functions.read_harm_coeff_from_file 
        %                   with the raw data from space agencies' servers 
        %                   to generate this matrix
        %
        %   max_deg (integer) - gravitational potential maximum degree 
        %
        %   SCONE2I (3 by 3 matrix) - rotation matrix from umbra cone to inertial frame
        %
        %   hCone (scalar, meters) - shadow cone's height
        %
        %   ConeApert (scalar, radians) - cone's aperture angle
        %
        %   mu_tb (scalar, meters cubed/second squared) - gravitational parameter 
        %                   of third body
        %
        %   A (scalar, meters squared) - spacecraft's reflective area
        %
        %   C (scalar, unitless) - spacecraft's reflective coefficient
        %
        %   m (scalar, kilograms) - spacecraft's mass
        %
        %   orb_elems (vector) - contains the orbital elements describing 
        %                   the orbit of the spacecraft in this order : 
        %                   (1) inclination (radians)
        %                   (2) longitude of ascending node (radians)
        %                   (3) argument of periapsis (radians)
        %                   (4) true anomaly (radians)
        %                   (5) semi-major axis (meters)
        %                   (6) eccentricity (unitless)
        %
        % Outputs :
        %   orb_elems_next (vector) - newly computed orbital elements vector, 
        %                   one time step forward. Elements remain in the 
        %                   same order as in the entry vector orb_elems.
        %
        %   Fg (vector, meters/second squared) - gravitational perturbation
        %
        %   Fsrp (vector, meters/second squared) - solar radiation pressure
        %
        %   Ftb (vector, meters/second squared) - third body gravity
        function [orb_elems_next, Fg, Fsrp, Ftb] = compute_next_orb_elems_RK4 ...
                (dt, mu, r0, rasc, dec, stime, OrbElemsMB, harms, max_deg, ...
                SCONE2I, hCone, ConeApert, mu_tb, A, C, m, orb_elems)
            
            [orb_elems_next, Fg, Fsrp, Ftb] = time_schemes.RK4(dt, ...
                @variational_eq_solvers.var_eq_RHS, mu, r0, rasc, dec, stime, ...
                OrbElemsMB, harms, max_deg, SCONE2I, hCone, ConeApert, ...
                mu_tb, A, C, m, orb_elems);
            
        end
        
        
        
        % COMPUTE RIGHT HAND SIDE OF VARIATIONAL EQUATIONS 
        % 
        % Inputs :
        %   mu (scalar, meters cubed/second squared) - gravitational parameter 
        %                   of orbited body
        %
        %   r0 (scalar, meters) - mean radius of orbited body
        %
        %   rasc (scalar, radians) - right ascension of body-fixed frame
        %                   with respect to inertial frame 
        %
        %   dec (scalar, radians) - declination of body-fixed frame with
        %                   respect to inertial frame 
        %
        %   stime (scalar, radians) - sidereal time of body-fixed frame
        %                   with respect to inertial frame
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
        %   harms - matrix containing the spherical harmonic coefficients
        %                   associated with the orbited body. Use the
        %                   function force_functions.read_harm_coeff_from_file 
        %                   with the raw data from space agencies' servers 
        %                   to generate this matrix
        %
        %   max_deg (integer) - gravitational potential maximum degree 
        %
        %   SCONE2I (3 by 3 matrix) - rotation matrix from umbra cone to inertial frame
        %
        %   hCone (scalar, meters) - shadow cone's height
        %
        %   ConeApert (scalar, radians) - cone's aperture angle
        %
        %   mu_tb (scalar, meters cubed/second squared) - gravitational parameter 
        %                   of third body
        %
        %   A (scalar, meters squared) - spacecraft's reflective area
        %
        %   C (scalar, unitless) - spacecraft's reflective coefficient
        %
        %   m (scalar, kilograms) - spacecraft's mass
        %
        %   orb_elems (vector) - contains the orbital elements describing 
        %                   the orbit of the spacecraft in this order : 
        %                   (1) inclination (radians)
        %                   (2) longitude of ascending node (radians)
        %                   (3) argument of periapsis (radians)
        %                   (4) true anomaly (radians)
        %                   (5) semi-major axis (meters)
        %                   (6) eccentricity (unitless)
        %
        % Outputs :
        %   RHS (vector) - contains the values of the right hand side of
        %                   the variational equations at the given orbital
        %                   elements and parameters
        %
        %   Fg (vector, meters/second squared) - gravitational perturbation
        %
        %   Fsrp (vector, meters/second squared) - solar radiation pressure
        %
        %   Ftb (vector, meters/second squared) - third body gravity
        function [RHS, Fg, Fsrp, Ftb] = var_eq_RHS(mu, r0, rasc, dec, ...
                stime, OrbElemsMB, harms, max_deg, SCONE2I, hCone, ...
                ConeApert, mu_tb, A, C, m, orb_elems)
            
            RHS = zeros(6, 1);
            
            i = orb_elems(1);
            O = orb_elems(2); 
            o = orb_elems(3);
            n = orb_elems(4);
            a = orb_elems(5);
            e = orb_elems(6);

            p = a*(1 - e^2);
            h = sqrt(mu*p);

            % radius vector from gravitating body's center to spacecraft
            rad_sc = spacecraft.setR(orb_elems);
            r = norm(rad_sc);
            
            % radius vector from spacecraft to Sun
            r_sc2sun = spacecraft.vec_sc2sun(orb_elems, OrbElemsMB);
            
            % unit vector from spacecraft to Sun in inertial frame
            e_sc2sun = r_sc2sun/norm(r_sc2sun); 
            
            % unit vector from main body to Sun in inertial frame
            e_mb2sun = (rad_sc + r_sc2sun)/norm(rad_sc + r_sc2sun);
            
            % compute gravitational perturbation from main body
            Fg = force_functions.GravAcc(orb_elems, mu, r0, rasc, dec, stime, ...
               harms, max_deg, 'oblateness');
            
            % compute solar radiation pressure 
            Fsrp = [0;0;0]; % default value
              
            if spacecraft.InUmbra(rad_sc, OrbElemsMB, SCONE2I, hCone, ...
                    ConeApert) == false
                
                Fsrp = force_functions.SRP ...
                    (norm(r_sc2sun), O, i, o, n, A, C, m, e_sc2sun);
                
            end 
            
            % compute third body gravity from Sun 
            Ftb = force_functions.third_body_gravity(norm(rad_sc), ... 
                mu_tb, norm(r_sc2sun), e_mb2sun, O, i, o, n); 
            
            % total perturbation
            f = Fg + Fsrp + Ftb;
            
            % right-hand side components
            RHS(1) = (r*cos(o + n)/h)*f(3);
            RHS(2) = (r*sin(o + n)/(h*sin(i)))*f(3);
            
            RHS(3) = -(p*cos(n)/(h*e))*f(1) ...
                + ((p + r)*sin(n)/(h*e))*f(2) ...
                - ((r*sin(o + n)*cot(i))/h)*f(3);
            
            RHS(4) = (h/(r^2)) ...
                + ((p*cos(n))/(h*e))*f(1) ...
                - ((p + r)*sin(n)/(h*e))*f(2);
            
            RHS(5) = (2*(a^2)*e*sin(n)/h)*f(1) + (2*(a^2)*p/(h*r))*f(2);
            RHS(6) = (p*sin(n)/h)*f(1) + (((p + r)*cos(n) + r*e)/h)*f(2);
            
        end 
        
    end 
    
end 