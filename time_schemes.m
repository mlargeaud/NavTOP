classdef time_schemes
    
    methods(Static)
        
        % SOLVING 1ST ORDER ODE USING RUNGE-KUTTA 4 METHOD 
        %
        % Inputs :
        %   dt (scalar, seconds) - time step
        %
        %   fun (function handle) - right-hand side function f in y' = f(t,y)
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
        %   y - state vector at initial time ti of size N
        %
        % Outputs :
        %   y_next - newly computed state vector at time ti+dt
        %
        %   Fg (vector, meters/second squared) - gravitational perturbation
        %
        %   Fsrp (vector, meters/second squared) - solar radiation pressure
        %
        %   Ftb (vector, meters/second squared) - third body gravity
        function [y_next, Fg, Fsrp, Ftb] = RK4(dt, fun, mu, r0, rasc, ...
                dec, stime, OrbElemsMB, harms, max_deg, SCONE2I, hCone, ...
                ConeApert, mu_tb, A, C, m, y)
           
            % compute right-hand side at y and collecting perturbing 
            % acceleration data at y
            [RHSy, Fg, Fsrp, Ftb] = feval(fun, mu, r0, rasc, dec, stime, ...
                OrbElemsMB, harms, max_deg, SCONE2I, hCone, ConeApert, ...
                mu_tb, A, C, m, y);
           
            k1 = dt*RHSy;
            
            k2 = dt*feval(fun, mu, r0, rasc, dec, stime, OrbElemsMB, ...
                harms, max_deg, SCONE2I, hCone, ConeApert, mu_tb, A, C, ...
                m, y + (k1/2));
            
            k3 = dt*feval(fun, mu, r0, rasc, dec, stime, OrbElemsMB, ...
                harms, max_deg, SCONE2I, hCone, ConeApert, mu_tb, A, C, ...
                m, y + (k2/2));
            
            k4 = dt*feval(fun, mu, r0, rasc, dec, stime, OrbElemsMB, ...
                harms, max_deg, SCONE2I, hCone, ConeApert, mu_tb, A, C, ...
                m, y + k3);
    
            y_next = y + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
            
            % angles modulo 2*pi
            y_next(1:4) = mod(y_next(1:4), 2*pi);
                       
        end
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
%         % SOLVING 1ST ORDER ODE USING RUNGE-KUTTA BUTCHER METHOD
%         %
%         % Inputs :
%         %   y - state vector at initial time ti of size N
%         %
%         %   dt - time step
%         %
%         %   fun - function handle corresponding to function f in y' = f(t, y)
%         %
%         %   mu (meters cubed/seconds squared) - gravitational parameter of
%         %                   the orbited body
%         %   
%         %   r0 (meters) - mean radius of orbited body
%         %
%         %   rasc (radians) - right ascension of orbited body with respect
%         %                   to inertial frame
%         %
%         %   dec (radians) - declination of orbited body with respect
%         %                   to inertial frame
%         %
%         %   stime (radians) - sidereal time of orbited body with respect to
%         %                   inertial frame
%         %
%         %   harms - matrix containing the spherical harmonic coefficients
%         %                   associated with the orbited body. Use the
%         %                   function force_functions.read_harm_coeff_from_file 
%         %                   with the raw data from space agencies' servers 
%         %                   to generate this matrix
%         %
%         % Outputs :
%         %   y_next - newly computed state vector at time ti+dt
%         function y_next = RKB(y, dt, fun, mu, f)
%             
%             k1 = dt*feval(fun, y, mu, f);
%             k2 = dt*feval(fun, y + (k1/3), mu, f);
%             k3 = dt*feval(fun, y + (2*k2/3), mu, f);
%             k4 = dt*feval(fun, y + (k1/12) + (k2/3) - (k3/12), mu, f);
%             k5 = dt*feval(fun, y - (k1/16) + (9*k2/8) - (3*k3/16) ...
%                 - (3*k4/8), mu, f);
%             k6 = dt*feval(fun, y + (9*k2/8) - (3*k3/8) - (3*k4/4) ...
%                 - (k5/2), mu, f);
%             k7 = dt*feval(fun, y + (9*k1/44) - (9*k2/11) + (63*k3/44) ...
%                 + (18*k4/11) - (16*k6/11), mu, f);
%             
%             y_next = y + (11/120)*k1 + (27/40)*k3 + (27/40)*k4 - (4/15)*k5 ...
%                 - (4/15)*k6 + (11/120)*k7;
%             
%         end 
        
%         %Modified Chebyshev-Picard Iteration Method
%         % y value of the function at time t
%         % dt time step
%         % fun is a function handle from y' = f(t, y)
%         function y_next = MCPI(y, t, ti, tf, dt, fun, mu, f)
%             
%             o1 = (tf + t0)/2;
%             o2 = (tf - t0)/2;
%             
%             % compute tau
%             tau = (t - o1)/o2;
%             
%         end 
%         
%         % test
%         function f = f(t, y)
%             f = [1/t; t];
%         end 
        
    end 
       
end 