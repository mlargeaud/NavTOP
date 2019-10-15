classdef rotation_functions
    
    methods(Static)
        
        % ROTATE 3D VECTOR FIELD FROM BODY FIXED FRAME TO SPHERICAL FRAME
        % 
        % Inputs : 
        %   lats (scalar, radians) - vector of size M containing the 
        %                   latitude coordinates of the locations associated 
        %                   with the input vector field
        %
        %   lons (scalar, radians) - vector of size N containing the longitude 
        %                   coordinates of the locations associated with
        %                   the input vector field
        %
        %   BFVecField (M by N by 3 matrix, any unit) -  contains the 
        %                   vectors to rotate
        %                   BFVecField(:,:,1) -> X comp. in body fixed frame
        %                   BFVecField(:,:,2) -> Y comp. in body fixed frame
        %                   BFVecField(:,:,3) -> Z comp. in body fixed frame
        %
        % Outputs :
        %   SvecField (M by N by 3 matrix, any unit) - contains the 
        %                   rotated vectors, now expressed in spherical
        %                   frame
        %                   SVecField(:,:,1) -> southward component
        %                   SVecField(:,:,2) -> eastward component
        %                   SVecField(:,:,3) -> radial component
        function SVecField = BF2Spherical(lats, lons, BFVecField)
            
            M = size(BFVecField, 1);
            N = size(BFVecField, 2);
            
            SVecField = zeros(M, N, 3);
            
            for i = 1:M
                
                for j = 1:N
                    f = [BFVecField(i, j, 1); BFVecField(i, j, 2); BFVecField(i, j, 3)];
                    R = rotation_functions.compute_rot_matrix(lons(j), ...
                        lats(i), 0, '323');
                    SVecField(i, j, :) = R*f;
                    
                end 
                
            end 
            
        end 
        

        
        % COMPUTE SPHERICAL COORDINATES FROM POSITION VECTOR EXPRESSED IN LOCAL FRAME
        % 
        % Inputs : 
        %   r_asc (scalar, radians) - right ascension of orbited body with 
        %                   respect to inertial frame
        %
        %   dec (scalar, radians) - declination of orbited body with respect 
        %                   to inertial frame
        %
        %   s_time (scalar, radians) - sidereal time of orbited body with 
        %                   respect to inertial frame
        %
        %   ecc (scalar, unitless) - eccentricity of spacecraft's orbit
        %
        %   semi_maj_axis (scalar, meters) - semi-major axis of spacecraft's 
        %                   orbit
        %
        %   lon_asc_node (scalar, radians) - longitude of ascending node of 
        %                   the spacecraft's orbit
        %
        %   inc (scalar, radians) - inclination of spacecraft's orbit
        %
        %   arg_per (scalar, radians) - argument of periapsis of spacecraft's 
        %                   orbit
        %
        %   tr_anom (scalar, radians) - spacecraft's true anomaly
        %
        % Outputs : 
        %   planetary_coord (vector) - contains the spherical coordinates
        %                   of the spacecraft, with respect to the body
        %                   fixed frame.
        %                   planetary_coord(1) -> planet-centric radius (meters)
        %                   planetary_coord(2) -> latitude (radians)
        %                   planetary_coord(3) -> longitude (radians)
        function planetary_coord = LocalToPlanetaryCoord(r_asc, dec, s_time, ...
            ecc, semi_maj_axis, lon_asc_node, inc, arg_per, tr_anom)

            planetary_coord = zeros(3, 1);

            r = semi_maj_axis*(1 - ecc^2)/(1 + ecc*cos(tr_anom));

            % compute rotation matrix from planetary frame to inertial frame
            R_PI = rotation_functions.compute_rot_matrix(r_asc, dec, s_time, '313');

            % compute rotation matrix from local orbit frame to inertial frame
            R_IO = rotation_functions.compute_rot_matrix(lon_asc_node, inc, ...
                arg_per + tr_anom, '313');

            % total rotation matrix from local orbit to planetary frame
            R_tot = transpose(R_PI)*R_IO;

            % position vector in planetary frame
            pos = R_tot*[r; 0.; 0.];

            % setting values close to the round-off error to 0
            for i = 1:3

                if (abs(pos(i)) < 10^(-9))

                    pos(i) = 0;

                end 

            end

            % planet-centric radius
            planetary_coord(1) = norm(pos);

            % latitude
            planetary_coord(2) = acos(pos(3)/planetary_coord(1));

            % longitude
            planetary_coord(3) = atan(pos(2)/pos(1));

            if (pos(1) < 0 && pos(2) < 0)
                planetary_coord(3) = atan(pos(2)/pos(1)) + pi;
            end 

            if (pos(1) > 0 && pos(2) < 0)
                planetary_coord(3) = 2*pi + atan(pos(2)/pos(1));
            end

            if (pos(1) < 0 && pos(2) > 0)
                planetary_coord(3) = pi + atan(pos(2)/pos(1));
            end 

            if (pos(1) > 0 && pos(2) > 0)
                planetary_coord(3) = atan(pos(2)/pos(1));
            end     

            if (pos(1) == 0 && pos(2) == 0)
                planetary_coord(3) = lon_asc_node;
            end 

            % if s/c passes through lon = 180 deg
            if (pos(2) == 0 && pos(1) < 0)
                planetary_coord(3) = pi;
            end 

            % particular case where inc = +-90 deg
            if (inc == pi/2 || inc == -pi/2)

                % if the s/c is past one pole
                if (tr_anom > pi/2 - arg_per && tr_anom < 3*(pi/2) - arg_per)

                    planetary_coord(3) = lon_asc_node + pi;

                end          

                % if the ascending node is 180 deg
                if (lon_asc_node == pi)

                    planetary_coord(3) = lon_asc_node;

                    % if the s/c is passed one pole
                    if (tr_anom > pi/2 - arg_per && tr_anom < 3*(pi/2) - arg_per)

                        planetary_coord(3) = lon_asc_node + pi;

                    end 

                end 

            end 

            planetary_coord(3) = mod(planetary_coord(3), 2*pi);

        end         
        
        
        
        % COMPUTE 3-1-3 OR 3-2-3 ROTATION MATRIX
        %
        % Inputs :
        %   a1, a2, a3 (scalars, radians) - first, second, and third rotation 
        %                   angles respectively
        %
        %   str (string) - parameter defining which rotation to perform :
        %                   str = '313' -> 3-1-3 rotation
        %                   str = '323' -> 3-2-3 rotation
        %
        % Outputs : 
        %   R (3 by 3 matrix) - rotation matrix 
        function R = compute_rot_matrix(a1, a2, a3, str)

            R = zeros(3, 3);

            if (str == '313') %#ok<BDSCA>
                
                R31 = [cos(a1) -sin(a1) 0; sin(a1) cos(a1) 0; 0 0 1];
                R12 = [1 0 0; 0 cos(a2) -sin(a2); 0 sin(a2) cos(a2)];
                R33 = [cos(a3) -sin(a3) 0; sin(a3) cos(a3) 0; 0 0 1];
                
                R = R31*R12*R33;

            end

            if (str == '323') %#ok<BDSCA>
                
                R31 = [cos(a1) -sin(a1) 0; sin(a1) cos(a1) 0; 0 0 1];
                R22 = [cos(a2) 0 sin(a2); 0 1 0; -sin(a2) 0 cos(a2)];
                R33 = [cos(a3) -sin(a3) 0; sin(a3) cos(a3) 0; 0 0 1];
                
                R = R31*R22*R33;
             
            end

        end 
        
    end 
    
end 