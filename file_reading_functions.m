% READ FROM FILE FUNCTIONS CLASS
classdef file_reading_functions
    
    methods(Static)
        
        % READING PROPAGATION PARAMETERS AND INITIAL CONDITIONS FROM INPUT FILE
        %
        % Inputs :
        %   filename (string) - input file 
        %
        % Outputs :
        %   case_name (string) - case name given by user
        %
        %   main_body (main_body object) - main gravitating body set by user 
        %
        %   sc_mass (scalar, kilograms) - spacecraft's mass
        %
        %   sc_refl_area (scalar, meters squared) - spacecraft's reflective
        %                   area
        %
        %   sc_refl_coeff (scalar, unitless) - spacecraft's reflective
        %                   coefficient
        %
        %   rasc, dec, stime (scalars, radians) - angles corresponding 
        %                   to body-fixed frame 3-1-3 rotation w.r.t. 
        %                   inertial frame
        %
        %   mb_init_OrbElems - vector containing the initial orbital 
        %                   elements describing the main body's orbit  
        %                   in this order : 
        %                   (1) inclination (scalar, radians)
        %                   (2) longitude of ascending node (scalar, radians)
        %                   (3) argument of periapsis (scalar, radians)
        %                   (4) true anomaly (scalar, radians)
        %                   (5) semi-major axis (scalar, meters)
        %                   (6) eccentricity (scalar, unitless)
        %
        %   sc_init_OrbElems - vector containing the initial orbital 
        %                   elements describing the orbit of the spacecraft 
        %                   in this order : 
        %                   (1) inclination (scalar, radians)
        %                   (2) longitude of ascending node (scalar, radians)
        %                   (3) argument of periapsis (scalar, radians)
        %                   (4) true anomaly (scalar, radians)
        %                   (5) semi-major axis (scalar, meters)
        %                   (6) eccentricity (scalar, unitless)
        %
        %   max_deg (integer) - gravitational potential maximum degree 
        %
        %   time_step (scalar, seconds) - propagation time step
        %
        %   final_time (scalar, seconds) - final time increment
        function [case_name, main_body, sc_mass, sc_refl_area, sc_refl_coeff, ...
                rasc, dec, stime, mb_init_OrbElems, sc_init_OrbElems, ...
                max_deg, time_step, final_time] = read_inputs(filename)
            
            fileID = fopen(filename);
            textLine = fgetl(fileID);
            while ischar(textLine)
                
                % find case name
                nameLoc = strfind(textLine, 'Case name');
                if nameLoc == 1
                    case_name = sscanf(textLine, 'Case name %100c');
                end
                
                % find main body
                mbLoc = strfind(textLine, 'Main gravitating body');
                if mbLoc == 1
                    main_body = sscanf(textLine, 'Main gravitating body %100c');
                end
                
                % find spacecraft mass
                scMassLoc = strfind(textLine, 'Spacecraft mass (kg)');
                if scMassLoc == 1
                    sc_mass = sscanf(textLine, 'Spacecraft mass (kg) %f');
                end
                
                % find spacecraft reflective area
                scReflAreaLoc = strfind ...
                    (textLine, 'Spacecraft reflective area (m2)');
                if scReflAreaLoc == 1
                    sc_refl_area = sscanf ...
                        (textLine, 'Spacecraft reflective area (m2) %f');
                end
                
                % find spacecraft reflective coefficient
                scReflCoeffLoc = strfind ...
                    (textLine, 'Spacecraft reflective coefficient');
                if scReflCoeffLoc == 1
                    sc_refl_coeff = sscanf ...
                        (textLine, 'Spacecraft reflective coefficient %f');
                end
                
                % find main body's initial orientation
                orientLoc = strfind ...
                    (textLine, 'Right ascension (deg)');
                if orientLoc == 1
                    rasc = (pi/180)*sscanf ...
                        (textLine, 'Right ascension (deg)  %f');
                    textLine = fgetl(fileID);
                    dec = (pi/180)*sscanf ...
                        (textLine, 'Declination (deg)  %f');
                    textLine = fgetl(fileID);
                    stime = (pi/180)*sscanf ...
                        (textLine, 'Sidereal time (deg)  %f');
                end
                
                % find main body's initial orbital elements
                mbInitOrbElems = strfind ...
                    (textLine, 'Main body''s orbital inclination (deg)');
                if mbInitOrbElems == 1
                    mb_init_OrbElems(1, 1) = (pi/180)*sscanf ...
                        (textLine, 'Main body''s orbital inclination (deg) %f');
                    textLine = fgetl(fileID);
                    mb_init_OrbElems(2, 1) = (pi/180)*sscanf ...
                        (textLine, 'Main body''s longitude of ascending node (deg) %f');
                    textLine = fgetl(fileID);
                    mb_init_OrbElems(3, 1) = (pi/180)*sscanf ...
                        (textLine, 'Main body''s argument of periapsis (deg) %f');
                    textLine = fgetl(fileID);
                    mb_init_OrbElems(4, 1) = (pi/180)*sscanf ...
                        (textLine, 'Main body''s true anomaly (deg) %f');
                    textLine = fgetl(fileID);
                    mb_init_OrbElems(5, 1) = sscanf ...
                        (textLine, 'Main body''s semi-major axis (m) %f');
                    textLine = fgetl(fileID);
                    mb_init_OrbElems(6, 1) = sscanf ...
                        (textLine, 'Main body''s eccentricity %f');
                end        
                
                % find spacecraft's initial orbital elements
                scInitOrbElems = strfind ...
                    (textLine, 'Inclination (deg)');
                if scInitOrbElems == 1
                    sc_init_OrbElems(1, 1) = (pi/180)*sscanf ...
                        (textLine, 'Inclination (deg) %f');
                    textLine = fgetl(fileID);
                    sc_init_OrbElems(2, 1) = (pi/180)*sscanf ...
                        (textLine, 'Longitude of ascending node (deg) %f');
                    textLine = fgetl(fileID);
                    sc_init_OrbElems(3, 1) = (pi/180)*sscanf ...
                        (textLine, 'Argument of periapsis (deg) %f');
                    textLine = fgetl(fileID);
                    sc_init_OrbElems(4, 1) = (pi/180)*sscanf ...
                        (textLine, 'True anomaly (deg) %f');
                    textLine = fgetl(fileID);
                    sc_init_OrbElems(5, 1) = sscanf(textLine, 'Semi-major axis (m) %f');
                    textLine = fgetl(fileID);
                    sc_init_OrbElems(6, 1) = sscanf(textLine, 'Eccentricity %f');
                end
                
                % find potential maximum degree
                maxDegLoc = strfind ...
                    (textLine, 'Potential maximum degree');
                if maxDegLoc == 1
                    max_deg = sscanf(textLine, 'Potential maximum degree %f');
                end
                
                % find propagation time step
                timeStepLoc = strfind ...
                    (textLine, 'Time step (sec)');
                if timeStepLoc == 1
                    time_step = sscanf(textLine, 'Time step (sec) %f');
                end
                
                % find final propagation time increment
                finalTimeLoc = strfind ...
                    (textLine, 'Final time increment (sec)');
                if finalTimeLoc == 1
                    final_time = sscanf(textLine, 'Final time increment (sec) %f');
                end
                
                textLine = fgetl(fileID);
                
            end 
            
            fclose(fileID);
            
        end 
        
        
        % READING MAIN GRAVITATING BODY'S PROPERTIES FROM FILE
        %
        % Inputs :
        %   filename (string) - input file
        %
        % Outputs :
        %   mu (scalar, meters cubed/second squared) - body's gravitational
        %                   parameter
        %
        %   r0 (scalar, meters) - body's mean radius
        %
        %   RotVecBF (vector of size 3, radians/second) - body's rotation 
        %                   vector expressed in body-fixed frame
        %
        %   InMomt (vector of size 3, kilograms.meter squared) - body's
        %                   principal moments of inertia, organized in 
        %                   accordance to the order of the principal axes 
        %                   X1, X2, X3 in PA2BF matrix (described below) :
        %                   InMomt(1) -> along X1 axis 
        %                   InMomt(2) -> along X2 axis 
        %                   InMomt(3) -> along X3 axis
        %
        %   PA2BF (3 by 3 matrix, unitless) - matrix containing the body's 
        %                   principal axes vectors (expressed in body fixed 
        %                   frame) within its columns, organized as a direct 
        %                   reference frame (X1, X2, X3) :
        %                   1st column -> X1
        %                   2nd column -> X2
        %                   3rd column -> X3 
        function [mu, r0, RotVecBF, InMomt, PA2BF] = ...
                read_body_phys_prop(filename)
      
            fileID = fopen(filename);
            textLine = fgetl(fileID);
            while ischar(textLine)
                
                % find gravitational parameter
                gravparamLoc = strfind ...
                    (textLine, 'Gravitational parameter (m3/s2)');
                if gravparamLoc == 1
                    mu = sscanf(textLine, 'Gravitational parameter (m3/s2) %f');
                end
                
                % find mean radius
                radiusLoc = strfind(textLine, 'Mean radius (m)');
                if radiusLoc == 1
                    r0 = sscanf(textLine, 'Mean radius (m) %f');
                end
                
                % find rotation vector
                omegaLoc = strfind(textLine, 'omega1');
                if omegaLoc == 1
                     RotVecBF(1, 1) = sscanf(textLine, 'omega1 %f');
                     textLine = fgetl(fileID);
                     RotVecBF(2, 1) = sscanf(textLine, 'omega2 %f');
                     textLine = fgetl(fileID);
                     RotVecBF(3, 1) = sscanf(textLine, 'omega3 %f');
                end
                
                % find principal moments of inertia
                princMomLoc = strfind(textLine, 'I1');
                if princMomLoc == 1
                     InMomt(1, 1) = sscanf(textLine, 'I1 %f');
                     textLine = fgetl(fileID);
                     InMomt(2, 1) = sscanf(textLine, 'I2 %f');
                     textLine = fgetl(fileID);
                     InMomt(3, 1) = sscanf(textLine, 'I3 %f');
                end
                
                % find rotation matrix from Principal Axes to Body Fixed frame
                PA2BFLoc = strfind(textLine, 'PA2BF11');
                if PA2BFLoc == 1
                     PA2BF(1, 1) = sscanf(textLine, 'PA2BF11 %f');
                     textLine = fgetl(fileID);
                     PA2BF(1, 2) = sscanf(textLine, 'PA2BF12 %f');
                     textLine = fgetl(fileID);
                     PA2BF(1, 3) = sscanf(textLine, 'PA2BF13 %f');
                     textLine = fgetl(fileID);
                     PA2BF(2, 1) = sscanf(textLine, 'PA2BF21 %f');
                     textLine = fgetl(fileID);
                     PA2BF(2, 2) = sscanf(textLine, 'PA2BF22 %f');
                     textLine = fgetl(fileID);
                     PA2BF(2, 3) = sscanf(textLine, 'PA2BF23 %f');
                     textLine = fgetl(fileID);
                     PA2BF(3, 1) = sscanf(textLine, 'PA2BF31 %f');
                     textLine = fgetl(fileID);
                     PA2BF(3, 2) = sscanf(textLine, 'PA2BF32 %f');
                     textLine = fgetl(fileID);
                     PA2BF(3, 3) = sscanf(textLine, 'PA2BF33 %f');
                end
                
                textLine = fgetl(fileID);
                
            end 
            
            fclose(fileID);
            
        end   
        
    end
    
end 