% WRITE IN FILE FUNCTION CLASS
classdef file_writing_functions
    
    methods(Static)
        
        % WRITE PROPAGATION RESULTS IN FILE 
        %
        % Inputs : 
        %   case_name (string) - case name given by user
        %
        %   prop_res - matrix containing the values taken by the spacrecraft's 
        %                   properties at each instant during propagation. 
        %                   Each row corresponds to a particular instant :
        %                   tab(:, 1) -> time increment
        %                   tab(:, 2:4) -> inertial radius
        %                   tab(:, 5:7) -> inertial velocity
        %                   tab(:, 8:10) -> radius in body-fixed frame
        %                   tab(:, 11:13) -> velocity in body-fixed frame
        %                   tab(:, 14) -> planet-centric radius
        %                   tab(:, 15) -> latitude
        %                   tab(:, 16) -> longitude
        %                   tab(:, 17) -> main body's right ascension
        %                   tab(:, 18) -> main body's declination
        %                   tab(:, 19) -> main body's sidereal time
        %                   tab(:, 20:22) -> gravitational perturbation
        %                   tab(:, 23:25) -> sol. rad. pressure 
        %                   tab(:, 26:28) -> third body gravity
        %                   tab(:, 29) -> inside/outside umbra index
        %                   tab(:, 30:35) -> orbital elements
        %                   tab(:, 36) -> RK4 computation error on orb. elems
        function PropData2File(case_name, prop_res)
            
            fileID = fopen(strcat('prop_data_', case_name, '.txt'), 'w');
            
            fprintf(fileID, '%-15s \n', strcat(case_name, ...
                ' NavTOP PROPAGATION RESULTS'));
            
            fprintf(fileID, '%-150s \n', ' ');
            
            fprintf(fileID, ...
                '%-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t\n', ...
                'Time (s)', 'x (ECI, m)', 'y (ECI, m)', 'z (ECI, m)', ...
                'vx (ECI, m/s)', 'vy (ECI, m/s)', 'vz (ECI, m/s)', ...
                'x (ECEF, m)', 'y (ECEF, m)', 'z (ECEF, m)', ...
                'vx (ECEF, m/s)', 'vy (ECEF, m/s)', 'vz (ECEF, m/s)', ...
                'Rad (m)', 'Lat (rad)', 'Lon (rad)', 'RA (rad)', ...
                'Dec (rad)', 'STime (rad)', 'Fg1 (m/s2)', 'Fg2 (m/s2)', ...
                'Fg3 (m/s2)', 'Fsrp1 (m/s2)', 'Fsrp2 (m/s2)', ...
                'Fsrp3 (m/s2)', 'Ftbg1 (m/s2)', ...
                'Ftbg2 (m/s2)', 'Ftbg3 (m/s2)', 'UIndex', ...
                'Inc (rad)', 'AscNode (rad)', 'ArgPer (rad)', ...
                'TrAnom (rad)', 'SMajAxis (m)', 'Ecc', 'Err');
            fprintf(fileID, ...
                '%-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t %-22s \t \t\n', ...
                '----------------------', '----------------------', '----------------------', ...
                '----------------------', '----------------------', '----------------------', ...
                '----------------------', '----------------------', '----------------------', ...
                '----------------------', '----------------------', '----------------------', ...
                '----------------------', '----------------------', '----------------------', ...
                '----------------------', '----------------------', '----------------------', ...
                '----------------------', '----------------------', '----------------------', ...
                '----------------------', '----------------------', '----------------------', ...
                '----------------------', '----------------------', '----------------------', ...
                '----------------------', '----------------------', '----------------------', ...
                '----------------------', '----------------------', '----------------------', ...
                '----------------------', '----------------------', '----------------------');
            fprintf(fileID, '%+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e \t \t %+.15e\n', prop_res');
            fclose(fileID);
                     
        end 
        
    end
    
end