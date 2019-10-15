classdef plotting_functions
    
    methods(Static)
        
        % PLOTTING SCALAR FIELD ON MAP 
        % Positions handled automatically, considering field generated on
        % whole map. Can choose to export the plot by specifying a file
        % name.
        %
        % Can be used as :
        %       plot_scalar_map(scalar_field, units, file_name)
        %   OR  plot_scalar_map(scalar_field, units)
        %
        % Inputs :
        %   scalar_field (NP by NT matrix) - scalar field to plot on map
        %                   NP -> # of latitudes
        %                   NT -> # of longitudes
        %
        %   units (string) - to display on plot giving the units of the field
        % 
        %   file_name (string, optional) - name of the output file containing 
        %                   the plot (exported in .png format)
        function plot_scalar_map(varargin)
            
            if nargin == 3
                
                scalar_field = varargin{1};
                units = varargin{2};
                file_name = varargin{3};  
                addpath ...
                ('/Users/maxime/Documents/ORBs project/Dynamic model/v79/export_fig');
                                
            else 
                
                scalar_field = varargin{1};
                units = varargin{2};
                
            end 

            load coastlines %#ok<LOAD>

            figure
            
            axesm ('eckert4', 'Frame', 'on', 'Grid', 'on');
            R = georasterref('RasterSize', size(scalar_field), ...
                'RasterInterpretation', 'cells', ...
                'ColumnsStartFrom', 'north', ...
                'Latlim',[-90 90], 'Lonlim', [0 360]);

            geoshow(scalar_field, R, 'DisplayType','texturemap')
            geoshow(coastlat, coastlon)

            sc = colorbar('southoutside');
            sc.Label.String = units;
            sc.FontSize = 12;
            sc.FontWeight = 'bold';
            
            if nargin == 3
                export_fig(file_name);
            end

        end
        
        
        
        % PLOTTING SCALAR FIELD ON MAP WITH SPECIFIED POSITIONS
        %
        % Can be used as :
        %       plot_scalar_map(positions, scalar_field, units, file_name)
        %   OR  plot_scalar_map(positions, scalar_field, units)
        % 
        % Inputs :
        %   positions (NP by NT by 2 matrix) - contains the coordinates 
        %                   of points considered 
        %                   NP # of latitudes
        %                   NT # of longitudes
        %                   positions(:,:,1) -> latitudes in range [0°, 180°]
        %                   positions(:,:,2) -> longitudes in range [0°, 360°] 
        %
        %   scalar_field (NP by NT by 2 matrix) - contains the values
        %                   to plot
        %
        %   units (string) - string to display on plot giving the units of 
        %                   the field
        %
        %   file_name (string, optional) - name of the output file containing 
        %                   the plot (exported in .png format)
        function plot_scalar_map2(varargin)
            
            if nargin == 4
                
                positions = varargin{1};
                scalar_field = varargin{2};
                units = varargin{3};
                file_name = varargin{4};  
                addpath ...
                ('/Users/maxime/Documents/ORBs project/Dynamic model/v81/export_fig');
                                
            else 
                
                positions = varargin{1};
                scalar_field = varargin{2};
                units = varargin{3};
                
            end 

            load coastlines %#ok<LOAD>
            
            figure
            axesm ('eckert4', 'Frame', 'on', 'Grid', 'on');
            
            geoshow(90 - positions(:, :, 1), positions(:, :, 2), ...
                scalar_field, 'DisplayType', 'texturemap')
            geoshow(coastlat, coastlon)
            
            sc = colorbar('southoutside');         
            sc.Label.String = units;
            sc.FontSize = 12;
            sc.FontWeight = 'bold';
            
            if nargin == 4
                export_fig(file_name);
            end

        end
        
        
        
        % PLOTTING VECTOR ON MAP
        %
        % Inputs :
        % X, Y, Z (NP by NT matrices, any unit) - respectively the southward, 
        %                   eastward and downward component matrices of the 
        %                   force
        %                   NP # of latitudes
        %                   NT # of longitudes
        %
        % r (scalar, meters) - altitude at which the field needs to be plotted
        function plot_grav_force_map(X, Y, Z, r)

            lat_step = 180/size(X, 1);
            lon_step = 360/size(X, 2);

            nb_lat = size(X, 1);
            nb_lon = size(X, 2);

            lats = zeros(nb_lat*nb_lon, 1);
            lons = zeros(nb_lat*nb_lon, 1);
            alts = zeros(nb_lat*nb_lon, 1);

            for i = 1:nb_lon

                for j = 1:nb_lat

                    lats((i-1)*nb_lat + j) = -90 + (j-1)*lat_step;
                    lons((i-1)*nb_lat + j) = (i-1)*lon_step; 
                    alts((i-1)*nb_lat + j) = r;

                end

            end

            X_q = reshape(X, [nb_lat*nb_lon, 1]);
            Y_q = reshape(Y, [nb_lat*nb_lon, 1]);
            Z_q = reshape(Z, [nb_lat*nb_lon, 1]);

            figure
            load coastlines %#ok<LOAD>
            axesm ('eckert4', 'Frame', 'on', 'Grid', 'on')
            plotm(coastlat,coastlon)
            quiver3m(lats, lons, alts, X_q, Y_q, Z_q)

        end
        
        
        
        % PLOTTING SUCCESSIVE ORBITS OF A SPACECRAFT ON MAP
        % Plots orbits with custom description 
        %
        % Inputs : 
        %   positions (NO by NPos by 2 matrix) - contains the successive 
        %                   coordinates of the s/c (each row of 'positions' 
        %                   describes a full orbit).
        %                   NO # of orbits 
        %                   Npos # of positions on each orbit coordinates 
        %                   positions(:,:,1) -> latitudes in range [0°, 180°]
        %                   positions(:,:,2) -> longitudes in range [0°, 360°] 
        % 
        %   step (scalar) - number of orbits to skip between each orbit plotted
        %   param (vector of size NO) - contains values to display on the orbits
        %   v (string) - gives the name of the values in param to display 
        function plot_orbit_on_map_earth(positions, step, param, v)

            figure
            load coastlines %#ok<LOAD>
            axesm('miller', 'Frame', 'on', 'Grid', 'on', 'MeridianLabel', 'on', ...
                'ParallelLabel', 'on')
            geoshow(coastlat, coastlon, 'Color', [0.45, 0.45, 0.45])

            if (size(positions, 1) == 1)

                [lat_max, id_lon] = max(positions(1, :, 2));

                scatterm(90 - positions(1, 1:size(positions, 2), 2), ...
                    positions(1, 1:size(positions, 2), 3), 10, [0 0 0], ...
                    'o', 'filled');

                textm(90 - lat_max, positions(1, id_lon, 2), ...
                    [v, ' = ', num2str(param(1)), '°'], ...
                    'FontSize', 14, 'VerticalAlignment', 'top', ...
                    'FontName', 'Times', 'FontWeight', 'Bold')

                return

            end 

            count = 0;
            for i = 1:step:(size(positions, 1))

                if (rem(count, 2) == 0)

                    scatterm(90 - positions(i, 1:size(positions, 2), 2), ...
                        positions(i, 1:size(positions, 2), 3), 15, [0 0 0], 'x');

                else

                    scatterm(90 - positions(i, 1:size(positions, 2), 2), ...
                        positions(i, 1:size(positions, 2), 3), 10, [0 0 0],  ...
                        'o', 'filled');

                end 

                [lat_max, id_lon] = max(positions(i, :, 2));

                % if all the latitudes are equal (inc = 0 deg)
                if (range(positions(i, :, 2)) < 10^(-9))
                    
                    id_lon = fix(size(positions, 2)/4);
                  
                end 

                % if polar orbit
                if (fix(range(positions(i, :, 3))) == 180)

                    lat_max = 170;

                end 

                textm(90 - lat_max, positions(i, id_lon, 3), ...
                    [v, ' = ', num2str(param(i)), '°'], ...
                    'FontSize', 12, 'VerticalAlignment', 'top', ...
                    'FontName', 'Times', 'FontWeight', 'Bold')

                hold on

                count = count + 1;

            end

        end 
        
        
        
        % PLOT SUCCESSIVE POSITIONS WITH ASSOCIATED SCALAR DATA
        %
        % Inputs :
        %   lat (vector of size N, degrees) - contains latitudes in [0; 90]° 
        %                   range
        %
        %   lon (vector of size N, degrees) - vector of size N containing 
        %                   longitudes in [0; 360]° range
        %
        %   r_asc (scalar, radians) - right ascension of orbited body with 
        %                   respect to inertial frame
        %
        %   dec (scalar, radians) - declination of orbited body with respect 
        %                   to inertial frame
        %
        %   s_time (scalar, radians) - sidereal time of orbited body with 
        %                   respect to inertial frame
        %
        %   orbelems (vector) contains the orbital elements describing 
        %                   the orbit of the spacecraft in this order : 
        %                   (1) inclination (scalar, radians)
        %                   (2) longitude of ascending node (scalar, radians)
        %                   (3) argument of periapsis (scalar, radians)
        %                   (4) true anomaly (scalar, radians)
        %                   (5) semi-major axis (scalar, meters)
        %                   (6) eccentricity (scalar, unitless)
        %
        %   units (string) - string to display on plot giving the units of 
        %                   the field
        %
        %   data (vector of size N) - contains scalar values associated 
        %                   with the coordinates from lat & lon
        %
        %   cond (vector of size N) - contains values indicating if the
        %                   positions to plot are inside or outside the 
        %                   body's umbra
        %                   0 -> outside the umbra
        %                   1 -> inside the umbra
        %
        %   figtitle (string) - title associated with data plots
        %
        %   exp (string) - name of the export file if needed
        %                   To prevent export -> set exp as '~'
        %                   To export -> enter file name
        function plot_wdata(lat, lon, r_asc, dec, s_time, orbelems, units, ...
                data, cond, figtitle, exp)
              
            load coastlines %#ok<LOAD>
            
            i = orbelems(1);
            O = orbelems(2);
            o = orbelems(3);
            a = orbelems(5);
            e = orbelems(6); 
            
            % sunlit orbital positions and time increments
            sunlitPos = [lat(cond == 0), lon(cond == 0)];
            sunlitTimes = data(cond == 0);
            
            % orbital positions in umbra and associated time increments
            InUmbraPos = [lat(cond == 1), lon(cond == 1)];
            InUmbraTimes = data(cond == 1);
            
            % Periapsis coordinates at final time
            per = rotation_functions.LocalToPlanetaryCoord...
                (r_asc, dec, s_time, e, a, O, i, o, 0);
            
            % Ascending node coordinates at final time
            nde = rotation_functions.LocalToPlanetaryCoord...
                (r_asc, dec, s_time, e, a, O, i, 0, 0);
            
            % Start point coordinates
            startPos = [lat(1); lon(1)];
            
            % End point coordinates
            endPos = [lat(size(lat,1)); lon(size(lon,1))];
            
            % Generating data to plot orbit at last time increment
            nu = 0;
            PresentOrb = zeros(360, 3);
            for index = 1:360
                
                PresentOrb(index, :) = rotation_functions.LocalToPlanetaryCoord...
                (r_asc, dec, s_time, e, a, O, i, o, nu*pi/180); 
            
                nu = nu + 1;
                
            end 

            figure 

            axesm('miller', 'Frame', 'on', 'Grid', 'on', ...
                'MeridianLabel', 'on', 'MlabelParallel', 'north', ...
                'MlabelLocation', [-180, -120, 120, 180], 'ParallelLabel', 'on')

            % plot landmarks
            geoshow(coastlat, coastlon, 'Color', ...
                [0.45, 0.45, 0.45], 'HandleVisibility', 'off')

            % plot sunlit orbital positions
            scatterm(90*ones(size(sunlitPos, 1), 1) - sunlitPos(:, 1), ...
                sunlitPos(:, 2), 7, sunlitTimes, 'HandleVisibility', 'off')
            
            % plot orbital positions in umbra
            scatterm(90*ones(size(InUmbraPos, 1), 1) - InUmbraPos(:, 1), ...
                InUmbraPos(:, 2), 0.9, InUmbraTimes, 'HandleVisibility', 'off')

            % plot orbit at last time increment
            orb = plotm(90*ones(360, 1) - PresentOrb(:, 2)*180/pi, ...
                PresentOrb(:, 3)*180/pi, '-k', 'LineWidth', 1.9);

            % plot periapsis location
            p = plotm(90 - per(2)*180/pi, per(3)*180/pi, '--g^',...
            'LineWidth', 2, 'MarkerSize', 10, 'MarkerEdgeColor', ...
            'k', 'MarkerFaceColor', 'r');

            % plot apoapsis location
            q = plotm(-90 + per(2)*180/pi, per(3)*180/pi + 180, ...
            '--g^', 'LineWidth', 2, 'MarkerSize', 10, ...
            'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g');

            % plot ascending node location
            r = plotm(-90 + nde(2)*180/pi, nde(3)*180/pi, '--gs',...
            'LineWidth', 2, 'MarkerSize', 10, 'MarkerEdgeColor', ...
            'k', 'MarkerFaceColor', 'w');

            % plot start point location
            s = plotm(90 - startPos(1), startPos(2), '*',...
            'LineWidth', 2, 'MarkerSize', 10, 'MarkerEdgeColor', ...
            'k', 'MarkerFaceColor', 'k');

            % plot end point location
            t = plotm(90 - endPos(1), endPos(2), 'x',...
            'LineWidth', 2, 'MarkerSize', 10, 'MarkerEdgeColor', ...
            'k', 'MarkerFaceColor', 'k');

            legend([orb, p, q, r, s, t], 'Final orbit', 'Periapsis', ...
                'Apoapsis', 'Ascending node', 'Start point', ...
                'End point', 'FontSize', 13, 'NumColumns', 3)
            
            textm(-85, 0, 'Thiner line : positions in umbra', 'FontSize', ...
                13, 'HorizontalAlignment', 'Center');

            sc = colorbar('southoutside');
            sc.Label.String = units;
            sc.FontSize = 12;
            sc.FontWeight = 'bold';

            title(figtitle, 'FontSize', 14);

            % if file export is required  
            if exp(1) ~= '~' 
                addpath ...
                ('/Users/maxime/Documents/ORBs project/Dynamic model/v82/export_fig');
                export_fig(exp);
            end
            
        end
        
        
        
        % PLOT ORBITAL ELEMENTS CURVES AFTER PROPAGATION
        %
        % Inputs : 
        %   times (vector of size N, hours) - contains the time increments
        %                   for which orbital elements have been computed
        %   
        %   tabOrbElems (N by 6 matrix) - contains the orbital elements
        %                   values through propagation :
        %                   (:, 1) -> inclinations (scalar, radians)
        %                   (:, 2) -> longitudes of ascending node (scalar, radians)
        %                   (:, 3) -> arguments of periapsis (scalar, radians)
        %                   (:, 4) -> true anomalies (scalar, radians)
        %                   (:, 5) -> semi-major axes (scalar, meters)
        %                   (:, 6) -> eccentricities (scalar, unitless)
        %
        %   tabAccs (N by 9 matrix) - contains perturbing acceleration
        %                   vectors values through propagation
        %                   (:, 1:3) -> gravitational perturbation
        %                   (:, 4,6) -> solar radiation pressure
        %                   (:, 7:9) -> third body gravity
        %
        % Outputs :
        %   fig1, fig2 (figure handles) - gives the figure handles of orbital 
        %                   elements and perturbing accelerations plots  
        function [fig1, fig2] = plotPropRes(times, tabOrbElems, tabAccs)
            
            is = tabOrbElems(:, 1);
            Os = tabOrbElems(:, 2);
            os = tabOrbElems(:, 3);
            ns = tabOrbElems(:, 4);
            as = tabOrbElems(:, 5);
            es = tabOrbElems(:, 6);
            
            Fgs = zeros(size(times, 1), 1);
            Fsrps = zeros(size(times, 1), 1);
            Ftbs = zeros(size(times, 1), 1);
            
            for i = 1:size(times, 1)
                
                Fgs(i) = norm(tabAccs(i, 1:3));
                Fsrps(i) = norm(tabAccs(i, 4:6));
                Ftbs(i) = norm(tabAccs(i, 7:9));
                
            end 
            
            % orbital elements evolution
            figure
            
            % inclination plot
            subplot(3, 2, 1)
            plot(times, is*180/pi)
            grid on
            title('Inclination (degrees)', 'FontSize', 12)
 
            % longitude of ascending node plot
            subplot(3, 2, 2)
            plot(times, Os*180/pi)
            grid on
            title('Longitude of asc. node (degrees)', 'FontSize', 12)
  
            % argument of periapsis plot
            subplot(3, 2, 3)
            plot(times, os*180/pi)
            grid on
            title('Argument of periapsis (degrees)', 'FontSize', 12)
            
            % true anomaly plot
            subplot(3, 2, 4)
            plot(times, ns*180/pi)
            grid on
            title('True anomaly (degrees)', 'FontSize', 12)
            
            % semi-major axis plot
            subplot(3, 2, 5)
            plot(times, as)
            grid on
            title('Semi-major axis (m)', 'FontSize', 12)
            xlabel('Time (hrs)')
            
            % eccentricity plot
            subplot(3, 2, 6)
            plot(times, es)
            grid on
            title('Eccentricity', 'FontSize', 12)
            xlabel('Time (hrs)')
            
            fig1 = gcf;
            
            % perturbing accelerations plot
            figure
                        
            % gravitational perturbation plot
            subplot(3, 2, [1,2])
            plot(times, Fgs); 
            grid on
            title('Gravitational perturbation (m/s^2)', 'FontSize', 12)
            
            % solar radiation pressure plot
            subplot(3, 2, [3,4])
            plot(times, Fsrps);
            grid on
            title('Solar radiation pressure perturbation (m/s^2)', 'FontSize', 12)
            
            % third body gravity plot
            subplot(3, 2, [5,6])
            semilogy(times, Ftbs);
            grid on
            title('Third body gravitational perturbation from Sun (m/s^2)', 'FontSize', 12)
            xlabel('Time (hrs)')
            
            fig2 = gcf;
        
        end 
        
    end  
    
end 