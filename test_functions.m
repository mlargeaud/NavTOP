classdef test_functions
    
    methods(Static)
        
        % PLOT MAIN BODY'S UMBRA CONE LOCATION AT A GIVEN POSITION IN INERTIAL FRAME
        %
        % Can be used as :
        %       PlotShadowCone(mb)
        %   OR  plot_scalar_map(mb, filename) if the plot is to be exported
        %
        % Inputs :
        %   mb (main body object) - celestial body to be studied
        %   
        %   filename (string, optional) - name of export file
        function PlotShadowCone(varargin)
            
            body = varargin{1};
            
            r0 = body.r0;
            hCone = body.hCone;
            
            % domain size around the cone
            length = hCone + 0.1*hCone;
            width = 2*(r0 + 0.1*r0);
            height = 2*(r0 + 0.1*r0);
            
            nb_wth = 101;
            nb_lth = 101;
            nb_hht = 101;
            
            % space steps 
            step_wth = width/(nb_wth-1);
            step_lth = length/(nb_lth-1);
            step_hht = height/(nb_hht-1);
            
            tab_pos_umbra = zeros(nb_wth, nb_lth, nb_hht, 3);
            tab_IsInUmbra = zeros(nb_wth, nb_lth, nb_hht);   

            % initial positions in umbra frame
            pos_wth = -width/2;

            for i = 1:nb_wth

                pos_lth = 0;

                for j = 1:nb_lth

                    pos_hht = -height/2;

                    for k = 1:nb_hht

                        % radius in inertial frame
                        radsc = body.SCONE2I*([pos_hht; pos_wth; pos_lth] ...
                            - hCone*[0;0;1]); 

                        sc = spacecraft(body, 0, radsc, [0;0;0]);            

                        % test whether the position coordinates verify the cone's equation
                        if spacecraft.InUmbra(sc) 

                            tab_IsInUmbra(i, j, k) = 1;

                        end 

                        tab_pos_umbra(i, j, k, :) = radsc;

                        pos_hht = pos_hht + step_hht;

                    end

                    pos_lth = pos_lth + step_lth;

                end 

                pos_wth = pos_wth + step_wth;

            end                 
            
            % plot results
            figure            
                
            % reshape position matrix 
            t1 = reshape(tab_pos_umbra(:, :, :, 1), 1, ...
                nb_wth*nb_lth*nb_hht);

            t2 = reshape(tab_pos_umbra(:, :, :, 2), 1, ...
                nb_wth*nb_lth*nb_hht);

            t3 = reshape(tab_pos_umbra(:, :, :, 3), 1, ...
                nb_wth*nb_lth*nb_hht);

            % reshape "in umbra" indicator matrix
            tu1 = reshape(tab_IsInUmbra(:, :, :), 1, ...
               nb_wth*nb_lth*nb_hht);

            scatter3(0, 0, 0, 'MarkerEdgeColor',[0 0 0],...
                'MarkerFaceColor',[0 0 0],...
                'LineWidth',2);

            hold on

            c = scatter3(t1(tu1 == 1)', t2(tu1 == 1)', ...
                t3(tu1 == 1)', '.');

            hold on                                 
            
            set(gca,'FontSize', 12)
            
            text(0, 0, 0, ...
                    'Earth center', ...
                    'FontSize', 13, 'VerticalAlignment', 'top', ...
                    'HorizontalAlignment', 'center')
                
            xlabel('Inertial X axis (vernal equinox)', 'FontSize', 14)
            ylabel('Inertial Y axis', 'FontSize', 14)
            zlabel('Inertial Z axis', 'FontSize', 14)
            
            legend(c, 'Positions inside the umbra', 'FontSize', 13)
            
            if nargin == 2
                addpath ...
                ('/Users/maxime/Documents/ORBs project/Dynamic model/v108/export_fig');
                export_fig(filename);
            end 
            
        end 
        
                 
        % COMPUTE ROTATION VECTOR OF ROTATING BODY OVER TIME
        % 
        % Inputs : 
        %   mb (main_body object) - celestial body for which the rotation
        %                   vector is computed over time 
        %
        %   dt (scalar, seconds) - time step
        %
        %   tf (scalar, seconds) - final time increment
        %
        % Outputs : 
        %   tab (N-by-13 matrix) - matrix containing output results :
        %                   tab(:, 1)     ->  sidereal time instant (seconds)
        %                   tab(:, 2:4)   ->  rotation vector in principal
        %                                     axes frame (radians/second)
        %                   tab(:, 5:7)   ->  principal frame rotation angles
        %                                     w.r.t. angular momentum-linked 
        %                                     inertial frame (radians)
        %                   tab(:, 8:10)  ->  body fixed rotation angles
        %                                     w.r.t. inertial frame (radians)
        %                   tab(:, 11:13) ->  first body fixed axis expressed 
        %                                     in inertial frame (unitless)
        function tab = RotVecCheck(mb, dt, tf)
            
            body = mb;
            N = (tf/dt) + 1;
            
            tab = zeros(N, 17); 
            t = 0;
            
            % Compute attitude through time
            for i = 1:N
                
                psi = body.psi;
                theta = body.theta;
                phi = body.phi;
                
                rasc = body.rasc;
                dec = body.dec;
                stime = body.stime;
                
                v = body.RotVecPA;
                v = v/norm(v);
                
                posTest = body.BF2I*[1;0;0];

                E = body.EccAnom;
                nu = body.OrbElemsMB(4);
                
                apert = body.ConeApert;
                hCone = body.hCone;
                
                tab(i, :) = [t/(86164.1), v', ...
                    psi, theta, phi, ...
                    rasc, dec, stime, ...
                    posTest', E*180/pi, nu*180/pi, hCone, apert];  

                body = main_body.update_angles(body, dt);
                
                t = t + dt;
                       
            end 
             
            % Plot results
            figure
            
            subplot(2, 1, 1)
            plot(tab(:, 1), tab(:, 2), 'LineWidth', 1);
            grid on
            hold on
            plot(tab(:, 1), tab(:, 3), 'LineWidth', 1);
            grid on
            hold on
            plot(tab(:, 1), sqrt(tab(:, 2).^2 + tab(:, 3).^2), '--', ...
                'LineWidth', 1, 'Color', 'k');
            hold on
            plot(tab(:, 1), -sqrt(tab(:, 2).^2 + tab(:, 3).^2), '--', ...
                'LineWidth', 1, 'Color', 'k');
            grid on
            title('X_1 and X_2 components', 'FontSize', 14);
            
            legend('$X_1$', '$X_2$', '$\pm\sqrt{X_1^2+X_2^2}$', ...
                'FontSize', 12, 'Location', 'southeast', ...
                'Interpreter', 'latex')
            
            subplot(2, 1, 2)
            plot(tab(:, 1), tab(:, 4));
            grid on
            xlabel('Sidereal days', 'FontSize', 12);
            title('X_3 component', 'FontSize', 14);
            
            figure
            
            subplot(3, 1, 1)
            plot(tab(:, 1), tab(:, 5).*180/pi, '.');
            grid on
            title('$\Psi$', 'FontSize', 17, 'Interpreter', 'latex');
            
            subplot(3, 1, 2)
            plot(tab(:, 1), tab(:, 6).*180/pi, '.');
            grid on
            title('$\theta$', 'FontSize', 17, 'Interpreter', 'latex');
            
            subplot(3, 1, 3)
            plot(tab(:, 1), tab(:, 7).*180/pi, '.');
            grid on
            xlabel('Sidereal days', 'FontSize', 12);
            title('$\Phi$', 'FontSize', 17, 'Interpreter', 'latex');          
            
        end 
        
        
        
        % COMPUTE BODY-FIXED COORDINATES OF FIXED POINTS IN INERTIAL FRAME
        % and plot on map for a given time span
        %
        % Inputs :
        %   mb (main body object) - celestial body
        %
        %   dt (scalar, seconds) - time step
        %
        %   tf (scalar, seconds) - final time increment
        %
        %   dnu (scalar, degrees) - true anomaly step between each point 
        %                   to be projected on the map
        %
        % Outputs :
        %   tab (NB_TIMES by 3 by NB_POINTS matrix) - contains the successive 
        %                   positions of the points (fixed in inertial frame) 
        %                   in body-fixed frame as the main body rotates
        %                   tab(:, 1, :) -> successive radii for each point
        %                   tab(:, 2, :) -> successive latitudes for each point
        %                   tab(:, 3, :) -> successive longitudes for each point
        function tab = BodyRotCheck(mb, dt, tf, dnu)
            
            body = mb;
            N = (tf/dt) + 1;
            
            dnu_rad = dnu*pi/180;
            nb_lat = (180/dlat) + 1;
            
            tab = zeros(N, 3, nb_lat);
            
            for i = 1:N
                
                rasc = body.rasc;
                dec = body.dec;
                stime = body.stime;
                
                for j = 1:nb_lat
                    
                    sph_coord = rotation_functions.LocalToPlanetaryCoord ...
                    (rasc, dec, stime, 0, 6371000, 0, pi/2, 0, (pi/2) - (j-1)*dnu_rad);
                
                    tab(i, :, j) = [sph_coord(1), 90 - sph_coord(2)*180/pi, ...
                        sph_coord(3)*180/pi];
                    
                end 
                
                body = main_body.update_angles(body, dt);
                
            end    
            
            % plot results on whole map
            colors = [0, 0.4470, 0.7410; 	          	
                        0.8500, 0.3250, 0.0980;	          	
                        0.9290, 0.6940, 0.1250;          	
                        0.4940, 0.1840, 0.5560;	          	
                        0.4660, 0.6740, 0.1880;	          	
                        0.3010, 0.7450, 0.9330;          	
                        0.6350, 0.0780, 0.1840;];
            
            figure
            load coastlines %#ok<LOAD>
            axesm('miller', 'Frame', 'on', 'Grid', 'on', 'MeridianLabel', ...
                'on', 'ParallelLabel', 'on', 'MlabelLocation', ...
                [-180, -120, -60, 0, 60, 120, 180])
            geoshow(coastlat, coastlon, 'Color', [0.45, 0.45, 0.45])
            
            for i = 1:nb_lat
                
                plotm((90-(i-1)*dnu)*ones(360, 1), (0:359)', '--k', ...
                    'LineWidth', 1.5);
                hold on
                p1 = plotm(tab(1, 2:3, i), '*', 'MarkerSize', 13, ...
                    'Color', [0, 0, 0]);
                hold on
                
                for j = 1:N-1
                    
                    tabPlots(j) = plotm(tab(j+1, 2:3, i), '.', 'MarkerSize', 20, ...
                     'Color', colors(j, :)); %#ok<AGROW>                
                    
                end 

            end
            
            legend([p1, tabPlots], 'Start point', '4 hrs', '8 hrs', ...
                '12 hrs', '16 hrs', '20 hrs', '24 hrs', 'NumColumns', 2, ...
                'Position', [0.67 0.745 0 0])
                        
            textm(-86, 0, ...
                    'Earth rotation \rightarrow \rightarrow', ...
                    'FontSize', 16, 'VerticalAlignment', 'bottom', ...
                    'HorizontalAlignment', 'center', 'FontName', 'Times', ... 
                    'FontWeight', 'Bold') 
                
            addpath ...
            ('/Users/maxime/Documents/ORBs project/Dynamic model/v97/export_fig');
            export_fig('rot_earth');
            
        end 
        
        
        
        % COMPUTE THIRD BODY GRAVITATIONAL PERTURBATION FOR A WHOLE ORBIT
        % Only circular orbits are considered
        %
        % Inputs :
        %   dnu (scalar, radians) - true anomaly step size
        %
        %   mu3 (scalar, meters cubed/second squared) - third body's
        %                   gravitational parameter
        %
        %   a (scalar, meters) - radius of spacecraft's orbit
        %
        %   d12 (scalar, meters) - distance between orbited and third body
        %
        % Outputs :
        %   TabVec - matrix containing the perturbation vectors, expressed 
        %                   in inertial frame. Each row is one location 
        %                   on orbit :
        %                   TabVec(:,1:3) -> X, Y, Z components of
        %                   acceleration, in Newtons/kilograms
        function TabVec = TBGCheck(dnu, mu3, a, d12)
            
            N = 2*pi/dnu;
            TabVec = zeros(N, 9);
            
            Omega = 0;
            inc = 0;
            omega = 0;
            nu = 0;
            
            rho_hat = [0; 0; 1];
            rho = d12*rho_hat;
            
            for i = 1:N
                
                % compute r in inertial frame
                R_IO = (rotation_functions.compute_rot_matrix ...
                    (Omega, inc, omega+nu, '313'));
                
                r = R_IO*[a; 0; 0];
                
                % vector from spacecraft to third body
                r3 = rho - r;
               
                % computing third body acceleration vector
                TabVec(i, :) = [r', (R_IO*force_functions. ...
                   third_body_gravity(a, mu3, d12, rho_hat, Omega, inc, ...
                   omega, nu))', r3'];
                
                nu = nu + dnu;
                
            end 
            
            figure
            quiver3(TabVec(:, 1), TabVec(:, 2), zeros(N, 1), TabVec(:, 4), TabVec(:, 5), zeros(N, 1), 'Marker', '*')
            hold on
            pos = zeros(360, 2);
            for i = 1:360
                pos(i, :) = [a*cos((i-1)*pi/180), a*sin((i-1)*pi/180)];
            end
            plot(pos(:, 1), pos(:, 2))
            plot(0, 0, 'o');
            alpha = 0;
            RotMat = [cos(alpha), -sin(alpha), 0; 
                        sin(alpha), cos(alpha), 0;
                        0, 0, 1];
            vec_moon = RotMat*[0; 0; 10^7];
            quiver3(0, 0, 4*10^6, vec_moon(1), vec_moon(2), vec_moon(3));
            
            txt = 'Earth''s center';
            t = text(0, -7.5*10^5, txt, 'HorizontalAlignment','center');
            t.FontSize = 12;
            
            txt2 = 'To Moon';
            t2 = text(10^6, 1.25*10^7, txt2, 'HorizontalAlignment','center');
            t2.FontSize = 12;
            
            set(gca,'FontSize', 12)
            legend('Orbital positions with \newline third body grav. pert. \newline (acc. vector, scaled)', ...
                    'Spacecraft''s orbit', 'FontSize', 13, 'Location', 'northwest');
                
            xlabel('Inertial axis X (m)')
            ylabel('Inertial axis Y (m)')
            zlabel('Inertial axis Z (m)')
            
        end
        
        
        
        
        % PLOTTING THIRD BODY AND RADIATION PRESSURE PERTURBATION MAGNITUDES
        % Using a varying planet-centric radius and condering 2 third
        % bodies acting on spacecraft. 
        %
        % Inputs : 
        %   r0 (scalar, meters) - mean radius of orbited body
        %
        %   ri (scalar, meters) - initial planet-centric radius to consider
        %
        %   rf (scalar, meters) - last planet-centric radius to consider
        %
        %   dr (scalar, meters) - radius step
        %
        %   mu1 (scalar, meters cubed/second squared) - gravitational
        %                   parameter of 1st third body to consider
        %
        %   mu2 (scalar, meters cubed/second squared) - gravitational
        %                   parameter of 2nd third body
        %
        %   d1 (scalar, meters) - distance between orbited and 1st third
        %                   body
        %
        %   d2 (scalar, meters) - distance between orbited and 2nd third
        %                   body
        %
        %   A (meters squared) - reflective area of the spacecraft
        %
        %   C (unitless) - coefficient associated with reflective 
        %                   property of the spacecraft
        %                   0 < C < 2
        %                   0 all radiation is absorbed 
        %                   2 all radiation is reflected
        %
        %   m (kilograms) - mass of the spacecraft
        %
        % Comments : 
        %   This function was designed to draw plots for Sun and Moon as
        %   third bodies, while orbiting the Earth. It is considered that
        %   the Earth, the spacecraft, the Moon and the Sun were all 
        %   aligned along the [1;0;0] inertial frame axis. If considering
        %   other bodies, change the labels on the plots.
        function plot_th_rp_acc(r0, ri, rf, dr, mu1, mu2, d1, d2, A, C, m)
            
            N = (rf - ri)/dr + 1;
            
            th_body_1 = zeros(N, 3);
            th_body_2 = zeros(N, 3);
            rad_press = zeros(N, 3);
            alts = zeros(N, 1);
            
            r = ri;
            e = [1; 0; 0];
            
            for i = 1:N
                
                alts(i) = r - r0;
                
                th_body_1(i, :) = force_functions.third_body_gravity(r, ... 
                mu1, d1 - r, e, 0, 0, 0, 0);
            
                th_body_2(i, :) = force_functions.third_body_gravity(r, ... 
                mu2, d2 - r, e, 0, 0, 0, 0);
            
                rad_press(i, :) = force_functions.SRP(d1 - r, ...
                    0, 0, 0, 0, A, C, m, e);
                
                r = r + dr;
      
            end 
            
            figure
            loglog(alts/1000, abs(th_body_1(:, 1)), 'LineWidth',2);
            hold on
            loglog(alts/1000, abs(th_body_2(:, 1)), 'LineWidth',2);
            hold on
            loglog(alts/1000, abs(rad_press(:, 1)), 'LineWidth',2);
            grid on
            set(gca,'FontSize', 12)
            xlabel('Altitude (km)')
            ylabel('Acceleration (m/s^2)')
            ylim([10^(-8), 10^2])
            legend('Sun', 'Moon', 'Radiation pressure', 'Location', 'northwest')
            
        end
        
        
           
        % COMPUTING GRAVITATIONAL PERTURBATION IN LOCAL ORBIT FRAME
        % From a given potential gradient expressed in spherical frame
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
        %   lon_asc_node (scalar, radians) - scalar giving the longitude of 
        %                   ascending node of the spacecraft's orbit
        %
        %   inc (scalar, radians) - scalar giving the inclination of the
        %                   spacecraft's orbit
        %
        %   arg_per (scalar, radians) - scalar giving the argument of 
        %                   periapsis of the spacecraft's orbit
        %
        %   tr_anom (scalar, radians) - scalar giving the true anomaly of 
        %                   the spacecraft
        %
        %   pot_grad (vector, Newtons/kilograms) - gravitational potential 
        %                   gradient, expressed in spherical frame
        %                   pot_grad(1) -> southward component
        %                   pot_grad(2) -> eastward component
        %                   pot_grad(3) -> upward component
        %
        % Output :
        %   f_grav (vector, Newtons/kilograms) - gravitational perturbation 
        %                    vector, expressed in local orbital frame
        %                   f_grav(1) -> radial component
        %                   f_grav(2) -> along-track component
        %                   f_grav(3) -> orbit-normal component
        function f_grav = compute_grav_force_loc_orbit_frame(r_asc, dec, s_time, ...
            ecc, semi_maj_axis, lon_asc_node, inc, arg_per, tr_anom, pot_grad)

            % compute the coordinates of the position vector in planetary frame
            plan_coord = rotation_functions.LocalToPlanetaryCoord(r_asc, dec, s_time, ...
            ecc, semi_maj_axis, lon_asc_node, inc, arg_per, tr_anom);

            lat = plan_coord(2);
            lon = plan_coord(3);

            % compute rotation matrix from spherical to planetary frame 
            R_VP = rotation_functions.compute_rot_matrix(lon, lat, 0, '323');

            % compute rotation matrix from planetary frame to inertial frame
            R_PI = rotation_functions.compute_rot_matrix(r_asc, dec, s_time, '313');

            % compute rotation matrix from local orbit frame to inertial frame
            R_IO = rotation_functions.compute_rot_matrix(lon_asc_node, inc, ...
                arg_per + tr_anom, '313');

            f_grav = (transpose(R_IO)*R_PI*R_VP)*pot_grad;

        end
        
        
        
        % COMPUTE GRAVITATIONAL PERTURBATION ON WHOLE MAP
        % Using polar orbits and spherical gradient
        %
        % Inputs :
        %   mu (meters cubed/seconds squared) - gravitational parameter 
        %                   of orbited body
        %
        %   r0 (meters) - mean radius of orbited body
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
        %   ecc (scalar, unitless) - eccentricity of spacecraft's orbit
        %
        %   semi_maj_axis (scalar, meters) - semi-major axis of spacecraft's 
        %                   orbit
        %
        %   arg_per (scalar, radians) - scalar giving the argument of 
        %                   periapsis of the spacecraft's orbit
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
        %   harm_parameter (string) - parameter used to consider J2 harmonic 
        %                   (or not) in the computations
        %                   "oblateness" -> with J2
        %                   "llumpiness" -> without J2
        %
        %   dnu (scalar, radians) - true anomaly step
        %
        %   dlon_asc_node (scalar, radians) - longitude of ascending node
        %                   step
        %
        %   dr (scalar, meters) - radius step size used in gradient
        %                   computations
        %
        %   dlat (scalar, radians) - latitude step size used in gradient
        %                   computations
        %
        %   dlon (scalar, radians) - longitude step size used in gradient
        %                   computations
        %
        % Outputs :
        %
        %   positions (NO by NA by 3 matrix) - spherical coordinates 
        %                   corresponding to the successive positions
        %                   occupied by the spacecraft
        %                   NO -> # of orbits (or latitudes)
        %                   NA -> # of true anomalies in one orbit (or
        %                   longitudes)
        %                   positions(1) -> southward component
        %                   positions(2) -> eastward component
        %                   positions(3) -> upward component
        %
        %   U_map (NO by NA matrix) - values of potential at successive
        %                   positions
        %
        %   U_map_tot (NO by NA by 3 matrix) - values of potential at 
        %                   successive positions for 3 values of
        %                   planet-centric radius
        %
        %   grad_U_map (NO by NA by 3 matrix) - potential gradient vectors 
        %                   at successive positions, expressed in spherical
        %                   frame
        %                   grad_U_map(1) -> southward component
        %                   grad_U_map(2) -> eastward component
        %                   grad_U_map(3) -> upward component
        %
        %   f_map (NO by NA by 3 matrix) - gravitational perturbation vectors
        %                   at successive positions, expressed in local
        %                   orbital frame
        %                   f_map(1) -> radial component
        %                   f_map(2) -> along-track component
        %                   f_map(3) -> orbit-normal component                
        function [positions, U_map, U_map_tot, grad_U_map, f_map] = ...
                f_grav_map(mu, r_0, r_asc, dec, s_time, ecc, semi_maj_axis, ...
                arg_per, raw_coeff, harm_parameter, dnu, dlon_asc_node, ...
                dr, dlat, dlon)
            
            inc = pi/2;
            
            nb_nu = (2*pi)/dnu;
            nb_orbs = pi/dlon_asc_node;

            nu_step = dnu;
            orbs_step = dlon_asc_node;

            U = zeros(3:3:3);
            harmonic_coeff = force_functions.read_harm_coeff_from_file(raw_coeff);

            positions = zeros(nb_orbs, nb_nu, 3); 
            U_map = zeros(nb_orbs, nb_nu);
            U_map_tot = zeros(nb_orbs, nb_nu, 3);
            grad_U_map = zeros(nb_orbs, nb_nu, 3);
            f_map = zeros(nb_orbs, nb_nu, 3);

            count = 0;
            total_it = (nb_orbs + 1)*(nb_nu + 1);

            for i = 1:nb_orbs + 1

                for j = 1:nb_nu + 1 

                    % computing the location of the s/c in spherical coordinates
                    sph_coord = rotation_functions.LocalToPlanetaryCoord(r_asc, dec, s_time, ...
                    ecc, semi_maj_axis, (i-1)*orbs_step, inc, ...
                    arg_per, (j-1)*nu_step);
                
                    positions(i, j, :) = [sph_coord(1); sph_coord(2:3)*(180/pi)];

                    % computing the potential at and around the location
                    for l = 1:3

                        for m = 1:3 

                            for n = 1:3
                           
                            U(l, m, n) = force_functions.compute_potential(mu, r_0, sph_coord(1) + dr*(n-2), ...
                                sph_coord(3) + dlon*(m-2), sph_coord(2) + dlat*(l-2), ...
                                harmonic_coeff, harm_parameter);

                            end

                        end 

                    end 

                    U_map(i, j) = U(2, 2, 2);
                    U_map_tot(i, j, :) = U(2, 2, :);

                    % computation of spherical gradient (J/kg/m)
                    r_start = sph_coord(1) - dr;
                    r_stop = sph_coord(1) + dr;

                    coord_r = [r_start; sph_coord(1); r_stop];
                    lats = [sph_coord(2) - dlat, sph_coord(2), sph_coord(2) + dlat];
                    [gradX_s, gradY_s, gradZ_s, ~] = ...
                        finite_difference_schemes.spherical_gradient_2 ...
                        ((10^3)*coord_r, lats*(180/pi), dlat*(180/pi), ...
                        dlon*(180/pi), (10^3)*dr, (10^6)*U);            

                    pot_grad = [gradX_s(2, 2, 2); gradY_s(2, 2, 2); gradZ_s(2, 2, 2)]; 
                    grad_U_map(i, j, :) = pot_grad;             

                    % compute rotation matrix from spherical to planetary frame 
                    R_VP = rotation_functions.compute_rot_matrix(sph_coord(3), sph_coord(2), 0, '323');

                    % compute rotation matrix from planetary frame to inertial frame
                    R_PI = rotation_functions.compute_rot_matrix(r_asc, dec, s_time, '313');

                    % compute rotation matrix from local orbit frame to inertial frame
                    R_IO = rotation_functions.compute_rot_matrix((i - 1)*orbs_step, inc, ...
                        arg_per + (j-1)*nu_step, '313');                 
                    
                    R = (transpose(R_IO)*R_PI*R_VP);
                    f_grav = R*pot_grad;

                    f_map(i, j, 1) = f_grav(1);
                    f_map(i, j, 2) = f_grav(2);
                    f_map(i, j, 3) = f_grav(3);

                    count = count + 1;
                    (count/total_it)*100 %#ok<NOPRT>
                    
                end 

            end  

        end
        
        
        
        % COMPUTE GRAVITATIONAL PERTURBATION ON WHOLE MAP
        % Using polar orbits and cartesian gradient 
        %
        % Inputs :
        %   mu (meters cubed/seconds squared) - gravitational parameter 
        %                   of orbited body
        %
        %   r0 (meters) - mean radius of orbited body
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
        %   ecc (scalar, unitless) - eccentricity of spacecraft's orbit
        %
        %   semi_maj_axis (scalar, meters) - semi-major axis of spacecraft's 
        %                   orbit
        %
        %   lon_asc_node (scalar, radians) - scalar giving the longitude of 
        %                   ascending node of the spacecraft's orbit
        %
        %   inc (scalar, radians) - scalar giving the inclination of the
        %                   spacecraft's orbit
        %
        %   arg_per (scalar, radians) - scalar giving the argument of 
        %                   periapsis of the spacecraft's orbit
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
        %   harm_parameter (string) - parameter used to consider J2 harmonic 
        %                   (or not) in the computations
        %                   "oblateness" -> with J2
        %                   "llumpiness" -> without J2
        %
        %   dnu (scalar, radians) - true anomaly step
        %
        %   dlon_asc_node (scalar, radians) - longitude of ascending node
        %                   step
        %
        %   dx (scalar, meters) - step size along 1st axis of inertial frame
        %
        %   dy (scalar, meters) - step size along 2nd axis of inertial frame
        %
        %   dz (scalar, meters) - step size along 3rd axis of inertial frame
        %
        % Outputs :
        %
        %   positions (NO by NA by 3 matrix) - spherical coordinates 
        %                   corresponding to the successive positions
        %                   occupied by the spacecraft
        %                   NO -> # of orbits (or latitudes)
        %                   NA -> # of true anomalies in one orbit (or
        %                   longitudes)
        %                   positions(1) -> southward component
        %                   positions(2) -> eastward component
        %                   positions(3) -> upward component
        %
        %   U_map (NO by NA matrix) - values of potential at successive
        %                   positions
        %
        %   U_map_tot (NO by NA by 3 matrix) - values of potential at 
        %                   successive positions for 3 values of
        %                   planet-centric radius
        %
        %   GradMap (NO by NA by 3 matrix) - potential gradient vectors 
        %                   at successive positions, expressed in spherical
        %                   frame
        %                   grad_U_map(1) -> southward component
        %                   grad_U_map(2) -> eastward component
        %                   grad_U_map(3) -> upward component
        %
        %   f_map (NO by NA by 3 matrix) - gravitational perturbation vectors
        %                   at successive positions, expressed in local
        %                   orbital frame
        %                   f_map(1) -> radial component
        %                   f_map(2) -> along-track component
        %                   f_map(3) -> orbit-normal component
        function [positions, U_map, U_map_tot, GradMap, f_map] = GravMap ...
                (mu, r_0, r_asc, dec, s_time, ecc, semi_maj_axis, arg_per, ...
                raw_coeff, harm_parameter, dnu, dlon_asc_node, dx, dy, dz)
            
            inc = pi/2;
            
            nb_nu = (2*pi)/dnu;
            nb_orbs = pi/dlon_asc_node;

            nu_step = dnu;
            orbs_step = dlon_asc_node;

            U = zeros(3:3:3);
            harmonic_coeff = force_functions.read_harm_coeff_from_file(raw_coeff);

            positions = zeros(nb_orbs, nb_nu, 3); 
            U_map = zeros(nb_orbs, nb_nu);
            U_map_tot = zeros(nb_orbs, nb_nu, 3);
            GradMap = zeros(nb_orbs, nb_nu, 3);
            f_map = zeros(nb_orbs, nb_nu, 3);

            count = 0;
            total_it = (nb_orbs + 1)*(nb_nu + 1);
            
            coordX = zeros(3, 3, 3);
            coordY = zeros(3, 3, 3);
            coordZ = zeros(3, 3, 3);

            for i = 1:nb_orbs + 1

                for j = 1:nb_nu + 1 

                    % computing the location of the s/c in spherical coordinates
                    sph_coord = rotation_functions.LocalToPlanetaryCoord...
                        (r_asc, dec, s_time, ecc, semi_maj_axis, ...
                        (i-1)*orbs_step, inc, arg_per, (j-1)*nu_step);
                
                    positions(i, j, :) = [sph_coord(1); sph_coord(2:3)*(180/pi)];
                    
                    r = sph_coord(1);
                    lat = sph_coord(2);
                    lon = sph_coord(3);                
                    
                    % Body-fixed frame coordinates
                    x = r*sin(lat)*cos(lon);
                    y = r*sin(lat)*sin(lon);
                    z = r*cos(lat);
                  
                    % computing the potential at and around the location
                    for n = 1:3

                        for l = 1:3 

                            for m = 1:3
                                
                            coordX(l, m, n) =  x + dx*(l-2);
                            coordY(l, m, n) =  y + dy*(m-2);
                            coordZ(l, m, n) =  z + dz*(n-2);
                            
                            U(l, m, n) = force_functions.PotentialFromBF...
                                (mu, r_0, x + dx*(l-2), ...
                                y + dy*(m-2), ...
                                z + dz*(n-2), ...
                                harmonic_coeff, harm_parameter);

                            end

                        end 

                    end 

                    U_map(i, j) = U(2, 2, 2);
                    U_map_tot(i, j, :) = U(2, 2, :);

                    % computation of gradient in planetary frame
                    [gradX, gradY, gradZ] = finite_difference_schemes. ...
                        gradient_o4((10^6)*U, dx*(10^3), dy*(10^3), dz*(10^3)); 

                    grad = [gradX(2, 2, 2); gradY(2, 2, 2); gradZ(2, 2, 2)]; 
                    
                    GradMap(i, j, :) = grad; 

                    % compute rotation matrix from planetary frame to inertial frame
                    R_PI = rotation_functions.compute_rot_matrix(r_asc, dec, s_time, '313');

                    % compute rotation matrix from local orbit frame to inertial frame
                    R_IO = rotation_functions.compute_rot_matrix((i - 1)*orbs_step, inc, ...
                        arg_per + (j-1)*nu_step, '313');                 
                    
                    R = (transpose(R_IO)*R_PI);
                    f_grav = R*grad;

                    f_map(i, j, 1) = f_grav(1);
                    f_map(i, j, 2) = f_grav(2);
                    f_map(i, j, 3) = f_grav(3);

                    count = count + 1;
                    (count/total_it)*100 %#ok<NOPRT>
                    
                end 

            end  

        end
        
        
        
        % COMPUTE POTENTIAL ON WHOLE MAP FOR RANGE OF RADII
        %
        % Inputs :
        %   r_ref (scalar, meters) - mean radius of orbited body
        %
        %   mu (meters cubed/seconds squared) - gravitational parameter 
        %                   of orbited body 
        %
        %   r_start(scalar, meters) - initial radius
        %
        %   r_stop(scalar, meters) - final radius
        %
        %   r_step(scalar, meters) - radius step size 
        %
        %   lat_step(scalar, degrees) - latitude step size
        %
        %   lon_step(scalar, degrees) - longitude step size
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
        %   harm_parameter - string parameter used to consider J2 harmonic (or not) 
        %                   in the computations
        %                   "oblateness" -> with J2
        %                   "llumpiness" -> without J2
        %
        % Outputs :
        %   U (NP by NT by NR matrix, joules/kilograms) - values of 
        %                   potential at each considered position and
        %                   radius
        %                   NP -> # of latitudes
        %                   NT -> # of longitudes
        %                   NR -> # of radii
        %
        %   coordX (NP by NT by NR matrix, meters) - 1st body-fixed axis
        %                   coordinate of considered positions
        %
        %   coordY (NP by NT by NR matrix, meters) - 2nd body-fixed axis
        %                   coordinate of considered positions
        %
        %   coordZ (NP by NT by NR matrix, meters) - 3rd body-fixed axis
        %                   coordinate of considered positions
        %
        %   lats (NP by NT by NR matrix, radians) - latitudes considered
        %
        %   lons (NP by NT by NR matrix, radians) - longitudes considered
        function [U, coordX, coordY, coordZ, lats, lons] = potential_3D ...
                (r_ref, mu, r_start, r_stop, r_step, lat_step, lon_step, ...
                raw_coeff, harm_parameter)

            nb_alt = ((r_stop - r_start)/r_step) + 1; 
            nb_lat = (180/lat_step) + 1;
            nb_lon = (360/lon_step) + 1;

            harmonic_coeff = force_functions.read_harm_coeff_from_file(raw_coeff);
            U = zeros(nb_lat, nb_lon, nb_alt);
            
            coordX = zeros(nb_lat, nb_lon, nb_alt);
            coordY = zeros(nb_lat, nb_lon, nb_alt);
            coordZ = zeros(nb_lat, nb_lon, nb_alt);
            
            lats = zeros(nb_lat, nb_lon, nb_alt);
            lons = zeros(nb_lat, nb_lon, nb_alt);
            rs = zeros(nb_lat, nb_lon, nb_alt);

            for z = 1:nb_alt

                for i = 1:nb_lat

                    for j = 1:nb_lon
                        
                        r = r_start + (z-1)*r_step;
                        lat = (i-1)*lat_step*(pi/180);
                        lon = (j-1)*lon_step*(pi/180);
                        
                        lats(i, j, z) = lat;
                        lons(i, j, z) = lon;
                        rs(i, j, z) = r;
                        
                        coordX(i, j, z) = r*sin(lat)*cos(lon);
                        coordY(i, j, z) = r*sin(lat)*sin(lon);
                        coordZ(i, j, z) = r*cos(lat);

                        U(i, j, z) = test_functions.compute_potential(mu, ...
                            r_ref, r_start + (z-1)*r_step, (j-1)*lon_step*(pi/180), ...
                            (i-1)*lat_step*(pi/180), harmonic_coeff, harm_parameter);  

                    end        

                end

            end

        end
        
        
        
        % COMPUTE POTENTIAL ON WHOLE MAP AT A GIVEN RADIUS
        %
        % Inputs :
        %   r_ref (scalar, meters) - mean radius of orbited body
        %
        %   mu (meters cubed/seconds squared) - gravitational parameter 
        %                   of orbited body 
        %
        %   r (scalar, meters) - planet-centric radius considered
        %
        %   lat_step(scalar, degrees) - latitude step size
        %
        %   lon_step(scalar, degrees) - longitude step size
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
        % Outputs :
        %   positions (NP by NT matrix, degrees) - spherical coordinates 
        %                   corresponding to the successive positions
        %                   occupied by the spacecraft
        %                   NP -> # of latitudes
        %                   NT -> # of longitudes
        %                   positions(:,1) -> longitudes
        %                   positions(:,2) -> latitudes
        %
        %   U_map (NP by NT matrix, joules/kilograms) - values of 
        %                   potential at successive positions for given
        %                   radius
        function [positions, U_map] = potential_2D_map(r_ref, mu, r, ...
                lat_step, lon_step, raw_coeff)

            nb_lat = (180/lat_step)+1;
            nb_lon = (360/lon_step)+1;
            
            positions = zeros(nb_lat, nb_lon, 2);

            %harmonic_coeff = compute_harm_coeff_earth;
            harmonic_coeff = force_functions.read_harm_coeff_from_file(raw_coeff);
            U_map = zeros(nb_lat, nb_lon);

            for i = 1:nb_lat

               for j = 1:nb_lon

                   U_map(i, j) = force_functions.compute_potential(mu, ...
                       r_ref, r, (j-1)*lon_step*(pi/180), (i-1)*lat_step*(pi/180), ...
                       harmonic_coeff, 'llumpiness');
                   
                   positions(i, j, 1) = (j-1)*lon_step;
                   positions(i, j, 2) = (i-1)*lat_step;

               end        

            end

        end
        
        
        
        % ANALYTICAL DERIVATION OF GRAVITATIONAL PERTURBATION
        %
        % Inputs :  
        %   mu (meters cubed/seconds squared) - gravitational parameter 
        %                   of orbited body 
        %
        %   r_0 (scalar, meters) - mean radius of orbited body
        %
        %   x, y, z (scalars, meters) - coordinates considered, in
        %                   body-fixed frame 
        %
        %   harm_coeff - matrix containing the spherical harmonic coefficients
        %                   associated with the orbited body. Use the
        %                   function force_functions.read_harm_coeff_from_file 
        %                   with the raw data from space agencies' servers 
        %                   to generate this matrix, but REMOVE the
        %                   denormalization by setting C = 1 in that
        %                   function.
        %
        %   harm_parameter - string parameter used to consider J2 harmonic (or not) 
        %                   in the computations
        %                   "oblateness" -> with J2
        %                   "llumpiness" -> without J2
        %
        % Outputs :
        %   gx, gy, gz (scalars, Newtons/kilograms) - components of 
        %                   gravitational perturbation vector, expressed
        %                   in body-fixed frame
        %
        % Comments :
        %   Spherical harmonics coefficients must be NORMALIZED in this
        %   function. 
        function [gx, gy, gz] = analyticalgrav(mu, r_0, x, y, z, ...
            harmonic_coeff, harm_parameter)

            r = sqrt(x^2 + y^2 + z^2);
            lambda = atan2(y, x);
            phi = (pi/2) - acos(z/r);

            drr = [x/r y/r z/r]';
            dpr = (1/sqrt(x^2+y^2))*[-x*z/(r^2) -y*z/(r^2) 1-(z^2/r^2)]';
            dlr = (1/(x^2+y^2))*[-y x 0]';

            inf_TS = size(harmonic_coeff, 1);
            
            dUr = 0;
            dUp = 0;
            dUl = 0;

            start_j = 2; %default value
            if (harm_parameter == 'oblateness') %#ok<BDSCA>
                start_j = 1;
            end 

            for i = 2:inf_TS
                
                % associated Legendre polynomials of degree i and order j = 0, ..., i  
                P = legendre(i, sin(phi), 'norm');

                N = i+1;
                for j = start_j:N

                    Cij = harmonic_coeff(i,j,2);
                    Sij = harmonic_coeff(i,j,3);
                    Pij = sqrt((i+j)*(i-j+1)*(2-eq(j,0))/2);
                    
                    dUr = dUr - (mu/r^2)*((r_0/r)^i)*(i+1)*P(j)* ...
                        (Cij*cos((j-1)*lambda) + Sij*sin((j-1)*lambda));
                    
                    if j < N
                        
                        dUp = dUp + (mu/r)*((r_0/r)^i)* ...
                            (P(j+1)*Pij - (j-1)*tan(phi)*P(j))* ...
                            (Cij*cos((j-1)*lambda) + Sij*sin((j-1)*(lambda)));
                    
                    end 

                    dUl = dUl + (mu/r)*((r_0/r)^i)*(j-1)*P(j)* ...
                        (Sij*cos((j-1)*lambda) - Cij*sin((j-1)*(lambda)));            

                end 

            end  

            g = dUr*drr + dUp*dpr + dUl*dlr;
            
            gx = g(1);
            gy = g(2);
            gz = g(3);

        end
        
        
        
        % COMPUTE POTENTIAL GRADIENT EXPRESSED IN SPHERICAL FRAME
        % At a given location
        %
        % Inputs :
        %   mu (meters cubed/seconds squared) - gravitational parameter 
        %                   of orbited body 
        %
        %   r_0 (scalar, meters) - mean radius of orbited body
        %
        %   sph_coord (3 by 1 vector) - contains the spherical coordinates
        %                   of location to consider
        %                   sph_coord(1) -> radius
        %                   sph_coord(2) -> longitude
        %                   sph_coord(3) -> latitude
        %
        %   dr (scalar, meters) - radius step size used in gradient
        %                   computations
        %
        %   dlat (scalar, radians) - latitude step size used in gradient
        %                   computations
        %
        %   dlon (scalar, radians) - longitude step size used in gradient
        %                   computations 
        %
        %   harmonic_coeff - matrix containing the spherical harmonic coefficients
        %                   associated with the orbited body. Use the
        %                   function force_functions.read_harm_coeff_from_file 
        %                   with the raw data from space agencies' servers 
        %                   to generate this matrix
        %
        %   harm_parameter (string) - parameter used to consider J2 harmonic 
        %                   (or not) in the computations
        %                   "oblateness" -> with J2
        %                   "llumpiness" -> without J2
        %
        % Outputs :
        %   potential (vector, Newtons/kilograms) - potential gradient
        %                   expressed in spherical frame 
        function potential = compute_potential_grad(mu, r_0, sph_coord, ...
                dr, dlat, dlon, harmonic_coeff, harm_parameter)
            
            U = zeros(3, 3, 3);
            for l = 1:3

                for m = 1:3 

                    for n = 1:3

                        U(l, m, n) = force_functions.compute_potential...
                            (mu, r_0, sph_coord(1) + dr*(n-2), ...
                            sph_coord(3) + dlon*(m-2), ...
                            sph_coord(2) + dlat*(l-2), ...
                            harmonic_coeff, harm_parameter);

                    end

                end 

            end 

            % computation of spherical gradient (J/kg/m)
            r_start = sph_coord(1) - dr;
            r_stop = sph_coord(1) + dr;

            coord_r = [r_start; sph_coord(1); r_stop];
            lats = [sph_coord(2) - dlat, sph_coord(2), sph_coord(2) + dlat];
            [gradX_s, gradY_s, gradZ_s, ~] = finite_difference_schemes. ...
                spherical_gradient_o4((10^3)*coord_r, ...
                lats*(180/pi), dlat, dlon, (10^3)*dr, (10^6)*U);            

            potential = [gradX_s(2, 2, 2); gradY_s(2, 2, 2); gradZ_s(2, 2, 2)];
            
        end
        
        
        
        % COMPUTE GRAVITATIONAL POTENTIAL
        % At given location
        %
        % Inputs :
        %   mu (meters cubed/seconds squared) - gravitational parameter 
        %                   of orbited body 
        %
        %   r_0 (scalar, meters) - mean radius of orbited body
        %
        %   r (scalar, meters) - planet-centric radius considered
        %
        %   lambda (scalar, radians) - spacecraft's longitude, in range [0;2*pi]
        %
        %   phi (scalar, radians) - spacecraft's latitude, in range [0;pi]
        %
        %   harm_coeff - matrix containing the spherical harmonic coefficients
        %                   associated with the orbited body. Use the
        %                   function force_functions.read_harm_coeff_from_file 
        %                   with the raw data from space agencies' servers 
        %                   to generate this matrix
        %
        %   harm_parameter (string) - parameter used to consider J2 harmonic 
        %                   (or not) in the computations
        %                   "oblateness" -> with J2
        %                   "llumpiness" -> without J2
        %
        % Outputs :
        %   U (scalar, Joules/kilograms) - potential value at given location
        function U = compute_potential(mu, r_0, r, lambda, phi, harmonic_coeff, harm_parameter) 

            inf_TS = size(harmonic_coeff, 1);
            temp_TS = 0.;

            start_j = 2; %default value
            if (harm_parameter == 'oblateness') %#ok<BDSCA>
                start_j = 1;
            end 

            for i = 2:inf_TS
                
                % associated Legendre polynomials of degree i and order j = 0, ..., i  
                P = legendre(i, cos(phi));

                N = i+1;
                for j = start_j:N
    
                    temp_TS = temp_TS + ((r_0/r)^i)*P(j)*(harmonic_coeff(i,j,2)*cos((j-1)*(-pi+lambda)) ...
                    + harmonic_coeff(i,j,3)*sin((j-1)*(-pi+lambda)));             

                end 

            end  

            % add 1 to temp_TS to have whole potential
            % default value : (mu/r)*(temp_TS), PERTURBATION potential  
            U = (mu/r)*(temp_TS);

        end
        
        
        
        % COMPUTE POSITIONS OF SPACECRAFT WITH INCLINATION CHANGE IN THE ORBITS
        % This function is aimed at testing the functions converting 
        % a position expressed in local frame to a position in spherical
        % frame
        %
        % Inputs :
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
        %   ecc (scalar, unitless) - eccentricity of spacecraft's orbit
        %
        %   semi_maj_axis (scalar, meters) - semi-major axis of spacecraft's 
        %                   orbit
        %
        %   step_inc (scalar, radians) - step size of inclination change
        %
        %   arg_per (scalar, radians) - scalar giving the argument of 
        %                   periapsis of the spacecraft's orbit
        %
        %   lon_asc_node (scalar, radians) - scalar giving the longitude of 
        %                   ascending node of the spacecraft's orbit
        %
        %   step_nu (scalar, radians) - true anomaly step size
        %
        % Outputs :
        %   positions (NO by NA by 3 matrix) - spherical coordinates
        %                   corresponding to positions occupied by spacecraft
        %                   NO -> # of orbits
        %                   NA -> # of true anomalies
        %                   positions(:,:,1) -> planet-centric radii
        %                   positions(:,:,2) -> latitudes
        %                   positions(:,:,3) -> longitudes
        function [positions] = map_orbits_inc(r_asc, dec, s_time, ecc, ...
            semi_maj_axis, step_inc, arg_per, lon_asc_node, step_nu)

            nb_orbs = 2*pi/step_inc;
            nb_nu_pos = 2*pi/step_nu;
    
            positions = zeros(nb_orbs, nb_nu_pos, 3);
    
            for i = 1:nb_orbs
        
                positions(i, :, :) = test_functions.compute_orbit_spherical_coord(r_asc, dec, ... 
                    s_time, ecc, semi_maj_axis, lon_asc_node, (i-1)*step_inc, arg_per, step_nu);
        
            end 

        end
        
        
        
        % COMPUTE POSITIONS OF SPACECRAFT WITH ASC. NODE CHANGE IN THE ORBITS
        % This function is aimed at testing the functions converting 
        % a position expressed in local frame to a position in spherical
        % frame
        %
        % Inputs :
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
        %   ecc (scalar, unitless) - eccentricity of spacecraft's orbit
        %
        %   semi_maj_axis (scalar, meters) - semi-major axis of spacecraft's 
        %                   orbit
        %
        %   inc (scalar, radians) - inclination of spacecraft's orbit
        %
        %   arg_per (scalar, radians) - scalar giving the argument of 
        %                   periapsis of the spacecraft's orbit
        %
        %   step_asc (scalar, radians) - longitude of asc. node step
        %
        %   step_nu (scalar, radians) - true anomaly step size
        %
        % Outputs :
        %   positions (NO by NA by 3 matrix) - spherical coordinates
        %                   corresponding to positions occupied by spacecraft
        %                   NO -> # of orbits
        %                   NA -> # of true anomalies
        %                   positions(:,:,1) -> planet-centric radii
        %                   positions(:,:,2) -> latitudes
        %                   positions(:,:,3) -> longitudes
        function [positions] = map_orbits_asc_node(r_asc, dec, s_time, ecc, ...
            semi_maj_axis, inc, arg_per, step_asc, step_nu)

            nb_orbs = 2*pi/step_asc;
            nb_nu_pos = 2*pi/step_nu;
    
            positions = zeros(nb_orbs, nb_nu_pos, 3);
    
            for i = 1:nb_orbs
        
                positions(i, :, :) = test_functions.compute_orbit_spherical_coord(r_asc, dec, ... 
                    s_time, ecc, semi_maj_axis, (i-1)*step_asc, inc, arg_per, step_nu);
        
            end 

        end
        
        
        
        % COMPUTE POSITIONS OF SPACECRAFT WITH RIGHT ASCENSION CHANGE 
        % This function is aimed at testing the functions converting 
        % a position expressed in local frame to a position in spherical
        % frame
        %
        % Inputs :
        %   
        %   step_rasc (scalar, radians) - right ascension step
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
        %   inc (scalar, radians) - inclination of spacecraft's orbit
        %
        %   arg_per (scalar, radians) - scalar giving the argument of 
        %                   periapsis of the spacecraft's orbit
        %
        %   asc_node (scalar, radians) - scalar giving the longitude of 
        %                   ascending node of the spacecraft's orbit
        %
        %   step_nu (scalar, radians) - true anomaly step size
        %
        % Outputs :
        %   positions (NO by NA by 3 matrix) - spherical coordinates
        %                   corresponding to positions occupied by spacecraft
        %                   NO -> # of orbits
        %                   NA -> # of true anomalies
        %                   positions(:,:,1) -> planet-centric radii
        %                   positions(:,:,2) -> latitudes
        %                   positions(:,:,3) -> longitudes
        function [positions] = map_orbits_rasc(step_rasc, dec, s_time, ecc, ...
            semi_maj_axis, inc, arg_per, asc_node, step_nu)

            nb_orbs = 2*pi/step_rasc;
            nb_nu_pos = 2*pi/step_nu;
    
            positions = zeros(nb_orbs, nb_nu_pos, 3);
    
            for i = 1:nb_orbs
        
                positions(i, :, :) = test_functions.compute_orbit_spherical_coord((i-1)*step_rasc, ...
                    dec, s_time, ecc, semi_maj_axis, asc_node, inc, arg_per, step_nu);
        
            end 

        end 
        
        
        
        % COMPUTE POSITIONS OF SPACECRAFT WITH DECLINATION CHANGE 
        % This function is aimed at testing the functions converting 
        % a position expressed in local frame to a position in spherical
        % frame
        %
        % Inputs :
        %   
        %   r_asc (scalar, radians) - right ascension of orbited body
        %                   with respect to inertial frame
        %
        %   step_dec (scalar, radians) - declination step
        %
        %   s_time (scalar, radians) - sidereal time of orbited body with 
        %                   respect to inertial frame
        %
        %   ecc (scalar, unitless) - eccentricity of spacecraft's orbit
        %
        %   semi_maj_axis (scalar, meters) - semi-major axis of spacecraft's 
        %                   orbit
        %
        %   inc (scalar, radians) - inclination of spacecraft's orbit
        %
        %   arg_per (scalar, radians) - scalar giving the argument of 
        %                   periapsis of the spacecraft's orbit
        %
        %   asc_node (scalar, radians) - scalar giving the longitude of 
        %                   ascending node of the spacecraft's orbit
        %
        %   step_nu (scalar, radians) - true anomaly step size
        %
        % Outputs :
        %   positions (NO by NA by 3 matrix) - spherical coordinates
        %                   corresponding to positions occupied by spacecraft
        %                   NO -> # of orbits
        %                   NA -> # of true anomalies
        %                   positions(:,:,1) -> planet-centric radii
        %                   positions(:,:,2) -> latitudes
        %                   positions(:,:,3) -> longitudes
        function [positions] = map_orbits_dec(r_asc, step_dec, s_time, ecc, ...
            semi_maj_axis, inc, arg_per, asc_node, step_nu)

            nb_orbs = 2*pi/step_dec;
            nb_nu_pos = 2*pi/step_nu;
    
            positions = zeros(nb_orbs, nb_nu_pos, 3);
    
            for i = 1:nb_orbs
        
                positions(i, :, :) = test_functions.compute_orbit_spherical_coord(r_asc, ...
                    (i-1)*step_dec, s_time, ecc, semi_maj_axis, asc_node, inc, arg_per, step_nu);
        
            end 

        end
        
        
        
        % COMPUTE POSITIONS OF SPACECRAFT WITH SIDEREAL TIME CHANGE 
        % This function is aimed at testing the functions converting 
        % a position expressed in local frame to a position in spherical
        % frame
        %
        % Inputs :
        %   
        %   r_asc (scalar, radians) - right ascension of orbited body
        %                   with respect to inertial frame
        %
        %   dec (scalar, radians) - declination of orbited body
        %                   with respect to inertial frame
        %
        %   step_s_time (scalar, radians) - sidereal time step
        %
        %   ecc (scalar, unitless) - eccentricity of spacecraft's orbit
        %
        %   semi_maj_axis (scalar, meters) - semi-major axis of spacecraft's 
        %                   orbit
        %
        %   inc (scalar, radians) - inclination of spacecraft's orbit
        %
        %   arg_per (scalar, radians) - scalar giving the argument of 
        %                   periapsis of the spacecraft's orbit
        %
        %   asc_node (scalar, radians) - scalar giving the longitude of 
        %                   ascending node of the spacecraft's orbit
        %
        %   step_nu (scalar, radians) - true anomaly step size
        %
        % Outputs :
        %   positions (NO by NA by 3 matrix) - spherical coordinates
        %                   corresponding to positions occupied by spacecraft
        %                   NO -> # of orbits
        %                   NA -> # of true anomalies
        %                   positions(:,:,1) -> planet-centric radii
        %                   positions(:,:,2) -> latitudes
        %                   positions(:,:,3) -> longitudes
        function [positions] = map_orbits_stime(r_asc, dec, step_s_time, ecc, ...
            semi_maj_axis, inc, arg_per, asc_node, step_nu)

            nb_orbs = 2*pi/step_s_time;
            nb_nu_pos = 2*pi/step_nu;
    
            positions = zeros(nb_orbs, nb_nu_pos, 3);
    
            for i = 1:nb_orbs
        
                positions(i, :, :) = test_functions.compute_orbit_spherical_coord(r_asc, ...
                    dec, (i-1)*step_s_time, ecc, semi_maj_axis, asc_node, inc, arg_per, step_nu);
        
            end 

        end
        
        
        
        % CONVERT SUCCESSIVE ORBIT POSITIONS TO SPHERICAL FRAME COORDINATES
        % This function is aimed at testing the functions converting 
        % a position expressed in local frame to a position in spherical
        % frame
        %
        % Inputs :
        %   
        %   r_asc (scalar, radians) - right ascension of orbited body
        %                   with respect to inertial frame
        %
        %   dec (scalar, radians) - declination of orbited body
        %                   with respect to inertial frame
        %
        %   s_time (scalar, radians) - sidereal time of orbited body
        %                   with respect to inertial frame
        %
        %   ecc (scalar, unitless) - eccentricity of spacecraft's orbit
        %
        %   semi_maj_axis (scalar, meters) - semi-major axis of spacecraft's 
        %                   orbit
        %
        %   lon_asc_node (scalar, radians) - scalar giving the longitude of 
        %                   ascending node of the spacecraft's orbit
        %
        %   inc (scalar, radians) - inclination of spacecraft's orbit
        %
        %   arg_per (scalar, radians) - scalar giving the argument of 
        %                   periapsis of the spacecraft's orbit
        %
        %   step (scalar, radians) - true anomaly step size
        %
        % Outputs :
        %   positions (NA by 3 matrix) - spherical coordinates
        %                   corresponding to positions occupied by spacecraft
        %                   NA -> # of true anomalies
        %                   positions(:,1) -> planet-centric radii
        %                   positions(:,2) -> latitudes
        %                   positions(:,3) -> longitudes
        function positions = compute_orbit_spherical_coord(r_asc, dec, s_time, ...
            ecc, semi_maj_axis, lon_asc_node, inc, arg_per, step)
    
            nb_pos = 2*pi/step;
            positions = zeros(nb_pos, 3);

            for i = 1:nb_pos

                % compute the coordinates of the position vector in spherical coordinates
                positions(i, :) = ...
                    transpose(rotation_functions.LocalToPlanetaryCoord ...
                    (r_asc, dec, s_time, ecc, semi_maj_axis, lon_asc_node, ...
                    inc, arg_per, step*(i-1)));
    
            end 
        
        end
        
    end
    
end 