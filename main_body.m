% MAIN ORBITED BODY CLASS
classdef main_body
    
    properties (Constant)
        
        muSun = 1.32712440042*10^20; % (meters cubed/second squared), gravitational parameter of the Sun
        r0Sun = 6.95700*10^8;        % (meters), mean radius of the Sun
        
    end 
    
    properties
        
        mu;                 % (scalar, meters cubed/second squared), gravitational parameter 
        r0;                 % (scalar, meters), mean radius
        
        RotVecIn;           % (3 by 1 vector, radians/second), rotation vector expressed in inertial frame
        RotVecBF;           % (3 by 1 vector, radians/second), rotation vector expressed in body-fixed frame 
        RotVecPA;           % (3 by 1 vector, radians/second), rotation vector expressed in principal axes frame
        
        InMomt;             % (3 by 1 vector, kilograms.meter squared), principal moments of inertia
        
        HIn;                % (3 by 1 vector, kilograms.meters squared/second), angular momentum expressed in inertial frame
        HAz;                % (scalar, radians), angular momentum azimuth angle w.r.t. inertial frame
        HEl;                % (scalar, radians), angular momentum elevation angle w.r.t. inertial frame
        
        % principal axes frame 3-1-3 rotation angles w.r.t. angular momentum-linked inertial frame)
        psi;                % (scalar, radians), first rotation angle 
        theta;              % (scalar, radians), second rotation angle
        phi;                % (scalar, radians), third rotation angle
        
        % body-fixed frame properties
        rasc;               % (scalar, radians), right ascension angle between body-fixed and inertial frame
        dec;                % (scalar, radians), declination angle between body-fixed and inertial frame
        stime;              % (scalar, radians), sidereal time angle between body-fixed and inertial frame
        
        % rotation matrices
        BF2I;               % (3 by 3 matrix), rotation matrix from body-fixed to inertial frame
        PA2BF;              % (3 by 3 matrix), rotation matrix from principal axes to body-fixed frame
        SCONE2I;            % (3 by 3 matrix), rotation matrix from shadow cone to inertial frame
        
        % Orbital properties (celestial body orbits the Sun)
        OrbElemsMB;         % (6 by 1 vector), contains orbital elements : 
                            % (1) inclination (scalar, radians)
                            % (2) longitude of ascending node (scalar, radians)
                            % (3) argument of periapsis (scalar, radians)
                            % (4) true anomaly (scalar, radians)
                            % (5) semi-major axis (scalar, meters)
                            % (6) eccentricity (scalar, unitless)                         
        EccAnom;            % (scalar, radians), eccentric anomaly 
        
        % Shadow cone properties
        hCone;              % (scalar, meters), shadow cone's height
        ConeApert;          % (scalar, radians), cone's aperture angle
        
    end 
    
    
    methods(Static)
    
        % CONSTRUCTOR
        % Creates a main_body object when called
        %
        % Inputs : 
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
        %
        %   rasc, dec, stime (scalars, radians) - angles corresponding 
        %                   to body-fixed frame 3-1-3 rotation w.r.t. 
        %                   inertial frame
        %
        %   OrbElemsMB (6 by 1 vector), contains orbital elements : 
        %                   (1) inclination (scalar, radians)
        %                   (2) longitude of ascending node (scalar, radians)
        %                   (3) argument of periapsis (scalar, radians)
        %                   (4) true anomaly (scalar, radians)
        %                   (5) semi-major axis (scalar, meters)
        %                   (6) eccentricity (scalar, unitless)
        %
        % Outputs :
        %   mb (main_body) - main_body object 
        function mb = main_body(mu, r0, RotVecBF, InMomt, PA2BF, rasc, ...
                dec, stime, OrbElemsMB)
        
            mb.mu = mu;
            mb.r0 = r0;
            mb.RotVecBF = RotVecBF;
            mb.InMomt = InMomt;
            mb.PA2BF = PA2BF;
            mb.rasc = rasc;
            mb.dec = dec;
            mb.stime = stime;
            mb.OrbElemsMB = OrbElemsMB;
            
            % rotation matrix from body-fixed to inertial frame
            mb.BF2I = rotation_functions.compute_rot_matrix ...
                (mb.rasc, mb.dec, mb.stime, '313');
            
            % rotation vector in principal axes and inertial frames
            mb.RotVecIn = mb.BF2I*mb.RotVecBF;
            mb.RotVecPA = mb.PA2BF'*mb.RotVecBF;
                        
            % angular momentum in principal axes frame
            HPA = diag(InMomt)*mb.RotVecPA;
            
            % angular momentum expressed in inertial frame (constant in our case)
            mb.HIn = mb.BF2I*mb.PA2BF*HPA;
            mb.HAz = mod(2*pi + atan2(mb.HIn(2), mb.HIn(1)), 2*pi);
            mb.HEl = (pi/2) - acos(mb.HIn(3)/norm(mb.HIn));

            % computing Euler angles defining rotation from angular
            % momentum-linked inertial frame to principal axes frame
            [psi, theta, phi] = main_body.BFx2Eul ...
                (mb.rasc, mb.dec, mb.stime, mb.PA2BF, mb.HAz, mb.HEl);
            mb.psi = psi;
            mb.theta = theta;
            mb.phi = phi;
            
            % compute initial eccentric anomaly
            n = mb.OrbElemsMB(4);
            a = mb.OrbElemsMB(5);
            e = mb.OrbElemsMB(6);
            r = (a*(1-e^2))/(1+e*cos(n));
            
            if n <= pi
                mb.EccAnom = acos((a-r)/(a*e));
            else
                mb.EccAnom = 2*pi - acos((a-r)/(a*e));
            end
            
            % compute rotation matrix from shadow cone to inertial frame
            i = mb.OrbElemsMB(1);
            O = mb.OrbElemsMB(2);
            o = mb.OrbElemsMB(3);
            R_IO = rotation_functions.compute_rot_matrix(O, i, o + n, '313'); % rot. mat. from local to inertial frame
            
            mb.SCONE2I(:, 3) = -R_IO(:, 1);
            mb.SCONE2I(:, 2) = R_IO(:, 2);
            mb.SCONE2I(:, 1) = R_IO(:, 3); 
            
            % compute shadow cone initial properties
            mb.hCone = (mb.r0*r)/(mb.r0Sun - mb.r0);
            mb.ConeApert = atan2(mb.r0, mb.hCone);
            
        end 
            
            
            
        % CONVERT EULER ANGLES TO BODY-FIXED FRAME ANGLES
        % Use rotation angles between PA and body-fixed frame to compute
        % rotation angles between body-fixed and inertial frame
        % 
        % Inputs : 
        %   psi, theta, phi (scalars, radians) - angles corresponding 
        %                   to principal axes 3-1-3 rotation w.r.t. angular
        %                   momentum-linked inertial frame
        %
        %   PA2BF (3 by 3 matrix) - rotation matrix from principal axes to 
        %                   body-fixed frame 
        %
        %   HAz, HEl (scalars, radians) - respectively the azimuth and
        %                   elevation of angular momentum vector w.r.t.
        %                   inertial frame
        %
        % Outputs :
        %   rasc, dec, stime (scalar, radians) - angles corresponding 
        %                   to body-fixed frame 3-1-3 rotation w.r.t. 
        %                   inertial frame
        %
        % Comments : this function computes 3-1-3 rotation angles from a 
        %                   rotation matrix using quaternions. When close
        %                   to singularity (corresponding to a rotation 
        %                   angle along the 1-axis close to 0° or 180°),
        %                   the results tend to be less precise. This could 
        %                   be improved in future work.
        function [rasc, dec, stime] = Eul2BFx(psi, theta, phi, PA2BF, HAz, HEl)

            % rotation matrix from body-fixed to principal axes frame
            R_BF2PA = PA2BF';

            % rotation matrix from principal axes to angular momentum-linked 
            % inertial frame
            R_PA2HIn = rotation_functions.compute_rot_matrix ... 
                (psi, theta, phi, '313');
            
            % rotation matrix from angular mom.-linked to inertial frame
            R_HIn2I = rotation_functions.compute_rot_matrix ...
                (HAz, (pi/2)-HEl, 0, '323');
            
            % total rotation matrix
            R = R_HIn2I*R_PA2HIn*R_BF2PA;
            
            % quaternion corresponding to rotation matrix
            q_comp = rotm2quat(R);              
            q = quaternion(q_comp);
            
            % computing 3-1-3 rotation angles from quaternion
            eul = euler(q, 'ZXZ', 'frame');

            % expressing the angles in [0; 2*pi] interval
            rasc = mod(2*pi + eul(1), 2*pi);
            dec = mod(2*pi + eul(2), 2*pi);
            stime = mod(2*pi + eul(3), 2*pi);

        end 
            
            
            
        % CONVERT BODY-FIXED FRAME ANGLES TO EULER ANGLES
        % Use rotation angles between body-fixed and inertial frame to compute
        % rotation angles between principal axes and body-fixed frame
        % 
        % Inputs : 
        %   rasc, dec, stime (scalar, radians) - angles corresponding 
        %                   to body-fixed frame 3-1-3 rotation w.r.t. 
        %                   inertial frame
        %
        %   PA2BF (3 by 3 matrix) - rotation matrix from principal axes to 
        %                   body-fixed frame 
        %
        %   HAz, HEl (scalars, radians) - respectively the azimuth and
        %                   elevation of angular momentum vector w.r.t.
        %                   inertial frame
        %
        % Outputs :
        %   psi, theta, phi (scalars, radians) - angles corresponding 
        %                   to principal axes 3-1-3 rotation w.r.t. angular
        %                   momentum-linked inertial frame
        %
        % Comments : this function computes 3-1-3 rotation angles from a 
        %                   rotation matrix using quaternions. When close
        %                   to singularity (corresponding to a rotation 
        %                   angle along the 1-axis close to 0° or 180°),
        %                   the results tend to be less precise. This could 
        %                   be improved in future work.
        function [psi, theta, phi] = BFx2Eul(rasc, dec, stime, PA2BF, HAz, HEl)

            % rotation matrix from PA to body-fixed frame
            R_PA2BF = PA2BF;

            % rotation matrix from body-fixed frame to inertial frame
            R_BF2I = rotation_functions.compute_rot_matrix(rasc, dec, stime, '313');
            
            % rotation matrix from inertial to angular momentum inertial frame
            R_I2HIn = (rotation_functions.compute_rot_matrix ...
                (HAz, (pi/2)-HEl, 0, '323'))';

            % total rotation matrix
            R = R_I2HIn*R_BF2I*R_PA2BF;
            
            % quaternion corresponding to rotation matrix
            q_comp = rotm2quat(R);
            q = quaternion(q_comp);
            
            % computing 3-1-3 rotation angles from quaternion
            eul = euler(q, 'ZXZ', 'frame');
            
            % expressing the angles in [0; 2*pi] interval
            psi = mod(2*pi + eul(1), 2*pi);
            theta = mod(2*pi + eul(2), 2*pi);
            phi = mod(2*pi + eul(3), 2*pi);
        
        end 
            
            
            
        % COMPUTE ORIENTATION ANGLES AND ROTATION VECTOR AT NEXT TIME INCREMENT
        % Initial conditions taken from principal axes Euler angles state 
        % in a  main_body given as an entry parameter. Also updates body-
        % fixed frame angles and rotation vector.
        %
        % Inputs : 
        %   mb (main_body) - main_body object
        %
        %   dt (seconds) - time step
        %
        % Outputs :
        %   mb (main_body) - modified main_body object given as an entry
        function mb = update_angles(mb, dt)

            % compute principal axes frame orientation angles w.r.t.
            % angular momentum-linked inertial frame at next time increment
            [psi, theta, phi] = main_body.SolveEulAngles ...
                ([mb.psi; mb.theta; mb.phi], mb.InMomt, dt, ...
                norm(mb.HIn));
            
            mb.psi = psi;
            mb.theta = theta;
            mb.phi = phi;

            % update body-fixed frame orientation angles w.r.t. inertial frame
            [r, d, s] = main_body.Eul2BFx ...
                (mb.psi, mb.theta, mb.phi, mb.PA2BF, mb.HAz, mb.HEl);
            
            mb.rasc = r;
            mb.dec = d;
            mb.stime = s;
            
            % update rotation matrix from body-fixed to inertial frame
            mb.BF2I = rotation_functions.compute_rot_matrix ...
                (mb.rasc, mb.dec, mb.stime, '313');
            
            % compute rotation vector at next time increment
            mb.RotVecPA = main_body.SolveEulEq ...
                (mb.InMomt, mb.RotVecPA, [0;0;0], [0, dt]);
            mb.RotVecBF = mb.PA2BF*mb.RotVecPA;
            mb.RotVecIn = mb.BF2I*mb.RotVecBF; 
            
            % compute eccentric anomaly around the Sun at next time increment
            a = mb.OrbElemsMB(5);
            e = mb.OrbElemsMB(6);
            fac = sqrt(mb.muSun/(a^3));
            E = mb.EccAnom;   
            mb.EccAnom = mod(E + (fac*dt)/(1 - e*cos(E)), 2*pi); 
            
            % compute true anomaly around the Sun at next time increment
            E = mb.EccAnom;
            CosNu = (cos(E) - e)/(1 - e*cos(E));
            SinNu = (sqrt(1 - e^2)*sin(E))/(1 - e*cos(E));
            mb.OrbElemsMB(4) = mod(2*pi + atan2(SinNu, CosNu), 2*pi);
            
            % update rotation matrix from shadow cone to inertial frame
            i = mb.OrbElemsMB(1);
            O = mb.OrbElemsMB(2);
            o = mb.OrbElemsMB(3);
            n = mb.OrbElemsMB(4);
            rad = (a*(1-e^2))/(1+e*cos(n)); 
            R_IO = rotation_functions.compute_rot_matrix(O, i, o + n, '313'); % rot. mat. from local to inertial frame
            
            mb.SCONE2I(:, 3) = -R_IO(:, 1);
            mb.SCONE2I(:, 2) = R_IO(:, 2);
            mb.SCONE2I(:, 1) = R_IO(:, 3);  
            
            % update shadow cone initial properties
            mb.hCone = (mb.r0*rad)/(mb.r0Sun - mb.r0);
            mb.ConeApert = atan2(mb.r0, mb.hCone);
            
        end
        
        
        
        % SOLVE FOR PRINCIPAL AXES ORIENTATION ANGLES AT NEXT TIME INCREMENT
        % 
        % Inputs : 
        %   EulAngles_init (3 by 1 vector, radians) - initial Euler 
        %                   (rotation) angles : 
        %                   (1) Psi (first rotation along 3-axis)
        %                   (2) Theta (second rotation along 1-axis)
        %                   (3) Phi (third rotation along 3-axis)
        %
        %   InMomt (vector of size 3, kilograms.meter squared) - body's
        %                   principal moments of inertia, organized in 
        %                   accordance to the order of the principal axes 
        %                   X1, X2, X3 in PA2BF matrix (described below) :
        %                   InMomt(1) -> along X1 axis 
        %                   InMomt(2) -> along X2 axis 
        %                   InMomt(3) -> along X3 axis  
        %
        %   dt (scalar, seconds) - time step
        %
        %   h (scalar, kilograms.meters squared/second) - body's angular 
        %                   momentum norm
        %
        % Outputs : 
        %   psi, theta, phi (scalars, radians) - angles corresponding 
        %                   to principal axes 3-1-3 rotation w.r.t. angular
        %                   momentum-linked inertial frame
        function [psi, theta, phi] = SolveEulAngles(EulAngles_init, InMomt, dt, h)
            
            sol = ode45...
                (@(t, EulAngles) main_body.RHS_EulAngles ...
                    (EulAngles, InMomt, h), [0 dt], EulAngles_init); 
            
            angles = sol.y(:, size(sol.y, 2));
            
            psi = mod(angles(1), 2*pi);
            theta = mod(angles(2), 2*pi);
            phi = mod(angles(3), 2*pi);
            
        end
        
        
        
        % COMPUTE RIGHT-HAND SIDE OF DIFFERENTIAL EQUATION DESCRIBING 
        % PRINCIPAL AXES ORIENTATION ANGLES W.R.T. ANGULAR MOMENTUM-LINKED 
        % INERTIAL FRAME EVOLUTION
        %
        % Inputs :
        %   EulAngles (3 by 1 vector, radians) - initial Euler 
        %                   (rotation) angles : 
        %                   (1) Psi (first rotation along 3-axis)
        %                   (2) Theta (second rotation along 1-axis)
        %                   (3) Phi (third rotation along 3-axis)
        %   
        %   InMomt (vector of size 3, kilograms.meter squared) - body's
        %                   principal moments of inertia, organized in 
        %                   accordance to the order of the principal axes 
        %                   X1, X2, X3 in PA2BF matrix (described below) :
        %                   InMomt(1) -> along X1 axis 
        %                   InMomt(2) -> along X2 axis 
        %                   InMomt(3) -> along X3 axis 
        %
        %   h (scalar, kilograms.meters squared/second) - body's angular 
        %                   momentum norm
        %
        % Outputs :
        %   f (3 by 1 vector, radians/second) - right-hand side of
        %                   differential equation value (corresponding to
        %                   instantaneous angles rates)
        function f = RHS_EulAngles(EulAngles, InMomt, h)
            
            f = zeros(3, 1);
            
            I1 = InMomt(1);
            I2 = InMomt(2);
            I3 = InMomt(3);
            
            psi = EulAngles(1);
            theta = EulAngles(2);
            phi = EulAngles(3); %#ok<NASGU>
            
            f(1) = h*(((cos(psi)^2)/I2) + ((sin(psi)^2)/I1));
            f(2) = h*((1/I1) - (1/I2))*sin(theta)*sin(psi)*cos(psi);
            f(3) = h*((1/I3) - ((cos(psi)^2)/I2) - ((sin(psi)^2)/I1))*cos(theta);
            
        end
        
        
        
        % SOLVE EULER'S EQUATIONS FOR ROTATION VECTOR
        %
        % Inputs :
        %   InMomt (vector of size 3, kilograms.meter squared) - body's
        %                   principal moments of inertia, organized in 
        %                   accordance to the order of the principal axes 
        %                   X1, X2, X3 in PA2BF matrix (described below) :
        %                   InMomt(1) -> along X1 axis 
        %                   InMomt(2) -> along X2 axis 
        %                   InMomt(3) -> along X3 axis 
        %
        %   omega_init (vector of size 3, radians/second) - initial
        %                   rotation vector expressed in principal 
        %                   axes frame
        %
        %   Torque (vector of size 3, kilograms.meters squared/second
        %                   squared) - torque vector acting on body
        %                   expressed in principal axes frame
        %
        %   dt (scalar, seconds) - time step
        %
        % Outputs :
        %   omega (vector of size 3, radians/second) - rotation vector at
        %                   next time increment, expressed in principal axes 
        %                   frame 
        function omega = SolveEulEq(InMomt, omega_init, Torque, dt)
            
            sol = ode45...
                (@(t, omega) main_body.RHS_EulEq(omega, InMomt, Torque), dt, ...
                omega_init); 
            
            omega = sol.y(:, size(sol.y, 2));
            
        end 
        
       
        
        % COMPUTE RIGHT-HAND SIDE OF EULER'S EQUATION FOR ROTATION VECTOR
        %
        % Inputs : 
        %   omega (vector of size 3, radians/second) - initial
        %                   rotation vector expressed in principal 
        %                   axes frame
        %
        %   InMomt (vector of size 3, kilograms.meter squared) - body's
        %                   principal moments of inertia, organized in 
        %                   accordance to the order of the principal axes 
        %                   X1, X2, X3 in PA2BF matrix (described below) :
        %                   InMomt(1) -> along X1 axis 
        %                   InMomt(2) -> along X2 axis 
        %                   InMomt(3) -> along X3 axis 
        %
        %   Torque (vector of size 3, kilograms.meters squared/second
        %                   squared) - torque vector acting on body
        %                   expressed in principal axes frame
        %
        % Outputs :
        %   f (3 by 1 vector, radians/second squared) - value of right-hand 
        %                   side of Euler's equation (corresponding to
        %                   instantaneous angles rates)   
        function f = RHS_EulEq(omega, InMomt, Torque)
            
            f = zeros(3, 1);
            
            I1 = InMomt(1);
            I2 = InMomt(2);
            I3 = InMomt(3);
            
            o1 = omega(1);
            o2 = omega(2);
            o3 = omega(3);
            
            M1 = Torque(1);
            M2 = Torque(2);
            M3 = Torque(3);
            
            f(1) = (-(I3-I2)*o2*o3 + M1)/I1;
            f(2) = (-(I1-I3)*o1*o3 + M2)/I2;
            f(3) = (-(I2-I1)*o1*o2 + M3)/I3;
            
        end 

    end 
    
end 