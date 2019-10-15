classdef finite_difference_schemes
    
    methods(Static)
        
        % COMPUTING GRADIENT IN CARTESIAN COORDINATES
        % Gradient is order 1 on the bounds, order 2 inside the domain
        %
        % Inputs : 
        %   scalar_field (NX by NY by NZ matrix) - contains the values of 
        %               the function to differentiate
        %               NX -> along first dimension (lines)
        %               NY -> along second dimension (columns)
        %               NZ -> along third dimension
        %                     
        %   dX, dY, dZ (scalars) - space steps (constants) along first, 
        %               second, and third dimension of the matrix 
        %
        % Outputs : 
        %   gradX, gradY, gradZ (NX by NY by NZ matrices) - gradients along 
        %               first, second, and third dimension of the entry matrix
        %               scalar_field
        function [gradX, gradY, gradZ] = gradient_o2(scalar_field, dX, dY, dZ)
        
            NX = size(scalar_field, 1);
            NY = size(scalar_field, 2);
            NZ = size(scalar_field, 3);
            
            % X direction
            gradX = zeros(NX, NY, NZ);

            % Y direction
            gradY = zeros(NX, NY, NZ);

            % Z direction
            gradZ = zeros(NX, NY, NZ);
            
            for i = 1:NX
                
                for j = 1:NY
                    
                    for k = 1:NZ
                        
                        % computing X component of the gradient 
                        if (i == 1)

                            gradX(i, j, k) = (scalar_field(i+1, j, k) - scalar_field(i, j, k)) ...
                                /(dX);

                        elseif (i == NX)

                            gradX(i, j, k) = (scalar_field(i, j, k) - scalar_field(i-1, j, k)) ...
                                /(dX);
                        else

                            gradX(i, j, k) = (scalar_field(i+1, j, k) - scalar_field(i-1, j, k)) ...
                                /(dX);

                        end 
                        
                        % computing Y component of the gradient 
                        if (j == 1)

                            gradY(i, j, k) = (scalar_field(i, j+1, k) - scalar_field(i, j, k)) ...
                                /(dY);

                        elseif (j == NY)

                            gradY(i, j, k) = (scalar_field(i, j, k) - scalar_field(i, j-1, k)) ...
                                /(dY);
                        else
                            
                            gradY(i, j, k) = (scalar_field(i, j+1, k) - scalar_field(i, j-1, k)) ...
                                /(dY);

                        end 
                        
                        % computing Z component of the gradient 
                        if (k == 1)

                            gradZ(i, j, k) = (scalar_field(i, j, k+1) - scalar_field(i, j, k)) ...
                                /(dZ);

                        elseif (k == NZ)

                            gradZ(i, j, k) = (scalar_field(i, j, k) - scalar_field(i, j, k-1)) ...
                                /(dZ);
                        else

                            gradZ(i, j, k) = (scalar_field(i, j, k+1) - scalar_field(i, j, k-1)) ...
                                /(dZ);

                        end                   
                        
                    end 
                    
                end 
                
            end 
                    
        end
        
        
        
        % COMPUTING GRADIENT IN CARTESIAN COORDINATES
        % Gradient is order 1 on the bounds, order 2 one step out of the 
        % bounds, order 4 inside the domain
        %
        % Inputs : 
        %   scalar_field (NX by NY by NZ matrix) - contains the values of 
        %               the function to differentiate
        %               NX -> along first dimension (lines)
        %               NY -> along second dimension (columns)
        %               NZ -> along third dimension
        %                     
        %   dX, dY, dZ (scalars) - space steps (constants) along first, 
        %               second, and third dimension of the matrix 
        %
        % Outputs : 
        %   gradX, gradY, gradZ (NX by NY by NZ matrices) - gradients along 
        %               first, second, and third dimension of the entry matrix
        %               scalar_field
        function [gradX, gradY, gradZ] = gradient_o4(scalar_field, dX, dY, dZ)
        
            NX = size(scalar_field, 1);
            NY = size(scalar_field, 2);
            NZ = size(scalar_field, 3);
            
            % X direction
            gradX = zeros(NX, NY, NZ);

            % Y direction
            gradY = zeros(NX, NY, NZ);

            % Z direction
            gradZ = zeros(NX, NY, NZ);
            
            for i = 1:NX
                
                for j = 1:NY
                    
                    for k = 1:NZ
                        
                        % computing X component of the gradient 
                        if (i == 1)

                            gradX(i, j, k) = (scalar_field(i+1, j, k) - scalar_field(i, j, k)) ...
                                /(dX);

                        elseif (i == NX)

                            gradX(i, j, k) = (scalar_field(i, j, k) - scalar_field(i-1, j, k)) ...
                                /(dX);
                        
                        elseif (i == 2 || i == NX-1)

                            gradX(i, j, k) = (scalar_field(i+1, j, k) - scalar_field(i-1, j, k)) ...
                                /(2*dX);
                            
                        else 
                            
                            gradX(i, j, k) = (-scalar_field(i+2, j, k) + 8*scalar_field(i+1, j, k) ...
                                - 8*scalar_field(i-1, j, k) + scalar_field(i-2, j, k))/(12*dX);

                        end 
                        
                        % computing Y component of the gradient 
                        if (j == 1)

                            gradY(i, j, k) = (scalar_field(i, j+1, k) - scalar_field(i, j, k)) ...
                                /(dY);

                        elseif (j == NY)

                            gradY(i, j, k) = (scalar_field(i, j, k) - scalar_field(i, j-1, k)) ...
                                /(dY);
                        
                        elseif (j == 2 || j == NY-1) 
                            
                            gradY(i, j, k) = (scalar_field(i, j+1, k) - scalar_field(i, j-1, k)) ...
                                /(2*dY);
                            
                        else 
                            
                            gradY(i, j, k) = (-scalar_field(i, j+2, k) + 8*scalar_field(i, j+1, k) ...
                                - 8*scalar_field(i, j-1, k) + scalar_field(i, j-2, k))/(12*dY);

                        end 
                        
                        % computing Z component of the gradient 
                        if (k == 1)

                            gradZ(i, j, k) = (scalar_field(i, j, k+1) - scalar_field(i, j, k)) ...
                                /(dZ);

                        elseif (k == NZ)

                            gradZ(i, j, k) = (scalar_field(i, j, k) - scalar_field(i, j, k-1)) ...
                                /(dZ);
                            
                        elseif (k == 2 || k == NZ-1) 
                            
                            gradZ(i, j, k) = (scalar_field(i, j, k+1) - scalar_field(i, j, k-1)) ...
                                /(2*dZ);
                            
                        else
                            
                            gradZ(i, j, k) = (-scalar_field(i, j, k+2) + 8*scalar_field(i, j, k+1) ...
                                - 8*scalar_field(i, j, k-1) + scalar_field(i, j, k-2))/(12*dZ);

                        end                   
                        
                    end 
                    
                end 
                
            end 
                    
        end
        
        
        
        % COMPUTING GRADIENT IN SPHERICAL COORDINATES
        % Gradient is order 1 on the bounds, order 2 inside the domain
        % This version considers a 3D function that is generated over the
        % entire map of a body (1st and 2nd dimensions of the matrix scalar_field), 
        % for several planet-centric radii (3rd dimension of the matrix scalar_field). 
        % Gradient is order 1 on the bounds, order 2 inside the domain
        %
        % Inputs : 
        %   coord_r (vector of size NR) - containing the values of 
        %               planet-centric radii used to generate the matrix 
        %               scalar_field 
        %   
        %   r_step (scalar) - defining the radius step used to generate the 
        %               matrix scalar_field
        %
        %   scalar_field (NP by NT by NR matrix) - contains the values 
        %               of the function to differentiate. 
        %               NP # of latitudes (along lines)
        %               NT # of longitudes (along columns)
        %               NR # of radii
        %
        % Outputs : 
        %   gradP, gradT, gradR (NP by NT by NR matrices) - gradients along 
        %               southward, eastward, and upward directions
        %
        % Comments : 
        %   This functions works fine for latitudes different than the
        %   poles. At these locations, the gradient along the eastward
        %   direction will go to infinity, since its computation involves a
        %   "1/sin(latitude)" term.
        function [gradP, gradT, gradR] = spherical_gradient(coord_r, r_step, scalar_field)

            NP = size(scalar_field, 1);
            NT = size(scalar_field, 2);
            NR = size(scalar_field, 3);

            % eastward direction
            gradT = zeros(NP, NT, NR);

            % southward direction
            gradP = zeros(NP, NT, NR);

            % upward direction
            gradR = zeros(NP, NT, NR);

            deltaR = r_step;
            deltaP = 180/(NP-1);
            deltaT = 360/(NT-1);

            for i = 1:NP

                for j = 1:NT

                    for k = 1:NR

                        % computing Phi (southward) component of the gradient 
                        if (i == 1)

                            gradP(i, j, k) = (scalar_field(i+1, j, k) - scalar_field(i, j, k)) ...
                                /(deltaP);

                        elseif (i == NP)

                            gradP(i, j, k) = (scalar_field(i, j, k) - scalar_field(i-1, j, k)) ...
                                /(deltaP);
                        else

                            gradP(i, j, k) = (scalar_field(i+1, j, k) - scalar_field(i-1, j, k)) ...
                                /(2*deltaP);

                        end 

                        gradP(i, j, k) =  (1./coord_r(k))*gradP(i, j, k);

                        % computing Theta (eastward) component of the gradient 
                        if (j == 1)

                            gradT(i, j, k) = (scalar_field(i, j+1, k) - scalar_field(i, j, k)) ...
                                /deltaT;

                        elseif (j == NT)

                            gradT(i, j, k) = (scalar_field(i, j, k) - scalar_field(i, j-1, k)) ...
                                /deltaT;

                        else

                            gradT(i, j, k) = (scalar_field(i, j+1, k) - scalar_field(i, j-1, k)) ...
                                /(2*deltaT);

                        end

                        parameter = 1/(coord_r(k)*sin((i-1)*deltaP*pi/180));
                        gradT(i, j, k) = parameter*gradT(i, j, k);

                        % computing R (upward) component of the gradient 
                        if (k == 1)

                            gradR(i, j, k) = (scalar_field(i, j, k+1) - scalar_field(i, j, k)) ...
                                /deltaR;

                        elseif (k == NR)

                            gradR(i, j, k) = (scalar_field(i, j, k) - scalar_field(i, j, k-1)) ...
                                /deltaR;

                        else

                            gradR(i, j, k) = (scalar_field(i, j, k+1) - scalar_field(i, j, k-1)) ...
                                /(2*deltaR);

                        end 

                    end 

                end

            end    

        end
        
        
        
        % COMPUTING GRADIENT IN SPHERICAL COORDINATES
        % Gradient is order 1 on the bounds, order 2 inside the domain
        % This version considers a 3D function that is generated over the
        % a range of latitudes, longitudes and radii (1st, 2nd and 3rd 
        % dimensions of the matrix scalar_field). 
        % Gradient is order 1 on the bounds, order 2 inside the domain
        %
        % Inputs : 
        %   coord_r (vector of size NR) - contains the values of 
        %               planet-centric radii used to generate the matrix 
        %               scalar_field 
        %
        %   lats (vector of size NP) - contains the values of latitudes
        %               used to generate the matrix scalar_field
        %
        %   lat_step (scalar) - defines the latitude step used to generate
        %               the matrix scalar_field
        %
        %   lon_step (scalar) - defines the longitude step used to generate
        %               the matrix scalar_field
        % 
        %   r_step (scalar) - defines the radius step used to generate the 
        %               matrix scalar_field
        %
        %   scalar_field (NP by NT by NR matrix) - contains the values 
        %               of the function to differentiate.
        %
        % Outputs : 
        %   gradP, gradT, gradR (NP by NT by NR matrices) - gradients along 
        %               eastward, southward, and upward directions
        %
        % Comments : 
        %   This functions works fine for latitudes different than the
        %   poles. At this location, the gradient along the eastward
        %   direction will go to infinity, since its computation involves a
        %   "1/sin(latitude)" term.
        function [gradP, gradT, gradR] = spherical_gradient_2(coord_r, ...
            lats, lat_step, lon_step, r_step, scalar_field)
         
            NP = size(scalar_field, 1);
            NT = size(scalar_field, 2);
            NR = size(scalar_field, 3);

            % eastward direction
            gradT = zeros(NP, NT, NR);

            % southward direction
            gradP = zeros(NP, NT, NR);

            % upward direction
            gradR = zeros(NP, NT, NR);

            deltaR = r_step;
            deltaP = lat_step;
            deltaT = lon_step;

            for i = 1:NP

                for j = 1:NT

                    for k = 1:NR

                        % computing Phi (southward) component of the gradient 
                        if (i == 1)

                            gradP(i, j, k) = (scalar_field(i+1, j, k) - scalar_field(i, j, k)) ...
                                /(deltaP);

                        elseif (i == NP)

                            gradP(i, j, k) = (scalar_field(i, j, k) - scalar_field(i-1, j, k)) ...
                                /(deltaP);
                        else

                            gradP(i, j, k) = (scalar_field(i+1, j, k) - scalar_field(i-1, j, k)) ...
                                /(2*deltaP);

                        end 

                        gradP(i, j, k) =  (1./coord_r(k))*gradP(i, j, k);

                        % computing Theta (eastward) component of the gradient 
                        if (j == 1)

                            gradT(i, j, k) = (scalar_field(i, j+1, k) - scalar_field(i, j, k)) ...
                                /deltaT;

                        elseif (j == NT)

                            gradT(i, j, k) = (scalar_field(i, j, k) - scalar_field(i, j-1, k)) ...
                                /deltaT;

                        else

                            gradT(i, j, k) = (scalar_field(i, j+1, k) - scalar_field(i, j-1, k)) ...
                                /(2*deltaT);

                        end
                        
                        if lats(i) == 0 || lats(i) == 180
                            parameter = NaN;
                        else 
                            parameter = 1/(coord_r(k)*sin(lats(i)*pi/180));
                        end 
                        
                        gradT(i, j, k) = parameter*gradT(i, j, k);

                        % computing R (upward) component of the gradient 
                        if (k == 1)

                            gradR(i, j, k) = (scalar_field(i, j, k+1) - scalar_field(i, j, k)) ...
                                /deltaR;

                        elseif (k == NR)

                            gradR(i, j, k) = (scalar_field(i, j, k) - scalar_field(i, j, k-1)) ...
                                /deltaR;

                        else

                            gradR(i, j, k) = (scalar_field(i, j, k+1) - scalar_field(i, j, k-1)) ...
                                /(2*deltaR);

                        end 

                    end 

                end

            end    

        end
        
        
        
        % COMPUTING GRADIENT IN SPHERICAL COORDINATES
        % Gradient is order 1 on the bounds, order 2 one step out of the
        % bounds, order 4 inside the domain
        % This version considers a 3D function that is generated over the
        % a range of latitudes, longitudes and radii (1st, 2nd and 3rd 
        % dimensions of the matrix scalar_field). 
        % Gradient is order 1 on the bounds, order 2 inside the domain
        %
        % Inputs : 
        %   coord_r (vector of size NR) - contains the values of 
        %               planet-centric radii used to generate the matrix 
        %               scalar_field 
        %
        %   lats (vector of size NP) - contains the values of latitudes
        %               used to generate the matrix scalar_field
        %
        %   lat_step (scalar) - defines the latitude step used to generate
        %               the matrix scalar_field
        %
        %   lon_step (scalar) - defines the longitude step used to generate
        %               the matrix scalar_field
        % 
        %   r_step (scalar) - defines the radius step used to generate the 
        %               matrix scalar_field
        %
        %   scalar_field (NP by NT by NR matrix) - contains the values 
        %               of the function to differentiate.
        %
        % Outputs : 
        %   gradP, gradT, gradR (NP by NT by NR matrices) - gradients along 
        %               eastward, southward, and upward directions
        %
        % Comments : 
        %   This functions works fine for latitudes different than the
        %   poles. At this location, the gradient along the eastward
        %   direction will go to infinity, since its computation involves a
        %   "1/sin(latitude)" term.
        function [gradP, gradT, gradR] = spherical_gradient_o4(coord_r, ...
            lats, lat_step, lon_step, r_step, scalar_field)
         
            NP = size(scalar_field, 1);
            NT = size(scalar_field, 2);
            NR = size(scalar_field, 3);

            % eastward direction
            gradT = zeros(NP, NT, NR);

            % southward direction
            gradP = zeros(NP, NT, NR);

            % upward direction
            gradR = zeros(NP, NT, NR);

            deltaR = r_step;
            deltaP = lat_step;
            deltaT = lon_step;

            for i = 1:NP

                for j = 1:NT

                    for k = 1:NR

                        % computing Phi (southward) component of the gradient 
                        if (i == 1)

                            gradP(i, j, k) = (scalar_field(i+1, j, k) - scalar_field(i, j, k)) ...
                                /(deltaP);

                        elseif (i == NP)

                            gradP(i, j, k) = (scalar_field(i, j, k) - scalar_field(i-1, j, k)) ...
                                /(deltaP);
                        
                        elseif (i == 2 || i == NP-1)

                            gradP(i, j, k) = (scalar_field(i+1, j, k) - scalar_field(i-1, j, k)) ...
                                /(2*deltaP);
                            
                        else 
                            
                            gradP(i, j, k) = (-scalar_field(i+2, j, k) + 8*scalar_field(i+1, j, k) ...
                                - 8*scalar_field(i-1, j, k) + scalar_field(i-2, j, k))/(12*deltaP);                            

                        end 

                        gradP(i, j, k) =  (1./coord_r(k))*gradP(i, j, k);

                        % computing Theta (eastward) component of the gradient 
                        if (j == 1)

                            gradT(i, j, k) = (scalar_field(i, j+1, k) - scalar_field(i, j, k)) ...
                                /deltaT;

                        elseif (j == NT)

                            gradT(i, j, k) = (scalar_field(i, j, k) - scalar_field(i, j-1, k)) ...
                                /deltaT;

                        elseif (j == 2 || j == NT-1)

                            gradT(i, j, k) = (scalar_field(i, j+1, k) - scalar_field(i, j-1, k)) ...
                                /(2*deltaT);
                            
                        else
                            
                            gradT(i, j, k) = (-scalar_field(i, j+2, k) + 8*scalar_field(i, j+1, k) ...
                                - 8*scalar_field(i, j-1, k) + scalar_field(i, j-2, k))/(12*deltaT);
                            

                        end

                        if lats(i) == 0 || lats(i) == 180
                            parameter = NaN;
                        else 
                            parameter = 1/(coord_r(k)*sin(lats(i)*pi/180));
                        end 
                        gradT(i, j, k) = parameter*gradT(i, j, k);

                        % computing R (upward) component of the gradient 
                        if (k == 1)

                            gradR(i, j, k) = (scalar_field(i, j, k+1) - scalar_field(i, j, k)) ...
                                /deltaR;

                        elseif (k == NR)

                            gradR(i, j, k) = (scalar_field(i, j, k) - scalar_field(i, j, k-1)) ...
                                /deltaR;

                        elseif (k == 2 || k == NR-1)

                            gradR(i, j, k) = (scalar_field(i, j, k+1) - scalar_field(i, j, k-1)) ...
                                /(2*deltaR);
                            
                        else 

                            gradR(i, j, k) = (-scalar_field(i, j, k+2) + 8*scalar_field(i, j, k+1) ...
                                - 8*scalar_field(i, j, k-1) + scalar_field(i, j, k-2))/(12*deltaR);

                        end                         

                    end 

                end

            end    

        end    

    end
    
end 
