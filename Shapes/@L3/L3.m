classdef L3
    %L3 3 node Lagrangian shape function
    
    properties
        ipcount1D % number of integration points per direction
        ipcount % total number of integration points
        
        rectangular % element is strictly rectangular (straight line in x or y)
        
        Nbase % basis functions
        Gbase % basis function derivatives
        
        wbase % integration point weights
        xbase % integration point coordiantes
    end
    
    methods
        function obj = L3(ipcount1D, rect, zeroWeight)
			% Initializations 
            [x1D, w1D] = obj.getIpscheme(ipcount1D, zeroWeight);
            
            if (zeroWeight)
                ipcount1D = ipcount1D+2;
            end
            obj.ipcount1D = ipcount1D;
            obj.ipcount = obj.ipcount1D;
            obj.rectangular = rect;
            
            for i=1:length(x1D)
                x = x1D(i);
                x1 = -1; x2 = 0; x3 = 1;
                Nx = [(x-x2)*(x-x3)/((x1-x2)*(x1-x3)); 
                      (x-x1)*(x-x3)/((x2-x1)*(x2-x3)); 
                      (x-x1)*(x-x2)/((x3-x1)*(x3-x2))];
                Dx = [1*(x-x3)/((x1-x2)*(x1-x3))+(x-x2)*1/((x1-x2)*(x1-x3)); 
                      1*(x-x3)/((x2-x1)*(x2-x3))+(x-x1)*1/((x2-x1)*(x2-x3)); 
                      1*(x-x2)/((x3-x1)*(x3-x2))+(x-x1)*1/((x3-x1)*(x3-x2))];
                
                Nbase(i,:) = Nx;
                Gbase(i,:,1) = Dx;
            end
            
            obj.xbase = x1D;
            obj.wbase = w1D;
            obj.Nbase = Nbase;
            obj.Gbase = Gbase;
        end
        
        function [N, G, w] = getVals(obj, X, Y)
			% Get shape function values, gradients, and accompanying
			% integration weights, providing the nodal coordinate vectors X
			% and Y

            N = obj.Nbase;
            
			if obj.rectangular 
                if (X(3)==X(1)) %vertical
                    G(:,:,1) = obj.Gbase(:,:,1)*2/abs((Y(3)-Y(1)));
                    w = obj.wbase * abs((Y(3)-Y(1)))/2;                    
                else %horizontal
                    G(:,:,1) = obj.Gbase(:,:,1)*abs(2/(X(3)-X(1)));
                    w = obj.wbase * abs((X(3)-X(1)))/2;
                end
            else
                dst = sqrt((X(3)-X(1))^2+(Y(3)-Y(1))^2);
                w = obj.wbase * dst/2;
				G(:,:,1) = obj.Gbase(:,:,1)*2/dst;
			end
        end
        
        function [n, t] = getNormals(obj, X, Y)
			% Get surface normal vector n and tangential vector t for all
			% integration points

            if obj.rectangular 
                dx = X(3)-X(1);
                dy = Y(3)-Y(1);
                
                n = [dy, -dx]./(dy^2+dx^2)^0.5;
                t = [-n(2), n(1)];
                
                n = ones(obj.ipcount1D,1)*n;
                t = ones(obj.ipcount1D,1)*t;
			else
                % not implemented
			end
        end
        
        function xy = getIPGlobal(obj, X,Y)
            % Get global coordinates for integration points

			if obj.rectangular 
                xy(1,:) = X(1) + (0.5*(obj.xbase+1))*(X(3)-X(1));
                xy(2,:) = Y(1) + (0.5*(obj.xbase+1))*(Y(3)-Y(1));
            else
                % not implemented
			end
        end
 
    end
    
    methods (Access = private)
        function [x1D, w1D] = getIpscheme(obj, ipcount1D, zeroWeight)
			% Gauss integration scheme

            if (ipcount1D == 1)
                x1D = 0;
                w1D = 2;
            elseif (ipcount1D == 2)
                x1D = [-1/sqrt(3); 1/sqrt(3)];
                w1D = [1; 1];                
            elseif (ipcount1D == 3)
                x1D = [-sqrt(3/5); 0; sqrt(3/5)];
                w1D = [5/9; 8/9; 5/9];                     
            elseif (ipcount1D ==4)
                x1D = [-sqrt(3/7+2/7*sqrt(6/5)); -sqrt(3/7-2/7*sqrt(6/5)); sqrt(3/7-2/7*sqrt(6/5)); sqrt(3/7+2/7*sqrt(6/5))];
                w1D = [(18-sqrt(30))/36; (18+sqrt(30))/36; (18+sqrt(30))/36; (18-sqrt(30))/36];                     
            elseif (ipcount1D == 5)
                x1D = [-1/3*sqrt(5+2*sqrt(10/7)); -1/3*sqrt(5-2*sqrt(10/7)); 0; 1/3*sqrt(5-2*sqrt(10/7)); 1/3*sqrt(5+2*sqrt(10/7))];
                w1D = [(322-13*sqrt(70))/900; (322+13*sqrt(70))/900; 128/225; (322+13*sqrt(70))/900; (322-13*sqrt(70))/900];                     
            else
                error("Higer order ip schemes not implemented in Shapes.Q9");
            end
            
            if (zeroWeight)
               x1D = [-1; x1D; 1];
               w1D = [0; w1D; 0];
            end
            
        end
    end
end