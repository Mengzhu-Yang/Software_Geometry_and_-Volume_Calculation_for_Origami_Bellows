function [solution] = calculatesolution(phi1Deg,phi2Deg,n)

%R = app.R.Value;  
 %d1overR = app.d1R.Value;  
 
%% Judge how many solutions 
% %convert angles
phi1 = phi1Deg*pi/180;
phi2 = phi2Deg*pi/180;

FOLD_count = 0;
noWrite = 0;
%Constraints on phi1 and phi2

      
      if phi1Deg <= 180/n   
        %  disp('0 solution');
            noWrite = 1;
            solution = 0;
      end
      
      if phi1Deg <= phi2Deg  
        %  disp('0 solution');
            noWrite = 1;
             solution = 0;
      end
      
      if phi1Deg <= 180/4+180/(2*n) 
          if phi2Deg>0 && phi2Deg<(phi1Deg-180/n)
%               x1 = tan(phi2)/tan(phi1);
%          x2= tan(pi/n);
%          beta1 = 2*atan(1/(2*x2)*(1-x1-(sqrt((x1-1)^2 - 4*x1*x2^2))));
%          beta2 = 2*atan(1/(2*x2)*(1-x1+(sqrt((x1-1)^2 - 4*x1*x2^2))));
%          beta = beta1;        
%             %1st solution (S2 more deployed)： beta = beta1
%             %2nd solution (S1 more compressed)： beta = beta2
%             % alpha is angular deflection from vertical (z-axis) 
%             %pi-alpha is the opening angle between adjacent bands of folds
%            alpha = asin(tan(beta/2)/tan(phi2));
%         
%           %%%The angle deflection from Z axis
%             if alpha > 80/180*pi
%       
%             noWrite = 1;
%             solution = 0;
%             else
         %   disp('1 solution');
         
            noWrite = 0;
             solution = 1;
         %   end 
          else
         %   disp('0 solution');          
            noWrite = 1;
             solution = 0;
          end
      end
      
      if phi1Deg>180/4+180/(2*n) 
        % Boundary of bistable region         
           limitangle1=(phi1Deg-180/n);
           limitangle2=(atan(tan(phi1)*(1-sin(pi/n))/(1+sin(pi/n))))/pi*180; 
         if phi2Deg>=phi1Deg-180/n && phi2Deg<(atan(tan(phi1)*(1-sin(pi/n))/(1+sin(pi/n))))/pi*180       
%          x1 = tan(phi2)/tan(phi1);
%          x2= tan(pi/n);
%          beta1 = 2*atan(1/(2*x2)*(1-x1-(sqrt((x1-1)^2 - 4*x1*x2^2))));
%          beta2 = 2*atan(1/(2*x2)*(1-x1+(sqrt((x1-1)^2 - 4*x1*x2^2))));
%          beta = beta2;        
%             %1st solution (S2 more deployed)： beta = beta1
%             %2nd solution (S1 more compressed)： beta = beta2
%             % alpha is angular deflection from vertical (z-axis) 
%             %pi-alpha is the opening angle between adjacent bands of folds
%            alpha = asin(tan(beta/2)/tan(phi2));
%         
%           %%%The angle deflection from Z axis
%             if alpha > 70/180*pi
%          %   disp('Flatfoldable');
%             noWrite = 1;
%             solution = 0;
%             else
          %      disp('2 solutions');
            solution = 2; 
%            end
            %1st solution (S2)： beta = beta1
            %2nd solution (S1)： beta = beta2
         elseif phi2>0 && phi2<phi1-pi/n
%                  x1 = tan(phi2)/tan(phi1);
%          x2= tan(pi/n);
%          beta1 = 2*atan(1/(2*x2)*(1-x1-(sqrt((x1-1)^2 - 4*x1*x2^2))));
%          beta2 = 2*atan(1/(2*x2)*(1-x1+(sqrt((x1-1)^2 - 4*x1*x2^2))));
%          beta = beta1;        
%             %1st solution (S2 more deployed)： beta = beta1
%             %2nd solution (S1 more compressed)： beta = beta2
%             % alpha is angular deflection from vertical (z-axis) 
%             %pi-alpha is the opening angle between adjacent bands of folds
%            alpha = asin(tan(beta/2)/tan(phi2));
%         
%           %%%The angle deflection from Z axis
%             if alpha > 70/180*pi
%          %   disp('Flatfoldable');
%             noWrite = 1;
%             solution = 0;
%             else
%           %  disp('1 solution');
            solution = 1;
        %    end
         elseif  phi2>atan(tan(phi1)*(1-sin(pi/n))/(1+sin(pi/n))) 
         %   disp('0 solution');
           solution = 0;
         else
             solution = 0;
         end
      end
end
