function [Rs] = MiuraFlatS1(phi1,phi2,r,n,m,d1overr)
%  n=4
%  m=3
%  phi1=75/180*pi
%  phi2=30/180*pi
%  d1overr=1
%  r=40
       
        unitcell = [];
        unitcell2 = [];
        faces = [];
        vertices = [];
        
         %% %calculate without scalefactor
         d1 = d1overr*r;
         d2 = 2*pi*r/n - d1;
         H = (2*pi*r/n-2*d2)/(1/tan(phi2) - 1/tan(phi1));
 
          if d1 >0 && d2>0 && H>0
          x1 = tan(phi2)/tan(phi1);
          x2= tan(pi/n);
             %1st solution (S2)： beta = beta1
            %2nd solution (S1)： beta = beta2
         beta1 = 2*atan(1/(2*x2)*(1-x1-(sqrt((x1-1)^2 - 4*x1*x2^2))));
         beta2 = 2*atan(1/(2*x2)*(1-x1+(sqrt((x1-1)^2 - 4*x1*x2^2))));
         beta = beta2;
         alpha =pi/2;
           %alpha = asin(tan(beta/2)/tan(phi2));
            %%%alpha will sometimes be imaginary at flatfoldable state
         % alpha = real(alpha)
        % nodes of single unit cell
        nodes(:,1) = [0 0 0]';
        nodes(:,2) = [d2 0 0]';
        nodes(:,3) = nodes(:,2)+ Rmat(-beta,[0 0 1])*[d1 0 0]';
        nodes(:,4) =  Rmat(alpha,[1 0 0])* Rmat(-phi1,[0 1 0])*(H*csc(phi1))*[1 0 0]';
        nodes(:,5) = nodes(:,2) + Rmat(alpha,[1 0 0])* Rmat(-phi2,[0 1 0])*(H*csc(phi2))*[1 0 0]';
        nodes(:,6) = nodes(:,5) + Rmat(-beta,[0 0 1])*Rmat(pi-alpha,[1 0 0])*[d2 0 0]';
        nodes(:,7) = nodes(:,4) -2 * [0 0 nodes(3,4)]';
        nodes(:,8) = nodes(:,5) -2 * [0 0 nodes(3,5)]';
        nodes(:,9) = nodes(:,6) -2 * [0 0 nodes(3,6)]'; 
        h=range(nodes(3,:))/2;
        unitcell(1).vertices = nodes';
        unitcell(1).faces = [1 2 5 4; 2 3 6 5; 2 8 9 3; 1 7 8 2];
        
        % 'unitcell.faces' is a reference for the first cell whereas
        % 'unitcell2.faces' is the actual node structure for the whole
        % bellows
        unitcell2(1).faces = [1 2 5 4; 2 3 6 5; 2 8 9 3; 1 7 8 2];
        
        bottomNodes = unitcell(1).vertices(7:9,:);
        % next unit cells in a layer: rotate and attach
        Node = unitcell(1).vertices;
        Panel = num2cell(unitcell2(1).faces,2);

        
        
        vertices.coords =  unitcell(1).vertices;
        faces.vertices =unitcell(1).faces;

%     figure()
%     PlotOri(Node,Panel,[],'PanelColor','g'); 
%     axis off;
%     light
% 
% % 
%       figure()
%       PlotOri(Node,Panel,[],'ShowNumber','on');
% %         
 
        % first unit cell
        unitcell(1).vertices = nodes';
        unitcell(1).faces = [1 2 5 4; 2 3 6 5; 2 8 9 3; 1 7 8 2];
        
        % 'unitcell.faces' is a reference for the first cell whereas
        % 'unitcell2.faces' is the actual node structure for the whole
        % bellows
        unitcell2(1).faces = [1 2 5 4; 2 3 6 5; 2 8 9 3; 1 7 8 2];
        
        bottomNodes = unitcell(1).vertices(7:9,:);
        % next unit cells in a layer: rotate and attach
        Node = unitcell(1).vertices;
        Panel = num2cell(unitcell2(1).faces,2);

        
        
%         vertices.coords =  unitcell(1).vertices;
%         faces.vertices =unitcell(1).faces;

%     figure()
%     PlotOri(Node,Panel,[],'PanelColor','g'); 
%     axis off;
%     light
% 
% 
%      figure()
%      PlotOri(Node,Panel,[],'ShowNumber','on');
      
      
      for ii = 2:n
            
            % Define the vertex locations (essentially a translation round
            % the bellows)
            prevUnitCellNodes = unitcell(ii-1).vertices';
          
           
            unitcell(ii).vertices = (repmat(prevUnitCellNodes(:,3),1,9) + Rmat((ii-1)*(2*pi/n),[0 0 1])*nodes)';
            
             
            % 'unitcell' is upadted so the structre can be plotted for
            % checking
            unitcell(ii).faces = [1 2 4 5; 2 3 6 5; 2 8 9 3; 1 7 8 2];
          
          
            % Run through all 4 faces of all 4 facets. If one of the sides
            % is on the previous cell then make sure the node numbers are
            % consistent. If not then increment for this cell
            for jj = 1:4
                for kk = 1:4
                    if unitcell(1).faces(jj,kk) == 1
                        unitcell2(ii).faces(jj,kk) = unitcell2(ii-1).faces(2,2);
                    elseif unitcell(1).faces(jj,kk) == 2
                        unitcell2(ii).faces(jj,kk) = unitcell2(ii-1).faces(1,2) + 8 - 2*heaviside(ii-2.5);
                    elseif unitcell(1).faces(jj,kk) == 3
                        if ii == n
                            unitcell2(ii).faces(jj,kk) = 1;
                        else
                            unitcell2(ii).faces(jj,kk) = unitcell2(ii-1).faces(1,2) + 9 - 2*heaviside(ii-2.5);
                        end
                        
                    elseif unitcell(1).faces(jj,kk) == 4
                        unitcell2(ii).faces(jj,kk) = unitcell2(ii-1).faces(2,3);
                    elseif unitcell(1).faces(jj,kk) == 5
                        unitcell2(ii).faces(jj,kk) = unitcell2(ii-1).faces(1,3) + 7 - heaviside(ii-2.5) - heaviside(ii-(n-0.5));
                    elseif unitcell(1).faces(jj,kk) == 6
                        if ii == n
                            unitcell2(ii).faces(jj,kk) = 4;
                        else
                            unitcell2(ii).faces(jj,kk) = unitcell2(ii-1).faces(1,3) + 8 - heaviside(ii-2.5);
                        end
                        
                    elseif unitcell(1).faces(jj,kk) == 7
                        unitcell2(ii).faces(jj,kk) = unitcell2(ii-1).faces(3,3);
                    elseif unitcell(1).faces(jj,kk) == 8
                        unitcell2(ii).faces(jj,kk) = unitcell2(ii-1).faces(3,2) + 6 - 2*heaviside(ii-(n-0.5));
                    elseif unitcell(1).faces(jj,kk) == 9
                        if ii == n
                            unitcell2(ii).faces(jj,kk) = 7;
                        else
                            unitcell2(ii).faces(jj,kk) = unitcell2(ii-1).faces(3,2) + 7;
                        end
                        
                    end
                end
            end
            
            bottomNodes = [bottomNodes; unitcell(ii).vertices(7:9,:)];
            
      end
   
        vec = bottomNodes - mean(bottomNodes);
        Rs = max(sqrt(vec(:,1).^2+vec(:,2).^2+vec(:,3).^2));


end