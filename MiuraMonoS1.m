function [Node, Panel, V1, Vs, Hs, Rs] = MiuraMonoS1(phi1,phi2,r,n,m,d1overr,Output)
%  n=4
%  m=5
%  phi1=75/180*pi
%  phi2=32/180*pi
%  d1overr=1
%  r=40
%  Output=1
%        
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
         beta = beta1;
          
           alpha = asin(tan(beta/2)/tan(phi2));
            %%%alpha will sometimes be imaginary at flatfoldable state
          %alpha = real(alpha)
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
        Hs=h;
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
      
      
%    
%         vec = bottomNodes - mean(bottomNodes);
%         R = max(sqrt(vec(:,1).^2+vec(:,2).^2+vec(:,3).^2));
%    
%Calculate V1 and Vs
% Write the final vertices and face node numbers in the format
% required by the rest of the code
for ii = 1:length(unitcell)
    vertices.coords((ii-1)*9+1:ii*9,:) = unitcell(ii).vertices;
    faces.vertices((ii-1)*4+1:ii*4,:) = unitcell2(ii).faces;
end
%         remove any duplicated vertices and maintain the same order
vertices.coords = remove_extra_vertices(vertices.coords);
Node1 = vertices.coords;
% Panel1 = num2cell(faces.vertices,2);
% set(0,'defaultfigurecolor','w')
% figure()
% %PlotOri(Node,Panel,[],'PanelColor','w');
%  PlotOri1(Node1,Panel1,[],'PanelColor','[0.96, 0.96, 0.96]/256')
%  axis equal;
%  axis off;
% %light

% %
% figure()
% PlotOri1(Node1,Panel1,[],'ShowNumber','on');

BX=mean(Node1((round(Node1(:,3),2)==min(round(Node1(:,3),2))),1));
BY= mean(Node1((round(Node1(:,3),2)==min(round(Node1(:,3),2))),2));
BZ=min(Node1(:,3));
TZ=max(Node1(:,3));
CTbottom=[BX BY BZ];
CTtop=[BX BY TZ];
CTmiddle=[BX BY 0];
Node1=[Node1;CTbottom;CTtop;CTmiddle];
[V1]=calvolume(Node1(2,:),Node1(3,:),Node1(6*n+3,:),Node1(6*n+2,:));
[V2]=calvolume(Node1(2,:),Node1(3,:),Node1(6,:),Node1(6*n+2,:));
[V3]=calvolume(Node1(2,:),Node1(5,:),Node1(6,:),Node1(6*n+2,:));
[V4]=calvolume(Node1(4,:),Node1(5,:),Node1(6*n+2,:),Node1(1,:));
[V5]=calvolume(Node1(6*n+2,:),Node1(5,:),Node1(1,:),Node1(2,:));
[V6]=calvolume(Node1(1,:),Node1(2,:),Node1(6*n+3,:),Node1(6*n+2,:));
V1=(V1+V2+V3+V4+V5+V6)*m*n;
Rs=norm(Node1(1,:)-CTmiddle);
Vs=m*h*pi*Rs^2;
        
        
% % %         % Update for layer 2
        ii=1;
        nbr = n + ii;
        
      %  Translate to the next layer
        unitcell(nbr).vertices = unitcell(ii).vertices + 2*repmat([0 0 nodes(3,4)],9,1);
        unitcell(nbr).faces = [1 2 5 4; 2 3 6 5; 2 8 9 3; 1 7 8 2];
        
        
        unitcell2(nbr).faces = unitcell2(nbr-n).faces + 12+6*(n-2) - (ii-1)*3;
        
        
        % Run through every face on every facet to see if any share nodes
        % with the layer below. If they do then make sure the nodes are
        % consistent
        for ii2 = 1:4
            for jj2 = 1:4
                if unitcell(1).faces(ii2,jj2) == 7
                    unitcell2(nbr).faces(ii2,jj2) = unitcell2(nbr - n).faces(1,4);
                elseif unitcell(1).faces(ii2,jj2) == 8
                    unitcell2(nbr).faces(ii2,jj2) = unitcell2(nbr - n).faces(1,3);
                elseif unitcell(1).faces(ii2,jj2) == 9
                    unitcell2(nbr).faces(ii2,jj2) = unitcell2(nbr - n).faces(2,3);
                end
            end
        end
        
        %                         p(nbr) = patch(unitcell(nbr),'FaceColor',0.9*[1 1 1]');
        
        % Run around the second layer defining the nodes and ensuring any
        % shared nodes are consistently numbered as before
        for ii = 2:n
            nbr = n + ii;
            unitcell(nbr).vertices = unitcell(ii).vertices + 2*repmat([0 0 nodes(3,4)],9,1);
            unitcell(nbr).faces = [1 2 5 4; 2 3 6 5; 2 8 9 3; 1 7 8 2];
            for ii2 = 1:4
                for jj2 = 1:4
                    if unitcell(1).faces(ii2,jj2) == 7
                        unitcell2(nbr).faces(ii2,jj2) = unitcell2(nbr - n).faces(1,4);
                    elseif unitcell(1).faces(ii2,jj2) == 8
                        unitcell2(nbr).faces(ii2,jj2) = unitcell2(nbr - n).faces(1,3);
                    elseif unitcell(1).faces(ii2,jj2) == 9
                        unitcell2(nbr).faces(ii2,jj2) = unitcell2(nbr - n).faces(2,3);
                    end
                    
                    if unitcell(1).faces(ii2,jj2) == 1
                        unitcell2(nbr).faces(ii2,jj2) = unitcell2(nbr-1).faces(2,2);
                    elseif unitcell(1).faces(ii2,jj2) == 2
                        unitcell2(nbr).faces(ii2,jj2) = unitcell2(nbr-1).faces(1,2) + 5 - heaviside(ii-2.5);
                    elseif unitcell(1).faces(ii2,jj2) == 3
                        if ii == n
                            unitcell2(nbr).faces(ii2,jj2) = unitcell2(nbr-n+1).faces(1,1);
                        else
                            unitcell2(nbr).faces(ii2,jj2) = unitcell2(nbr-1).faces(1,2) + 6 - heaviside(ii-2.5);
                        end
                        
                    elseif unitcell(1).faces(ii2,jj2) == 4
                        unitcell2(nbr).faces(ii2,jj2) = unitcell2(nbr-1).faces(2,3);
                    elseif unitcell(1).faces(ii2,jj2) == 5
                        unitcell2(nbr).faces(ii2,jj2) = unitcell2(nbr-1).faces(1,3) + 4 - heaviside(ii-(n-0.5));
                    elseif unitcell(1).faces(ii2,jj2) == 6
                        if ii == n
                            unitcell2(nbr).faces(ii2,jj2) = unitcell2(nbr - n + 1).faces(1,4);
                        else
                            unitcell2(nbr).faces(ii2,jj2) = unitcell2(nbr-1).faces(1,3) + 5;% - heaviside(ii-2.5);
                        end
                    end
                end
            end
            %                             p(nbr) = patch(unitcell(nbr),'FaceColor',0.9*[1 1 1]');
        end
        
        % repeat for m layers this time the node numbers are 8+4*(n-2)
        % greater than before unless they are shared with beneath so
        % therefore should be consistent
        for jj = 3:m
            for ii = 1:n
                nbr = (jj-1)*n + ii;
                unitcell(nbr).vertices = unitcell(ii).vertices + (jj-1)*2*repmat([0 0 nodes(3,4)],9,1);
                unitcell(nbr).faces = [1 2 5 4; 2 3 6 5; 2 8 9 3; 1 7 8 2];
                unitcell2(nbr).faces = unitcell2(nbr-n).faces + 8+4*(n-2);
                for ii2 = 1:4
                    for jj2 = 1:4
                        if unitcell(1).faces(ii2,jj2) == 7
                            unitcell2(nbr).faces(ii2,jj2) = unitcell2(nbr - n).faces(1,4);
                        elseif unitcell(1).faces(ii2,jj2) == 8
                            unitcell2(nbr).faces(ii2,jj2) = unitcell2(nbr - n).faces(1,3);
                        elseif unitcell(1).faces(ii2,jj2) == 9
                            unitcell2(nbr).faces(ii2,jj2) = unitcell2(nbr - n).faces(2,3);
                        end
                    end
                end
                %                                 p(nbr) = patch(unitcell(nbr),'FaceColor',0.9*[1 1 1]');
            end
       end
% %         
% %         
        % Write the final vertices and face node numbers in the format
        % required by the rest of the code
        for ii = 1:length(unitcell)
            vertices.coords((ii-1)*9+1:ii*9,:) = unitcell(ii).vertices;            
            faces.vertices((ii-1)*4+1:ii*4,:) = unitcell2(ii).faces;
        end
        
       % remove any duplicated vertices and maintain the same order
       vertices.coords = remove_extra_vertices(vertices.coords);     
       Node = vertices.coords;
       Panel = num2cell(faces.vertices,2);  % This is the Node and Panel of Miura unit containing 2 stories
       [Node,index]=sortrows(Node,3);                       
       Panel=cell2mat(Panel);
       Panelstore=zeros(length(Panel),4);
       for j = 1:length(index)
           Panelstore=Panelstore+j.*(Panel==index(j,1));
       end
       Panel=Panelstore;      
       deletenodes=find(Node(:,3)>(m-1)*h);
       Node(any(Node(:,3)>(m-1)*h,3),:)=[];          %Miura unit containing 1 story
          end
          
          for j = 1:length(deletenodes)
              Panel(any(Panel == deletenodes(j,1),2),:)=[];
          end
          Panel=num2cell(Panel,2);
          V=CTbottom-[0 0 0];
          Node=Node-V.*ones(length(Node),3);
          vertices.coords = Node;            
          faces.vertices= cell2mat(Panel);
          
%           figure()
% 
%           PlotOri1(Node,Panel,[],'PanelColor','g');
%           axis on;
%           light
%           
%           figure()
%           PlotOri1(Node,Panel,[],'ShowNumber','on');
%           hold on


%%%Output.Fold
%%
if Output==1
% .FOLD Version
file.spec = 2.0;

% Made in MATLAB - Keep the quotes
file.creator = '"MATLAB Code"';

% Author - Keep the quotes
file.author = '"Mengzhu Yang"';

% Class:
% - "singleModel": A single origami model
file.classes = '"singleModel"';

frame.title = '"Miura_S1"';
% This code is geared towards only doing the folded form here
frame.classes = '"foldedForm"';

% Again this code is only aimed towards 3D here
frame.attributes = '"3D"';

% NOTE: all units are in mm throughout
frame.unit = '"mm"';

FOLD_count = 0;

%Write Output

filename = sprintf('Miura_S1.FOLD');
fid = fopen(filename,'w');
fprintf(fid,'{\n');
FOLD_count = FOLD_count + 1;
filenameList{FOLD_count} = filename;
%RStore(FOLD_count) = R;
%save('BellowsLandscapeValidOptions.mat','filenameList','RStore');

file_names = fieldnames(file);
for ii = 1:length(file_names)
    command_str = sprintf('var_contents = num2str(file.%s);',file_names{ii});
    eval(command_str)
    fprintf(fid,'  "file_%s": %s,\n',file_names{ii},var_contents);
end

% Write the contents of the 'frame' structure
frame_names = fieldnames(frame);
for ii = 1:length(frame_names)
    command_str = sprintf('var_contents = num2str(frame.%s);',frame_names{ii});
    eval(command_str)
    fprintf(fid,'  "frame_%s": %s,\n',frame_names{ii},var_contents);
end

% Write the vertex locations
fprintf(fid,'  "vertices_coords": [\n');
for ii = 1:size(vertices.coords,1)-1
    fprintf(fid,'    [%g,%g,%g],\n',vertices.coords(ii,1),vertices.coords(ii,2),vertices.coords(ii,3));
end
fprintf(fid,'    [%g,%g,%g]\n',vertices.coords(ii+1,1),vertices.coords(ii+1,2),vertices.coords(ii+1,3));
fprintf(fid,'  ],\n');

% Write the faces
fprintf(fid,'  "faces_vertices": [\n');
for ii = 1:size(faces.vertices,1)-1
    fprintf(fid,'    [%g,%g,%g,%g],\n',faces.vertices(ii,1),faces.vertices(ii,2),faces.vertices(ii,3),faces.vertices(ii,4));
end
fprintf(fid,'    [%g,%g,%g,%g]\n',faces.vertices(ii+1,1),faces.vertices(ii+1,2),faces.vertices(ii+1,3),faces.vertices(ii+1,4));
fprintf(fid,'  ],\n');


fprintf(fid,'}');
[~] = fclose(fid);
end

end