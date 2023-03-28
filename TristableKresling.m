%input n,m,H,phi1,phi2
%clear all; close all; clc;


%%  Input geometric parameters 
% Flat state, according to Reid
% https://arxiv.org/pdf/1609.01354.pdf

 n=4;
 m=3;
 H = 45;        
%phi1 and phi2 at state 1&2
 phi1Deg=77;
 phi2Deg=35;
%corresponding phi1 and phi2 at state 3
%         phi1Deg=134.0214;
%         phi2Deg=77;
% 
% %convert angles
phi1 = phi1Deg*pi/180;
phi2 = phi2Deg*pi/180;



%%
% .FOLD Version
file.spec = 2.0;

% Made in MATLAB - Keep the quotes
file.creator = '"MATLAB Code"';

% Author - Keep the quotes
file.author = '"Mengzhu Yang"';

% Class:
% - "singleModel": A single origami model
file.classes = '"singleModel"';

frame.title = '"Tristable_Kresling"';
% This code is geared towards only doing the folded form here
frame.classes = '"foldedForm"';

% Again this code is only aimed towards 3D here
frame.attributes = '"3D"';

% NOTE: all units are in mm throughout
frame.unit = '"mm"';

FOLD_count = 0;
noWrite = 0;


%Constraints on phi1 and phi2     
      if phi1Deg <= 180/n   
          disp('0 solution');
            noWrite = 1;
            solution = 0;
      end
      
      if phi1Deg <= phi2Deg  
          disp('0 solution');
            noWrite = 1;
             solution = 0;
      end
      
      if phi1Deg>180/4+180/(2*n) 
        % Boundary of bistable region         
           limitangle1=(phi1Deg-180/n);
           limitangle2=(atan(tan(phi1)*(1-sin(pi/n))/(1+sin(pi/n))))/pi*180; 
         if phi2Deg>phi1Deg-180/n && phi2Deg<(atan(tan(phi1)*(1-sin(pi/n))/(1+sin(pi/n))))/pi*180               
            disp('2 solutions');
            solution = 2;             
            %1st solution (S2)： beta = beta1
            %2nd solution (S1)： beta = beta2
         elseif phi2>0 && phi2<phi1-pi/n
            disp('1 solution');
            solution = 1; 
         elseif  phi2>atan(tan(phi1)*(1-sin(pi/n))/(1+sin(pi/n))) 
            disp('0 solution');
            solution = 0;
         else
             solution = 0;
         end
      end

      
      %%%%preclude flat-foldable
      if phi1Deg == phi2Deg +180/n 
          disp('Flatfoldable');
            noWrite = 1;
            solution = 0;          
      end
      
      %% 0 solution 
      if  solution == 0
        error('No solution')  
      end
      
      
      %% 1 solution 
      if  solution == 1   %S1 && S3
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
      %S1
      x1 = tan(phi2)/tan(phi1);
      x2= tan(pi/n);
      beta1 = 2*atan(1/(2*x2)*(1-x1-(sqrt((x1-1)^2 - 4*x1*x2^2))));
      beta2 = 2*atan(1/(2*x2)*(1-x1+(sqrt((x1-1)^2 - 4*x1*x2^2))));
      beta = beta1;        
            % alpha is angular deflection from vertical (z-axis) 
            %pi-alpha is the opening angle between adjacent bands of folds
          alpha = asin(tan(beta/2)/tan(phi2));
          alphadeg=90-alpha/pi*180;
        %constraints on alpha beta
          tf = isreal(beta) ;
        if tf==0
            noWrite = 1;     
        end
        
        tf2 = isreal(alpha) ;
        if tf2==0
            noWrite = 1;     
        end
        
         %Calculate the length d
          d=H/(tan(phi2))-H/(tan(phi1));  

        
        %An unit cell in 3D   
        nodes(:,1) = [0 0 0]';      
        nodes(:,2) =Rmat(-beta,[0 0 1])*[d 0 0]';  
        nodes(:,3) = Rmat(alpha,[1 0 0])* Rmat(-phi1,[0 1 0])*(H*csc(phi1))*[1 0 0]';
        nodes(:,4) = Rmat(alpha,[1 0 0])* Rmat(-phi2,[0 1 0])*(H*csc(phi2))*[1 0 0]'; 
        %height after folding
        h=range(nodes(3,:)); % a layer
        % first unit cell
        unitcell(1).vertices = nodes';       
        unitcell(1).faces = [1 2 4;4 1 3];
        
     
        % 'unitcell.faces' is a reference for the first cell whereas
        % 'unitcell2.faces' is the actual node structure for the whole
        % bellows        
        unitcell2(1).faces = [1 2 4;4 1 3];
        bottomNodes = unitcell(1).vertices(1:2,:);
        Node = unitcell(1).vertices;
        Panel = num2cell(unitcell2(1).faces,2);

        % next unit cells in a layer: rotate and attach
        for ii = 2:n           
            % Define the vertex locations (essentially a translation round
            % the bellows)
            prevUnitCellNodes = unitcell(ii-1).vertices';   
            unitcell(ii).vertices = (repmat(prevUnitCellNodes(:,2),1,4) +Rmat((ii-1)*(2*pi/n),[0 0 1])*nodes)'; 
            % 'unitcell' is updated 
            unitcell(ii).faces = [1 2 4;4 1 3];
         
            Node = unitcell(ii).vertices;
            Panel = num2cell(unitcell(ii).faces,2);
%     PlotOri(Node,Panel,[],'PanelColor','g'); 
%     axis off;
%     light
%     hold on
     end 
        
 % Covert into the calculated model 
 %https://journals.aps.org/pre/abstract/10.1103/PhysRevE.100.063001
        R = d/(2*sind(180/n));   %radius
        layer = m+1;  %number of layer
        a=R;b=R;
        c=sqrt((nodes(1,1)-nodes(1,3))^2+(nodes(2,1)-nodes(2,3))^2+(nodes(3,3)-h-nodes(3,1)).^2); 
        if phi1>pi/2    
        c=sqrt((nodes(1,1)-nodes(1,4))^2+(nodes(2,1)-nodes(2,4))^2+(nodes(3,4)-h-nodes(3,1)).^2);
        theta=acos((a^2+b^2-c^2)/(2*a*b)); % rotation angle between layer-1 and layer       
        thetadegree = theta/pi*180;
        else 
        theta=acos((a^2+b^2-c^2)/(2*a*b)); % rotation angle between layer-1 and layer 
        thetadegree = theta/pi*180;
        end
        
        
%%%calculation
angle = zeros(layer,n);
zori = zeros(layer,n);

  for i = 1:n 
   for t = 1:layer
        if mod(t,2) == 0
        angle(t,:) = [theta: 2*pi/n : theta+2*pi/n*(n-1)];
        else
        angle(t,:) = [0: 2*pi/n : 2*pi/n*(n-1)];
        end
   zori(t,:) = [(t-1)*h];
   end 
end   

angleres = reshape(angle',layer*n,1);
%%%coordinate of each point
x = R*cos(angleres);
y = R*sin(angleres);
z = reshape(zori',layer*n,1);
Node = [x,y,z];
%Rotate for an angle to twist
 Node = Rmat(-theta,[0 0 1])*Node'; %%%S2 n4
 Node = Node';
vertices.coords = Node;

Nop = zeros(2*n*(layer-1),3);  %panel:connection of  nodes
if phi1>pi/2
for t = 1:layer-1
for i = 1:n
  if mod(t,2)~= 0
   if i == 1
     Nop(2*(t-1)*n+i,:) = [i+(t-1)*n;i+n+(t-1)*n; n+t*n;];
    Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+1+(t-1)*n; i+n+(t-1)*n];
  elseif mod(i,n)~= 0
    Nop(2*(t-1)*n+i,:) = [i+(t-1)*n;i+n+(t-1)*n-1; i+n+(t-1)*n;];
    Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+1+(t-1)*n; i+n+(t-1)*n];
    else
    Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+n+(t-1)*n; i+n+(t-1)*n-1];
    Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+(t-1)*n-(n-1); i+n+(t-1)*n];
    end  
  else
    if i == 1
      Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+(n-1)+(t-1)*n; i+(2*n-1)+(t-1)*n];
      Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+n+(t-1)*n; i+(2*n-1)+(t-1)*n];
    else
     Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+(n-1)+(t-1)*n; i+n+(t-1)*n];
     Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i-1+(t-1)*n; i+(n-1)+(t-1)*n];
    end 
  end
end
end   
else
for t = 1:layer-1
for i = 1:n
  if mod(t,2)~= 0
    if mod(i,n)~= 0
    Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+n+(t-1)*n; i+n+1+(t-1)*n];
    Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+1+(t-1)*n; i+1+n+(t-1)*n];
    else
    Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+n+(t-1)*n; i+1+(t-1)*n];
    Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+(t-1)*n-(n-1); i+1+(t-1)*n];
    end  
  else
    if i == 1
      Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+(n-1)+(t-1)*n; i+(2*n-1)+(t-1)*n];
      Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+n+(t-1)*n; i+(2*n-1)+(t-1)*n];
    else
     Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+(n-1)+(t-1)*n; i+n+(t-1)*n];
     Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i-1+(t-1)*n; i+(n-1)+(t-1)*n];
    end 
  end
end
end
end


faces.vertices=Nop;      
Panel=num2cell(Nop,2);
Node = vertices.coords;
Panel = num2cell(faces.vertices,2);


set(0,'defaultfigurecolor','w')
figure(1)
PlotOri(Node,Panel,[],'PanelColor','[0.67, 0.83, 0.9]/256'); 
title('S1')
axis off;
axis equal
light
   

%Write Output

filename = sprintf('Kresling_S1.FOLD');
fid = fopen(filename,'w');
fprintf(fid,'{\n');
FOLD_count = FOLD_count + 1;
filenameList{FOLD_count} = filename;
RStore(FOLD_count) = R;
save('BellowsLandscapeValidOptions.mat','filenameList','RStore')

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
    fprintf(fid,'    [%g,%g,%g],\n',faces.vertices(ii,1),faces.vertices(ii,2),faces.vertices(ii,3));
end
fprintf(fid,'    [%g,%g,%g]\n',faces.vertices(ii+1,1),faces.vertices(ii+1,2),faces.vertices(ii+1,3));
fprintf(fid,'  ],\n');

fprintf(fid,'}');
[~] = fclose(fid);
        
        
      
        
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
 %S3
  diagnalphi1=180+atan((tan(phi1)*tan(phi2))/(2*tan(phi2)-tan(phi1)))/pi*180;
  phi2=phi1;
  phi1=diagnalphi1/180*pi;

 
        unitcell = [];
        unitcell2 = [];
        faces = [];
        vertices = [];
        
   
        x1 = tan(phi2)/tan(phi1);
          x2= tan(pi/n);
          beta1 = 2*atan(1/(2*x2)*(1-x1-(sqrt((x1-1)^2 - 4*x1*x2^2))));
          beta2 = 2*atan(1/(2*x2)*(1-x1+(sqrt((x1-1)^2 - 4*x1*x2^2))));
          beta = beta1;    
          alpha = asin(tan(beta/2)/tan(phi2));
          alphadeg=90-alpha/pi*180;
        nodes(:,1) = [0 0 0]';      
        nodes(:,2) =Rmat(-beta,[0 0 1])*[d 0 0]';  
        nodes(:,3) = Rmat(alpha,[1 0 0])* Rmat(-phi1,[0 1 0])*(H*csc(phi1))*[1 0 0]';
        nodes(:,4) = Rmat(alpha,[1 0 0])* Rmat(-phi2,[0 1 0])*(H*csc(phi2))*[1 0 0]';
        %nodes(:,5) =  Rmat(beta2,[0 0 1])*[d 0 0]'; 
        %height after folding
        hdeploy=range(nodes(3,:));    
         % first unit cell
        unitcell(1).vertices = nodes';       
        unitcell(1).faces = [1 2 4;4 1 3];
        % 'unitcell.faces' is a reference for the first cell whereas
        % 'unitcell2.faces' is the actual node structure for the whole
        % bellows        
        unitcell2(1).faces = [1 2 4;4 1 3];
        bottomNodes = unitcell(1).vertices(1:2,:);
        % next unit cells in a layer: rotate and attach
        Node = unitcell(1).vertices;
        Panel = num2cell(unitcell2(1).faces,2);
   
 % Covert into the calculated model 
 %https://journals.aps.org/pre/abstract/10.1103/PhysRevE.100.063001
        layer = m+1;  %number of layer
        a=R;b=R;
        c=sqrt((nodes(1,1)-nodes(1,4))^2+(nodes(2,1)-nodes(2,4))^2+(nodes(3,4)-hdeploy-nodes(3,1)).^2);
        theta=acos((a^2+b^2-c^2)/(2*a*b)); % rotation angle between layer-1 and layer            
       
%%%calculation
angle = zeros(layer,n);
zori = zeros(layer,n);


  for i = 1:n 
   for t = 1:layer
        if mod(t,2) == 0
        angle(t,:) = [theta: 2*pi/n : theta+2*pi/n*(n-1)]; 
        %angle(t,:) = [(-1)*(t-1)*theta: 2*pi/n : (-1)*(t-1)*theta+2*pi/n*(n-1)]; 
        else
        angle(t,:) = [0: 2*pi/n : 2*pi/n*(n-1)];
       % angle(t,:) = [(-1)*(t-1)*theta: 2*pi/n : (-1)*(t-1)*theta+2*pi/n*(n-1)]; 
        end
   zori(t,:) = [(t-1)*hdeploy];
   end 
end   

angleres = reshape(angle',layer*n,1);
%%%coordinate of each point
x = R*cos(angleres);
y = R*sin(angleres);
z = reshape(zori',layer*n,1);
Node = [x,y,z];
%Rotate for an angle to twist
 Node = Rmat(-theta,[0 0 1])*Node'; %%%S2 n4
 Node = Node';
vertices.coords = Node;

Nop = zeros(2*n*(layer-1),3);  %panel:connection of  nodes
for t = 1:layer-1
for i = 1:n
  if mod(t,2)~= 0
   if i == 1
     Nop(2*(t-1)*n+i,:) = [i+(t-1)*n;i+n+(t-1)*n; n+t*n;];
    Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+1+(t-1)*n; i+n+(t-1)*n];
    elseif mod(i,n)~= 0
    Nop(2*(t-1)*n+i,:) = [i+(t-1)*n;i+n+(t-1)*n-1; i+n+(t-1)*n;];
    Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+1+(t-1)*n; i+n+(t-1)*n];
    else
    Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+n+(t-1)*n; i+n+(t-1)*n-1];
    Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+(t-1)*n-(n-1); i+n+(t-1)*n];
    end  
  else
    if i == 1
    Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+(n-1)+(t-1)*n; i+n+(t-1)*n];
    Nop(2*(t-1)*n+i+n,:)= [i+(n-1)+(t-1)*n; i+n+(t-1)*n; i+(2*n-1)+(t-1)*n];
    else
     Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i-1+(t-1)*n; i+n+(t-1)*n];
     Nop(2*(t-1)*n+i+n,:)= [i+n+(t-1)*n; i-1+(t-1)*n; i+(n-1)+(t-1)*n];
    end 
  end
end
end   

faces.vertices=Nop;      
Panel=num2cell(Nop,2);

%Plot to check
Node = vertices.coords;
Panel = num2cell(faces.vertices,2);

% 
figure(2)
PlotOri(Node,Panel,[],'PanelColor','[0.67, 0.83, 0.9]/256'); 
title('S3')
axis off;
axis equal
light  
 
 
 %Write Output

filename = sprintf('Kresling_S3.FOLD');
fid = fopen(filename,'w');
fprintf(fid,'{\n');
FOLD_count = FOLD_count + 1;
filenameList{FOLD_count} = filename;
RStore(FOLD_count) = R;
save('BellowsLandscapeValidOptions.mat','filenameList','RStore')

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
    fprintf(fid,'    [%g,%g,%g],\n',faces.vertices(ii,1),faces.vertices(ii,2),faces.vertices(ii,3));
end
fprintf(fid,'    [%g,%g,%g]\n',faces.vertices(ii+1,1),faces.vertices(ii+1,2),faces.vertices(ii+1,3));
fprintf(fid,'  ],\n');

fprintf(fid,'}');
[~] = fclose(fid);

    
end   
     

      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
 %% 2 solutions 
 %S1 S2 & S3
      if  solution == 2   
          
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
      %S1
      x1 = tan(phi2)/tan(phi1);
      x2= tan(pi/n);
      beta1 = 2*atan(1/(2*x2)*(1-x1-(sqrt((x1-1)^2 - 4*x1*x2^2))));
      beta2 = 2*atan(1/(2*x2)*(1-x1+(sqrt((x1-1)^2 - 4*x1*x2^2))));
      beta = beta2;        
            % alpha is angular deflection from vertical (z-axis) 
            %pi-alpha is the opening angle between adjacent bands of folds
          alpha = asin(tan(beta/2)/tan(phi2));
          alphadeg=90-alpha/pi*180;
        %constraints on alpha beta
          tf = isreal(beta) ;
        if tf==0
            noWrite = 1;     
        end
        
        tf2 = isreal(alpha) ;
        if tf2==0
            noWrite = 1;     
        end
        
         %Calculate the length d
          d=H/(tan(phi2))-H/(tan(phi1));  

        
        %An unit cell in 3D   
        nodes(:,1) = [0 0 0]';      
        nodes(:,2) =Rmat(-beta,[0 0 1])*[d 0 0]';  
        nodes(:,3) = Rmat(alpha,[1 0 0])* Rmat(-phi1,[0 1 0])*(H*csc(phi1))*[1 0 0]';
        nodes(:,4) = Rmat(alpha,[1 0 0])* Rmat(-phi2,[0 1 0])*(H*csc(phi2))*[1 0 0]'; 
        %height after folding
        h=range(nodes(3,:)); % a layer
        % first unit cell
        unitcell(1).vertices = nodes';       
        unitcell(1).faces = [1 2 4;4 1 3];
        
     
        % 'unitcell.faces' is a reference for the first cell whereas
        % 'unitcell2.faces' is the actual node structure for the whole
        % bellows        
        unitcell2(1).faces = [1 2 4;4 1 3];
        bottomNodes = unitcell(1).vertices(1:2,:);
        Node = unitcell(1).vertices;
        Panel = num2cell(unitcell2(1).faces,2);

        % next unit cells in a layer: rotate and attach
        for ii = 2:n           
            % Define the vertex locations (essentially a translation round
            % the bellows)
            prevUnitCellNodes = unitcell(ii-1).vertices';   
            unitcell(ii).vertices = (repmat(prevUnitCellNodes(:,2),1,4) +Rmat((ii-1)*(2*pi/n),[0 0 1])*nodes)'; 
            % 'unitcell' is updated 
            unitcell(ii).faces = [1 2 4;4 1 3];
         
            Node = unitcell(ii).vertices;
            Panel = num2cell(unitcell(ii).faces,2);
%     PlotOri(Node,Panel,[],'PanelColor','g'); 
%     axis off;
%     light
%     hold on
     end 
        
 % Covert into the calculated model 
 %https://journals.aps.org/pre/abstract/10.1103/PhysRevE.100.063001
        R = d/(2*sind(180/n));   %radius
        layer = m+1;  %number of layer
        a=R;b=R;
        c=sqrt((nodes(1,1)-nodes(1,3))^2+(nodes(2,1)-nodes(2,3))^2+(nodes(3,3)-h-nodes(3,1)).^2); 
        if phi1>pi/2    
        c=sqrt((nodes(1,1)-nodes(1,4))^2+(nodes(2,1)-nodes(2,4))^2+(nodes(3,4)-h-nodes(3,1)).^2);
        theta=acos((a^2+b^2-c^2)/(2*a*b)); % rotation angle between layer-1 and layer       
        thetadegree = theta/pi*180;
        else 
        theta=acos((a^2+b^2-c^2)/(2*a*b)); % rotation angle between layer-1 and layer 
        thetadegree = theta/pi*180;
        end
        
        
%%%calculation
angle = zeros(layer,n);
zori = zeros(layer,n);

  for i = 1:n 
   for t = 1:layer
        if mod(t,2) == 0
        angle(t,:) = [theta: 2*pi/n : theta+2*pi/n*(n-1)];
        else
        angle(t,:) = [0: 2*pi/n : 2*pi/n*(n-1)];
        end
   zori(t,:) = [(t-1)*h];
   end 
end   

angleres = reshape(angle',layer*n,1);
%%%coordinate of each point
x = R*cos(angleres);
y = R*sin(angleres);
z = reshape(zori',layer*n,1);
Node = [x,y,z];
%Rotate for an angle to twist
 Node = Rmat(-theta,[0 0 1])*Node'; %%%S2 n4
 Node = Node';
vertices.coords = Node;

Nop = zeros(2*n*(layer-1),3);  %panel:connection of  nodes
if phi1>pi/2
for t = 1:layer-1
for i = 1:n
  if mod(t,2)~= 0
   if i == 1
     Nop(2*(t-1)*n+i,:) = [i+(t-1)*n;i+n+(t-1)*n; n+t*n;];
    Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+1+(t-1)*n; i+n+(t-1)*n];
  elseif mod(i,n)~= 0
    Nop(2*(t-1)*n+i,:) = [i+(t-1)*n;i+n+(t-1)*n-1; i+n+(t-1)*n;];
    Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+1+(t-1)*n; i+n+(t-1)*n];
    else
    Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+n+(t-1)*n; i+n+(t-1)*n-1];
    Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+(t-1)*n-(n-1); i+n+(t-1)*n];
    end  
  else
    if i == 1
      Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+(n-1)+(t-1)*n; i+(2*n-1)+(t-1)*n];
      Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+n+(t-1)*n; i+(2*n-1)+(t-1)*n];
    else
     Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+(n-1)+(t-1)*n; i+n+(t-1)*n];
     Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i-1+(t-1)*n; i+(n-1)+(t-1)*n];
    end 
  end
end
end   
else
for t = 1:layer-1
for i = 1:n
  if mod(t,2)~= 0
    if mod(i,n)~= 0
    Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+n+(t-1)*n; i+n+1+(t-1)*n];
    Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+1+(t-1)*n; i+1+n+(t-1)*n];
    else
    Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+n+(t-1)*n; i+1+(t-1)*n];
    Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+(t-1)*n-(n-1); i+1+(t-1)*n];
    end  
  else
    if i == 1
      Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+(n-1)+(t-1)*n; i+(2*n-1)+(t-1)*n];
      Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+n+(t-1)*n; i+(2*n-1)+(t-1)*n];
    else
     Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+(n-1)+(t-1)*n; i+n+(t-1)*n];
     Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i-1+(t-1)*n; i+(n-1)+(t-1)*n];
    end 
  end
end
end
end


faces.vertices=Nop;      
Panel=num2cell(Nop,2);
Node = vertices.coords;
Panel = num2cell(faces.vertices,2);


set(0,'defaultfigurecolor','w')
figure(1)
PlotOri(Node,Panel,[],'PanelColor','[0.67, 0.83, 0.9]/256'); 
title('S1')
axis off;
axis equal
light
   

%Write Output

filename = sprintf('Kresling_S1.FOLD');
fid = fopen(filename,'w');
fprintf(fid,'{\n');
FOLD_count = FOLD_count + 1;
filenameList{FOLD_count} = filename;
RStore(FOLD_count) = R;
save('BellowsLandscapeValidOptions.mat','filenameList','RStore')

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
    fprintf(fid,'    [%g,%g,%g],\n',faces.vertices(ii,1),faces.vertices(ii,2),faces.vertices(ii,3));
end
fprintf(fid,'    [%g,%g,%g]\n',faces.vertices(ii+1,1),faces.vertices(ii+1,2),faces.vertices(ii+1,3));
fprintf(fid,'  ],\n');

fprintf(fid,'}');
[~] = fclose(fid);
 

 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %S2     
      x1 = tan(phi2)/tan(phi1);
      x2= tan(pi/n);
      beta1 = 2*atan(1/(2*x2)*(1-x1-(sqrt((x1-1)^2 - 4*x1*x2^2))));
      beta2 = 2*atan(1/(2*x2)*(1-x1+(sqrt((x1-1)^2 - 4*x1*x2^2))));
      beta = beta1;        
            % alpha is angular deflection from vertical (z-axis) 
            %pi-alpha is the opening angle between adjacent bands of folds
          alpha = asin(tan(beta/2)/tan(phi2));
          alphadeg=90-alpha/pi*180;
        %constraints on alpha beta
          tf = isreal(beta) ;
        if tf==0
            noWrite = 1;     
        end
        
        tf2 = isreal(alpha) ;
        if tf2==0
            noWrite = 1;     
        end
        
         %Calculate the length d
          d=H/(tan(phi2))-H/(tan(phi1));  

        
        %An unit cell in 3D   
        nodes(:,1) = [0 0 0]';      
        nodes(:,2) =Rmat(-beta,[0 0 1])*[d 0 0]';  
        nodes(:,3) = Rmat(alpha,[1 0 0])* Rmat(-phi1,[0 1 0])*(H*csc(phi1))*[1 0 0]';
        nodes(:,4) = Rmat(alpha,[1 0 0])* Rmat(-phi2,[0 1 0])*(H*csc(phi2))*[1 0 0]'; 
        %height after folding
        h=range(nodes(3,:)); % a layer
        % first unit cell
        unitcell(1).vertices = nodes';       
        unitcell(1).faces = [1 2 4;4 1 3];
        
     
        % 'unitcell.faces' is a reference for the first cell whereas
        % 'unitcell2.faces' is the actual node structure for the whole
        % bellows        
        unitcell2(1).faces = [1 2 4;4 1 3];
        bottomNodes = unitcell(1).vertices(1:2,:);
        Node = unitcell(1).vertices;
        Panel = num2cell(unitcell2(1).faces,2);

        % next unit cells in a layer: rotate and attach
        for ii = 2:n           
            % Define the vertex locations (essentially a translation round
            % the bellows)
            prevUnitCellNodes = unitcell(ii-1).vertices';   
            unitcell(ii).vertices = (repmat(prevUnitCellNodes(:,2),1,4) +Rmat((ii-1)*(2*pi/n),[0 0 1])*nodes)'; 
            % 'unitcell' is updated 
            unitcell(ii).faces = [1 2 4;4 1 3];
         
            Node = unitcell(ii).vertices;
            Panel = num2cell(unitcell(ii).faces,2);
%     PlotOri(Node,Panel,[],'PanelColor','g'); 
%     axis off;
%     light
%     hold on
     end 
        
 % Covert into the calculated model 
 %https://journals.aps.org/pre/abstract/10.1103/PhysRevE.100.063001
        R = d/(2*sind(180/n));   %radius
        layer = m+1;  %number of layer
        a=R;b=R;
        c=sqrt((nodes(1,1)-nodes(1,3))^2+(nodes(2,1)-nodes(2,3))^2+(nodes(3,3)-h-nodes(3,1)).^2); 
        if phi1>pi/2    
        c=sqrt((nodes(1,1)-nodes(1,4))^2+(nodes(2,1)-nodes(2,4))^2+(nodes(3,4)-h-nodes(3,1)).^2);
        theta=acos((a^2+b^2-c^2)/(2*a*b)); % rotation angle between layer-1 and layer       
        thetadegree = theta/pi*180;
        else 
        theta=acos((a^2+b^2-c^2)/(2*a*b)); % rotation angle between layer-1 and layer 
        thetadegree = theta/pi*180;
        end
        
        
%%%calculation
angle = zeros(layer,n);
zori = zeros(layer,n);

  for i = 1:n 
   for t = 1:layer
        if mod(t,2) == 0
        angle(t,:) = [theta: 2*pi/n : theta+2*pi/n*(n-1)];
        else
        angle(t,:) = [0: 2*pi/n : 2*pi/n*(n-1)];
        end
   zori(t,:) = [(t-1)*h];
   end 
end   

angleres = reshape(angle',layer*n,1);
%%%coordinate of each point
x = R*cos(angleres);
y = R*sin(angleres);
z = reshape(zori',layer*n,1);
Node = [x,y,z];
%Rotate for an angle to twist
 Node = Rmat(-theta,[0 0 1])*Node'; %%%S2 n4
 Node = Node';
vertices.coords = Node;

Nop = zeros(2*n*(layer-1),3);  %panel:connection of  nodes
if phi1>pi/2
for t = 1:layer-1
for i = 1:n
  if mod(t,2)~= 0
   if i == 1
     Nop(2*(t-1)*n+i,:) = [i+(t-1)*n;i+n+(t-1)*n; n+t*n;];
    Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+1+(t-1)*n; i+n+(t-1)*n];
  elseif mod(i,n)~= 0
    Nop(2*(t-1)*n+i,:) = [i+(t-1)*n;i+n+(t-1)*n-1; i+n+(t-1)*n;];
    Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+1+(t-1)*n; i+n+(t-1)*n];
    else
    Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+n+(t-1)*n; i+n+(t-1)*n-1];
    Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+(t-1)*n-(n-1); i+n+(t-1)*n];
    end  
  else
    if i == 1
      Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+(n-1)+(t-1)*n; i+(2*n-1)+(t-1)*n];
      Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+n+(t-1)*n; i+(2*n-1)+(t-1)*n];
    else
     Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+(n-1)+(t-1)*n; i+n+(t-1)*n];
     Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i-1+(t-1)*n; i+(n-1)+(t-1)*n];
    end 
  end
end
end   
else
for t = 1:layer-1
for i = 1:n
  if mod(t,2)~= 0
    if mod(i,n)~= 0
    Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+n+(t-1)*n; i+n+1+(t-1)*n];
    Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+1+(t-1)*n; i+1+n+(t-1)*n];
    else
    Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+n+(t-1)*n; i+1+(t-1)*n];
    Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+(t-1)*n-(n-1); i+1+(t-1)*n];
    end  
  else
    if i == 1
      Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+(n-1)+(t-1)*n; i+(2*n-1)+(t-1)*n];
      Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+n+(t-1)*n; i+(2*n-1)+(t-1)*n];
    else
     Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+(n-1)+(t-1)*n; i+n+(t-1)*n];
     Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i-1+(t-1)*n; i+(n-1)+(t-1)*n];
    end 
  end
end
end
end


faces.vertices=Nop;      
Panel=num2cell(Nop,2);
Node = vertices.coords;
Panel = num2cell(faces.vertices,2);


set(0,'defaultfigurecolor','w')
figure(2)
PlotOri(Node,Panel,[],'PanelColor','[0.67, 0.83, 0.9]/256'); 
title('S2')
axis off;
axis equal
light
   

%Write Output

filename = sprintf('Kresling_S2.FOLD');
fid = fopen(filename,'w');
fprintf(fid,'{\n');
FOLD_count = FOLD_count + 1;
filenameList{FOLD_count} = filename;
RStore(FOLD_count) = R;
save('BellowsLandscapeValidOptions.mat','filenameList','RStore')

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
    fprintf(fid,'    [%g,%g,%g],\n',faces.vertices(ii,1),faces.vertices(ii,2),faces.vertices(ii,3));
end
fprintf(fid,'    [%g,%g,%g]\n',faces.vertices(ii+1,1),faces.vertices(ii+1,2),faces.vertices(ii+1,3));
fprintf(fid,'  ],\n');

fprintf(fid,'}');
[~] = fclose(fid);


        
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
 %S3
  diagnalphi1=180+atan((tan(phi1)*tan(phi2))/(2*tan(phi2)-tan(phi1)))/pi*180;
  phi2=phi1;
  phi1=diagnalphi1/180*pi;

 
        unitcell = [];
        unitcell2 = [];
        faces = [];
        vertices = [];
        
   
        x1 = tan(phi2)/tan(phi1);
          x2= tan(pi/n);
          beta1 = 2*atan(1/(2*x2)*(1-x1-(sqrt((x1-1)^2 - 4*x1*x2^2))));
          beta2 = 2*atan(1/(2*x2)*(1-x1+(sqrt((x1-1)^2 - 4*x1*x2^2))));
          beta = beta1;    
          alpha = asin(tan(beta/2)/tan(phi2));
          alphadeg=90-alpha/pi*180;
        nodes(:,1) = [0 0 0]';      
        nodes(:,2) =Rmat(-beta,[0 0 1])*[d 0 0]';  
        nodes(:,3) = Rmat(alpha,[1 0 0])* Rmat(-phi1,[0 1 0])*(H*csc(phi1))*[1 0 0]';
        nodes(:,4) = Rmat(alpha,[1 0 0])* Rmat(-phi2,[0 1 0])*(H*csc(phi2))*[1 0 0]';
        %nodes(:,5) =  Rmat(beta2,[0 0 1])*[d 0 0]'; 
        %height after folding
        hdeploy=range(nodes(3,:));    
         % first unit cell
        unitcell(1).vertices = nodes';       
        unitcell(1).faces = [1 2 4;4 1 3];
        % 'unitcell.faces' is a reference for the first cell whereas
        % 'unitcell2.faces' is the actual node structure for the whole
        % bellows        
        unitcell2(1).faces = [1 2 4;4 1 3];
        bottomNodes = unitcell(1).vertices(1:2,:);
        % next unit cells in a layer: rotate and attach
        Node = unitcell(1).vertices;
        Panel = num2cell(unitcell2(1).faces,2);
   
 % Covert into the calculated model 
 %https://journals.aps.org/pre/abstract/10.1103/PhysRevE.100.063001
        layer = m+1;  %number of layer
        a=R;b=R;
        c=sqrt((nodes(1,1)-nodes(1,4))^2+(nodes(2,1)-nodes(2,4))^2+(nodes(3,4)-hdeploy-nodes(3,1)).^2);
        theta=acos((a^2+b^2-c^2)/(2*a*b)); % rotation angle between layer-1 and layer            
       
%%%calculation
angle = zeros(layer,n);
zori = zeros(layer,n);


  for i = 1:n 
   for t = 1:layer
        if mod(t,2) == 0
        angle(t,:) = [theta: 2*pi/n : theta+2*pi/n*(n-1)]; 
        %angle(t,:) = [(-1)*(t-1)*theta: 2*pi/n : (-1)*(t-1)*theta+2*pi/n*(n-1)]; 
        else
        angle(t,:) = [0: 2*pi/n : 2*pi/n*(n-1)];
       % angle(t,:) = [(-1)*(t-1)*theta: 2*pi/n : (-1)*(t-1)*theta+2*pi/n*(n-1)]; 
        end
   zori(t,:) = [(t-1)*hdeploy];
   end 
end   

angleres = reshape(angle',layer*n,1);
%%%coordinate of each point
x = R*cos(angleres);
y = R*sin(angleres);
z = reshape(zori',layer*n,1);
Node = [x,y,z];
%Rotate for an angle to twist
 Node = Rmat(-theta,[0 0 1])*Node'; %%%S2 n4
 Node = Node';
vertices.coords = Node;

Nop = zeros(2*n*(layer-1),3);  %panel:connection of  nodes
for t = 1:layer-1
for i = 1:n
  if mod(t,2)~= 0
   if i == 1
     Nop(2*(t-1)*n+i,:) = [i+(t-1)*n;i+n+(t-1)*n; n+t*n;];
    Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+1+(t-1)*n; i+n+(t-1)*n];
    elseif mod(i,n)~= 0
    Nop(2*(t-1)*n+i,:) = [i+(t-1)*n;i+n+(t-1)*n-1; i+n+(t-1)*n;];
    Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+1+(t-1)*n; i+n+(t-1)*n];
    else
    Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+n+(t-1)*n; i+n+(t-1)*n-1];
    Nop(2*(t-1)*n+i+n,:)= [i+(t-1)*n; i+(t-1)*n-(n-1); i+n+(t-1)*n];
    end  
  else
    if i == 1
    Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i+(n-1)+(t-1)*n; i+n+(t-1)*n];
    Nop(2*(t-1)*n+i+n,:)= [i+(n-1)+(t-1)*n; i+n+(t-1)*n; i+(2*n-1)+(t-1)*n];
    else
     Nop(2*(t-1)*n+i,:) = [i+(t-1)*n; i-1+(t-1)*n; i+n+(t-1)*n];
     Nop(2*(t-1)*n+i+n,:)= [i+n+(t-1)*n; i-1+(t-1)*n; i+(n-1)+(t-1)*n];
    end 
  end
end
end   

faces.vertices=Nop;      
Panel=num2cell(Nop,2);

%Plot to check
Node = vertices.coords;
Panel = num2cell(faces.vertices,2);

% 
figure(3)
PlotOri(Node,Panel,[],'PanelColor','[0.67, 0.83, 0.9]/256'); 
title('S3')
axis off;
axis equal
light  
 
 
 %Write Output

filename = sprintf('Kresling_S3.FOLD');
fid = fopen(filename,'w');
fprintf(fid,'{\n');
FOLD_count = FOLD_count + 1;
filenameList{FOLD_count} = filename;
RStore(FOLD_count) = R;
save('BellowsLandscapeValidOptions.mat','filenameList','RStore')

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
    fprintf(fid,'    [%g,%g,%g],\n',faces.vertices(ii,1),faces.vertices(ii,2),faces.vertices(ii,3));
end
fprintf(fid,'    [%g,%g,%g]\n',faces.vertices(ii+1,1),faces.vertices(ii+1,2),faces.vertices(ii+1,3));
fprintf(fid,'  ],\n');

fprintf(fid,'}');
[~] = fclose(fid);

    
end   