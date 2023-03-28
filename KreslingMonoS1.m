function [Node, Panel,V1,Vs,Hs,Rs] = KreslingMonoS1(phi1,phi2,r,n,m,Output)
%   n=4
%  m=3
%  phi1=75/180*pi
%  phi2=30/180*pi
%  d1overr=1
%  r=40     
%% 1 solution 
      x1 = tan(phi2)/tan(phi1);
      x2= tan(pi/n);
      beta1 = 2*atan(1/(2*x2)*(1-x1-(sqrt((x1-1)^2 - 4*x1*x2^2))));
      beta2 = 2*atan(1/(2*x2)*(1-x1+(sqrt((x1-1)^2 - 4*x1*x2^2))));
      beta = beta1;        
            %1st solution (S2 more deployed)： beta = beta1
            %2nd solution (S1 more compressed)： beta = beta2
            % alpha is angular deflection from vertical (z-axis) 
            %pi-alpha is the opening angle between adjacent bands of folds
          alpha = asin(tan(beta/2)/tan(phi2));
%            %%%alpha will sometimes be imaginary at flatfoldable state
%           alpha = real(alpha)
         
%         %constraints on alpha beta
%           tf = isreal(beta) ;
%         if tf==0
%             noWrite = 1;     
%         end
%         
%         tf2 = isreal(alpha) ;
%         if tf2==0
%             noWrite = 1;     
%         end
%          
%         if noWrite == 1
%             Objective=[100000];
%         end
            
         %Calculate the length d
        % r=n*R*sin(pi/n)/pi;
        % R=r*pi/(n*sin(pi/n));
         d=2*pi*r/n;
         H=d/(1/(tan(phi2))-1/(tan(phi1)));
        % calculate Vstowed
        
        %An unit cell in 3D   
        nodes(:,1) = [0 0 0]';      
        nodes(:,2) =Rmat(-beta,[0 0 1])*[d 0 0]';  
        nodes(:,3) = Rmat(alpha,[1 0 0])* Rmat(-phi1,[0 1 0])*(H*csc(phi1))*[1 0 0]';
        nodes(:,4) = Rmat(alpha,[1 0 0])* Rmat(-phi2,[0 1 0])*(H*csc(phi2))*[1 0 0]'; 
        %height after folding
        h=range(nodes(3,:)); % height of a layer
         Hs=h;
   %Vstowed=m*h*pi*R^2;
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
   
%         figure()
%         PlotOri(Node,Panel,[],'PanelColor','g'); 
%         axis off;
%         light
%         PlotOri(Node,Panel,[],'ShowNumber','on');
%         light
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
        Rs=R;
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

%%%Calculate V1 and Vs
Node1=Node(1:2*n,:);
Nop = zeros(2*n*(2-1),3);  %panel:connection of  nodes
if phi1>pi/2
for t = 1:2-1
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
for t = 1:2-1
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


%faces.vertices=Nop;      
Panel1=num2cell(Nop,2);
%Panel1 = num2cell(faces.vertices,2);
Vs=h*pi*R^2*m;

%% Add the middle points
     Nn = size(Node1,1);
    pseudofold = nan(6*length(Panel1),4);
    panelctr = nan(length(Panel1),1);
    cntb = 0;   cntc = 0;

 for i = 1:length(Panel1)
    cntc = cntc+1;

      
           if i>n
          % if i<n
            L1 = norm((Node1(Panel1{i}(1),:)-Node1(Panel1{i}(3),:)));
            L2 = norm((Node1(Panel1{i}(1),:)-Node1(Panel1{i}(2),:)));
            if L1>L2            
             mm = [Node1(Panel1{i}(3),:)-Node1(Panel1{i}(1),:)]'/L1; %vector1
             coeff = ([mm]\[Node1(Panel1{i}(3),:)-Node1(Panel1{i}(1),:)]')/2;
              ctrnode = Node1(Panel1{i}(1),:)+coeff(1)*mm';
            else
             mm = [Node1(Panel1{i}(2),:)-Node1(Panel1{i}(1),:)]'/L2; %vector2
             coeff = ([mm]\[Node1(Panel1{i}(2),:)-Node1(Panel1{i}(1),:)]')/2;
             ctrnode = Node1(Panel1{i}(1),:)+coeff(1)*mm';   
            end 
         else
             L1 = norm((Node1(Panel1{i}(1),:)-Node1(Panel1{i}(3),:)));
            L2 = norm((Node1(Panel1{i}(3),:)-Node1(Panel1{i}(2),:)));
            if L1>L2            
             mm = [Node1(Panel1{i}(3),:)-Node1(Panel1{i}(1),:)]'/L1; %vector1
             coeff = ([mm]\[Node1(Panel1{i}(3),:)-Node1(Panel1{i}(1),:)]')/2;
             ctrnode = Node1(Panel1{i}(1),:)+coeff(1)*mm';
            else
             mm = [Node1(Panel1{i}(2),:)-Node1(Panel1{i}(3),:)]'/L2; %vector2
             coeff = ([mm]\[Node1(Panel1{i}(2),:)-Node1(Panel1{i}(3),:)]')/2;
             ctrnode = Node1(Panel1{i}(3),:)+coeff(1)*mm';   
            end 
           end
        


            %judge if there are repeated nodes 
                judge=sum((Node1 - repmat(ctrnode,length(Node1),1)),2);
               %if ismember(0,judge)==1
               if min(abs(judge))<1e-8
                  %  [m,n]=find(abs(judge)<1e-8);
                  [mm,nn]=find(abs(judge)<1e-8);
                    Panelctr(i)=mm;  
                    cntc=cntc-1;
                else
                   Node1(Nn+cntc,:) = ctrnode;
                    Panelctr(i)=Nn+cntc;  
                end
            
            for k=1:numel(Panel1{i})         
                cntb = cntb+1;
                ref1 = mod(k-2,numel(Panel1{i}))+1;
                ref2 = mod(k,numel(Panel1{i}))+1;               
               pseudofold(cntb,:) = [panelctr(i),Panel1{i}(k),Panel1{i}(ref1),Panel1{i}(ref2)];
            end
end
  pseudofold(isnan(pseudofold(:,1)),:)=[];  
    panelctr = [];
    pseudofold(isnan(pseudofold(:,1)),:)=[];
    pseudofold1=pseudofold;
    pseudofold=zeros(length(pseudofold1)/3*2,4);
    for ii= 0:(length(pseudofold1)/3-1)
        pseudofold(2*ii+1,:)=pseudofold1(3*ii+2,:);
        pseudofold(2*ii+2,:)=[pseudofold1(3*ii+3,1),pseudofold1(3*ii+3,3),pseudofold1(3*ii+3,2),pseudofold1(3*ii+3,4)];
    end   
               
Panel1=num2cell(pseudofold,2);
Triglraw = [pseudofold(:,[1,2,3])];
Panel1 = num2cell(sort(Triglraw,2),2);   




CTbottom=[0 0 0];
CTtop=[0 0 h];
CTmiddle=[0 0 h/2];
Node1=[Node1;CTbottom;CTtop;CTmiddle];
% set(0,'defaultfigurecolor','w')
% figure()
% PlotOri(Node1,Panel1,[],'PanelColor','w');
% axis equal; axis off;
% light
% % Inspect nodal index assignment
% figure()
% PlotOri(Node,Panel,[],'ShowNumber','on');
% axis equal
% hold on

% Vmath=n*R^2*h*sin(pi/n)*cos(pi/n)*2/3;
% [Vadd]=calvolume(Node1(1,:),Node1(n+1,:),Node1(n+2,:),Node1(3*n+1,:));
% Vaddmath=1/6*R^2*h*abs((sin(theta)-sin(theta+2*pi/n)));
V0M=[Node1(3*n+2,:),1; Node1(n+1,:),1; Node1(n+2,:),1; Node1(3*n+1,:),1;];
V0=1/6*det(V0M);
V1M=[Node1(2*n+1,:),1; Node1(n+1,:),1; Node1(n+2,:),1; Node1(3*n+1,:),1;];
V11=1/6*det(V1M);
V2M=[Node1(3*n+2,:),1; Node1(2*n+1,:),1; Node1(n+2,:),1; Node1(3*n+1,:),1;];
V22=1/6*det(V2M);
V3M=[Node1(3*n+2,:),1; Node1(n+1,:),1; Node1(2*n+1,:),1; Node1(3*n+1,:),1;];
V33=1/6*det(V3M);
V4M=[Node1(3*n+2,:),1; Node1(n+1,:),1; Node1(n+2,:),1; Node1(2*n+1,:),1;];
V44=1/6*det(V4M);

if V0*V11>0 & V0*V22>0 & V0*V33>0  & V0*V44>0             %Inside
[Vadd1]=calvolume(Node1(1,:),Node1(2,:),Node1(3*n+1,:),Node1(3*n,:));
[Vadd2]=calvolume(Node1(2,:),Node1(2*n+1,:),Node1(3*n,:),Node1(3*n+1,:));
[Vadd3]=calvolume(Node1(2*n+1,:),Node1(3*n+1,:),Node1(3*n,:),Node1(3*n+3,:));
[Vadd4]=calvolume(Node1(1,:),Node1(2*n+1,:),Node1(3*n,:),Node1(n+1,:));
V1=((Vadd1+Vadd2+Vadd3)*2+Vadd4)*n*m;
else
    
Vmath=n*R^2*h*sin(pi/n)*cos(pi/n)*2/3;
%[Vadd]=calvolume(Node1(1,:),Node1(n+1,:),Node1(n+2,:),Node1(3*n+1,:));
Vaddmath=1/6*R^2*h*abs((sin(theta)-sin(theta+2*pi/n)));
V1=m*(Vmath+n*Vaddmath);
end

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

frame.title = '"Kresling_S1"';
% This code is geared towards only doing the folded form here
frame.classes = '"foldedForm"';

% Again this code is only aimed towards 3D here
frame.attributes = '"3D"';

% NOTE: all units are in mm throughout
frame.unit = '"mm"';

FOLD_count = 0;

%Write Output

filename = sprintf('Kresling_S1.FOLD');
fid = fopen(filename,'w');
fprintf(fid,'{\n');
FOLD_count = FOLD_count + 1;
filenameList{FOLD_count} = filename;
RStore(FOLD_count) = R;
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
    fprintf(fid,'    [%g,%g,%g],\n',faces.vertices(ii,1),faces.vertices(ii,2),faces.vertices(ii,3));
end
fprintf(fid,'    [%g,%g,%g]\n',faces.vertices(ii+1,1),faces.vertices(ii+1,2),faces.vertices(ii+1,3));
fprintf(fid,'  ],\n');

fprintf(fid,'}');
[~] = fclose(fid);
end
end
