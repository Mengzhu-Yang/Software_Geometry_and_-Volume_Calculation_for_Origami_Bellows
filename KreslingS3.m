function [Node, Panel,V3,Vd, Hd, Rd] = KreslingS3(phi1,phi2,r,n,m,Output)
%   n=4
%  m=3
%  phi1=75/180*pi
%  phi2=30/180*pi
%  d1overr=1
%  r=40     
% Output=0

%S3
diagnalphi1=180+atan((tan(phi1)*tan(phi2))/(2*tan(phi2)-tan(phi1)))/pi*180;
phi2=phi1;
phi1=diagnalphi1/180*pi;

d=2*pi*r/n;
H=d/(1/(tan(phi2))-1/(tan(phi1)));

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
Hd=hdeploy;
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
R = d/(2*sind(180/n));   %radius
Rd=R;
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
 


%%%Calculate V3 and Vd
Node1=Node(1:2*n,:);
CTbottom=[0 0 0];
CTtop=[0 0 hdeploy];
Node1=[Node1;CTbottom;CTtop];
  A0=[0 0];
 %calculate C1 :the intersection point of the projection of A1A2 and line B1B2
   X1=Node1(1,1:2);
   Y1=Node1(2,1:2);
   X2=Node1(n+1,1:2);
   Y2=Node1(n+2,1:2);

if X1(1)==Y1(1)
   C(1,1)=X1(1);
   k2=(Y2(2)-X2(2))/(Y2(1)-X2(1));
   b2=X2(2)-k2*X2(1); 
   C(1,2)=k2*C(1,1)+b2;
end
if X2(1)==Y2(1)
   C(1,1)=X2(1);
   k1=(Y1(2)-X1(2))/(Y1(1)-X1(1));
   b1=X1(2)-k1*X1(1);
   C(1,2)=k1*C(1,1)+b1;
end
if X1(1)~=Y1(1)& X2(1)~=Y2(1)
   k1=(Y1(2)-X1(2))/(Y1(1)-X1(1));
   k2=(Y2(2)-X2(2))/(Y2(1)-X2(1));
   b1=X1(2)-k1*X1(1);
   b2=X2(2)-k2*X2(1);  
   C(1,1)=(b2-b1)/(k1-k2);
   C(1,2)=k1*C(1,1)+b1;
end  


%calculate C2 :the intersection point of the projection of A2A3 and line B1B2
   X1=Node1(3,1:2);
   Y1=Node1(2,1:2);
   X2=Node1(n+1,1:2);
   Y2=Node1(n+2,1:2);


if X1(1)==Y1(1)
   C2(1,1)=X1(1);
   k2=(Y2(2)-X2(2))/(Y2(1)-X2(1));
   b2=X2(2)-k2*X2(1); 
   C2(1,2)=k2*C2(1,1)+b2;
end
if X2(1)==Y2(1)
   C2(1,1)=X2(1);
   k1=(Y1(2)-X1(2))/(Y1(1)-X1(1));
   b1=X1(2)-k1*X1(1);
   C2(1,2)=k1*C2(1,1)+b1;
end
if X1(1)~=Y1(1)& X2(1)~=Y2(1)
   k1=(Y1(2)-X1(2))/(Y1(1)-X1(1));
   k2=(Y2(2)-X2(2))/(Y2(1)-X2(1));
   b1=X1(2)-k1*X1(1);
   b2=X2(2)-k2*X2(1);  
   C2(1,1)=(b2-b1)/(k1-k2);
   C2(1,2)=k1*C2(1,1)+b1;
end  
  p=(norm(A0-C)+norm(A0-C2)+norm(C2-C))/2;
  s=sqrt(p*(p-norm(A0-C))*(p-norm(A0-C2))*(p-norm(C2-C)));
   Vd=s*2*n*hdeploy*m;

%    Vds3=2*n*h*(norm(A0-C))^2*sin(pi/(2*n))*cos(pi/(2*n));

Vmath=n*R^2*hdeploy*sin(pi/n)*cos(pi/n)*2/3;
[Vadd]=calvolume(Node1(2,:),Node1(n+1,:),Node1(n+2,:),Node1(2*n+1,:));
Vaddmath=1/6*R^2*hdeploy*abs((sin(theta)+sin(-theta+2*pi/n)));
V3=(Vmath+n*Vaddmath)*m;



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

frame.title = '"Kresling_S3"';
% This code is geared towards only doing the folded form here
frame.classes = '"foldedForm"';

% Again this code is only aimed towards 3D here
frame.attributes = '"3D"';

% NOTE: all units are in mm throughout
frame.unit = '"mm"';

FOLD_count = 0;

%Write Output

filename = sprintf('Kresling_S3.FOLD');
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