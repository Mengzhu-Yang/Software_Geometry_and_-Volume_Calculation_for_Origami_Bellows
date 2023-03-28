function [Rs] =KreslingFlatS1(phi1,phi2,r,n,m,d1overr)
d=2*pi*r/n;
% Covert into the calculated model
%https://journals.aps.org/pre/abstract/10.1103/PhysRevE.100.063001
Rs = d/(2*sind(180/n));   %radius

end