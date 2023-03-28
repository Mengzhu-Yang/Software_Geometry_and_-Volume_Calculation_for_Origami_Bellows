function R = Rmat(phi,n)
% ROTATION_MATRIX
%
%     R = Rmat(phi,n)
%
%     Creates rotation matrix R using Euler Parameters, see:
%     http://mathworld.wolfram.com/EulerParameters.html
%
%     3D rotation axis is given by n, and rotation angle by phi
%
%     Author:  Mark Schenk (ms652)
%     Created: 28 January 2008

n = n/norm(n);
e(1) = cos(phi/2);
e(2:4) = sin(phi/2)*n;

e0 = e(1);
e1 = e(2);
e2 = e(3);
e3 = e(4);

R(1,1) = e0^2+e1^2-e2^2-e3^2;
R(1,2) = 2*(e1*e2-e0*e3);
R(1,3) = 2*(e1*e3+e0*e2);
R(2,1) = 2*(e2*e1+e0*e3);
R(2,2) = e0^2-e1^2+e2^2-e3^2;
R(2,3) = 2*(e2*e3-e0*e1);
R(3,1) = 2*(e3*e1-e0*e2);
R(3,2) = 2*(e3*e2+e0*e1);
R(3,3) = e0^2-e1^2-e2^2+e3^2;
end
