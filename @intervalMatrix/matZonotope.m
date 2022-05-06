function matZ = matZonotope(matI)
% matZonotopes - converts an interval matrix to a matrix zonotope
%
% Syntax:  
%    matZ = matZonotope(matI)
%
% Inputs:
%    matI - interval matrix 
%
% Outputs:
%    matZ - zonotope matrix
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Matthias Althoff
% Written:      21-June-2010 
% Last update:  25-July-2016 (intervalhull replaced by interval)
% Last revision:---

%------------- BEGIN CODE --------------
% obtain center
% c = center(matI);
% 
% % construct generator matrix G
% G = diag(rad(matI));
% 
% % instantiate zonotope
% Z = matZonotope(Z,6);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %convert to interval
temp = matI.int;
col = size(temp.sup,2);
I=interval(matI);

%convert to zonotope
Z=zonotope(I);
Z=deleteZeros(Z); %delete zero generators

%convert to matrix zonotope
matZ=matZonotopeConv(Z,col);%6

%------------- END OF CODE --------------