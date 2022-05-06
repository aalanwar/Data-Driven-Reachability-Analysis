function res = matZonotope(obj,varargin)
% zonotope - over-approximates a constrained matrix zonotope with a matrix zonotope
%
% Syntax:
%    res = matZonotope(obj)
%    res = matZonotope(obj,alg)
%
% Inputs:
%    obj - conZonotope object
%    alg - algorithm used to compute enclosure ('nullSpace' or 'reduce')
%
% Outputs:
%    res - zonotope object
%
% Example:
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: interval

% Author:
% Written:
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if nargin >= 2
    alg = varargin{1};
else
    alg = 'nullSpace';
end

if(strcmp(alg,'liftup'))
    % Z_up = zonotope([X_data_cmz{i+1,1}.Z; -X_data_cmz{i+1,1}.b, X_data_cmz{i+1,1}.A]);
    
    newcen = [obj.center;-obj.B{1}];
    for i=1:length(obj.generator)
      newgen{i} = [obj.generator{i};obj.A{1}{i} ] ;
    end
   res =  matZonotope(newcen,newgen);
else
    numOfColc = size(obj.center,2);
    CZ=conZonotope(obj);
    alg = 'nullSpace';
    if nargin >= 3
        alg = varargin{2};
    end
    %convert to zonotope
    CZred = zonotope(CZ, alg);
    
    %convert back to matrix zonotope
    res = matZonotopeConv(CZred,numOfColc);
end
% parse input arguments
end

%------------- END OF CODE --------------