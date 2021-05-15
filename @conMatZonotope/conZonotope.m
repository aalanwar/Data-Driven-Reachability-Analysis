function Z = conZonotope(CmatZ)
% zonotope - Converts a constrained matrix zonotope into a constranied zonotope 
%
% Syntax:  
%    Z = zonotope(matZ)
%
% Inputs:
%    matZ - matrix zonotope
%
% Outputs:
%    Z - zonotope
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Amr Alanwar
% Written:      18-December-2020 
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

%concert center
center=mat2vec(CmatZ.center);

%convert generators
for i=1:CmatZ.gens
    generatorMatrix(:,i) = mat2vec(CmatZ.generator{i});
end

for i=1:length(CmatZ.A)
    AMatrix(:,i) = mat2vec(CmatZ.A{i}); %TODO remove the {1}
end
Bvec=  mat2vec(CmatZ.B);
%instantiate constrained zonotope
Z=conZonotope([center,generatorMatrix],AMatrix,Bvec);

%------------- END OF CODE --------------