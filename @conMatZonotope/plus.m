function R = plus(summand1,summand2)
% plus - Overloaded '+' operator for the Minkowski addition of a
%        constrained matrix zonotope with other set representations
%
% Syntax:
%    cZ = plus(summand1,summand2)
%
% Inputs:
%    summand1 - conMatZonotope object or numerical vector
%    summand2 - conMatZonotope object or numerical vector
%
% Outputs:
%    R - constrained matrix Zonotope after Minkowski addition
%
% References:
%    [1] Amr Alanwar, Anne Koch, Frank Allgöwer,Karl Henrik Johansson
%   "Data-Driven Reachability Analysis from Noisy Data"
%
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: mtimes
%
% References:
%
% Author:       Amr Alanwar
% Written:  12-December 2020
% Last update: 9-May 2021
%
% Last revision:---

%------------- BEGIN CODE --------------

% find a conMatZonotope object
if isa(summand1, 'conMatZonotope')
    matCZ=summand1;
    summand=summand2;
elseif isa(summand2, 'conMatZonotope')
    matCZ=summand2;
    summand=summand1;
end

% handle different classes of the second summand
if isa(summand, 'conMatZonotope')
    
    %Calculate minkowski sum
    %Add the centers
    newcenter = matCZ.center + summand.center;
    
    if isempty(matCZ.generator)
        newgenerator = summand.generator;
        newA = summand.A;
        newB = summand.B;
    elseif isempty(summand.generator)
        newgenerator = matCZ.generator;
        newA = matCZ.A;
        newB = matCZ.B;
    else
        %concatenate matrix generators
        newgenerator = matCZ.generator;
        newgenerator((end+1):(end+summand.gens)) = summand.generator;
        for i=1:length(matCZ.A)
            newA{i} =blkdiag(matCZ.A{i},zeros(size(summand.A{1})));
        end
        index=1;
        for i=length(matCZ.A)+1:length(matCZ.A)+length(summand.A)
            newA{i} =blkdiag(zeros(size(matCZ.A{1})),summand.A{index});
            index = index +1;
        end
        newB =blkdiag(matCZ.B,summand.B);
    end
    
    
    R= conMatZonotope(newcenter,newgenerator,newA,newB);

elseif isa(summand, 'matZonotope')
    newcenter = matCZ.center + summand.center;
    
    if isempty(matCZ.generator)
        newgenerator = summand.generator;
        R= conMatZonotope(newcenter,newgenerator,matCZ.A,matCZ.B);
    elseif isempty(summand.generator)
        newgenerator = matCZ.generator;
        R= conMatZonotope(newcenter,newgenerator,matCZ.A,matCZ.B);
    else        
        %concatenate matrix generators
        newgenerator = matCZ.generator;
        newA = matCZ.A;
        newgenerator((end+1):(end+summand.gens)) = summand.generator;
        
        %add zeros to A to have the same number of matrices 
        numOldA = length(matCZ.A);
        for i=1:summand.gens
            newA{numOldA+i} = zeros(size(summand.generator{1}));
        end
        
        
        R= conMatZonotope(newcenter,newgenerator,newA,matCZ.B);
    end

    
elseif isnumeric(summand)
    R= conMatZonotope(matCZ.center+summand,matCZ.generator,matCZ.A,matCZ.B);
end

%------------- END OF CODE --------------