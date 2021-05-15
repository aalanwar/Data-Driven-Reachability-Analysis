classdef conMatZonotope < matZonotope
    % conZonotope - object constructor for constrained matrix zonotopes [1]
    %
    % Syntax:
    %    C = [0 0; 0 0];
    %    G{1} = [1 3; -1 2];
    %    G{2} = [2 0; 1 -1];
    %    A{1} = [1 1; 0 1];
    %    A{2} = [0.3 0.2;-0.2 0.1];
    %    B = [0.5 0.4; 0.3 0.2];
    %
    %    cmz = conMatZonotope(C,G,A,B);
    %
    % Outputs:
    %    cmz - Generated Object
    %
    % Example:
    %
    % References:
    %    [1] Amr Alanwar, Anne Koch, Frank Allgöwer,Karl Henrik Johansson
    %   "Data-Driven Reachability Analysis from Noisy Data"
    %
    % Other m-files required: none
    % Subfunctions: none
    % MAT-files required: none
    %
    % See also: 
    %
    % Author:        Amr Alanwar
    % Written:       03-December-2020
    % Last update:   ---
    % Last revision:
    
    %------------- BEGIN CODE --------------
    
    properties (SetAccess = private, GetAccess = public)
        A = [];
        B = [];
    end
    
    methods
        
        % class constructor
        function obj = conMatZonotope(varargin)
            
            if nargin==3
                if isa(varargin{1},'conZonotope')
                    CZ = varargin{1};
                    numofcolc = varargin{2};
                    numofcolA = varargin{3};
                    %extract center
                    matrixCenter=CZ.Z(:,1);
                    %extract generator matrix
                    oldgenerator=CZ.Z(:,2:end);
                    %obtain matrix center
                    center = vec2mat(matrixCenter,numofcolc);
                    %obtain matrix generators
                    if ~isempty(oldgenerator)
                        for i=1:length(oldgenerator(1,:))
                            generator{i}=vec2mat(oldgenerator(:,i),numofcolc);
                        end
                    else
                        generator{1} = zeros(size(matrixCenter));
                    end
                    
                    
                    if ~isempty(CZ.A)
                        for i=1:length(CZ.A(1,:))
                            A{i}=vec2mat(CZ.A(:,i),numofcolA);
                        end
                    else
                        A{1} = zeros(size(matrixCenter));
                    end
                    B = vec2mat(CZ.b,numofcolA);                   
                end
            elseif nargin == 4
                center = varargin{1};
                generator = varargin{2};
                A = varargin{3};
                B = varargin{4};
            else
                error('This class takes 4 inputs.')
            end
            
            % create a zonotope
            obj@matZonotope(center,generator);
            

            obj.A = A;
            obj.B = B;
        end
        
        
        % methods in seperate files
        cZ = plus(summand1,summand2);
        cZ = mtimes(factor1,factor2);
        res = reduce(obj,method,orderG,varargin);
        res = conZonotope(obj,varargin);
        %     res = and(obj,S);
        %     cZ = cartProd(cZ1,cZ2);
        %     res = center(obj);
        %     cZ = convHull(cZ1,varargin);
        %     d = dim(obj);
        %     display(obj);
        %     cZ = enclose(varargin);
        %     res = in(obj1,obj2,varargin);
        %     res = interval(obj);
        %     res = intervalMultiplication(obj,I);
        %     res = isempty(obj);
        %     res = isFullDim(obj);
        %     res = isIntersecting(obj1,obj2,varargin);
        %     res = mptPolytope(obj);
        
        %     cZ = or(cZ1, varargin);
        %     handle = plot(obj,varargin);
        %     handle = plotZono(obj,varargin);
        %     pZ = polyZonotope(obj);
        %     obj = project(obj,dims);
        %     cZquad = quadMap(varargin);
        
        %     res = rescale(obj,varargin);
        %     cZsplit = split(obj,varargin);
        %     [val,x,ksi] = supportFunc(obj,dir,varargin);
        %     V = vertices(obj);
        %     res = zonoBundle(obj);
        
        
    end
    
    
end

%------------- END OF CODE -------