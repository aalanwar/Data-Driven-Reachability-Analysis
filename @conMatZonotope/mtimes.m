function cMZ = mtimes(factor1,factor2)
% mtimes - Overloaded '.*' operator for the multiplication of a matrix or an
%          interval matrix with a constrained matrix zonotope
%
% Syntax:
%    cZ = times(cmz,cZ)
%
% Inputs:
%    cmz - constrained matrix zonotope
%    cz - conZonotope or zonotope object
%
% Outputs:
%    cZ - constrained zonotpe after multiplication
%
% Example:
%
%    G = [0 1 0 1;0 1 2 -1];
%    A = [-2 1 -1];
%    b = 2;
%    cz = conZonotope(Z,A,b);
%    cMul = [3 1;2 4] * cz;
%
%    hold on
%    plot(cZono,[1,2],'r');
%    plot(cMul,[1,2],'b');
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: plus

% Author:       Amr Alanwar
% Written:
% Last update:  ---
% Last revision:---

%------------- BEGIN CODE --------------

if isa(factor2,'conZonotope')
    if isempty(factor2.A)
        % The multiplication consists of three terms as shown in [1]
        matrixczono = factor1;
        czono = factor2;
        % Znew has the center which is (cmz.center*cz.center)
        % add to Znew the cmz.center * cz.generator
        Znew=matrixczono.center*czono.Z;
        for i=1:matrixczono.gens
            % add to Znew the cmz.generator * cz.center
            Zadd=matrixczono.generator{i}*czono.Z(:,1);
            Znew(:,(end+1):(end+length(czono.Z(1,1))))=Zadd;
        end
        %find constrain from cmz
        Amat =[];
        for j=1:length(factor1.A)
            Amat(:,j) = mat2vec(factor1.A{j}) ;
        end
        bmat = mat2vec(factor1.B);
        %noOfGen = matrixczono.gens + size(czono.Z(:,2:end),2) +matrixczono.gens * size(czono.Z(:,2:end),2);
        %adding constraint in the following way
        % 0(for cz.A)  cmz.A  0(remaining)
        sizLeftZeros = size(czono.Z(:,2:end),2);
        sizeRightZeros = matrixczono.gens * size(czono.Z(:,2:end),2);
        Anew= [ zeros(size(Amat,1),sizLeftZeros),Amat,...
            zeros(size(Amat,1),sizeRightZeros)];
        %%
        Beta_matzono=minBeta(conZonotope(ones(2,size(Amat,2)+1),Amat,bmat));
        % find max beta
        index=1;
        for i =1:size(Amat,2)
            for j=1:size(czono.Z(:,2:end),2)
                Beta_t(index) = Beta_matzono(i)*interval(-1,1);%noneed for multiplication
                temp = Beta_t(index);
                maxfac(index) = max(abs(temp.inf),abs(temp.sup));
                index = index+1;
            end
        end
        %multiply betas and find the factor f
        index=1;
        for i=1:matrixczono.gens
            %generators * generators
            for j=1:size(czono.Z(:,2:end),2)
                Zadd=matrixczono.generator{i}*czono.Z(:,j+1)*maxfac(index);
                Znew(:,end+1)=Zadd;
                index = index+1;
            end
        end
        
        cMZ = conZonotope( Znew,...
            Anew,bmat);
    else
        matrixczono = factor1;
        czono = factor2;
        % The multiplication consists of three terms as shown in [1]
        % Znew has the center which is (cmz.center*cz.center)
        % add to Znew the cmz.center * cz.generator
        Znew=matrixczono.center*czono.Z;
        for i=1:matrixczono.gens
            % add to Znew the cmz.generator * cz.center
            Zadd=matrixczono.generator{i}*czono.Z(:,1);
            Znew(:,(end+1):(end+length(czono.Z(1,1))))=Zadd;
        end
        
        %Total number of generators = #of cmz.generators + #of cz.generator
        %+ #of cmz.generators * #of cz.generator
        noOfGen = matrixczono.gens + length(czono.Z(:,2:end)) +matrixczono.gens * length(czono.Z(:,2:end));
        % concatinate zeros to cz.A to have the same number of generator
        diffczA = noOfGen - size(czono.A,2);
        %the constraints from the cz
        Anew =[ czono.A, zeros(size(czono.A,1),diffczA)];
        Bnew =[ czono.b];
        %%
        %add the constraints from cmz
        Amat =[];
        for i=1:matrixczono.gens
            Amat(:,i) = mat2vec(matrixczono.A{i});
        end
        Bmat =  mat2vec(matrixczono.B);
        Bnew = [ Bnew;Bmat];
        
        %adding constraint in the following way
        % 0(for cz.A)  cmz.A  0(remaining)
        numOfRemainZeors = noOfGen  - size(Amat,2)-size(czono.A,2);
        Anew =[Anew;  zeros(size(Amat,1),size(czono.A,2)),Amat,...
            zeros(size(Amat,1),numOfRemainZeors)];
        %%
        % find beta for cz
        Beta_czono=minBeta(czono);
        % find beta for cmz
        Beta_matzono=minBeta(conZonotope(ones(2,size(Amat,2)+1),Amat,Bmat));
        
        %multiply betas and find the factor f
        index=1;
        for i =1:size(Amat,2)
            for j=1:size(czono.A,2)
                Beta_t(index) = Beta_matzono(i)*Beta_czono(j);
                temp = Beta_t(index);
                maxfac(index) = max(abs(temp.inf),abs(temp.sup));
                index = index+1;
            end
        end
        %add the remaining generators multiplied by a factor
        index=1;
        for i=1:matrixczono.gens
            %generators * generators
            for j=1:size(czono.A,2)
                Zadd=matrixczono.generator{i}*czono.Z(:,j+1)*maxfac(index);
                Znew(:,end+1)=Zadd;
                index = index+1;
            end
        end
        
        cMZ = conZonotope( Znew,...
            Anew,Bnew);
    end
elseif isa(factor2,'zonotope')
    %same code as constrained zonotope with empty A
    matrixczono = factor1;
    czono = factor2;
    Znew=matrixczono.center*czono.Z;
   
    for i=1:matrixczono.gens
        %mat generators * center
        Zadd=matrixczono.generator{i}*czono.Z(:,1);
        Znew(:,(end+1):(end+length(czono.Z(1,1))))=Zadd;
    end
    Amat =[];
    for j=1:length(factor1.A)
        Amat(:,j) = mat2vec(factor1.A{j}) ;
    end
    bmat = mat2vec(factor1.B); %*factor2.Z(:,1)];
    %noOfGen = matrixczono.gens + size(czono.Z(:,2:end),2) +matrixczono.gens * size(czono.Z(:,2:end),2);
    sizLeftZeros = size(czono.Z(:,2:end),2);
    sizeRightZeros = matrixczono.gens * size(czono.Z(:,2:end),2);
    Anew= [ zeros(size(Amat,1),sizLeftZeros),Amat,...
        zeros(size(Amat,1),sizeRightZeros)];
    %%
    Beta_matzono=minBeta(conZonotope(ones(2,size(Amat,2)+1),Amat,bmat));
    
    index=1;
    for i =1:size(Amat,2)
        for j=1:size(czono.Z(:,2:end),2)
            Beta_t(index) = Beta_matzono(i)*interval(-1,1);%noneed for multiplication
            temp = Beta_t(index);
            maxfac(index) = max(abs(temp.inf),abs(temp.sup));
            index = index+1;
        end
    end
    index=1;
    for i=1:matrixczono.gens
        %generators * generators
        for j=1:size(czono.Z(:,2:end),2)
            Zadd=matrixczono.generator{i}*czono.Z(:,j+1)*maxfac(index);
            Znew(:,end+1)=Zadd;
            index = index+1;
        end
    end

    cMZ = conZonotope( Znew,...
        Anew,bmat);
elseif  isnumeric(factor2)
    cMZCenGen = mtimes@matZonotope(factor1,factor2);
    cMZ = conMatZonotope(cMZCenGen.center,cMZCenGen.generator,factor1.A,factor1.B);
else
    cMZ = mtimes@matZonotope(factor1,factor2);
end

%------------- END OF CODE --------------