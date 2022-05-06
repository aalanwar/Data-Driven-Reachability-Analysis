function E = logRemainder(intMat,maxOrder,maxabs)
%there is a problem in this file which is 
%compute absolute value bound
% M = abs(intMat);
% %next line problem
% M = logm(M);
M = maxabs;
%compute exponential matrix
eM =expm(M);

% no value is infinity
if ~any(any(isnan(eM)))

    %compute first Taylor terms
    Mpow = eye(intMat.dim);
    eMpartial = eye(intMat.dim);
    for i=1:maxOrder
        Mpow = M*Mpow;
        eMpartial = eMpartial + Mpow/factorial(i);
    end

    W = eM-eMpartial;

    %instantiate remainder
    E = intervalMatrix(zeros(intMat.dim),W);
else
    %instantiate remainder
    E = intervalMatrix(zeros(intMat.dim),inf*ones(intMat.dim));
end

%------------- END OF CODE --------------