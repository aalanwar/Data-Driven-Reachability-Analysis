function [Rtp ,Rtp_data] = linReach_DT(obj,Rinit,R_data,options)
% linReach - computes the reachable set after linearization
%
% Syntax:  
%    [Rtp] = linReach_DT(obj,Rinit,options)
%
% Inputs:
%    obj - nonlinearSysDT system object
%    Rinit - initial reachable set
%    options - options struct
%
% Outputs:
%    Rtp - resulting reachable set
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 

% Author:       Amr Alanwar, Matthias Althoff, Niklas Kochdumper 
% Written:      21-August-2012
% Last update:  29-January-2018 (NK)
%               29-October-2020 (Amr) add data driven reachability
% Last revision:---

%------------- BEGIN CODE --------------

% linearize nonlinear system
[obj,A_lin,U] = linearize_DT(obj,Rinit,options); 

%translate Rinit by linearization point
Rdelta = Rinit + (-obj.linError.p.x);

% compute reachable set of linearized system
Rtp = A_lin*Rdelta + U;

% obtain linearization error
if options.tensorOrder > 2
    Verror = linError_thirdOrder(obj, options, Rdelta); 
else
    Verror = linError_mixed_noInt_DT(obj, options, Rdelta);   
end


%add interval of actual error
Rtp=Rtp+Verror+options.W;
% %%%%%%%%-------------------Data driven reachability-----------
options.Uorig= options.U +  options.uTrans;
 xStar = R_data.center;
 uStar =options.Uorig.center;
 xStarMat = repmat(xStar,1,size(options.X_0T,2));
 uStarMat = repmat(uStar,1,size(options.U_full,2));
 oneMat = repmat([1],1,size(options.U_full,2));
 IAB = (options.X_1T )*pinv([oneMat; options.X_0T+(-1*xStarMat);options.U_full+-1*uStarMat]);

V =  options.X_1T + -1*(IAB*[oneMat; options.X_0T+(-1*xStarMat);options.U_full+-1*uStarMat] + options.Wmatzono);
 VInt = intervalMatrix(V);
 leftLimit = VInt.Inf;
 rightLimit = VInt.Sup;
 
 V_one= zonotope(interval(min(leftLimit')',max(rightLimit')'));
 
 
 
 Rtp_data = IAB*cartProd([1],cartProd(R_data+(-1*xStar),options.Uorig+(-1*uStar))) +V_one+ options.W ;



%------------- END OF CODE --------------