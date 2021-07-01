function [R ,R_data]= reach_DT(obj,params,options,varargin)
% reach - computes the reachable sets of the discrete time system
%
% Syntax:  
%    R = reach_DT(obj,params,options)
%    [R,res] = reach_DT(obj,params,options,spec)
%
% Inputs:
%    obj - nonlinearSysDT object
%    params - parameter defining the reachability problem
%    options - options for the computation of the reachable set
%    spec - object of class specification 
%
% Outputs:
%    R - object of class reachSet storing the reachable set
%    res  - 1 if specifications are satisfied, 0 if not
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: nonlinearSysDT

% Author:       Matthias Althoff, Niklas Kochdumper, Amr Alanwar
% Written:      21-August-2012
% Last update:  29-January-2018
% Last revision:---

%------------- BEGIN CODE --------------

    % options preprocessing
    options = params2options(params,options);
    options = checkOptionsReach(obj,options,0);
    
    spec = [];
    if nargin >= 4
       spec = varargin{1}; 
    end

    % compute symbolic derivatives
    derivatives(obj,options);

    % initialize cell array that stores the reachable sets
    t = options.tStart:obj.dt:options.tFinal;

    steps = length(t)-1;
    R = cell(steps+1,1);
    R_data = cell(steps+1,1);
    R{1} = params.R0;
    R_data{1} = params.R0;
    % loop over all reachablity steps
    for i = 1:steps

        % if a trajectory should be tracked
        if isfield(options,'uTransVec')
            options.uTrans = options.uTransVec(:,i);
        end  
        %reduce
        R{i} = reduce(R{i},'girard',20);
        R_data{i} = reduce(R_data{i},'girard',100);
        % compute next reachable set
        [R{i+1},R_data{i+1}] = linReach_DT(obj,R{i},R_data{i},options);

        if isfield(options,'verbose') && options.verbose 
            disp(t(i));
        end
        
        % check specification
        if ~isempty(spec)
           if ~check(spec,R{i+1})
               timePoint.set = R(2:i+1);
               timePoint.time = num2cell(t(2:i+1)');
               R = reachSet(timePoint);
               return;
           end
        end
    end

    % create reachable set object
    timePoint.set = R(2:end);
    timePoint.time = num2cell(t(2:end)');
    
    timePoint_data.set = R_data(2:end);
    timePoint_data.time = num2cell(t(2:end)');
    
    R = reachSet(timePoint);
    R_data = reachSet(timePoint_data);
end

%------------- END OF CODE --------------