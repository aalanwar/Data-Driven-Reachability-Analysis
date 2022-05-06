% t_nonlinearDT - computes the data driven reachable set of discrete time systems
% x(k+1) = f(x(k),u(k)) + w(k)  
% The approach is based on [1]. The nonlinear system is found in [2]
% 
% 
%
% Syntax:  
%    t_nonlinearDT
%
% Inputs:
%    no
%
% Outputs:
%    no
%
% Example: 
%
% References:
%    [1] Amr Alanwar, Anne Koch, Frank Allgöwer, Karl Johansson 
%       "Data Driven Reachability Analysis Using Matrix Zonotopes"
%    [2] J.M. Bravo, Robust MPC of constrained discrete-time
%        nonlinear systems based on approximated reachable sets, 2006
%    
% 
% Author:       Amr Alanwar
% Written:      29-October-2020
% Last update:  
% Last revision:---


%------------- BEGIN CODE --------------

clear all
close all
%addpath('@nonlinearDT')
rand('seed',1);
dt =0.015;
NN=5;
params.tFinal = dt*NN;

%input set
params.U = zonotope([[0.01;0.01],diag([0.1;0.2])]);  

%initial set
params.R0 = zonotope([[-1.9;-20],diag([0.005;0.3])]);
% dimension of x
options.dim_x=2;

%Number of trajectories
initpoints=30;
%Number of time steps
steps=20;

%Totoal number of samples
totalsamples = steps*initpoints;

%noise zonotope
wfac = 1e-4;
options.W = zonotope(zeros(options.dim_x,1),wfac*ones(options.dim_x,1)); % disturbance

%noise matrix zonotope
for i=1:size(options.W.generators,2)
    vec=options.W.Z(:,i+1);
    for j=0:totalsamples-1
        GW{j+i}= [ zeros(options.dim_x,j),vec,zeros(options.dim_x,totalsamples-j-1)];
    end
end
options.Wmatzono= matZonotope(zeros(options.dim_x,totalsamples),GW);

% Reachability Settings 

options.zonotopeOrder = 100;
options.tensorOrder = 2;
options.errorOrder = 5;


% System Dynamics  
fun = @(x,u) cstrDiscr(x,u,dt);

%input random sample points
for i=1:totalsamples
    u(:,i) = randPointExtreme(params.U);
end

%get state trajectories
x(:,1) = randPoint(params.R0);
index=1;
for j=1:options.dim_x:initpoints*options.dim_x
    x(j:j+options.dim_x-1,1) = randPoint(params.R0);
    x_free(j:j+options.dim_x-1,1) = x(j:j+options.dim_x-1,1);
    for i=1:steps
        x_free(j:j+options.dim_x-1,i+1) = fun(x(j:j+options.dim_x-1,i),u(:,index));
        x(j:j+options.dim_x-1,i+1) = fun(x(j:j+options.dim_x-1,i),u(:,index)) +randPoint(options.W);
        index=index+1;
    end
end


%combine trajectories
index_0 =1;
index_1 =1;
for j=1:options.dim_x:initpoints*options.dim_x
    for i=2:steps+1        
        x_meas_vec_1(:,index_1) = x(j:j+options.dim_x-1,i);
        x_free_vec_1(:,index_1) = x_free(j:j+options.dim_x-1,i);
        index_1 = index_1 +1;
    end
    for i=1:steps
        x_free_vec_0(:,index_0) = x_free(j:j+options.dim_x-1,i);
        x_meas_vec_0(:,index_0) = x(j:j+options.dim_x-1,i);
        index_0 = index_0 +1;
    end
end

stepsLip=1;
initpointsLip=50;
 [gamma,L]= compLipConst(fun,params.U,params.R0,stepsLip,initpointsLip,options.dim_x);

 eps(1)= L(1) .* gamma(1)/2;
 eps(2)= L(2) .* gamma(2)/2;
%options.Zeps = zonotope([zeros(2,1),eps*diag(ones(2,1))]);
options.Zeps = zonotope([zeros(2,1),diag(eps)]);
% flag to add Zeps 0 to skip, 1 to add the Zeps
options.ZepsFlag = 1;

% X_+ is X_1T
% X_- is X_0T
options.U_full = u(:,1:totalsamples);
options.X_0T = x_meas_vec_0(:,1:totalsamples);
options.X_1T = x_meas_vec_1(:,1:totalsamples);

% define system 
sysDisc = nonlinearDT('stirredTankReactor',fun,dt,2,2);


% Reachability Analysis ---------------------------------------------------
% compute model based reachability (R) and data driven one (R_data)
tic
[R ,R_data]= reach_DT(sysDisc,params,options);
tComp = toc;
disp("Computation time: " + tComp);

if options.ZepsFlag
    for i = 1:NN
        R_data.timePoint.set{i} = R_data.timePoint.set{i} + options.Zeps;
    end
end


% Visualization -----------------------------------------------------------

figure; hold on; box on;

% plot initial set
handleX0=plot(params.R0,[1,2],'k-','LineWidth',2);

% plot model based reachable set
handleModel=plot(R,[1 2],'b','Filled',true,'FaceColor',[.8 .8 .8],'EdgeColor','b');

% plot data driven reachable set
handleData=plot(R_data,[1 2],'r','Filled',false);


% formatting
xlabel('x_1');
ylabel('x_2');

% skip warning for extra legend entries
warOrig = warning; warning('off','all');
legend([handleX0,handleModel,handleData],...
    'Initial set $\mathcal{X}_0$','Set from model $\mathcal{R}_k$','Set from data $\mathcal{R}_k^{\prime}$','Location','northwest','Interpreter','latex');
warning(warOrig);
ax = gca;
ax.FontSize = 16;
%set(gcf, 'Position',  [50, 50, 800, 400])
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = outerpos(3) - ti(1) - ti(3)-.01;
ax_height = outerpos(4) - ti(2) - ti(4);
ax.Position = [left bottom ax_width ax_height];
%set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
% example completed
completed = 1;

%------------- END OF CODE --------------