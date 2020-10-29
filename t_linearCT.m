% t_linearCT - computes the data driven reachable set of continuous  time
% systems using x dot and x data samples in discrete time
%
%
% Example:
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: Amr Alanwar, Anne Koch, Frank Allgöwer, Karl Johansson "Data Driven Reachability Analysis from Noisy Data with Unknown
% System Model"
%
%
%
% Author:       Amr Alanwar
% Written:      28-October-2020
% Last update:
% Last revision:---

%------------- BEGIN CODE --------------

close all
clear all
% set model parameters ----------------------------------------------------
%time step size for reachable set computation
options.timeStep = 0.02;
% Number of steps
nrOfSteps = 5;
%final time
params.tFinal = options.timeStep * nrOfSteps;
%dimension of the state
dim_x=5;
% initial set X_0
params.R0 = zonotope([ones(dim_x,1),0.1*diag(ones(dim_x,1))]);
%input set
U = zonotope([1,0.1]);
% noise zonotope
W = zonotope(zeros(dim_x,1),0.0001*ones(dim_x,1));
%Number of data samples
totalsamples = 5000;
% concatenate W to get matrix zonotope Wmatzono
for i=1:size(W.generators,2)
    vec=W.Z(:,i+1);
    for j=0:totalsamples-1
        GW{j+i}= [ zeros(dim_x,j),vec,zeros(dim_x,totalsamples-j-1)];
    end
end
Wmatzono= matZonotope(zeros(dim_x,totalsamples),GW);

%noise zonotope for xdot data
G = zonotope(zeros(dim_x,1),0.025*ones(dim_x,1)); % disturbance
for i=1:size(G.generators,2)
    vec=G.Z(:,i+1);
    for j=0:totalsamples-1
        GG{j+i}= [ zeros(dim_x,j),vec,zeros(dim_x,totalsamples-j-1)];
    end
end
Gmatzono= matZonotope(zeros(dim_x,totalsamples),GG);


% set reachability options ------------------------------------------------
options.taylorTerms = 5; %number of taylor terms for reachable sets
options.zonotopeOrder = 50; %zonotope order
options.samplingtime = options.timeStep;
options.intermediateOrder = 2;
% specify continuous dynamics----------------------------------------------
A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B = ones(dim_x,1);

%--------------------------------------------------------------------------

%% system dynamics

C = [1,zeros(1,dim_x-1)];
D = 0;
% define continuous time system
sys_c = ss(A,B,C,D);
% convert to discrete system
sys_d = c2d(sys_c,options.samplingtime);



%% simulate the discrete system
%getting random inputs
for i=1:totalsamples+1
    u(i) = randPoint(U);%pinv(B_ss)*randPoint(params.U);
end

x0 = ones(dim_x,1);
x2(:,1) = x0;
xdot(:,1) = zeros(dim_x,1);
for i=1:totalsamples+1
    x2(:,i+1) = sys_d.A*x2(:,i) + sys_d.B*u(i)+ randPoint(W);
end
% using lsim is the same as using discrete A and B, but you can't process
% noise W to lsim.

x0=ones(dim_x,1);
t= 0:options.timeStep:(totalsamples)*options.timeStep;
[y,t,x] = lsim(sys_c,u,t,x0);
x=x';
for i=1:totalsamples+1
    xdot(:,i) = A * x(:,i) + B*u(i)+ randPoint(G);
    % measument noise W also
    x_noisy(:,i) = x2(:,i) + randPoint(W);
end


U_full = u(1:totalsamples);
X_0T = x_noisy(:,1:totalsamples);
X_1T = x_noisy(:,2:totalsamples+1);
Xdot_1T = xdot(:,1:totalsamples);


%fit approx model
AB_con = (Xdot_1T  - Wmatzono.center - Gmatzono.center)*pinv([X_0T;U_full]);
%find upper and lower bound on AV and model mismatch
AW = Xdot_1T + -1* AB_con*[X_0T;U_full] + -1* Gmatzono   + -1* Wmatzono;

AWInt = intervalMatrix(AW);
leftLimit = AWInt.Inf;
rightLimit = AWInt.Sup;

AW_one= zonotope(interval(min(leftLimit')',max(rightLimit')'));

A_con = AB_con(:,1:dim_x);
B_con = AB_con(:,end);
fiveDimSys_zono2 = linearSys('fiveDimSys',A_con,1);; %instantiate system

params_data = params;
%integrate the the AW with the inputs
params_data.U =   B_con*U +W  + AW_one + G;
params_data.R0 = params.R0 + G;
Rdata = reach(fiveDimSys_zono2, params_data, options);


fiveDimSys=linearSys('fiveDimSys',A,1); %initialize system
params.U =   B*U +W;
Rcont = reach(fiveDimSys, params, options);
for i=1:length(Rdata.timeInterval.set)
    % to get the x from tilde(x)
    Rdata.timeInterval.set{i} = Rdata.timeInterval.set{i} + -1* G;
end


% plot results-------------------------------------------------------------
for plotRun=1:3
    % plot different projections
    if plotRun==1
        projectedDims=[1 2];
    elseif plotRun==2
        projectedDims=[3 4];
    elseif plotRun==3
        projectedDims=[4 5];
    end
    
    
    % plot reachable sets
    %subplot(1,3,plotRun); hold on; box on;
    figure('Renderer', 'painters', 'Position', [10 10 700 900])
    % plot initial set
    handleX0= plot(params.R0,projectedDims,'k-','lineWidth',2);
    
    hold on
    %model based
    handleModel = plot(Rcont,projectedDims,'b','Filled',true,'FaceColor',[.8 .8 .8],'EdgeColor','b');
    
    
    % data-driven
    handleData = plot(Rdata,projectedDims,'r','Filled',false);
    
    
    handleX0= plot(params.R0,projectedDims,'k-','lineWidth',2);
    
    
    % label plot
    xlabel(['x_{',num2str(projectedDims(1)),'}']);
    ylabel(['x_{',num2str(projectedDims(2)),'}']);
    
    % skip warning for extra legend entries
    warOrig = warning; warning('off','all');
    legend([handleX0,handleModel,handleData],...
        'Initial Set','Set from Model','Set from Data','Location','northwest');
    warning(warOrig);
    ax = gca;
    ax.FontSize = 16;
    %set(gcf, 'Position',  [50, 50, 800, 400])
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    
end
