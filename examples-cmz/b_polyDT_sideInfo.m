% example_poly_sideInfo - example of reachability analysis for polynomial systems with side
% information. This example generates the reachable regoins using matrix
% zonotope, constrained matrix zonotope with exact noise description, and
% side information.
%
% This example can be found in [1].
%
% Syntax:
%    example_sideInfo
%
% Inputs:
%    no
%
% Outputs:
%    no
%
% References:
%    [1] Amr Alanwar, Anne Koch, Frank Allgöwer,Karl Henrik Johansson
%   "Data-Driven Reachability Analysis from Noisy Data"
% Author:       Amr Alanwar
% Written:      17-Jan-2021
% Last update:  21-May-2021
% Last revision:---

%------------- BEGIN CODE --------------

clear all
close all
rand('seed',50);
%time steps
N=3;
dt = 0.015;
%final time for model based reachability
params.tFinal = dt*N;

%input set
params.U = zonotope([[0.2;0.3],diag([0.01,0.02])]); 
%initial set
params.R0 = zonotope([[1;2],diag([0.05;0.3])]);


% dimension of x
options.dim_x=2;

%Number of trajectories
initpoints=25;
%Number of time steps
steps=7;

%Totoal number of samples
totalsamples = steps*initpoints;

%noise zonotope
wfac = 0.7e-4;
options.W = zonotope(zeros(options.dim_x,1),wfac*ones(options.dim_x,1)); % disturbance

%noise matrix zonotope
for i=1:size(options.W.generators,2)
    vec=options.W.Z(:,i+1);
    for j=0:totalsamples-1
        GW{j+i}= [ zeros(options.dim_x,j),vec,zeros(options.dim_x,totalsamples-j-1)];
    end
end
Wmatzono= matZonotope(zeros(options.dim_x,totalsamples),GW);

% Reachability Settings
options.zonotopeOrder = 200;
options.tensorOrder = 2;
options.errorOrder = 5;


% True system dynamics
fun = @(x,u) polyFun(x,u,dt);

%input random sample points
for i=1:totalsamples
    u(:,i) = randPoint(params.U);
end

%get state trajectories
x(:,1) = randPoint(params.R0);
index=1;
for j=1:options.dim_x:initpoints*options.dim_x
    x(j:j+options.dim_x-1,1) = randPoint(params.R0);
    for i=1:steps
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
        index_1 = index_1 +1;
    end
    for i=1:steps
        x_meas_vec_0(:,index_0) = x(j:j+options.dim_x-1,i);
        index_0 = index_0 +1;
    end
end




% X_+ is X_1T
% X_- is X_0T
% U_- is U_full
X_1T = x_meas_vec_1(:,1:totalsamples);
U_full = u(:,1:totalsamples);
X_0T = x_meas_vec_0(:,1:totalsamples);

%compute monomials of trajectories
X_2 = X_0T .*X_0T;
X1X2 = X_0T(1,:) .* X_0T(2,:);
U_2 = U_full .*U_full;
U1U2 = U_full(1,:) .* U_full(2,:);
XU = X_0T .*U_full;
X1U2X2U1 = X_0T .* U_full([2 1],:);
data=[ones(1,totalsamples);X_0T;X_2;X1X2;U_full;U_2;U1U2;XU;X1U2X2U1];


%%%%%%%%%%%%%%%
dim_x =2;
nn=null(data);
basis = nn;%*nn';%*pinv(nn);%null([X_0T;U_full])*pinv(null([X_0T;U_full]));
for i=1:length(GW)
    Acon{i}=GW{i}*basis;
end
Bcon = (X_1T-zeros(dim_x,totalsamples))*basis;

Wconmat = conMatZonotope(zeros(dim_x,totalsamples),GW,Acon,Bcon);
%%
RANK = rank(data)
% set of A B that is consistent with the data
AB = (X_1T + -1* Wmatzono)*pinv(data);
%[    ones  ; X_0T   ;X_2    ;X1X2;U_full;U_2 ;U1U2;XU ;X1U2X2U1]);
truedyn=[ 0,0.7 0   ,0.32 0 ,  0 , 1 0  ,0 0 , 0  ,0 0,  0 0; ...
    0,0.09 0  ,0  0.4 ,  0 , 0 0  ,0 0 , 0  ,0 0   , .320 0];
% check if the true dynamics is within AB
intAB11 = intervalMatrix(AB);
intAB1 = intAB11.int;
intAB1.sup >= truedyn
intAB1.inf <= truedyn

% constrain matrix zonotpe that is consitent with data and exact noise
% description
X1W_cen =  X_1T - Wconmat.center;
X1W_cmz = conMatZonotope(X1W_cen,Wconmat.generator,Wconmat.A,Wconmat.B);
AB_cmz = X1W_cmz * pinv(data);

%side information
[r_c c_c]= size(truedyn);
Q = eye(r_c);
Y = zeros(r_c,c_c);
rfac = 0.1;
     %[ one  ; X_0T ;X_2    ;X1X2   ;U_full      ;U_2       ;U1U2   ;XU       ;X1U2X2U1]);
Rmat =[ rfac,1 rfac ,1 rfac ,  rfac , 1.2 rfac   ,rfac rfac , rfac  ,rfac rfac   , rfac  rfac; ...
        rfac,1  rfac,rfac 1 ,  rfac , rfac rfac  ,rfac rfac , rfac  ,rfac rfac   ,1  rfac];
index=1;
for i=1:r_c
    for j=1:c_c
        %  if  R(i,j) ~=rfac
        RList{index} =zeros(r_c,c_c);
        RList{index}(i,j) = Rmat(i,j);
        index = index+1;
        %  end
    end
end


%add constraint from side information 
for i =1:length(AB_cmz.generator)
    Aside{i} = [Q*AB_cmz.generator{i}];
    [row_side,col_side]=size(Aside{i});
    [row_A, col_A]= size(AB_cmz.A{i});
    AB_cmz_new_A{i} = [ AB_cmz.A{i};Aside{i},zeros(row_side,col_A-col_side)];
end
AB_cmz_gen_new = AB_cmz.generator;
for i=1:length(RList)
    AB_cmz_new_A{length(AB_cmz.generator)+i} = [zeros(row_A,col_A);-RList{i},zeros(row_side,col_A-col_side) ];
    AB_cmz_gen_new{length(AB_cmz.generator)+i} = zeros(size(AB_cmz.generator{1}));
end

Bside =Y- Q * AB_cmz.center ;
AB_cmz_new_B = [ AB_cmz.B;Bside,zeros(row_side,col_A-col_side)] ;

% constrain matrix zonotpe that is consitent with data, exact noise
% description and side information
AB_cmz_side = conMatZonotope(AB_cmz.center,AB_cmz_gen_new,AB_cmz_new_A,AB_cmz_new_B);
%%
% Model based reachability
sysDisc = nonlinearSysDT('stirredTankReactor',fun,0.015,2,2);
%reachable set using matrix zonotope
X_data = cell(N+1,1);
%reachable set using constrained matrix zonotope and exact noise
%description 
X_data_cmz = cell(N+1,1);
%reachable set using constrained matrix zonotope, exact noise
%description and side information
X_data_cmz_side = cell(N+1,1);
X_data{1} = params.R0;
X_data_cmz{1} = conZonotope(params.R0);
X_data_cmz_side{1} = conZonotope(params.R0);
execTimeModel = 0;
execTimeData = 0; 
execTimeCmz = 0;
execTimeSide = 0;
for i = 1:N
    i
    %compute monomial of each reachable set
    X_z1 = interval(X_data{i});
    U_int = interval(params.U);    
    tic()
    card =zonotope([interval([1]);X_z1;X_z1.*X_z1;X_z1(1)*X_z1(2);U_int;U_int.*U_int;U_int(1)*U_int(2);X_z1.*U_int;X_z1(1)*U_int(2);X_z1(2)*U_int(1)]);
    X_data{i+1} = AB *card +options.W ;
    execTimeData=execTimeData + toc()/60;
    X_data{i+1,1}=reduce(X_data{i+1,1},'girard',options.zonotopeOrder); 
    
    X_z1 = interval(X_data_cmz{i});
    tic()
    card_cmz =zonotope([interval([1]);X_z1;X_z1.*X_z1;X_z1(1)*X_z1(2);U_int;U_int.*U_int;U_int(1)*U_int(2);X_z1.*U_int;X_z1(1)*U_int(2);X_z1(2)*U_int(1)]);
    X_data_cmz{i+1} = AB_cmz *card_cmz+options.W ;
    execTimeCmz=execTimeCmz + toc()/60;
    X_data_cmz{i+1,1}=reduce(X_data_cmz{i+1,1},'girard',options.zonotopeOrder);

    X_z1 = interval(X_data_cmz_side{i});
    tic()
    card_cmz_side=zonotope([interval([1]);X_z1;X_z1.*X_z1;X_z1(1)*X_z1(2);U_int;U_int.*U_int;U_int(1)*U_int(2);X_z1.*U_int;X_z1(1)*U_int(2);X_z1(2)*U_int(1)]);
    X_data_cmz_side{i+1} = AB_cmz_side *card_cmz_side+options.W ;
    execTimeSide=execTimeSide + toc()/60;
    X_data_cmz_side{i+1}=reduce(X_data_cmz_side{i+1},'girard',options.zonotopeOrder);
end
execTimeData=execTimeData/N
execTimeCmz=execTimeCmz/N
execTimeSide=execTimeSide/N
% Reachability Analysis ---------------------------------------------------

tic
R = reach(sysDisc,params,options);
tComp = toc/60;
disp("Computation time: " + tComp);

for i=1:length(R.timePoint.set)
    R.timePoint.set{i} = R.timePoint.set{i} +options.W ;
end


% Simulation --------------------------------------------------------------

simOpt.points = 500;
simOpt.fracVert = 0.5;
simOpt.fracInpVert = 0.5;
simOpt.inpChanges = 3;

simRes = simulateRandom(sysDisc, params, simOpt);


% Visualization -----------------------------------------------------------



% plot initial set
figure('Renderer', 'painters', 'Position', [10 10 700 900])
hold on; box on;
handleX0 = plot(params.R0,[1,2],'k-','LineWidth',2);


pts =[];
for i =1:length(simRes.x)
    for j = 1:length(simRes.x{i})
        noisept = randPoint(options.W);
        pts = [ pts ; simRes.x{i}(j,:) + noisept'];
    end
end
handelpts = plot(pts(:,1),pts(:,2),'.k');
% formatting
xlabel('x_1');
ylabel('x_2');
%
projectedDims = {[1 2]};%,[3 4],[4 5]};
axx{1} = [0.75,1.5,0.5,4]; axx{2} = [0.75,3,0.8,2.2];axx{3} = [0.75,2.3,0.75,2.8];
index=1;
numberofplots = N+1;%length(X_model)
for plotRun=1:length(projectedDims)
    
    plotRun
    
    % plot reachable sets starting from index 2, since index 1 = X0
    for iSet=2:numberofplots
        handleModel=  plot( R.timePoint.set{iSet-1},projectedDims{plotRun},'b');
    end
    
    for iSet=2:numberofplots
        handleDatazono  =  plot(X_data{iSet},projectedDims{plotRun},'k');
    end
    for iSet=2:numberofplots
        handleDataconside   =  plot(X_data_cmz_side{iSet},projectedDims{plotRun},'r*-','Template',200);
    end
    for iSet=2:numberofplots
        handleDatacon   =  plot(X_data_cmz{iSet},projectedDims{plotRun},'r','Template',200);
    end
    
    
    xlabel(['x_{',num2str(projectedDims{plotRun}(1)),'}']);
    ylabel(['x_{',num2str(projectedDims{plotRun}(2)),'}']);
    %axis(axx{plotRun});
    % skip warning for extra legend entries
    warOrig = warning; warning('off','all');
    legend([handleX0,handleModel,handelpts,handleDatazono,handleDatacon,handleDataconside],...
        'Initial Set','Set from model','Simulation points','Set from data (Zono)','Set from data (ConZono)','Set from data (Side-info)','Location','northwest');
    warning(warOrig);
    ax = gca;
    ax.FontSize = 20;
    %set(gcf, 'Position',  [50, 50, 800, 400])
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3)-0.01;
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
end

% example completed
completed = 1;

%------------- END OF CODE --------------