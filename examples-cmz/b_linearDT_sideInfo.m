
% example_sideInfo - example of linear reachability analysis with side
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
rand('seed',1);

clear all
close all
clock
%% system dynamics
dim_x = 5;
A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B_ss = ones(5,1);
C = [1,0,0,0,0];
D = 0;
% define continuous time system
sys_c = ss(A,B_ss,C,D);
% convert to discrete system
samplingtime = 0.05;
sys_d = c2d(sys_c,samplingtime);

% number of trajectories
initpoints =10;
% number of step per trajectory
steps = 3;

% total number of data points
totalsamples = initpoints*steps;
%% initial set and input
X0 = zonotope(ones(dim_x,1),0.1*diag(ones(dim_x,1)));
U = zonotope(10,0.25);

%noise zontope W
W = zonotope(zeros(dim_x,1),0.005*ones(dim_x,1));
GW={};
index=1;
for i=1:size(W.generators,2)
    vec=W.Z(:,i+1);
    GW{index}= [ vec,zeros(dim_x,totalsamples-1)];
    for j=1:totalsamples-1
        GW{j+index}= [GW{index+j-1}(:,2:end) GW{index+j-1}(:,1)];
    end
    index = j+index+1;
end

% concatinate W to get matrix zonotope Wmatzono
Wmatzono= matZonotope(zeros(dim_x,totalsamples),GW);




% randomly choose constant inputs for each step / sampling time
for i=1:totalsamples
    u(i) = randPoint(U);
end


%  simulate the discrete system starting from rand point inside X0
index=1;
for j=1:dim_x:initpoints*dim_x
    x(j:j+dim_x-1,1) = randPoint(X0);
    for i=1:steps
        utraj(j,i) = u(index);
        x(j:j+dim_x-1,i+1) = sys_d.A*x(j:j+dim_x-1,i) + sys_d.B*u(index) + randPoint(W);      
        index=index+1;
    end
end



index_0 =1;
index_1 =1;
for j=1:dim_x:initpoints*dim_x
    for i=2:steps+1
        x_meas_vec_1(:,index_1) = x(j:j+dim_x-1,i);
        index_1 = index_1 +1;
    end
    for i=1:steps
        u_mean_vec_0(:,index_0) = utraj(j,i);
        x_meas_vec_0(:,index_0) = x(j:j+dim_x-1,i);
        index_0 = index_0 +1;
    end
end

%U_-
U_full = u_mean_vec_0(:,1:totalsamples); %same as u 
%X_-
X_0T = x_meas_vec_0(:,1:totalsamples);
%X_+
X_1T = x_meas_vec_1(:,1:totalsamples);


%compute A_Nw and B_Nw
 basis=null([X_0T;U_full]);
for i=1:length(GW)
    Acon{i}=GW{i}*basis;
end
Bcon = (X_1T-zeros(dim_x,totalsamples))*basis  ;
%N_w
Wconmat = conMatZonotope(zeros(dim_x,totalsamples),GW,Acon,Bcon);

% plot simulated trajectory
figure;
subplot(1,2,1); hold on; box on; plot(x(1,:),x(2,:),'b'); xlabel('x_1'); ylabel('x_2');
subplot(1,2,2); hold on; box on; plot(x(3,:),x(4,:),'b'); xlabel('x_3'); ylabel('x_4');
close;


X1W_cen =  X_1T - Wconmat.center;
X1W_cmz = conMatZonotope(X1W_cen,Wconmat.generator,Wconmat.A,Wconmat.B);
X1W = matZonotope(X1W_cen,Wconmat.generator);

%compute N_sigma
AB_cmz = X1W_cmz * pinv([X_0T;U_full]);
%compute M_sigma
AB = X1W  *pinv([X_0T;U_full]);

%double check that the true system dyn are inside AB
 intAB11 = intervalMatrix(AB);
 intAB1 = intAB11.int;
 intAB1.sup >= [sys_d.A,sys_d.B]
 intAB1.inf <= [sys_d.A,sys_d.B]

%specify side informaiton about the A B
 Y = zeros(5,6);
 Q = eye(5);
  rfac = 0.005;
  R = [ 1 1 rfac rfac rfac 1;...
       1 1 rfac rfac rfac 1;...
       rfac rfac 1 1 rfac 1;...
       rfac rfac 1 1 rfac 1;...
       rfac rfac rfac rfac 1 1];


 %list of R (
 index=1;
 for i=1:5
     for j=1:6
        RList{index} =zeros(5,6);
        RList{index}(i,j) = R(i,j);
        index = index+1;
     end
 end

   
 % compute N_s =  <CNΣ, G_Ns, A_Ns, BNs>
 for i =1:length(AB_cmz.generator)
   Aside{i} = [Q*AB_cmz.generator{i}];
   [row_side,col_side]=size(Aside{i});
   [row_A, col_A]= size(AB_cmz.A{i});
   AB_cmz_new_A{i} = [ AB_cmz.A{i};Aside{i},zeros(row_side,col_A-col_side)];
 end
 %add one more beta for A and generators
 AB_cmz_gen_new = AB_cmz.generator;
 for i=1:length(RList)
    AB_cmz_new_A{length(AB_cmz.generator)+i} = [zeros(row_A,col_A);-RList{i},zeros(row_side,col_A-col_side) ];
 AB_cmz_gen_new{length(AB_cmz.generator)+i} = zeros(size(AB_cmz.generator{1}));
 end
 
 Bside =Y- Q * AB_cmz.center ;
 AB_cmz_new_B = [ AB_cmz.B;Bside,zeros(row_side,col_A-col_side)] ;
 
 %N_s
 AB_cmz_side = conMatZonotope(AB_cmz.center,AB_cmz_gen_new,AB_cmz_new_A,AB_cmz_new_B);




%% compute next step sets from model / data

% set number of steps in analysis
totalsteps = 3;
% \mathcal{R}
X_model = cell(totalsteps+1,1);
% \hat{\mathcal{R}}
X_data = cell(totalsteps+1,1);
% \bar{\mathcal{R}}
X_data_cmz = cell(totalsteps+1,1);
% \bar{\mathcal{R}}^s
X_data_cmz_side = cell(totalsteps+1,1);

% init sets for loop
X_model{1} = X0; 
X_data{1} = X0;
X_data_cmz{1} = conZonotope(X0);
X_data_cmz_side{1} = conZonotope(X0);

delFlag = 0;
redFact = 200;
execTimeModel = 0;
execTimeData = 0; 
execTimeCmz = 0;
execTimeSide = 0;
for i=1:totalsteps
    i
    % 1) model-based computation 
    tic()
    X_model{i+1,1} = sys_d.A * X_model{i} + sys_d.B * U+W;
    execTimeModel=execTimeModel + toc()/60;
    X_model{i+1,1}=reduce(X_model{i+1,1},'girard',redFact);
    
    %2) data driven matrix zonotope
    tic()
    X_data{i+1,1} = AB*cartProd(X_data{i},U) +W;
    execTimeData=execTimeData + toc()/60;
    X_data{i+1,1} = reduce(X_data{i+1,1},'girard',redFact);

    
    % 3) Data Driven approach using exact noise description 
    tic()
    cart_cmz = cartProd(X_data_cmz{i},conZonotope(U));
    X_data_cmz{i+1,1} = AB_cmz*cart_cmz +W;
    execTimeCmz=execTimeCmz + toc()/60;
    X_data_cmz{i+1,1}=reduce(X_data_cmz{i+1,1},'girard',redFact); %390 
    
    % 4) Data Driven approach using side information 
    tic()
    X_data_cmz_side{i+1,1} = AB_cmz_side*cartProd(X_data_cmz_side{i},conZonotope(U)) +W;
    execTimeSide=execTimeSide + toc()/60;
    X_data_cmz_side{i+1,1}=reduce(X_data_cmz_side{i+1,1},'girard',redFact);
    
end
execTimeModel=execTimeModel/totalsteps
execTimeData=execTimeData/totalsteps
execTimeCmz=execTimeCmz/totalsteps
execTimeSide=execTimeSide/totalsteps

%% visualization

projectedDims = {[1 2],[3 4],[4 5]};
axx{1} = [0.75,1.5,0.5,4]; axx{2} = [0.75,3,0.8,2.2];axx{3} = [0.75,2.3,0.75,2.8];
index=1;
numberofplots = totalsteps+1;%length(X_model)
for plotRun=1:length(projectedDims)
    plotRun
    figure('Renderer', 'painters', 'Position', [10 10 700 900])
    
    % set axis
      
    index=index+1;
    % plot initial set
    handleX0 = plot(X0,projectedDims{plotRun},'k-','LineWidth',2);
    hold on;
    
    
    % plot reachable sets starting from index 2, since index 1 = X0
    
    % plot reachable sets from model
    for iSet=2:numberofplots
        handleModel=  plot(X_model{iSet},projectedDims{plotRun},'b','Filled',true,'FaceColor',[.8 .8 .8],'EdgeColor','b');
    end
    
    % plot reachable sets from data
    for iSet=2:numberofplots
     handleDatacon   =   plot(X_data_cmz{iSet},projectedDims{plotRun},'r','Template',128);%,'Splits',0);
    end
    
    % plot reachable sets from data
    for iSet=2:numberofplots
     handleDataconside   =   plot(X_data_cmz_side{iSet},projectedDims{plotRun},'r*-','Template',128);%,'Splits',0);
    end
    % plot reachable sets from data
    for iSet=2:numberofplots
      handleDatazono  =   plot(X_data{iSet},projectedDims{plotRun},'k');
    end
    % label plot
    xlabel(['x_{',num2str(projectedDims{plotRun}(1)),'}']);
    ylabel(['x_{',num2str(projectedDims{plotRun}(2)),'}']);
    %axis(axx{plotRun});
    % skip warning for extra legend entries
    warOrig = warning; warning('off','all');
    legend([handleX0,handleModel,handleDatazono,handleDatacon,handleDataconside],...
        'Initial set $\mathcal{X}_0$','Set from model $\mathcal{R}_k$','Set from data $\hat{\mathcal{R}}_k$','Set from data $\bar{\mathcal{R}}_k$','Set from data $\bar{\mathcal{R}}_k^{s}$','Interpreter','latex');
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
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
end




%------------- END OF CODE --------------