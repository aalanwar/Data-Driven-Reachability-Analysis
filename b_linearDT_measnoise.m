% t_linear_measnoise - computes the data driven reachable set of discrete
% time systems with modeling and measurement noise
% x(k+1) = Ax(k) + Bu(k) + w(k)
% \tilde(x)(k) = x(k) + v(k)
%
% This example can be found in [1].
%
% Syntax:  
%    example_measnoise
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
%% system dynamics
dim_x = 5;
A = [-1 -4 0 0 0; 4 -1 0 0 0; 0 0 -3 1 0; 0 0 -1 -3 0; 0 0 0 0 -2];
B_ss = ones(dim_x,1);
C = [1,0,0,0,0];
D = 0;
% define continuous time system
sys_c = ss(A,B_ss,C,D);
% convert to discrete system
samplingtime = 0.05;
sys_d = c2d(sys_c,samplingtime);
% Number of trajectories
initpoints =10;
% Number of time steps
steps = 2;
totalsamples = initpoints*steps;
% initial set and input
X0 = zonotope(ones(dim_x,1),0.1*diag(ones(dim_x,1)));
U = zonotope(10,0.25);

%noise zontope W
W = zonotope(zeros(dim_x,1),0.005*ones(dim_x,1));
%Construct matrix zonotpe \mathcal{M}_w
for i=1:size(W.generators,2)
    vec=W.Z(:,i+1);
     GW{i}= [ vec,zeros(dim_x,totalsamples-1)];
    for j=1:totalsamples-1
        GW{j+i}= [GW{i+j-1}(:,2:end) GW{i+j-1}(:,1)];
    end
end
Wmatzono= matZonotope(zeros(dim_x,totalsamples),GW);

%measument noise v
V = zonotope(zeros(dim_x,1),0.002*ones(dim_x,1));
%take care to change the center of Vmatzono if you change V
%Construct matrix zonotpe \mathcal{M}_v
index=1;
for i=1:size(V.generators,2)
    vec=V.Z(:,i+1);
    GV{index}= [ vec,zeros(dim_x,totalsamples-1)];
    for j=1:totalsamples-1
        GV{j+index}= [GV{index+j-1}(:,2:end) GV{index+j-1}(:,1)];
    end
    index = j+index+1;
end
Vmatzono= matZonotope(zeros(dim_x,totalsamples),GV);

%Construct matrix zonotpe \mathcal{M}_Av
AVmatzono = sys_d.A * Vmatzono;


% randomly choose constant inputs for each step / sampling time
for i=1:totalsamples
    u(i) = randPoint(U);
end


% simulate the discrete system starting from random initial points
x(:,1) = randPoint(X0);
index=1;
for j=1:dim_x:initpoints*dim_x
    x(j:j+dim_x-1,1) = randPoint(X0);
    x_v(j:j+dim_x-1,1) = x(j:j+dim_x-1,1);
    for i=1:steps
        utraj(j,i) = u(index);
        x(j:j+dim_x-1,i+1) = sys_d.A*x(j:j+dim_x-1,i) + sys_d.B*u(index) + randPoint(W);  
        x_v(j:j+dim_x-1,i+1) =  x(j:j+dim_x-1,i+1) + randPoint(V);
        index=index+1;
    end
end


% concatenate the data trajectories 
index_0 =1;
index_1 =1;
for j=1:dim_x:initpoints*dim_x
    for i=2:steps+1
        x_meas_vec_1_v(:,index_1) = x_v(j:j+dim_x-1,i);
        index_1 = index_1 +1;
    end
    for i=1:steps
        u_mean_vec_0(:,index_0) = utraj(j,i);
        x_meas_vec_0_v(:,index_0) = x_v(j:j+dim_x-1,i);
        x_meas_vec_0(:,index_0) = x(j:j+dim_x-1,i);
        index_0 = index_0 +1;
    end
end

% X_+ is X_1T
% X_- is X_0T
% U_- is U_full
U_full = u_mean_vec_0(:,1:totalsamples); %same as u 
X_0T = x_meas_vec_0_v(:,1:totalsamples);
X_1T = x_meas_vec_1_v(:,1:totalsamples);
X_0T_pure = x_meas_vec_0(:,1:totalsamples);

mink= -1*Wmatzono +  -1*Vmatzono + AVmatzono;
%compute A_tildeNsigma and B_tildeNsigma
basis=null([X_0T;U_full]);
for i=1:length(mink.generator)
    Acon{i}=(mink.generator{i})*basis;
end
Bcon = (X_1T-zeros(dim_x,totalsamples))*basis  ;

%% AB with AV assumption (tildeMsigma)
AB_av = (X_1T + -1*Wmatzono + -1* Vmatzono + AVmatzono)*pinv([X_0T;U_full]);


% constrained matrix zonotope AB with AV assumption (tildeNsigma)
AB_cmz = conMatZonotope(AB_av.center,AB_av.generator,Acon,Bcon);


% plot simulated trajectory
figure;
subplot(1,2,1); hold on; box on; plot(x(1,:),x(2,:),'b'); xlabel('x_1'); ylabel('x_2');
subplot(1,2,2); hold on; box on; plot(x(3,:),x(4,:),'b'); xlabel('x_3'); ylabel('x_4');
close;



%% AB wihtout AV assumption (Algorithm 4 LTIMeasReachability) 
X1W_cen =  X_1T - Wmatzono.center- Vmatzono.center;
AB = X1W_cen  *pinv([X_0T;U_full]);
minTerm = min((X_1T - AB * ([X_0T;U_full]))')';
maxTerm = max((X_1T - AB * ([X_0T;U_full]))')';
AV_oneterm = zonotope(interval(minTerm,maxTerm))+ -1*W + -1*V;
%%%%%%%%%%%%%%%%%%%%
%%
%checking if true A B are within the set
AV_minus_matzono = X_1T - AB * ([X_0T;U_full])+ -1*Wmatzono + -1*Vmatzono;
AB_f= (AB * ([X_0T;U_full]) + AV_minus_matzono)*pinv([X_0T_pure;U_full]);
intAB11 = intervalMatrix(AB_f);
intAB1 = intAB11.int;
intAB1.sup >= [sys_d.A,sys_d.B]
intAB1.inf <= [sys_d.A,sys_d.B]
VInt = intervalMatrix(AV_minus_matzono);
leftLimit = VInt.Inf;
rightLimit = VInt.Sup;
AV_minus= zonotope(interval(min(leftLimit')',max(rightLimit')'));
figure;plot(AV_minus,[ 1 2],'r');hold on;plot(AV_oneterm,[ 1 2],'k');
close;

%%

%% compute next step sets from model / data

% set number of steps in analysis
totalsteps = 3;
%model based reachability 
X_model = cell(totalsteps+1,1);
X_model{1} = X0; 

%Algorithm 4
X_data = cell(totalsteps+1,1);
X_data_tilde = cell(totalsteps+1,1);
X_data{1} = X0; 
X_data_tilde{1} = X0;



%using matrix zonotope
X_data_av = cell(totalsteps+1,1);
X_data_av{1} = X0;

%reach using constrained matrix zonotope 
X_data_cmz = cell(totalsteps+1,1);
X_data_cmz{1} = conZonotope(X0);

tic
redOrder = 390;
for i=1:totalsteps
    i
    % 1) model-based computation
    
    X_model{i+1,1} = sys_d.A * X_model{i} + sys_d.B * U+W;
    X_model{i+1,1}=reduce(X_model{i+1,1},'girard',redOrder);

    % 2) Data Driven approach with Av assumption
    X_data_av{i+1,1} = AB_av * (cartProd(X_data_av{i},U)) +W;
    X_data_av{i+1,1}=reduce(X_data_av{i+1,1},'girard',redOrder);
    
    
    % 3) Data Driven approach cmz with Av assumption
    cart_cmz = cartProd(X_data_cmz{i},conZonotope(U));
    %cart_cmz = cartProd(zonotope(X_data_cmz{i}.Z),U);
    X_data_cmz{i+1,1} = AB_cmz*cart_cmz +W;
    X_data_cmz{i+1,1}=reduce(X_data_cmz{i+1,1},'girard',redOrder);  
    
    % 4) Data Driven approach without Av assumption (algorithm 4)
    if i==1
        X_data{i+1,1} = AB * (cartProd(X_data{i},U))+AV_oneterm +W;
        X_data_tilde{i+1} = X_data{i+1}+V;
    else
        X_data{i+1,1} = AB * (cartProd(X_data{i}+ V,U))+AV_oneterm +W;
        %or
        %X_data_tilde{i+1,1} = AB * (cartProd(X_data_tilde{i},U))+AV_minus +W +V;
        %X_data{i+1,1}=X_data_tilde{i+1,1}+-1*V;
    end
    X_data{i+1,1}=reduce(X_data{i+1,1},'girard',redOrder);
end
exec_time = toc/60



%% visualization

projectedDims = {[1 2],[3 4],[4 5]};
axx{1} = [0.75,1.5,0.5,4]; axx{2} = [0.75,3,0.8,2.2];axx{3} = [0.75,2.3,0.75,2.8];
index=1;
numberofplots = totalsteps+1;%
for plotRun=1:length(projectedDims)
    plotRun
    figure('Renderer', 'painters', 'Position', [10 10 700 900])
    
      
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
        handleData=   plot(X_data{iSet},projectedDims{plotRun},'k-+');
    end
    
    % plot reachable sets from data
    for iSet=2:numberofplots
        handleData_av=   plot(X_data_av{iSet},projectedDims{plotRun},'k');
    end
    
        % plot reachable sets from data
    for iSet=2:numberofplots
        handleData_cmz=   plot(X_data_cmz{iSet},projectedDims{plotRun},'r','Template',128);
    end
    % label plot
    xlabel(['x_{',num2str(projectedDims{plotRun}(1)),'}']);
    ylabel(['x_{',num2str(projectedDims{plotRun}(2)),'}']);
    %axis(axx{plotRun});
    % skip warning for extra legend entries
    warOrig = warning; warning('off','all');
    legend([handleX0,handleModel,handleData,handleData_av,handleData_cmz],...
        'Initial Set','Set from Model','Set from Data','Set From Data (Av-Zono)','Set From Data (Av-ConZono)','Location','northwest');
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