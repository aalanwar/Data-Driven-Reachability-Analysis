% t_linearDT - computes the data driven reachable set of discrete time systems
% x(k+1) = Ax(k) + Bu(k) 
%
% Example: 
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: 
% [1] Amr Alanwar, Mark Wetzlinger , Matthias Althoff, Florian Dörfler, Karl Johansson 
% "Data-Enabled Reachability Analysis using the Fundamental Lemma"
%
%
%
% Author:       Amr Alanwar, Mark Wetzlinger
% Written:      28-October-2020
% Last update:  
% Last revision:---

%------------- BEGIN CODE --------------


res = 1;

close all; clear all;
rand('seed',1);


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


%% initial set and input
X0 = zonotope(ones(dim_x,1),0.1*diag(ones(dim_x,1)));
U = zonotope(1,0.25);


%% simulate the discrete system 
totalsamples = 520;
% random sample input
for i=1:totalsamples
    u(i) = randPoint(U);
end

x(:,1) = randPoint(X0);
for i=1:totalsamples
	x(:,i+1) = sys_d.A*x(:,i) + sys_d.B*u(i);
end

% plot simulated trajectory
figure;
subplot(1,2,1); hold on; box on; plot(x(1,:),x(2,:),'b'); xlabel('x_1'); ylabel('x_2');
subplot(1,2,2); hold on; box on; plot(x(3,:),x(4,:),'b'); xlabel('x_3'); ylabel('x_4');
close;


%%  get function x_k+1 = G(x,u) x_k|u_k
numofsamples = totalsamples - 1;
t = 1;
T = numofsamples - 1;
U_full = u(1:T);
X_0T = x(:,1:T);
X_1T = x(:,2:T+1);
Gxu = X_1T * pinv([U_full;X_0T]);



%% compute next step sets from model / data

% set number of steps in analysis
totalsteps = 50;
X_model = cell(totalsteps+1,1); X_data = cell(totalsteps+1,1);
% init sets for loop
X_model{1} = X0; X_data{1} = X0;

for i=1:totalsteps

    % 1) model-based computation
    X_model{i+1,1} = sys_d.A * X_model{i} + sys_d.B * U;
    
    % 2) data-driven computation
    X_data{i+1,1} = Gxu * cartProd(U,X_data{i,1});
    
end



%% visualization

projectedDims = {[1 2],[3 4],[4 5]};
axx{1} = [-1.5,1.5,-1,2]; axx{2} = [0,1.2,0,1.2];axx{3} = [0,2,0,2];

numberofplots = 15;%length(X_model)
for plotRun=1:length(projectedDims)
    
    figure('Renderer', 'painters', 'Position', [10 10 700 900])
    % choose subplot
   % plot(1,length(projectedDims),plotRun); hold on; box on;
 
    % set axis
    axis(axx{plotRun});
    
    % plot initial set
    handleX0 = plot(X0,projectedDims{plotRun},'k-','LineWidth',2);
    hold on; 
    % plot reachable sets starting from index 2, since index 1 = X0
    
    % plot reachable sets from model
    for iSet=2:numberofplots
        %handleModel = plotFilled(X_model{iSet},projectedDims{plotRun},...
         %   [.8 .8 .8],'EdgeColor','b');
        handleModel = plot(X_model{iSet},projectedDims{plotRun},'b','Filled',true,'FaceColor',[.8 .8 .8],'EdgeColor','b');
    end
    
    % plot reachable sets from data
    for iSet=2:numberofplots
        handleData = plot(X_data{iSet},projectedDims{plotRun},'r');
    end

    % label plot
    xlabel(['x_{',num2str(projectedDims{plotRun}(1)),'}']);
    ylabel(['x_{',num2str(projectedDims{plotRun}(2)),'}']);
    
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
    %set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
end



%------------- END OF CODE --------------