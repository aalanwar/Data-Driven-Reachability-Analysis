function res = o_linearDT_noise

% t_linearDT - computes the data driven reachable set of discrete time systems
% x(k+1) = Ax(k) + Bu(k) 
% z(k) = x(k) + w(k)
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
% [2] De Persis, Tesi: Formulas for Data-Driven Control: Stabilization,
%     Optimality, and Robustness, IEEE Transaction on Automatic Control, 2020.
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


%% set definitions
X0 = zonotope(ones(dim_x,1),0.1*diag(ones(dim_x,1))); % initial set
U = zonotope(1,0.25); % input set
W = zonotope(zeros(dim_x,1),0.005*diag(ones(dim_x,1))); % disturbance


%% simulate the discrete system
% requirement (acc. to [1], (6)):
% T >= (m+1)*n+m ... where n = dim_x, m = dim_u
steps = 100;
numofinit =100;
T = steps;
totalsamples = steps*numofinit;
% randomly choose constant inputs for each step / sampling time
for i=1:totalsamples
    u(i) = randPoint(U);
end
index =1;
for j=1:5:numofinit*5
    %x0 = rand(5,1);
    x0 = randPointExtreme(X0);
    x_noiseless(j:j+4,1) = x0; % without noise
    x_measured(j:j+4,1) = x0; % without noise
    for i=1:steps
        x_noiseless(j:j+4,i+1) = sys_d.A*x_noiseless(j:j+4,i) + sys_d.B*u(index); % without noise
        index = index+1;
        x_measured(j:j+4,i+1) = x_noiseless(j:j+4,i+1) + randPoint(W);
    end
end

% x_meas_vec = x_measured(1:5,:);
index =1;
for j=1:5:numofinit*5
    for i=1:steps
         x_meas_vec(:,index) = x_measured(j:j+4,i);
         index = index +1;
    end
end

% concatenate the data trajectories 
index_0 =1;
index_1 =1;
for j=1:dim_x:numofinit*dim_x
    for i=2:steps+1
        x_meas_vec_1(:,index_1) = x_measured(j:j+dim_x-1,i);
        index_1 = index_1 +1;
    end
    for i=1:steps
        x_meas_vec_0(:,index_0) = x_measured(j:j+dim_x-1,i);
        index_0 = index_0 +1;
    end
end

%% get function x_k+1 = f(x,u) x_k|u_k
U_01T = u(:,1:totalsamples);  % U_0,1,T-L+1   [1] page 3, column 1
X_0T = x_meas_vec_0(:,1:totalsamples);             % X_0,T-L+1     [1] page 2, column 2
X_1T = x_meas_vec_1(:,1:totalsamples);               % X_1,T-L+1     [1] page 3, column 2
fxu = X_1T * pinv([U_01T;X_0T]);        %               [1] page 3, (7)



%% compute next step sets from model / data

% set number of steps in analysis
totalsteps = steps;
%X_model_true = cell(totalsteps+1,1);
X_model_noise = cell(totalsteps+1,1);
X_data = cell(totalsteps+1,1);
% init sets for loop
%X_model_true{1}  = X0;
X_model_noise{1} = X0;
X_data{1} = X0;

for i=1:steps
    % 1) model-based computation 
    %X_model_true{i+1,1} = sys_d.A * X_model_true{i} + sys_d.B * U;
    X_model_noise{i+1,1} = sys_d.A * X_model_noise{i} + sys_d.B * U ;
    
    % 2) data-driven computation
    X_data{i+1,1} = fxu * cartProd(U,X_data{i});
    
end

for i=1:steps
    X_data{i} = reduce(X_data{i},'girard',2);
end


% obtain W by comparing measurement points with computed set
A = [];
b = [];
Aeq = [];
beq = [];

noisept =[];
index = 1;
for j=1:5:numofinit*5
    for i=1:steps
        %check if the point is inside the zonotope
        if ~containsPoint(X_data{i},x_measured(j:j+4,i))
            %if not compute the distance 
            point = x_measured(j:j+4,i);
            x0 = zeros(size(X_data{i}.generators,2),1);
            %alpha from -1 to 1
            lb=-ones(size(X_data{i}.generators,2),1);
            ub=ones(size(X_data{i}.generators,2),1);
            nonlcon = @circlecon;
            alpha = fmincon(@fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);
            noisept(:,index) = point - X_data{i}.center - X_data{i}.generators*alpha;
            index = index + 1;
        end
    end
end
% to make sure that we do not shift by adding a zonotope
noisept(:,index) = zeros(5,1);
if(length(noisept)>1)
     noisezono = zonotope.enclosePoints(noisept);
    X_data{1} = X0;
    for i=1:steps
        X_data{i} = X_data{i} + noisezono;
    end
end

% plotting ----------------------------------------------------------------
projectedDims = {[1 2],[3 4],[4 5]};
axx{1} = [-1.5,1.5,-1,2]; axx{2} = [0,1.2,0,1.2];axx{3} = [0,2,0,2];
numofplots = 15;
for plotRun=1:length(projectedDims)

    figure
    % set axis
    axis(axx{plotRun});
    
    % plot initial set
    handleX0 = plot(X0,projectedDims{plotRun},'k-','LineWidth',2);
    hold on

    
    % plot reachable sets from model (no noise)
    for iSet=1:numofplots
         handleModel=  plot(X_model_noise{iSet},projectedDims{plotRun},'b','Filled',true,'FaceColor',[.8 .8 .8],'EdgeColor','b');

    end
    
    % plot reachable sets from data (without noise)
    for iSet=1:numofplots
       handleData =  plot(X_data{iSet},projectedDims{plotRun},'r');
    end

    % label plot
    xlabel(['x_{',num2str(projectedDims{plotRun}(1)),'}']);
    ylabel(['x_{',num2str(projectedDims{plotRun}(2)),'}']);
    
    % skip warning for extra legend entries
    warOrig = warning; warning('off','all');
    legend([handleX0,handleModel,handleData], 'Initial Set','Set from Model','Set from Data','Location','northwest');
    %legend('Initial Set','Set from Model','Set from Data','Location','northwest');

    warning(warOrig);
    ax = gca;
    ax.FontSize = 16;
    set(gcf, 'Position',  [50, 50, 800, 400])
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left bottom ax_width ax_height];
end
%close;


    function nfro = fun(x)
        
        gen = X_data{i}.generators;
        nfro=  norm(X_data{i}.center + (gen*x) - point );
    end


    function [c,ceq] = circlecon(x)
        c = [];%abs(x) -ones(length(x),1);
        ceq = [] ;
    end

end
%------------- END OF CODE --------------