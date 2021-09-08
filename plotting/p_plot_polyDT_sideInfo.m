clear all
%close all

load('workspace\paper\run_poly2.mat')
% plot initial set
figure('Renderer', 'painters', 'Position', [10 10 800 900])
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
        'Initial Set $\mathcal{X}_0$','Set from model','Simulation points','Set from data $\hat{\mathcal{R}}_k^p$','Set from data $\bar{\mathcal{R}}_k^p$','Set from data $\bar{\mathcal{R}}_k^{s,p}$','Location','northwest','Interpreter','latex');
  
 
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
