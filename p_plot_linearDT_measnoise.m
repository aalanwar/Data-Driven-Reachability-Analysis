clear all
close all

load('AMRTest\workspace\good\run_3appIn10S2_0_005w_0_002v_red390.mat')

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
 %   legend([handleX0,handleModel,handleData,handleData_av,handleData_cmz],...
 %       'Initial Set','Set from Model','Set from Data (practical)','Set From Data (Zono-Asm. 2)','Set From Data (ConZono-Asm. 2)','Location','northwest');
    legend([handleX0,handleModel,handleData,handleData_av,handleData_cmz],...
        'Initial set $\mathcal{X}_0$','Set from model $\mathcal{R}_k$','Set from data $\tilde{\mathcal{R}}_k^{m}$','Set from data $\hat{\mathcal{R}}_k^{m}$','Set from data $\bar{\mathcal{R}}_k^{m}$','Interpreter','latex');

    warning(warOrig);
    ax = gca;
    ax.FontSize = 22;
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

