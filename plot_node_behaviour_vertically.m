figure(6)
clf
tl = tiledlayout(num_models,2,"TileIndexing","rowmajor","TileSpacing","compact");

models_short = ["Empirical","ER","SBM","Deg-corr SBM","Config","Comm-config","Analytic"];

for model_setup=1:num_models
xlabel(tl,'Time','FontSize',12)
    
    % Plot the average node state in C0
    nexttile
    plot1a = plot(times,ave_node_state_C0_store(:,:,model_setup),'Color',[0.00 0.45 0.74],'LineWidth',0.1);
    hold on
    plot1c = plot(times,expected_state_C0(:,model_setup),'k--','LineWidth',1.5);
    plot1b = plot(times,observed_mean_state_C0(:,model_setup),'k','LineWidth',1.5);
    ylim([0 1])
    yyaxis left
    ylabel('State','FontSize',10)
    yyaxis right
    set(gca,'YTick',[]);

    % Plot the average node state in C1
    nexttile
    plot1d = plot(times,ave_node_state_C1_store(:,:,model_setup),'Color',[0.8500 0.3250 0.098],'LineWidth',0.1);
    hold on
    plot(times,expected_state_C1(:,model_setup),'k--','LineWidth',1.5)
    plot(times,observed_mean_state_C1(:,model_setup),'k','LineWidth',1.5)
    ylim([0 1])
    yyaxis right
    set(gca,'YTick',[]);
    ylabel({'',models_short(model_setup)},'FontSize',12,'Color',[0 0 0],'FontWeight','bold')
end

nexttile(1)
title('C0','FontSize',12)
nexttile(2)
title('C1','FontSize',12)