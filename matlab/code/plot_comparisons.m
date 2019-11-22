addpath(genpath('code/'))
addpath('../../matlab_additions/immoptibox/')
addpath(('../../matlab_additions/gridLegend_v1.4/'))
addpath(('../../matlab_additions/export_fig/'))

run('code/colsandlinestyles.m')

commitid = ''; %'cebe70d0becef29b49bf9b53a5a9c1b03c3b1a6d'; %'d544d947b1e29f1648f3685e0a1de5903fec08c3';
load(['results/', commitid, '/aucs.mat'])
savefig = true;
%% check all data is there
legend_position = [0.8 0.59 0.095 0.3];
aucs = [];
aucs{1} = auc_lda;
aucs{2} = auc_dgtda;
aucs{3} = auc_dater;
aucs{4} = auc_datereig;
aucs{5} = auc_cmda;
aucs{6} = auc_hoda;
aucs{7} = auc_ManPDA;
aucs{8} = auc_ManTDA;
aucs{9} = auc_ManPDA_normsratio;
aucs{10} = auc_ManTDA_normsratio;
aucs{11} = auc_tucker;
aucs{12} = auc_parafac;
aucs{13} = auc_tucker2;
LDAcol = [0, 0, 1];
HODAcol = [0.5, 0, 0];
method_cols = {LDAcol, DGTDAcol, DATERcol, DATEReigcol, CMDAcol, ...
    HODAcol, ManPDAcol, ManTDAcol, ManPDA_normsratiocol,...
    ManTDA_normsratiocol, Tuckercol, Parafaccol, Tucker2col};
method_markers = {'d', 'o', '*', '*', '*', 'o', ...
    'v', 'v', 'v', 'v', 'p', 'p', 'p'};

method_names = {'LDA', 'DGTDA', 'DATER', 'DATEReig', 'CMDA', 'HODA', ...
    'ManPDA', 'ManTDA', 'ManPDA\_sr', 'ManTDA\_sr', 'Tucker', 'PARAFAC',...
    'Tucker2'};

% check all iterations were actually run
nmethods = length(aucs); % if the last method has values, all methods have been run
cur_aucs = aucs{nmethods};
for iconfig = 1:size(cur_aucs, 3)
    for i = 1:length(trainobs)
        for iit = 1:size(cur_aucs, 1)
            if any(isnan(cur_aucs(iit, i, iconfig)))
                display(['iconfig: ', num2str(iconfig), ...
                    ', i: ', num2str(i), ...
                    ', iit: ', num2str(iit)])
            end
        end
    end
end


%% log x axis plots of mean
for iconfig = 1:3
    figh = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    
    leghs = [];
    for imethod = 1:nmethods
        cur_aucs = aucs{imethod};
        mean_to_plot = mean(cur_aucs(:,:,iconfig), 1);
        leghs(imethod) = semilogx(trainobs, mean_to_plot, ...
            'Color', method_cols{imethod}, ...
            'Marker', method_markers{imethod}, 'MarkerSize', 10);
        
        hold on
    end
    
    legend1 = legend(leghs, 'LDA', 'DGTDA', 'DATER', 'DATEReig', 'CMDA', 'HODA',...
        'ManPDA', 'ManTDA', 'ManPDA\_sr', 'ManTDA\_sr', ...
        'Tucker', 'PARAFAC', 'Tucker2');
    set(legend1, 'Position', legend_position, 'FontSize',23.4);
    
    title('Mean AUC vs number of training samples', 'FontSize', titlefontsize)
    set(gca, 'FontSize', gcafontsize)
    xlabel('Number of training samples')
    ylabel('Mean Area Under ROC Curve')
    set(gcf,'color','w');
    ylim([0.4, 1])
    if savefig
        export_fig(figh, ['figures/', commitid, '/auc_avgss_signal_strengh_', configs{iconfig}, '.pdf'], '-pdf')
    end
end


%% log x axis plots of std

for iconfig = 1:3
    figh = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    leghs = [];
    for imethod = 1:nmethods
        cur_aucs = aucs{imethod};
        mean_to_plot = mean(cur_aucs(:,:,iconfig), 1);
        std_of_mean = std(cur_aucs(:,:,iconfig), [], 1)/...
            sqrt(size(cur_aucs, 1));
        leghs(imethod) = semilogx(trainobs, std_of_mean, ...
            'Color', method_cols{imethod}, ...
            'Marker', method_markers{imethod}, 'MarkerSize', 10);
        
        hold on
    end
    ylim([0, 0.05])
    
    legend1 = legend(leghs, 'LDA', 'DGTDA', 'DATER', 'DATEReig', 'CMDA', 'HODA',...
        'ManPDA', 'ManTDA', 'ManPDA\_sr', 'ManTDA\_sr', ...
        'Tucker', 'PARAFAC', 'Tucker2');
    set(legend1, 'Position', legend_position, 'FontSize',23.4);
    set(gca, 'FontSize', gcafontsize)
    xlabel('Number of training samples')
    
    title('Standard deviation of AUC vs number of training samples', 'FontSize', titlefontsize)
    set(gca, 'FontSize', gcafontsize)
    xlabel('Number of training samples')
    ylabel('Standard deviation of Area Under ROC Curve')
    set(gcf,'color','w');
    
    if savefig
        export_fig(figh, ['figures/', commitid, '/auc_stds_signal_strengh_', configs{iconfig}, '.pdf'], '-pdf')
    end
end

%% plot means and stds of means on log x axis
%savefig = false
for iconfig = 1:3
    figh = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    leghs = [];
    for imethod = 1:nmethods
        cur_aucs = aucs{imethod};
        mean_to_plot = mean(cur_aucs(:,:,iconfig), 1);
        std_of_mean = std(cur_aucs(:,:,iconfig), [], 1)/...
            sqrt(size(cur_aucs, 1));
        leghs(imethod) = semilogx(trainobs, mean_to_plot, ...
            'Color', method_cols{imethod}, ...
            'Marker', method_markers{imethod}, 'MarkerSize', 10);
        
        hold on
        mean_m_std = mean_to_plot-1*std_of_mean;
        mean_p_std = mean_to_plot+1*std_of_mean;
        fill([trainobs;flipud(trainobs)],...
            [mean_p_std;flipud(mean_m_std)], method_cols{imethod})
    end
    
    legend1 = legend(leghs, 'LDA', 'DGTDA', 'DATER', 'DATEReig', 'CMDA', 'HODA',...
        'ManPDA', 'ManTDA', 'ManPDA\_sr', 'ManTDA\_sr', ...
        'Tucker', 'PARAFAC', 'Tucker2');
    set(legend1, 'Position', legend_position, 'FontSize',23.4);
    set(gca, 'FontSize', gcafontsize)
    xlabel('Number of training samples')
    ylabel('Area Under ROC Curve')
    title('Mean plus/minus standard deviation over 25 simulations')
    set(gcf,'color','w');
    
    if savefig
        export_fig(figh, ['figures/', commitid, ...
            '/auc_mean_p_m_stds_signal_strengh_', configs{iconfig}, '.pdf'], '-pdf')
    end
end

%% plot barplot for chosen trainingsample sizes
% (40, 160, 640, 1280, 2560)
indices = [3, 5, 7, 8, 9];
trainobs(indices)
for iconfig=1:3
    aucs_in_bar_plot = NaN(length(trainobs), nmethods);
    std_aucs_in_bar_plot = NaN(length(trainobs), nmethods);
    for imethod = 1:nmethods
        cur_aucs = aucs{imethod};
        aucs_in_bar_plot(:, imethod) = mean(cur_aucs(:, :, iconfig), 1);
        std_aucs_in_bar_plot(:, imethod) = std(cur_aucs(:, :, iconfig),...
            [], 1)/sqrt(size(cur_aucs, 1));
        
        
    end
    figh = figure('units', 'normalized', 'outerposition', [0 0 1 0.7]);%, ...
    hold on
    colormap(reshape(cell2mat(method_cols), 3, [])')
    hbar = bar(1:5, aucs_in_bar_plot(indices, :));
    errorbar(repmat(linspace(0.63, 1.37, nmethods), 5, 1) + repmat(0:4, nmethods, 1)', aucs_in_bar_plot(indices, :), ...
        std_aucs_in_bar_plot(indices, :),'k.')
    ylim([0.4 1.1])
    
    gridLegend(hbar, 7, method_names, 'Location', 'North', 'FontSize', legfontsize);
    
    set(gcf,'color','w');
    set(gca, 'FontSize', gcafontsize,...
        'XTickLabel', {'' trainobs(indices(1)) '' trainobs(indices(2)) '' ,...
        trainobs(indices(3)) '' trainobs(indices(4)) '' trainobs(indices(5))}, ...
        'YTickLabel', {0.4:0.1:1.1 '' ''})
    
    ylabel('Area under ROC curve on test data')
    xlabel('Training observations')
    
    if savefig
        export_fig(figh, ['figures/', commitid, '/barplot_', configs{iconfig}, '.pdf'], '-pdf')
    end
end