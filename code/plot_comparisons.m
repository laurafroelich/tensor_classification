addpath(genpath('code/'))
addpath('../matlab_additions/immoptibox/')
run('colsandlinestyles.m')
savefig = false;

load('results/aucs.mat')
%% check all data is there
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
LDAcol = 'b';
HODAcol = 'r';
method_cols = {LDAcol, DGTDAcol, DATERcol, DATEReigcol, CMDAcol, ...
    HODAcol, ManPDAcol, ManTDAcol, ManPDA_normsratiocol,...
    ManTDA_normsratiocol, Tuckercol, Parafaccol, Tucker2col};
method_markers = {'d', 'o', '*', '*', '*', '*', 'o', ...
    's', 's', 's', 's', 'x', 'x', 'x'};

% check all iterations were actually run
nmethods = length(aucs); % if the last method has values, all methods have been run
cur_aucs = aucs{nmethods};
for iconfig = 1:size(cur_aucs, 3)
    for i = 1:length(trainobs)
        for iit = 1:size(cur_aucs, 1)
            if any(cur_aucs(iit, i, iconfig)==0)
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
    semilogx(trainobs, mean(auc_lda(:,:,iconfig), 1))
    hold on;
    semilogx(trainobs, mean(auc_dgtda(:,:,iconfig), 1), 'o-')
    semilogx(trainobs, mean(auc_dater(:,:,iconfig), 1))
    semilogx(trainobs, mean(auc_datereig(:,:,iconfig), 1))
    semilogx(trainobs, mean(auc_cmda(:,:,iconfig), 1))
    semilogx(trainobs, mean(auc_hoda(:,:,iconfig), 1), 'o-')
    
    
    semilogx(trainobs, mean(auc_ManPDA(:,:,iconfig), 1), 'r-', 'LineWidth', 2)
    semilogx(trainobs, mean(auc_ManTDA(:,:,iconfig), 1), 'g-', 'LineWidth', 2)
    semilogx(trainobs, mean(auc_ManPDA_normsratio(:,:,iconfig), 1), 'b-', 'LineWidth', 2)
    semilogx(trainobs, mean(auc_ManTDA_normsratio(:,:,iconfig), 1), 'k-', 'LineWidth', 2)
    semilogx(trainobs, mean(auc_tucker(:,:,iconfig), 1), '-x')
    semilogx(trainobs, mean(auc_parafac(:,:,iconfig), 1), '-x')
    semilogx(trainobs, mean(auc_tucker2(:,:,iconfig), 1), '-x')
    %plot(trainobs, auc_parafac2, '-x')
    legend('LDA', 'DGTDA', 'DATER', 'DATEReig', 'CMDA', 'HODA',...
        'ManPDA', 'ManTDA', 'ManPDA\_sr', 'ManTDA\_sr', ...
        'Tucker', 'PARAFAC', 'Tucker2',... % 'parafac2', ...
        'Location', 'NorthEast')
    title('Mean AUC vs number of training samples', 'FontSize', titlefontsize)
    set(gca, 'FontSize', gcafontsize)
    xlabel('Number of training samples')
    ylabel('Mean Area Under ROC Curve')
    set(gcf,'color','w');
    ylim([0.4, 1])
    if savefig
        export_fig(figh, ['figures/auc_avgss_signal_strengh_', configs{iconfig}, '.pdf'], '-pdf')
    end
end


%% log x axis plots of std

for iconfig = 1:3
    
    figh = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    semilogx(trainobs, std(auc_lda(:,:,iconfig), [], 1))
    hold on;
    semilogx(trainobs, std(auc_dgtda(:,:,iconfig), [], 1), '-o')
    semilogx(trainobs, std(auc_dater(:,:,iconfig), [], 1))
    semilogx(trainobs, std(auc_datereig(:,:,iconfig), [], 1))
    semilogx(trainobs, std(auc_cmda(:,:,iconfig), [], 1))
    semilogx(trainobs, std(auc_hoda(:,:,iconfig), [], 1), '-o')
    
    
    semilogx(trainobs, std(auc_ManPDA(:,:,iconfig), [], 1), 'r-', 'LineWidth', 2)
    semilogx(trainobs, std(auc_ManTDA(:,:,iconfig), [], 1), 'g-', 'LineWidth', 2)
    semilogx(trainobs, std(auc_ManPDA_normsratio(:,:,iconfig), [], 1), 'b-', 'LineWidth', 2)
    semilogx(trainobs, std(auc_ManTDA_normsratio(:,:,iconfig), [], 1), 'k-', 'LineWidth', 2)
    semilogx(trainobs, std(auc_tucker(:,:,iconfig), [], 1), '-x')
    semilogx(trainobs, std(auc_parafac(:,:,iconfig), [], 1), '-x')
    semilogx(trainobs, std(auc_tucker2(:,:,iconfig), [], 1), '-x')
    %plot(trainobs, auc_parafac2, '-x')
    legend('LDA', 'DGTDA', 'DATER', 'DATEReig', 'CMDA', 'HODA',...
        'ManPDA', 'ManTDA', 'ManPDA\_sr', 'ManTDA\_sr', ...
        'Tucker', 'PARAFAC', 'Tucker2',... % 'parafac2', ...
        'Location', 'NorthEast')
    title('Standard deviation of AUC vs number of training samples', 'FontSize', titlefontsize)
    set(gca, 'FontSize', gcafontsize)
    xlabel('Number of training samples')
    ylabel('Standard deviation of Area Under ROC Curve')
    set(gcf,'color','w');
    ylim([0, 0.3])
    
    if savefig
        export_fig(figh, ['figures/auc_stds_signal_strengh_', configs{iconfig}, '.pdf'], '-pdf')
    end
end

%% plot means and stds of means on log x axis
for iconfig = 1:3
    figh = figure('units', 'normalized', 'outerposition', [0 0 1 1]);
    leghs = [];
    hold on
    for imethod = 1:nmethods
        cur_aucs = aucs{imethod};
        mean_to_plot = mean(cur_aucs(:,:,iconfig), 1);
        std_of_mean = std(cur_aucs(:,:,iconfig), [], 1)/...
            sqrt(size(cur_aucs, 1));
        leghs(imethod) = semilogx(trainobs, mean_to_plot, ...
            'Color', method_cols{imethod}, ...
            'Marker', method_markers{imethod});
        mean_m_std = mean_to_plot-2*std_of_mean;
        mean_p_std = mean_to_plot+2*std_of_mean;
        fill([trainobs;flipud(trainobs)],...
        [mean_p_std;flipud(mean_m_std)], method_cols{imethod})

        
        
    end
    ylim([0.4, 1])
    
    legend(leghs, 'LDA', 'DGTDA', 'DATER', 'DATEReig', 'CMDA', 'HODA',...
        'ManPDA', 'ManTDA', 'ManPDA\_sr', 'ManTDA\_sr', ...
        'Tucker', 'PARAFAC', 'Tucker2',...
        'Location', 'NorthEast')
    set(gca, 'FontSize', gcafontsize)
    xlabel('Number of training samples')
    ylabel('Area Under ROC Curve')
    set(gcf,'color','w');
end



