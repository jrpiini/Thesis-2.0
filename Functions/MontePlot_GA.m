function MontePlot_GA(MonteCarloResults, mean_or_median)

close all
% 1: diplay mean
% 2: display median

coe_true = curtis_coe(MonteCarloResults.r0_truth, MonteCarloResults.v0_truth, 398600.44189);
a_truth = coe_true(end);
ecc_truth = coe_true(2);
inc_truth = coe_true(4);
RAAN_truth = coe_true(3);
AoP_truth = coe_true(5);
TA_truth = coe_true(6);


% row1 = [a_truth     median(MonteCarloResults.a_t0) abs(median(MonteCarloResults.a_t0) - a_truth)/a_truth*100   medianian(MonteCarloResults.a_t0)   abs(medianian(MonteCarloResults.a_t0) - a_truth)/a_truth*100];
% row2 = [ecc_truth median(MonteCarloResults.ecc_t0) abs(median(MonteCarloResults.ecc_t0) - ecc_truth)/ecc_truth*100 medianian(MonteCarloResults.ecc_t0) abs(medianian(MonteCarloResults.ecc_t0) - ecc_truth)/ecc_truth*100];
% row3 = [inc_truth median(MonteCarloResults.inc_t0) abs(median(MonteCarloResults.inc_t0) - inc_truth)/inc_truth*100 medianian(MonteCarloResults.inc_t0) abs(medianian(MonteCarloResults.inc_t0) - inc_truth)/inc_truth*100];
% row4 = [RAAN_truth median(MonteCarloResults.RAAN_t0) abs(median(MonteCarloResults.RAAN_t0) - RAAN_truth)/RAAN_truth*100 medianian(MonteCarloResults.RAAN_t0) abs(medianian(MonteCarloResults.RAAN_t0) - RAAN_truth)/RAAN_truth*100];
% row5 = [AoP_truth median(MonteCarloResults.AoP_t0) abs(median(MonteCarloResults.AoP_t0) - AoP_truth)/AoP_truth*100   medianian(MonteCarloResults.AoP_t0)   abs(medianian(MonteCarloResults.AoP_t0) - AoP_truth)/AoP_truth*100];
% row6 = [TA_truth median(MonteCarloResults.TA_t0) abs(median(MonteCarloResults.TA_t0) - TA_truth)/TA_truth*100    medianian(MonteCarloResults.TA_t0)    abs(medianian(MonteCarloResults.TA_t0) - TA_truth)/TA_truth*100];

% GA_coes = [MonteCarloResults.a_t0_GA MonteCarloResults.ecc_t0_GA MonteCarloResults.inc_t0_GA MonteCarloResults.RAAN_t0_GA MonteCarloResults.AoP_t0_GA MonteCarloResults.TA_t0_GA];
% 
% GA_a_t0_percent_error = abs(median(MonteCarloResults.a_t0_GA) - a_truth)/a_truth*100;

if mean_or_median == 1
    row1 = [a_truth     mean(MonteCarloResults.a_t0_GA)         abs(mean(MonteCarloResults.a_t0_GA) - a_truth)/a_truth*100                 mean(MonteCarloResults.a_t0_GA_ACO)         abs(mean(MonteCarloResults.a_t0_GA_ACO) - a_truth)/a_truth*100               mean(MonteCarloResults.a_t0_GA_3pts)         abs(mean(MonteCarloResults.a_t0_GA_3pts) - a_truth)/a_truth*100                       mean(MonteCarloResults.a_t0_GA_3pts_ACO)         abs(mean(MonteCarloResults.a_t0_GA_3pts_ACO) - a_truth)/a_truth*100                 mean(MonteCarloResults.a_t0_GE) abs(mean(MonteCarloResults.a_t0_GE) - a_truth)/a_truth*100                   mean(MonteCarloResults.a_t0_GE_ACO) abs(mean(MonteCarloResults.a_t0_GE_ACO) - a_truth)/a_truth*100                    ];
    row2 = [ecc_truth   mean(MonteCarloResults.ecc_t0_GA)       abs(mean(MonteCarloResults.ecc_t0_GA) - ecc_truth)/ecc_truth*100           mean(MonteCarloResults.ecc_t0_GA_ACO)       abs(mean(MonteCarloResults.ecc_t0_GA_ACO) - ecc_truth)/ecc_truth*100         mean(MonteCarloResults.ecc_t0_GA_3pts)       abs(mean(MonteCarloResults.ecc_t0_GA_3pts) - ecc_truth)/ecc_truth*100                 mean(MonteCarloResults.ecc_t0_GA_3pts_ACO)       abs(mean(MonteCarloResults.ecc_t0_GA_3pts_ACO) - ecc_truth)/ecc_truth*100           mean(MonteCarloResults.ecc_t0_GE) abs(mean(MonteCarloResults.ecc_t0_GE) - ecc_truth)/ecc_truth*100           mean(MonteCarloResults.ecc_t0_GE_ACO) abs(mean(MonteCarloResults.ecc_t0_GE_ACO) - ecc_truth)/ecc_truth*100            ];
    row3 = [inc_truth   mean(MonteCarloResults.inc_t0_GA)       abs(mean(MonteCarloResults.inc_t0_GA) - inc_truth)/inc_truth*100           mean(MonteCarloResults.inc_t0_GA_ACO)       abs(mean(MonteCarloResults.inc_t0_GA_ACO) - inc_truth)/inc_truth*100         mean(MonteCarloResults.inc_t0_GA_3pts)       abs(mean(MonteCarloResults.inc_t0_GA_3pts) - inc_truth)/inc_truth*100                 mean(MonteCarloResults.inc_t0_GA_3pts_ACO)       abs(mean(MonteCarloResults.inc_t0_GA_3pts_ACO) - inc_truth)/inc_truth*100           mean(MonteCarloResults.inc_t0_GE) abs(mean(MonteCarloResults.inc_t0_GE) - inc_truth)/inc_truth*100           mean(MonteCarloResults.inc_t0_GE_ACO) abs(mean(MonteCarloResults.inc_t0_GE_ACO) - inc_truth)/inc_truth*100            ];
    row4 = [RAAN_truth  mean(MonteCarloResults.RAAN_t0_GA)      abs(mean(MonteCarloResults.RAAN_t0_GA) - RAAN_truth)/RAAN_truth*100        mean(MonteCarloResults.RAAN_t0_GA_ACO)      abs(mean(MonteCarloResults.RAAN_t0_GA_ACO) - RAAN_truth)/RAAN_truth*100      mean(MonteCarloResults.RAAN_t0_GA_3pts)      abs(mean(MonteCarloResults.RAAN_t0_GA_3pts) - RAAN_truth)/RAAN_truth*100              mean(MonteCarloResults.RAAN_t0_GA_3pts_ACO)      abs(mean(MonteCarloResults.RAAN_t0_GA_3pts_ACO) - RAAN_truth)/RAAN_truth*100        mean(MonteCarloResults.RAAN_t0_GE) abs(mean(MonteCarloResults.RAAN_t0_GE) - RAAN_truth)/RAAN_truth*100       mean(MonteCarloResults.RAAN_t0_GE_ACO) abs(mean(MonteCarloResults.RAAN_t0_GE_ACO) - RAAN_truth)/RAAN_truth*100        ];
    row5 = [AoP_truth   mean(MonteCarloResults.AoP_t0_GA)       abs(mean(MonteCarloResults.AoP_t0_GA) - AoP_truth)/AoP_truth*100           mean(MonteCarloResults.AoP_t0_GA_ACO)       abs(mean(MonteCarloResults.AoP_t0_GA_ACO) - AoP_truth)/AoP_truth*100         mean(MonteCarloResults.AoP_t0_GA_3pts)       abs(mean(MonteCarloResults.AoP_t0_GA_3pts) - AoP_truth)/AoP_truth*100                 mean(MonteCarloResults.AoP_t0_GA_3pts_ACO)       abs(mean(MonteCarloResults.AoP_t0_GA_3pts_ACO) - AoP_truth)/AoP_truth*100           mean(MonteCarloResults.AoP_t0_GE) abs(mean(MonteCarloResults.AoP_t0_GE) - AoP_truth)/AoP_truth*100           mean(MonteCarloResults.AoP_t0_GE_ACO) abs(mean(MonteCarloResults.AoP_t0_GE_ACO) - AoP_truth)/AoP_truth*100            ];
    row6 = [TA_truth    mean(MonteCarloResults.TA_t0_GA)        abs(mean(MonteCarloResults.TA_t0_GA) - TA_truth)/TA_truth*100              mean(MonteCarloResults.TA_t0_GA_ACO)        abs(mean(MonteCarloResults.TA_t0_GA_ACO) - TA_truth)/TA_truth*100            mean(MonteCarloResults.TA_t0_GA_3pts)        abs(mean(MonteCarloResults.TA_t0_GA_3pts) - TA_truth)/TA_truth*100                    mean(MonteCarloResults.TA_t0_GA_3pts_ACO)        abs(mean(MonteCarloResults.TA_t0_GA_3pts_ACO) - TA_truth)/TA_truth*100              mean(MonteCarloResults.TA_t0_GE) abs(mean(MonteCarloResults.TA_t0_GE) - TA_truth)/TA_truth*100               mean(MonteCarloResults.TA_t0_GE_ACO) abs(mean(MonteCarloResults.TA_t0_GE_ACO) - TA_truth)/TA_truth*100               ];
    
    T = array2table([row1;row2;row3;row4;row5;row6],'VariableNames',{'Truth','GA Mean','GA % Error','ACO GA Mean','ACO MC Mean % Error','3 pts. GA Mean','3 pts. GA Mean % Error','3 pts. ACO GA Mean','3 pts. ACO GA Mean % Error','Gauss Ext. Mean','Gauss Ext. Mean % Error', 'Gauss Ext. ACO Mean', 'Gauss Ext. ACO Mean % Error'},'RowName',{'a','ecc','inc','RAAN','AoP','TA'}); 
    T{:,:} = round(T{:,:}, 5); 
    format longG
    disp(T) 

else
    row1 = [a_truth     median(MonteCarloResults.a_t0_GA)         abs(median(MonteCarloResults.a_t0_GA) - a_truth)/a_truth*100                 median(MonteCarloResults.a_t0_GA_ACO)         abs(median(MonteCarloResults.a_t0_GA_ACO) - a_truth)/a_truth*100               median(MonteCarloResults.a_t0_GA_3pts)         abs(median(MonteCarloResults.a_t0_GA_3pts) - a_truth)/a_truth*100                       median(MonteCarloResults.a_t0_GA_3pts_ACO)         abs(median(MonteCarloResults.a_t0_GA_3pts_ACO) - a_truth)/a_truth*100                 median(MonteCarloResults.a_t0_GE) abs(median(MonteCarloResults.a_t0_GE) - a_truth)/a_truth*100                   median(MonteCarloResults.a_t0_GE_ACO) abs(median(MonteCarloResults.a_t0_GE_ACO) - a_truth)/a_truth*100                    ];
    row2 = [ecc_truth   median(MonteCarloResults.ecc_t0_GA)       abs(median(MonteCarloResults.ecc_t0_GA) - ecc_truth)/ecc_truth*100           median(MonteCarloResults.ecc_t0_GA_ACO)       abs(median(MonteCarloResults.ecc_t0_GA_ACO) - ecc_truth)/ecc_truth*100         median(MonteCarloResults.ecc_t0_GA_3pts)       abs(median(MonteCarloResults.ecc_t0_GA_3pts) - ecc_truth)/ecc_truth*100                 median(MonteCarloResults.ecc_t0_GA_3pts_ACO)       abs(median(MonteCarloResults.ecc_t0_GA_3pts_ACO) - ecc_truth)/ecc_truth*100           median(MonteCarloResults.ecc_t0_GE) abs(median(MonteCarloResults.ecc_t0_GE) - ecc_truth)/ecc_truth*100           median(MonteCarloResults.ecc_t0_GE_ACO) abs(median(MonteCarloResults.ecc_t0_GE_ACO) - ecc_truth)/ecc_truth*100            ];
    row3 = [inc_truth   median(MonteCarloResults.inc_t0_GA)       abs(median(MonteCarloResults.inc_t0_GA) - inc_truth)/inc_truth*100           median(MonteCarloResults.inc_t0_GA_ACO)       abs(median(MonteCarloResults.inc_t0_GA_ACO) - inc_truth)/inc_truth*100         median(MonteCarloResults.inc_t0_GA_3pts)       abs(median(MonteCarloResults.inc_t0_GA_3pts) - inc_truth)/inc_truth*100                 median(MonteCarloResults.inc_t0_GA_3pts_ACO)       abs(median(MonteCarloResults.inc_t0_GA_3pts_ACO) - inc_truth)/inc_truth*100           median(MonteCarloResults.inc_t0_GE) abs(median(MonteCarloResults.inc_t0_GE) - inc_truth)/inc_truth*100           median(MonteCarloResults.inc_t0_GE_ACO) abs(median(MonteCarloResults.inc_t0_GE_ACO) - inc_truth)/inc_truth*100            ];
    row4 = [RAAN_truth  median(MonteCarloResults.RAAN_t0_GA)      abs(median(MonteCarloResults.RAAN_t0_GA) - RAAN_truth)/RAAN_truth*100        median(MonteCarloResults.RAAN_t0_GA_ACO)      abs(median(MonteCarloResults.RAAN_t0_GA_ACO) - RAAN_truth)/RAAN_truth*100      median(MonteCarloResults.RAAN_t0_GA_3pts)      abs(median(MonteCarloResults.RAAN_t0_GA_3pts) - RAAN_truth)/RAAN_truth*100              median(MonteCarloResults.RAAN_t0_GA_3pts_ACO)      abs(median(MonteCarloResults.RAAN_t0_GA_3pts_ACO) - RAAN_truth)/RAAN_truth*100        median(MonteCarloResults.RAAN_t0_GE) abs(median(MonteCarloResults.RAAN_t0_GE) - RAAN_truth)/RAAN_truth*100       median(MonteCarloResults.RAAN_t0_GE_ACO) abs(median(MonteCarloResults.RAAN_t0_GE_ACO) - RAAN_truth)/RAAN_truth*100        ];
    row5 = [AoP_truth   median(MonteCarloResults.AoP_t0_GA)       abs(median(MonteCarloResults.AoP_t0_GA) - AoP_truth)/AoP_truth*100           median(MonteCarloResults.AoP_t0_GA_ACO)       abs(median(MonteCarloResults.AoP_t0_GA_ACO) - AoP_truth)/AoP_truth*100         median(MonteCarloResults.AoP_t0_GA_3pts)       abs(median(MonteCarloResults.AoP_t0_GA_3pts) - AoP_truth)/AoP_truth*100                 median(MonteCarloResults.AoP_t0_GA_3pts_ACO)       abs(median(MonteCarloResults.AoP_t0_GA_3pts_ACO) - AoP_truth)/AoP_truth*100           median(MonteCarloResults.AoP_t0_GE) abs(median(MonteCarloResults.AoP_t0_GE) - AoP_truth)/AoP_truth*100           median(MonteCarloResults.AoP_t0_GE_ACO) abs(median(MonteCarloResults.AoP_t0_GE_ACO) - AoP_truth)/AoP_truth*100            ];
    row6 = [TA_truth    median(MonteCarloResults.TA_t0_GA)        abs(median(MonteCarloResults.TA_t0_GA) - TA_truth)/TA_truth*100              median(MonteCarloResults.TA_t0_GA_ACO)        abs(median(MonteCarloResults.TA_t0_GA_ACO) - TA_truth)/TA_truth*100            median(MonteCarloResults.TA_t0_GA_3pts)        abs(median(MonteCarloResults.TA_t0_GA_3pts) - TA_truth)/TA_truth*100                    median(MonteCarloResults.TA_t0_GA_3pts_ACO)        abs(median(MonteCarloResults.TA_t0_GA_3pts_ACO) - TA_truth)/TA_truth*100              median(MonteCarloResults.TA_t0_GE) abs(median(MonteCarloResults.TA_t0_GE) - TA_truth)/TA_truth*100               median(MonteCarloResults.TA_t0_GE_ACO) abs(median(MonteCarloResults.TA_t0_GE_ACO) - TA_truth)/TA_truth*100               ];
    
    T = array2table([row1;row2;row3;row4;row5;row6],'VariableNames',{'Truth','GA Median','GA % Error','ACO GA Median','ACO GA Median % Error','3 pts. GA Median','3 pts. GA Median % Error','3 pts. ACO GA Median','3 pts. ACO GA Median % Error','Gauss Ext. Median','Gauss Ext. Median % Error', 'Gauss Ext. ACO Median', 'Gauss Ext. ACO Median % Error'},'RowName',{'a','ecc','inc','RAAN','AoP','TA'});     T{:,:} = round(T{:,:}, 5); 
    format longG
    disp(T) 
end


[sorta, sortedIndices] = sort(MonteCarloResults.a_t0_GA, 'descend');
GA_50th = sortedIndices(floor(height(sortedIndices)/2));
ecc_50th = MonteCarloResults.ecc_t0_GA(sortedIndices(floor(height(sortedIndices)/2)));
[sorta, sortedIndices] = sort(MonteCarloResults.a_t0_GE, 'descend');
GE_50th = sortedIndices(floor(height(sortedIndices)/2));
[sorta, sortedIndices] = sort(MonteCarloResults.a_t0_GA_3pts, 'descend');
GA_3pts_50th = sortedIndices(floor(height(sortedIndices)/2));


figure
subplot(1,3,1)
scatter(a_truth, ecc_truth, 70, 'LineWidth',2)
hold on
scatter(MonteCarloResults.a_t0_GA, MonteCarloResults.ecc_t0_GA, 2, 'filled')
scatter(MonteCarloResults.a_t0_GA(GA_50th), MonteCarloResults.ecc_t0_GA(GA_50th), 'green','filled')
xlabel('a [km]')
ylabel('ecc')
legend('Truth', 'MC Results', 'MC Median', 'Location','northwest')
grid on
grid minor
title('Lambert GA')
subplot(1,3,2)
scatter(a_truth, ecc_truth, 70, 'LineWidth',2)
hold on
scatter(MonteCarloResults.a_t0_GA_3pts, MonteCarloResults.ecc_t0_GA_3pts, 2, 'filled')
scatter(MonteCarloResults.a_t0_GA_3pts(GA_3pts_50th), MonteCarloResults.ecc_t0_GA_3pts(GA_3pts_50th), 'green','filled')
xlabel('a [km]')
ylabel('ecc')
legend('Truth', 'MC Results', 'MC Median', 'Location','northwest')
grid on
grid minor
title('3 pts. Lambert GA')
subplot(1,3,3)
scatter(a_truth, ecc_truth, 70, 'LineWidth',2)
hold on
scatter(MonteCarloResults.a_t0_GE, MonteCarloResults.ecc_t0_GE, 2, 'filled')
scatter(MonteCarloResults.a_t0_GE(GE_50th), MonteCarloResults.ecc_t0_GE(GE_50th), 'green','filled')
xlabel('a [km]')
ylabel('ecc')
legend('Truth', 'MC Results', 'MC Median', 'Location','northwest')
grid on
grid minor
title('Gauss Method')
sgtitle("Monte Carlo ecc-a for 15 second arc (3'' 1-\sigma noise)")


end