% draw pain detection threshold event
% author: sile Hu
% date: 2017-3-16
function lgd = drawPValueThresh(options,edges,lgd,subplotIdx)
subplot(subplotIdx);
hold on;
p_thr165 = plot([edges(1) edges(end)],[options.threshold options.threshold],[options.plotOpt.threshColor '--'],'linewidth',1,'DisplayName',['p-value threshold:' num2str(options.threshold)]);
lgd = [lgd p_thr165];