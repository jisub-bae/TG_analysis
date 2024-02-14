%% plot SDI

cd("C:\Users\user\Desktop\20231018 50T_5_shift decoding\results\SDI\Results4")
list=dir("*.mat");
list={list.name}';

for a=1:length(list);
    
    load(list{a})
    
    figure
    
    set(gcf,'color','w')
    set(gca,'linewidth',2);
    
    h = boundedline(SDI.time,  SDI.di, SDI.di_se, 'alpha', 'cmap', [1 0 0], 'transparency', 0.3);
%     h = boundedline(SDI.time,  SDI.di, SDI.di_se, 'alpha', 'cmap', [1 0.3 0], 'transparency', 0.3);
    h.LineWidth = 1;
    hold on;
    h = boundedline(SDI.time,  SDI.si, SDI.si_se, 'alpha', 'cmap', [0 0 1], 'transparency', 0.3);
%     h = boundedline(SDI.time,  SDI.si, SDI.si_se, 'alpha', 'cmap', [0 0.7 0], 'transparency', 0.3);
    h.LineWidth = 1;
    
    h=plot(SDI.time, zeros(1,length(SDI.time)), 'k');
    h=plot([0 0],[-0.2 1],'k')
    
    xticks(-0.5:0.1:1)
    xlim([-0.2 1]); ylim([-0.2 1.2]);
    
    set(gcf,'OuterPosition', [3, 270, 840, 420])
    
end

