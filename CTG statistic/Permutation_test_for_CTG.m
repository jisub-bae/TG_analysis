clear all; close all; clc;
%% directory setting for fieldtrip
addpath('Z:\toolbox\ssPreprocessTool\fieldtrip-20161023' );
ft_defaults;

%% data loading..
path="C:\Users\user\Desktop\20231018 50T_5_shift decoding\results\CTG statistic";
cd(path)

load CTG_DATA_frontal;

%% set experiment time
% saveMode=1;
saveMode=0;
trainAxis = 1:1:16;

twin = [-0.500  1.000];
Ts=500;
winsize=10;
% winsize = 0.04; %0.02
% stpsize = 0.01; %0.002
% stw=twin(1):stpsize:(twin(2)-winsize); % Starting time point of each window
stw=linspace(twin(1), twin(2), size(CTG_DATA,3)); % Time point
enum=-0.2;
[~,edts]=min(abs(stw-enum));


%% number of subjects
NoS = size(CTG_DATA,1);
%% define the first data to be compared 
for i = 1:size(CTG_DATA,2)
    M_HO{i,1} = [];
    M_HO{i,1}.label = {'1'};
    M_HO{i,1}.freq = 1:1:size(CTG_DATA,3);
    M_HO{i,1}.time = stw;
    M_HO{i,1}.dimord = 'subj_chan_freq_time';
    M_HO{i,1}.powspctrm = CTG_DATA(:,i,:,:);
end; clear i;
%% define the second data to be compared / if you are going to compare the first data with 0, creat a data with a value of 0
bmat = zeros(size(CTG_DATA,1),size(CTG_DATA,2),size(CTG_DATA,3),size(CTG_DATA,4))+0.5;
for i = 1:size(bmat,2)
    BL{i,1} = [];
    BL{i,1}.label = {'1'};
    BL{i,1}.freq = 1:1:size(CTG_DATA,3);
    BL{i,1}.time = stw;
    BL{i,1}.dimord = 'subj_chan_freq_time';
    BL{i,1}.powspctrm = bmat(:,i,:,:);
end; clear i;
%% Let's configure cfg file for permutation test
St= -0.5; Et = 1.0;
cfg = [];
cfg.channel     = 1;
cfg.latency     = [St+0.3 Et]; % Time of interest
cfg.neighbours = [];
cfg.method      = 'montecarlo'; % 'montecarlo' / 'analytic' / 'stats' / 'crossvalidate'
cfg.statistic   = 'ft_statfun_depsamplesT'; %'ft_statfun_depsamplesT' %  if you want to compare more than two experimental conditions, you should choose an F-statistic ('ft_statfun_indepsamplesF' or 'ft_statfun_depsamplesFmultivariate') // 'ft_statfun_depsamplesT' //'depsamplesT' // 'ttest'
cfg.correctm    = 'cluster';
cfg.clusteralpha= 0.05;
cfg.clusterstatistic = 'maxsum'; % test statistic that will be evaluated under the  permutation distribution.
cfg.clusterthreshold = 'nonparametric_individual';
cfg.tail        = 1; % -1, 1 or 0 (default = 0); one-sided or two-sided test
cfg.clustertail = 1; % 0, 1
cfg.alpha = 0.05;
% cfg.correcttail = 'prob'; % correct p-values or alpha-values when doing a two-sided test, 'alpha','prob' or 'no' (default = 'no')
cfg.numrandomization = 1000; % number of draws from the permutation distribution
% cfg.frequency        = 'all';
% cfg.frequency        = [240 740];
cfg.frequency        = [150 740];
%% cfg.design (the independent variance, the subject number)
Nsub = NoS;
cfg.design(1,1:2*Nsub)  = [ones(1,Nsub) 2*ones(1,Nsub)];
cfg.design(2,1:2*Nsub)  = [1:Nsub 1:Nsub];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
%% permutation test!
for i = 1:size(CTG_DATA,2)
    [stat] = ft_freqstatistics(cfg, M_HO{i,1},BL{i,1})
    stat_CTG{i,1} = stat;
end; clear i;
%% set stat.label 1
for i = 1:size(CTG_DATA,2)
    if  isfield(stat_CTG{i,1},'posclusters') ~= 1
        stat_CTG{i,1}.posclusters = [];
    end
    
    if  isfield(stat_CTG{i,1}.posclusters,'prob') ~= 1
        stat_CTG{i,1}.posclusters.prob = [];
    end
    
    if  isfield(stat_CTG{i,1},'posclusterslabelmat') ~= 1
        stat_CTG{i,1}.posclusterslabelmat = zeros(1,size(CTG_DATA,3),size(CTG_DATA,4));
    end
end; clear i;

%% set stat.label 2

sign = 0.05;

for i = 1:size(CTG_DATA,2)
    M0sp = find([stat_CTG{i,1}.posclusters.prob]<sign);
    M0_ps = squeeze(stat_CTG{i,1}.posclusterslabelmat);
    stat_S_ps{1,i} = M0sp;
    stat_S_ps{2,i} = M0_ps;
    clear M0sp M0_ps;
end; clear i M0sp M0_ps;

for i = 1:size(CTG_DATA,2)
    prew1 = stat_S_ps{1,i};
    prew2 = stat_S_ps{2,i};
    for i1 = 1:size(prew2,1)
        for i2 = 1:size(prew2,2)
            for i3 = 1:size(prew1,2)
                if prew2(i1,i2) == prew1(i3)
                    prew2(i1,i2) = 999;
                end
            end
        end
    end; clear i1 i2 i3;
    stat_S_ps2{i,1} = prew2;
end; clear i prew1 prew2;

for i = 1:size(CTG_DATA,2)
    prew1 = stat_S_ps2{i,1};
    for i1 = 1:size(prew1,1)
        for i2 = 1:size(prew1,2)
            if prew1(i1,i2) ~= 999
                prew1(i1,i2) = 0;
            end
        end
    end; clear i1 i2;
    stat_S_ps3{i,1} = [zeros(size(CTG_DATA,3)-size(prew1,1),size(CTG_DATA,4));...
        [zeros(size(prew1,1),size(CTG_DATA,4)-size(prew1,2)) prew1]];
end; clear i prew1;

for i = 1:size(CTG_DATA,2)
    smat{i,1}=squeeze(stat_CTG{i,1}.stat);
end
%%  make the plot (stat)

% cmin = min(min([smat{1,1} smat{1,2} smat{1,3}]));
% cmax = max(max([smat{1,1} smat{1,2} smat{1,3}]));

mdata=squeeze(mean(CTG_DATA,1));

for kl=1:size(CTG_DATA,2),
    figure,
    imagesc(stw(edts:end),stw(edts:end),squeeze(mdata(kl,(edts:end),(edts:end))));colormap('jet');hh=colorbar;ylabel(hh,'AUC');hold on;
    contour(stw(edts:end),stw(edts:end), stat_S_ps3{kl,1}(edts:end,edts:end),[900 1000], 'k', 'LineWidth', 1);
    xlabel('Testing time(s)');
    ylabel('Training time(s)');
    set(gca,'YDir','normal');
    plot([0 0],[-0.2 1],'color', 'k', 'LineWidth', 1);
    plot([-0.2 1],[0 0],'color', 'k', 'LineWidth', 1);
%     plot([0 0],[0 0.15],'k', 'LineWidth', 5);
%     plot([0 0.15],[0 0],'k', 'LineWidth', 5);
    plot([-0.2 1],[-0.2 1],'color', [0.3 0.3 0.3], 'LineWidth', 1, 'LineStyle', '--');
    axis square;
    %     set(gca,'CLim',[cmin cmax]);
    
    set(gca,'CLim',[0.44 0.56]);
%     colormap redblue;
%     
    %     if saveMode,
    %         cPlotDir = fullfile('xx:\xxx\xxx\xxx\xxxxxx');
    %         if exist(cPlotDir,'dir')==0, mkdir(cPlotDir); end
    %         print([cPlotDir,strcat('\xxx_',num2str(kl))],'-dpng','-r0');
    %         savefig([cPlotDir,strcat('\xxx_',num2str(kl))]);
    %         close gcf
    %     end
end; clear kl;
%%  make the plot (slope)
% CMin = squeeze(min(min(min(min(CTG_DATA)))));
% CMax = squeeze(max(max(max(max(CTG_DATA)))));
% 
% for kl=1:size(CTG_DATA,2),
%     figure,
%     imagesc(stw(1:edts),stw(1:edts),squeeze(mean(CTG_DATA(:,kl,1:edts,1:edts),1)));colormap('jet');hh=colorbar;ylabel(hh,'');hold on;
%     contour(stw(1:edts),stw(1:edts), stat_S_ps3{kl,1}(1:edts,1:edts),[900 1000], 'k', 'LineWidth', 4);
%     xlabel('Generalization time');
%     ylabel('Training time');
%     set(gca,'YDir','normal');
%     plot([0 0],[-0.2 0.6],'color', [0.3 0.3 0.3], 'LineWidth', 3, 'LineStyle', '--');
%     plot([-0.2 0.6],[0 0],'color', [0.3 0.3 0.3], 'LineWidth', 3, 'LineStyle', '--');
%     plot([0 0],[0 0.15],'k', 'LineWidth', 5);
%     plot([0 0.15],[0 0],'k', 'LineWidth', 5);
%     plot([-0.2 0.6],[-0.2 0.6],'color', [0.3 0.3 0.3], 'LineWidth', 3);
%     axis square;
%     %          set(gca,'CLim',[CMin CMax]);
%     set(gca,'CLim',[-0.3 0.3]);
% %     colormap redblue;
%     
%     %     if saveMode,
%     %         cPlotDir = fullfile('xx:\xxx\xxx\xxx\xxxxxx');
%     %         if exist(cPlotDir,'dir')==0, mkdir(cPlotDir); end
%     %         print([cPlotDir,strcat('\xxx_',num2str(kl))],'-dpng','-r0');
%     %         savefig([cPlotDir,strcat('\xxx_',num2str(kl))]);
%     %         close gcf
%     %     end
% end; clear kl;
% %%  make the plot (diagonal)
% for kl=1:size(CTG_DATA,2),
%     figure,
%     
%     prew1 = squeeze(CTG_DATA(:,kl,1:edts,1:edts));
%     for yy1 = 1:size(prew1,1)
%         for yy2 = 1:size(prew1,2)
%             prew2(yy1,yy2) = prew1(yy1,yy2,yy2);
%         end
%     end; clear yy1 yy2;
%     
%     for yy = 1:size(prew1,2)
%         preSts(1,yy) = stat_S_ps3{kl,1}(yy,yy);
%     end; clear yy;
%     
%     H1=plot(stw(1:edts), mean(prew2,1), 'Color', 'k', 'LineWidth', 8);hold on;
%     H(1) = shadedErrorBar(stw(1:edts), prew2(:,1:edts), {@mean, @(x) std(x)/sqrt(NoS)}, {'-k', 'LineWidth', 8}, 0);
%     hfig = gcf; set(hfig,'color',[1 1 1]);
%     xlims    = [-0.2 stw(edts)]; colmap   = [1 1 1]*0.70;
%     ylims    = [-0.15 0.3];
%     plot(xlims,[0 0],'k','LineWidth', 6); hold on;
%     plot([0 0],ylims,'k','LineWidth', 6);
%     plot([-0.05 0],[ylims(2) ylims(2)],'k','LineWidth', 6);
%     plot([-0.05 0],[ylims(1) ylims(1)],'k','LineWidth', 6);
%     plot([xlims(1) xlims(1)],[-0.025 0.025],'k','LineWidth', 6);
%     plot([0.6 0.6],[-0.025 0.025],'k','LineWidth', 6);
%     plot([0.2 0.2],[-0.025 0.025],'k','LineWidth', 6);
%     plot([0.4 0.4],[-0.025 0.025],'k','LineWidth', 6);
%     
%     if sum(find(preSts >= 999))
%         for iz = 1:size(stat_CTG{kl,1}.prob,3)
%             stastat.prob_B(1,iz)=squeeze(stat_CTG{kl,1}.prob(1,iz,iz));
%             stastat.psmat(1,iz)=squeeze(stat_CTG{kl,1}.posclusterslabelmat(1,iz,iz));
%         end; clear iz;
%         for ii = 1:size(stat_CTG{kl,1}.posclusters,2)
%             if stat_CTG{kl,1}.posclusters(ii).prob < 0.05 && sum(find(stastat.psmat==ii))~= 0
%                 ptp=(find(stastat.psmat==ii));
%                 stat_etm_B = stat_CTG{kl,1}.time(ptp);
%                 tat = find(stastat.prob_B<0.05 & stastat.psmat == ii);
%                 stat_Bm=zeros(size(tat,2),size(tat,2));
%                 for hj = 1:size(tat,2)
%                     if hj == 1
%                         yuu = 1;
%                     else
%                         yuu = inu;
%                     end
%                     if hj == 1
%                         stat_Bm(1,hj) = stat_etm_B(hj);
%                         inu = hj;
%                     elseif hj < size(tat,2) && tat(hj) - 1 == tat(hj-1)
%                         stat_Bm(yuu,hj) = stat_etm_B(hj);
%                     elseif hj < size(tat,2) && tat(hj) - 1 > tat(hj-1) && tat(hj) + 1 == tat(hj+1)
%                         inu = 1+hj;
%                         yuu = inu;
%                         stat_Bm(yuu,hj) = stat_etm_B(hj);
%                     elseif hj < size(tat,2) && tat(hj) - 1 > tat(hj-1) && tat(hj) + 1 < tat(hj+1)
%                         inu = 1+hj;
%                         yuu = inu;
%                         stat_Bm(yuu,hj) = stat_etm_B(hj);
%                     elseif hj == size(tat,2) && tat(hj) - 1 == tat(hj-1)
%                         stat_Bm(yuu,hj) = stat_etm_B(hj);
%                     elseif hj == size(tat,2) && tat(hj) - 1 > tat(hj-1)
%                         yuu = hj;
%                         stat_Bm(yuu,hj) = stat_etm_B(hj);
%                     end
%                 end
%                 for iy = 1:size(stat_Bm,1)
%                     if sum(stat_Bm(iy,:) ~=0)
%                         prewu1 = stat_Bm(iy,:); prewu2 = find(prewu1==0); prewu1(prewu2) = [];
%                         plot([prewu1(1) prewu1(end)+0.001],[0.002 0.002],'r', 'LineWidth', 10);  % marking of the stimulus duration
%                         hold on;
%                         clear prewu1 prewu2;
%                     end
%                 end
%                 clear stat_Bm hj inu yuu iy;
%             end
%         end
%     end; clear ii;
%     
%     %     if saveMode,
%     %         cPlotDir = fullfile('xx:\xxx\xxx\xxx\xxxxxx');
%     %         if exist(cPlotDir,'dir')==0, mkdir(cPlotDir); end
%     %         print([cPlotDir,strcat('\xxx_',num2str(kl))],'-dpng','-r0');
%     %         savefig([cPlotDir,strcat('\xxx_',num2str(kl))]);
%     %         close gcf
%     %     end
% end; clear kl;
% 
% %% plot the decoding weight values of each channel
% 
% load CTG_DATA_WEIGHT;
% 
% % z_mat=[1 17 32 43 44 48 49 56 107 113 114 119 120 125 126 127 128];
% for ii = 1:size(CTG_DATA_WEIGHT,1)
%     for ii2 = 1:size(CTG_DATA_WEIGHT,2)
%         for ii3 = 1:size(CTG_DATA_WEIGHT,3)
%             tp_mat{ii,ii2,ii3} = zeros(128,16,size(CTG_DATA_WEIGHT{1,1},3));
%         end
%     end
% end
% 
% for ii = 1:size(CTG_DATA_WEIGHT,1)
%     for ii2 = 1:size(CTG_DATA_WEIGHT,2)
%         tp_mat{ii,ii2}([2:16 18:31 33:42 45:47 50:55 57:106 108:112 115:118 121:124],:,:) = CTG_DATA_WEIGHT{ii,ii2}(:,:,:);
%     end
% end; clear ii ii2;
% 
% if ~true
%     V = mrC_readELPfile(elpFile,true,true,[-2 1 3]);
%     V = mrC_flattenELPdata(V);
% else
%     tmp = load('defaultFlatNet.mat');
%     V = [ tmp.xy, zeros(128,1) ];
% end
% 
% F = mrC_EGInetFaces(false);
% 
% wls=[5 8 11 6 9 12]; % set 1 : 5 8 11 / set2 : 6 9 12
% tls=[6 21 36 51 66 81];
% 
% for ii = 1:size(tp_mat,1)
%     for ii2 = 1:size(tp_mat,2)
%         for ii3 = 1:(size(tls,2)-1)
%             tp_M{ii,ii2,ii3} = mean(squeeze(tp_mat{ii,ii2}(:,wls((ii*ii2)+((ii-1)*(3-ii2))),((6+(15*ii3)-15)):(6+(15*ii3)))),2);
%         end
%     end
% end; clear ii ii2 ii3;
% 
% for a1 = 1:size(tp_M,2)
%     for a2 = 1:size(tp_M,3)
%         %         tp_G{a1,a2} = mean([tp_M{1,a1,a2} tp_M{2,a1,a2}],2) - mean([tp_M{1,a1,1} tp_M{2,a1,1}],2);
%         tp_G{a1,a2} = mean([tp_M{1,a1,a2} tp_M{2,a1,a2}],2);
%     end
% end
% 
% figure;
% for i = 1:size(tp_G,1)% 1 2 3
%     for k = 1:size(tp_G,2) % 1 2 3 4
%         subplot(size(tp_G,1),size(tp_G,2),(((i-1)*5)+k));
%         
%         colormap jet;
%         H = patch('Vertices',V(1:128,:),'Faces',F,'FaceColor','interp');
%         set(gca,'dataaspectratio',[1 1 1])
%         set(H,'facevertexcdata',tp_G{i,k}(:,1));
%         hold on;
%         
%         if i == 1 && k == 5
%             title('450-600 ms / first tone');
%         elseif i == 2  && k == 5
%             title('450-600 ms / second tone');
%         elseif i == 3  && k == 5
%             title('450-600 ms / third tone');
%         end
%         
%         if k == 1
%             title('-150-0 ms');
%         elseif k == 2
%             title('0-150 ms');
%         elseif k == 3
%             title('150-300 ms');
%         elseif k == 4
%             title('300-450 ms');
%         elseif k == 5
%             title('450-600 ms');
%         end
%         
%         axis off
%         colorbar
%         set(gca,'CLim',[-6e+06 6e+06]);
%     end
% end
% 
% % if saveMode,
% %     cPlotDir = fullfile('xx:\xxx\xxx\xxx\xxxxxx');
% %     if exist(cPlotDir,'dir')==0, mkdir(cPlotDir); end
% %     print([cPlotDir,'\xxxxx_xxxx'],'-dpng','-r0');
% %     savefig([cPlotDir,'\xxxxx_xxxx']);
% %     close gcf
% % end
% 
% tp_GG{1,1} = mean([tp_G{1,1} tp_G{2,1} tp_G{3,1}],2);
% tp_GG{1,2} = mean([tp_G{1,2} tp_G{2,2} tp_G{3,2}],2);
% tp_GG{1,3} = mean([tp_G{1,3} tp_G{2,3} tp_G{3,3}],2);
% tp_GG{1,4} = mean([tp_G{1,4} tp_G{2,4} tp_G{3,4}],2);
% tp_GG{1,5} = mean([tp_G{1,5} tp_G{2,5} tp_G{3,5}],2);
% 
% figure;
% for jj = 1:size(tp_GG,2)
%     
%     subplot(1,5,jj);
%     
%     colormap jet;
%     H = patch('Vertices',V(1:128,:),'Faces',F,'FaceColor','interp');
%     set(gca,'dataaspectratio',[1 1 1])
%     set(H,'facevertexcdata',tp_GG{1,jj}(:,1));
%     hold on;
%     
%     if jj == 1
%         title('-150-0 ms');
%     elseif jj == 2
%         title('0-150 ms');
%     elseif jj == 3
%         title('150-300 ms');
%     elseif jj == 4
%         title('300-450 ms');
%     elseif jj == 5
%         title('450-600 ms');
%     end
%     
%     axis off
%     colorbar
%     set(gca,'CLim',[-4e+06 4e+06])
% end
% 
% % if saveMode,
% %     cPlotDir = fullfile('xx:\xxx\xxx\xxx\xxxxxx');
% %     if exist(cPlotDir,'dir')==0, mkdir(cPlotDir); end
% %     print([cPlotDir,'\xxxxx_xxxx'],'-dpng','-r0');
% %     savefig([cPlotDir,'\xxxxx_xxxx']);
% %     close gcf
% % end
% 
% 
% 
% 
% 
