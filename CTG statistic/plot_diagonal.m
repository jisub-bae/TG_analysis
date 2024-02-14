type=3 %1: occipital, 2: parietal, 3: frontal

sigWind=xlsread('sigWind.xlsx');
list={'CTG_DATA_occipital','CTG_DATA_parietal','CTG_DATA_frontal'};

load(list{type})

data=CTG_DATA(:,5,:,:);
data=squeeze(data);

ff=[];
for a=1:740;
    
    dd=data(:,a,a);
    ff=[ff dd];
end

win=10; %timewindow (ms), if sampling rate is 500hz, win=10 is 20ms
Ts=500; %sampling rate

final=[];
for k=1:length(ff(:,1));
    
    temp=ff(k,:);
    s_temp=smoothdata(temp,'gaussian',win);
    final=[final;s_temp];
end

data1=final;
data2=final;
data3=final;

num_timepoints=740;
ss=1:20;

figure
set(gcf,'color','w')
% x_scale=linspace(-0.5+1/Ts,1-1/Ts,num_timepoints)*1000;
x_scale=linspace(-0.5+1/Ts,1-1/Ts,num_timepoints);
set(gca,'linewidth',2);

%     h_cat_grad=plot(x_scale,m_data, 'color',[0,0,1])

e_cat_grad=std(data1,0,1)/sqrt(length(ss));
m_data=mean(data1);

h_cat_grad=boundedline(x_scale,m_data(1:num_timepoints),...
    e_cat_grad,'cmap',[0,0,1],'transparency', 0.2,'alpha'); %plot category grad
set(h_cat_grad,'linewidth',2);

hold on

%
plot(x_scale, zeros(1,length(x_scale))+0.5, 'k');
xline(0,'k','linewidth',2)

hvsig=sigWind(type,:);
hvsig=hvsig(~isnan(hvsig));
% 
% if length(hvsig)~=0
%     for k = 1:length(hvsig)
%         plot(hvsig(k), 0.49, 'b.','MarkerSize',20);
%     end
% end

file_Name=sprintf('dia_%d.fig',type)
savefig(file_Name)