%% load data

path="C:\Users\user\Desktop\20230118 CTG statistic";

DataFile=fullfile(path,"CTG_DATA.mat");
load(DataFile)   %participants X experiment condition X training time X testing time

%% Analyze SDI

path2="C:\Users\user\Desktop\20230119 SDI calculate\Results"


for a=1:size(CTG_DATA,2);
    
    raw=squeeze(CTG_DATA(:,a,:,:));
    strtemp=sprintf("SDI_%d",a);

    [SDI.si, SDI.di, SDI.time] = calculateSDI(raw,[-0.5 1],zeropoint);
    [SDI.si_se, SDI.di_se] = bootstrapSDI(raw, [-0.5 1],zeropoint);
    
    fileDir=fullfile(path2,strtemp);
    
    save(fileDir, 'SDI');
    
end