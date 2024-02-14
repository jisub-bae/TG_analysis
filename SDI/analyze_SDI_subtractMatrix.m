% addpath("C:\Users\user\Desktop\20230823 shift decoding\frontal\SDI")

%% load data

path="C:\Users\user\Desktop\20231018 50T_5_shift decoding\results\SDI";

DataFile=fullfile(path,"CTG_DATA_all.mat");
load(DataFile)   %participants X experiment condition X training time X testing time

%% Analyze SDI

path2="C:\Users\user\Desktop\20231018 50T_5_shift decoding\results\SDI\Results4";
zeropoint=0.5;

for a=1:size(CTG_DATA,2);
    
    raw=squeeze(CTG_DATA(:,a,:,:));
    strtemp=sprintf("SDI_%d",a);

    [SDI.si, SDI.di, SDI.time] = calculateSDI(raw,[-0.5 1],zeropoint);
    [SDI.si_se, SDI.di_se, SDI.si_data, SDI.di_data] = bootstrapSDI(raw, [-0.5 1],zeropoint);
    
    fileDir=fullfile(path2,strtemp);
    
    save(fileDir, 'SDI');
    
end