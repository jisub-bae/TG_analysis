%% Load data

path="C:\Users\user\Desktop\20231018 50T_5_shift decoding\data";
path2='C:\Users\user\Desktop\20231018 50T_5_shift decoding\results\CTG statistic';

for ds={'occipital','parietal','frontal'};
    
    cd(path)
    total_data=[];
    for a=0:4;

        name=sprintf("ctg_%s_%d.mat",ds{1},a);
        load(name)
        total_data=cat(4,total_data,ctg);

    end 
    
    CTG_DATA=permute(total_data,[1,4,2,3]);
    name2=sprintf("CTG_DATA_%s.mat",ds{1});
    
    cd(path2)
    save(name2,'CTG_DATA','-v7.3')
    
end