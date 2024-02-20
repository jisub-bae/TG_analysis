load GSN-HydroCel-128

lay.pos(1:3,:)=[]
lay.width(1:3)=[]
lay.height(1:3)=[]
lay.label(1:3)=[]

cfg.layout=lay;
cfg.method        = 'triangulation'
cfg.neighbourdis=10

[neighbours, cfg] = ft_prepare_neighbours(cfg)

ch={neighbours.neighblabel};

save('target_ch.mat','ch')