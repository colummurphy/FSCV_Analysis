sessionnum=69;

ncschannels={'p1','p2','p1-p2','p5','p2-p5','p1-p5','pl1','pl1-p1',...
    'cl1','cl1-cl4','cl3','cl4','cl3-cl4'...
    'cl6','cl4-cl6','cl3-cl6','s6','s5','s6-s5','s4','s5-s4','s3','s4-s3',...
    's2','s3-s2',...
    's1','s2-s1','eyed','lickx','pulse'};     

ncschannels={'p1','p2','p1-p2','p3','p1-p3','p2-p3','pl2','pl2-p1','pl3',...
    'pl2-pl3','p1-pl3',...
    'cl1','cl4','cl1-cl4',...
    'cl5','cl4-cl5','cl1-cl5','s5','s4','s5-s4','s3','s4-s3',...
    'eyed','lickx','pulse'};    %reduced midbrain
ncschannels={'p1-p3','pl2-p1','pl2-pl3','p1-pl3',...
    'cl1-cl4','cl4-cl5','cl1-cl5',...
    'eyex','eyed','lickx','pulse'};    %reduced midbrain

letterdrive='Y:';
fscvdir='patra_fscv';

paths{1}=fullfile(letterdrive,'data_MIT',fscvdir,'patra_chronic69_05152018','1dr','cvtotxt',filesep);

%paths{1}=fullfile(homedir,graserver,fscvdir,'patra_chronic69_05152018','1dr','cvtotxt',filesep);
%paths{2}=fullfile(homedir,graserver,'patra2','2018-05-15_11-18-17',filesep);
