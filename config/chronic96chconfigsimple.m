sessionnum=96;

ncschannels={'pl1-p5','p1-p5','p2-p5','cl1-cl4','cl3-cl4','cl4-cl6','cl3-cl6',...
    'eyex','eyed','lickx','pulse'};   %chronic65/67

letterdrive='Y:';
fscvdir='patra_fscv';
paths{1}=fullfile(letterdrive,'data_MIT',fscvdir,'patra_chronic96_06272018','1dr','cvtotxt',filesep); %COPYHERE

%paths{1}=fullfile(homedir,graserver,fscvdir,'patra_chronic96_06272018','1dr','cvtotxt',filesep);
%paths{2}=fullfile(homedir,graserver,'patra2','2018-06-27_11-05-56',filesep);