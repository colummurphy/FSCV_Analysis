%chronic 38 channels

ncschannels={'p2','p3','p2-p3','p5','p2-p5','p3-p5','pl1','pl2','pl1-pl2',...
    'cl1','cl4','cl1-cl4','cl5','cl4-cl5','cl1-cl5',...
    's6','s5','s6-s5','s4','s5-s4','s4-s3','s3','s2','s3-s2','s1','s2-s1',...
    'eyed','eyex','lickx','pulse'};          %chronic 38

ncschannels={'p2','p3','p2-p3','p5','p2-p5','p3-p5','pl1','pl2','pl1-pl2',...
    'pl2-p5','pl1-p3','pl2-p3','cl1','cl4','cl1-cl4','cl5','cl4-cl5','cl1-cl5',...
    's5','s4','s5-s4','s4-s3','s3','s2','s3-s2',...
    'eyex','eyey','eyed','lickx','pulse'};  %reduced midbrain

%get info on whether pc or some other system (ie. chunky) to determine dir
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
homedir='Z:';
if ~ispc
    homedir=[filesep 'smbshare'];
end
graserver='inj-monkey';
fscvdir='patra_fscv';

paths{1}=fullfile(homedir,graserver,'raw',fscvdir,'patra_chronic38_03162018','1dr','cvtotxt',filesep);
paths{2}=fullfile(homedir,graserver,'raw','patra','2018-03-16_11-42-37',filesep);