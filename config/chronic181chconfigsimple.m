%chronic 174 and beyond p6 & pm3 connected
sessionnum=181;
ncschannels={'p1','p2','p1-p2','p3','p5','p6','pm3','p2-p3','p1-p3','p1-p5','p3-p5',...
    'pl1','pl2','pl2-p1','pl2-p3','pl1-p5','p5-pm3','p3-pm3','pl2-pm3','p5-p6',...
    'pl1-p1','pl2-p2','pl1-p3','pl1-pl2',...
    'cl1','cl1-cl4','cl3','cl4','cl3-cl4'...
    'cl5','cl4-cl5','cl3-cl5','cl1-cl5','s5','s4','s5-s4','s3','s4-s3',...
    'eyex','eyed','lickx','pulse'};    %reduced midbrain

ncschannels={'p1-p3','p3-p5','pl2-p3','p5-pm3',...
    'p5-p6','pl1-pl2','pl1-p3',...
    'cl1-cl4','cl3-cl4','cl1-cl5',...
    'eyex','eyed','lickx','pulse'};    %reduced midbrain


%get info on whether pc or some other system (ie. chunky) to determine dir
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
homedir='Z:';
if ~ispc
    homedir=[filesep 'smbshare'];
end
graserver='inj-monkey2';
fscvdir='patra_fscv2';

paths{1}=fullfile(homedir,graserver,fscvdir,'patra_chronic181_11072018','1dr','cvtotxt',filesep);
paths{2}=fullfile(homedir,graserver,'patra2','2018-11-07_08-47-19',filesep);

