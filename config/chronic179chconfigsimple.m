%chronic 179 channels
%fscv chs p5, pl3, cl6
sessionnum=179;

ncschannels={'p1','p2','p1-p2','p3','p2-p3','p1-p3','pl1','pl2','pl2-p1','pl2-p3',...
    'pl1-p1','pl2-p2','pl1-p3','pl1-pl2',...
    'cl1','cl1-cl4','cl3','cl4','cl3-cl4'...
    'cl5','cl4-cl5','cl3-cl5','s5','s4','s5-s4','s3','s4-s3',...
    'eyed','lickx','pulse'};    %reduced midbrain

ncschannels={'p2-p3','p1-p3','pl2-p1','pl2-p3',...
    'pl1-p1','pl1-pl2',...
    'cl1-cl4','cl3-cl4'...
    'cl4-cl5','cl3-cl5',...
    'eyex','eyed','lickx','pulse'};    %reduced midbrain

%pl1-p1 good appears modulated pre-fix
%get info on whether pc or some other system (ie. chunky) to determine dir
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
homedir='Z:';
if ~ispc
    homedir=[filesep 'smbshare'];
end
graserver='inj-monkey2';
fscvdir='patra_fscv2';

paths{1}=fullfile(homedir,graserver,fscvdir,'patra_chronic179_11052018','1dr','cvtotxt',filesep);
paths{2}=fullfile(homedir,graserver,'patra2','2018-11-05_10-39-49',filesep);
