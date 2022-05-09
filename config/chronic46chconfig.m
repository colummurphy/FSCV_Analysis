%chronic 46 channels

ncschannels={'p1','p2','p1-p2','pl1','pl3','pl1-pl3','pl1-p1','p1-pl3',...
    'cl1','cl4','cl1-cl4',...    
    'cl5','cl1-cl5','cl4-cl5','cl6','cl5-cl6','cl4-cl6',...
    's6','s5','s6-s5','s4','s5-s4','s4-s3','s3','s2','s3-s2','s1','s2-s1',...
    'eyex','eyed','lickx','pulse'};          

%get info on whether pc or some other system (ie. chunky) to determine dir
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
homedir='Z:';
if ~ispc
    homedir=[filesep 'smbshare'];
end
graserver='inj-monkey2';
fscvdir='patra_fscv2';

paths{1}=fullfile(homedir,graserver,fscvdir,'patra_chronic46_03302018','1dr','cvtotxt',filesep);
paths{2}=fullfile(homedir,graserver,'patra2','2018-03-30_11-11-16',filesep);