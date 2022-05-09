%chronic 58 channels

ncschannels={'p1','p2','p3','p1-p2','p2-p3','p1-p3','pl2','pl2-p1','pl3','pl2-pl3',...
    'cl1','cl1-cl4','cl4',...
    'cl5','cl4-cl5','s6','s5','s6-s5','s4','s5-s4','s3','s4-s3','s2','s3-s2',...
    's1','s2-s1','eyed','lickx','pulse'};          %plot ncs channels chronic58 all

ncschannels={'p1-p3','p1-p5','p3-p5',...
    'pl2-p1','cl1-cl4','cl3-cl4','cl1-cl3',...
    'cl4-cl5','cl1-cl5','cl3-cl5',...
    'eyed','lickx','pulse'};    %reduced midbrain  

%ncschannels={'p3','p1-p3','pl3','cl1','cl1-cl4','cl5','cl4-cl5',...
%    's1','s2-s1','eyed','lickx','pulse'}; 

%p3 no good
%p1-p3 no good
%p1 maybe ok, need to check more carefully
%pl2 maybe ok
%cl5 no good, cl4-cl5 no good
%s1 no good
%s2 maybe ok
%s3 good
%s4,s5,s6 no good

%get info on whether pc or some other system (ie. chunky) to determine dir
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
homedir='Z:';
if ~ispc
    homedir=[filesep 'smbshare'];
end
graserver='inj-monkey2';
fscvdir='patra_fscv2';

paths{1}=fullfile(homedir,graserver,fscvdir,'patra_chronic192_12112018','1dr','cvtotxt',filesep);
paths{2}=fullfile(homedir,graserver,'patra2','2018-12-11_14-11-16',filesep);

