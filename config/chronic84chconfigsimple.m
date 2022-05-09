sessionnum=84;

ncschannels={'p1','p2','p1-p2','p5','p2-p5','p1-p5','pl1','pl1-p1',...
    'cl1','cl1-cl4','cl3','cl4','cl3-cl4'...
    'cl6','cl4-cl6','cl3-cl6','s6','s5','s6-s5','s4','s5-s4','s3','s4-s3',...
    's2','s3-s2',...
    's1','s2-s1','eyed','lickx','pulse'};     

ncschannels={'p1','p2','p1-p2','p5','p2-p5','p1-p5','pl1','pl1-p5',...
    'cl1','cl1-cl4','cl3','cl4','cl3-cl4'...
    'cl6','cl4-cl6','cl3-cl6','s5','s4','s5-s4','s3','s4-s3',...
    'eyed','lickx','pulse'};    %reduced midbrain
ncschannels={'pl1-p5','p1-p5','p2-p5','cl1-cl4','cl3-cl4','cl4-cl6','cl3-cl6',...
    'eyex','eyed','lickx','pulse'};   %chronic65/67


%get info on whether pc or some other system (ie. chunky) to determine dir
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
homedir='Z:';
if ~ispc
    homedir=[filesep 'smbshare'];
end
graserver='inj-monkey2';
fscvdir='patra_fscv2';

paths{1}=fullfile(homedir,graserver,fscvdir,'patra_chronic84_06112018','1dr','cvtotxt',filesep);
paths{2}=fullfile(homedir,graserver,'patra2','2018-06-11_12-22-29',filesep);