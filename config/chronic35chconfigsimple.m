%chronic 35 channels
%CL5 is BAD only CL1-CL4 good CN
ncschannels={'p1-p2','p1','p2','p5','pl1-pl2','pl1','pl2',...
    'cl1','cl1-cl4','cl4','cl5','cl4-cl5',...
    's6','s5','s6-s5','s4','s5-s4','s3','s2','s3-s2','s1','s2-s1',...
    'eyed','eyex','lickx'};          %chronic 31

% no pulse meter installed at this time
%fscv channels are p3, cl6, cl3, pl3
%P1 is recorded by FSCV

ncschannels={'p1','p2','p1-p2','p5','p2-p5','p1-p5','pl1','pl2','pl1-pl2','pl1-p1',...
    'pl1-p5','pl2-p5','pl2-p1'...
    'cl1','cl1-cl4','cl4','cl5','cl4-cl5','cl1-cl5',...
    's5','s4','s5-s4','s3','s4-s3','s2','s3-s2',...
    'eyed','lickx'};         %chronic 31 reduced midbrain


ncschannels={'p2-p5','pl1-pl2',...
    'pl1-p5','pl2-p5',...
    'cl1-cl4','cl4-cl5','cl1-cl5',...
    'eyex','eyed','lickx'};         


%get info on whether pc or some other system (ie. chunky) to determine dir
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
homedir='Z:';
if ~ispc
    homedir=[filesep 'smbshare'];
end
graserver='inj-monkey';
fscvdir='patra_fscv';

paths{1}=fullfile(homedir,graserver,'raw',fscvdir,'patra_chronic35_03082018','1dr','cvtotxt',filesep);
paths{2}=fullfile(homedir,graserver,'raw','patra','2018-03-08_10-36-24',filesep);