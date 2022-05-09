sessionnum=109;

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

paths{1}=fullfile(homedir,graserver,fscvdir,'patra_chronic109_07182018','1dr','cvtotxt',filesep);
paths{2}=fullfile(homedir,graserver,'patra2','2018-07-18_10-43-03',filesep);