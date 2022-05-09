ncschannels={'p3-p5','p2-p5','p1-p5','p1-p3',...
    'pl1-p5','pl1-p1','pl1-p3','pl1-pl2','pl2-p3','pl2-p1','pl2-p5',...
    'pl2-pl3','p1-pl3','p5-pl3','p3-pl3',...
    'cl1-cl4','cl1-cl3','cl3-cl4','cl4-cl5','cl3-cl5',...
    'cl4-cl6','cl3-cl6','cl5-cl6',...
    'eyex','eyey','eyed','lickx','pulse'};    %reduced midbrain

%get info on whether pc or some other system (ie. chunky) to determine dir
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
homedir='Z:';
if ~ispc
    homedir=[filesep 'smbshare'];
end
graserver='inj-monkey2';
fscvdir='patra_fscv2';
paths{1}=fullfile(homedir,graserver,fscvdir,'patra_chronic134c_08292018','1dr','cvtotxt',filesep);
paths{2}=fullfile(homedir,graserver,'patra2','2018-08-29_09-19-19',filesep);