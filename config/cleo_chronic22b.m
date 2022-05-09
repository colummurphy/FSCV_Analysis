%chronic 22b channels

ncschannels={'cl6','cl7','cl6-cl7', 'c5' ,  'c4', 'c5-c4','c3','c3-c5', 'c2',...
   'c2-c4' , 'c1', 'c2-c1', 'c3-c1', 'cl5', 'cl6-cl5','cl5-cl7','p6',...
     'p1','p5-p1','m1','m3', 'm4','m2','eyed'  , 'eyex' };

%get info on whether pc or some other system (ie. chunky) to determine dir
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
homedir='Z:';
if ~ispc
    homedir=[filesep 'smbshare'];
end
graserver='inj-monkey';
fscvdir='cleo_fscv';

paths{1}=fullfile(homedir,graserver,'raw',fscvdir,'cleo_chronic22b_05082017','1dr','cvtotxt',filesep);
paths{2}=fullfile(homedir,graserver,'raw','cleo','2017-05-08_13-35-21',filesep);