%chronic 19b channels

ncschannels={'c5','c4','c3','c2','c1','cl6','cl5','cl4','cl3','cl2','p6','p3',...
    'cl6-c5', 'c4-c5', 'c2-c4',  'c1-c3', 'cl6-cl5', 'cl4-cl5', ...
    'cl3-cl4', 'cl2-cl3', 'c2-cl6', 'p6-p3', 'm1', 'm3', 'm4', 'm2', 'eyed','eyex' };

%get info on whether pc or some other system (ie. chunky) to determine dir
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
homedir='Z:';
if ~ispc
    homedir=[filesep 'smbshare'];
end
graserver='inj-monkey';
fscvdir='cleo_fscv';

paths{1}=fullfile(homedir,graserver,'raw',fscvdir,'cleo_chronic19b_04202017','1dr','cvtotxt',filesep);
paths{2}=fullfile(homedir,graserver,'raw','cleo','2017-04-20_08-54-07',filesep);