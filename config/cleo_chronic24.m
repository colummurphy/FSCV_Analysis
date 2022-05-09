ncschannels={};

%get info on whether pc or some other system (ie. chunky) to determine dir
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
homedir='Z:';
if ~ispc
    homedir=[filesep 'smbshare'];
end
graserver='inj-monkey';
fscvdir='cleo_fscv';

paths{1}=fullfile(homedir,graserver,'raw',fscvdir,'cleo_chronic24_4ch_05122017','1dr','cvtotxt',filesep);
paths{2}=fullfile(homedir,graserver,'raw','cleo','2017-05-12_10-02-49',filesep);