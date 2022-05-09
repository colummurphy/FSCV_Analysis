%chronic 46 channels

ncschannels={     'p6-p1'...
    'pl2-p4'...
    'cl7-p6'...    
    'p4-p6'...
      'eyed'   ...
    'eyex' ...
    'lickx' };

%get info on whether pc or some other system (ie. chunky) to determine dir
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
homedir='Z:';
if ~ispc
    homedir=[filesep 'smbshare'];
end
graserver='inj-monkey';
fscvdir='cleo_fscv';

paths{1}=fullfile(homedir,graserver,'raw',fscvdir,'cleo_chronic25_4ch_05202017','1dr','cvtotxt',filesep);
paths{2}=fullfile(homedir,graserver,'raw','cleo','2017-05-19_11-23-28',filesep);