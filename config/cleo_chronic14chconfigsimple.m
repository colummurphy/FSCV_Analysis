%chronic 46 channels

ncschannels={   
    'c5-p4'...
     'c5-cl3'...
     'p2-p0' ...
       'p4-p2' ...
      'cl4-cl3'...
    'cl4-c5' ...       
     'eyed'   ...
    'eyex'  
    };

%get info on whether pc or some other system (ie. chunky) to determine dir
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
homedir='Z:';
if ~ispc
    homedir=[filesep 'smbshare'];
end
graserver='inj-monkey';
fscvdir='cleo_fscv';

paths{1}=fullfile(homedir,graserver,'raw',fscvdir,'cleo_chronic14_4ch_04052017','1dr3','cvtotxt',filesep);
paths{2}=fullfile(homedir,graserver,'raw','cleo','2017-04-05_12-02-18',filesep);