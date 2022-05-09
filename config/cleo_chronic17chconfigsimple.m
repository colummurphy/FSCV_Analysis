%chronic 46 channels



ncschannels={  'cl4-cl5'...
       'cl2-cl3'...     
  'p3-p4'...
  'p1-p3'...
     'p2-p1'...     
       'eyed'   ...
    'eyex' };

%get info on whether pc or some other system (ie. chunky) to determine dir
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
homedir='Z:';
if ~ispc
    homedir=[filesep 'smbshare'];
end
graserver='inj-monkey';
fscvdir='cleo_fscv';

paths{1}=fullfile(homedir,graserver,'raw',fscvdir,'cleo_chronic17_4ch_04122017','1dr','cvtotxt',filesep);
paths{2}=fullfile(homedir,graserver,'raw','cleo','2017-04-12_11-31-14',filesep);