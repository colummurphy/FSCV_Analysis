function [chs,cscmap,paths]=getconfig(subject,sessnum)
%get channel map configuration based on subject and session day recording
sessid=num2str(sessnum);
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
%default on putamen pc in lab
configdir='A:\mit\injectrode\experiments\fscv\matlab\analysis\analysis\config\';
if ~ispc
    %chunky dir
    configdir=fullfile(filesep,'home','schwerdt','matlab','analysis','analysis','config',filesep);
end
%2022 PITT
analysisdir=fullfile(userpath,'fscv','analysis',filesep);    %default analysis dir
configdir=fullfile(analysisdir,'config',filesep);   %path with config files for data directory and channels to be analyzed

d=dir(configdir);
filenames={d.name};
targfiles=strfind(filenames,['chronic' sessid 'chconfigsimple.m']);
if strcmp(subject, 'cleo')
    targfiles=strfind(filenames,['cleo_chronic' sessid '.m']);
end
processfiles=find(~cellfun(@isempty,targfiles));
targconfigname=filenames{processfiles};
run([configdir targconfigname]);  %get ncschannels
run([configdir 'patra_map_bipolar.m']); %get csc_map based on ncschannels
if strcmp(subject,'cleo')
    switch sessid
        case '17'
            run([configdir 'cleo_map_bipolar17']);  %get csc_map based on ncschannels
        case '16'
            run([configdir 'cleo_map_bipolar16']);
        case '25'
                run([configdir 'cleo_map_bipolar25']);
        case '14'
                run([configdir 'cleo_map_bipolar14']);
    end 
end
chs=ncschannels;
cscmap=csc_map;
end