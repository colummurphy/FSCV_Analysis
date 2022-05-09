function extractsessioncleo(sessnum,varargin)
%script to run all functions to get compiled trial data lfp/da
%After running syncsigs to get "raw data" to make individual bigreward/smallreward/etc.
%folders
sessid='';
%assume patra recording
if ~isnumeric(sessnum)
    %user gave character also ie. 19b for ephys recording
    sessid=sessnum;
else
    sessid=num2str(sessnum);
end
doreconvert=1;
doauto=1;
docompile=1;
argnum=1;
alltypes=1;
nlxsel=0;           %specify selective chs or get from chronic..config
nlxchsel=[];
chs=[1 2 3 4];
types={'big','small','targetbreak','fixbreak','fixedintervals'};   
types={'big','small','targetbreak','fixbreak'};      

numpcs=3;
nofscv=0;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'reconvertauto'
            %no compile
            doauto=1;
            docompile=0;
            doreconvert=1;
        case 'compile'
            %only run compiletrials
            doreconvert=0;
            doauto=0;
            docompile=1;
        case 'autocompile'
            %run auto and compile
            doauto=1;
            docompile=1;
            doreconvert=0;
        case 'types'
            %user provides trial types to process
            argnum=argnum+1;
            alltypes=0;
            types=varargin{argnum};       
        case 'fscvchs'
            argnum=argnum+1;    %user provides fscv chs selected
            chs=varargin{argnum};
        case 'nlxsel'
            %only reconvert specified nl channels, not reconverting da
            %again
            argnum=argnum+1;
            nlxchsel=varargin{argnum};
        case 'numpcs'
            %# pcs to use in autotrialdir
            argnum=argnum+1;
            numpcs=varargin{argnum};
       case 'nofscv'
            nofscv=1;
    end
    argnum=argnum+1;
end

%get dir with config files
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
%default on putamen pc in lab
configdir='A:\mit\injectrode\experiments\fscv\matlab\analysis\analysis\config\';
if ~ispc
    %chunky dir
    configdir=fullfile(filesep,'home','schwerdt','matlab','analysis','analysis','config',filesep);
end

d=dir(configdir);
cd(configdir);
filenames={d.name};
targfiles=strfind(filenames,['cleo_chronic' sessid]);
processfiles=find(~cellfun(@isempty,targfiles));
targconfigname=filenames{processfiles};
run([configdir targconfigname ]);            %MUST BE RUN FIRST FOR SESSIONNUM
run([configdir 'cleo_map_bipolar' sessid ]);   %run patra_map_bipolar for ch settings
%run config file to get ncschannels & paths

targpath={};
if ~isempty(types)
    %set path names for different trial types folders
    for itype=1:length(types)
        name=types{itype};
        if strcmp(name,'big') || strcmp(name,'small')
            name=[name 'reward'];
        end
        targpath=setfield(targpath,types{itype},...
            fullfile(paths{1}, 'matlab',name,filesep));
    end
end

assignin('base','fcnStatus',targpath)   %store targpath in workspace

trialtypes=fieldnames(targpath);
numtypes=length(trialtypes);

%auto & compile trials
for ii=1:numtypes
    currpath=getfield(targpath,trialtypes{ii});
    disp(currpath)
    if ~isdir(currpath)
        %files/folder not created by sync sigs, skip
        disp(['folder not created, skipping : ' char(10) currpath]);
        continue
    end
    %reconvert directory
    if doreconvert
        disp(['reconvert ' trialtypes{ii}])
        if ~nlxsel
            if ~nofscv
                if sessnum>24
                   reconvertfscv(currpath,ncschannels,'map',csc_map,'split','cleomoves','cleolick37')
                else
                  reconvertfscv(currpath,ncschannels,'map',csc_map,'split','cleomoves')
                end
            else
                if sessnum>24
                   reconvertfscv(currpath,ncschannels,'map',csc_map,'split','cleomoves','cleolick37','nofscv')
                else
                  reconvertfscv(currpath,ncschannels,'map',csc_map,'split','cleomoves','nofscv')
                end
            end
        else
           if ~nofscv
            %only reconvert select nlx chs, no da
                %already previously reconverted everything lese
                %already ran syncsigs again (by default merges everything..
                %so have to run all chs in this step)
                if sessnum>24
                    reconvertfscv(currpath,nlxchsel,'map',csc_map,'split','nlxchsel','cleomoves','cleolick37')
                else
                    reconvertfscv(currpath,nlxchsel,'map',csc_map,'split','nlxchsel','cleomoves')
                end
           else
                if sessnum>24
                    reconvertfscv(currpath,nlxchsel,'map',csc_map,'split','nlxchsel','cleomoves','cleolick37','nofscv')
                else
                    reconvertfscv(currpath,nlxchsel,'map',csc_map,'split','nlxchsel','cleomoves','nofscv')
                end
           end
        end
    end
    pathauto=[currpath(1:end-1) '_pro' filesep];
    if doauto
        if strcmp(trialtypes{ii},'fixedintervals') 
            %no autotrial for fixed intervals
            continue
        end 
        disp(['auto trial ' trialtypes{ii}])
        if ~nofscv
        if sessnum~=14
            autotrialdir(pathauto,'fscvchs',chs,'numpcs',numpcs,'config',['cleo' sessid],'settings','cleo')
        else
           autotrialdir(pathauto,'fscvchs',chs,'numpcs',numpcs,'config',['cleo' sessid],'settings','cleo14')
        end
        else
                      autotrialdir(pathauto,'fscvchs',chs,'numpcs',numpcs,'config',['cleo' sessid],'settings','cleo','nofscv')
        end 
    end
    pathcomp=[pathauto 'analyzed' filesep];
    if docompile
        if strcmp(trialtypes{ii},'fixedintervals')
            %no autotrial for fixed intervals
            continue
        end 
        disp(['compile trials ' trialtypes{ii}])
        if ~nofscv
            compiletrialsdir(pathcomp)
        else
                   compiletrialsdir(pathcomp,'nofscv')
        end 
    end
end

end
