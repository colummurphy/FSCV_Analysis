function plotmultiple(sessnums, varargin)
%plot xinfo multiple sessions/days
%get dir with config files
argnum=1;
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'dasites'
            argnum=argnum+1;
            dasites=varargin{argnum};
    end
    argnum=argnum+1;
end
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
%default on putamen pc in lab
configdir='A:\mit\injectrode\experiments\fscv\matlab\analysis\analysis\config\';
if ~ispc
    %chunky dir
    configdir=fullfile(filesep,'home','schwerdt','matlab','analysis','analysis','config',filesep);
end
d=dir(configdir);
filenames={d.name};

sesspaths={};
xinfos={};
xinfos_burstgrps={};
nlxchs={};
commonnlx={};           %get lfp chs shared across sessions
for isess=1:length(sessnums)
    sessnum=sessnums(isess);
    sessid=num2str(sessnum);
    targfiles=strfind(filenames,['chronic' sessid 'chconfigsimple.m']);
    processfiles=find(~cellfun(@isempty,targfiles));
    targconfigname=filenames{processfiles};
    run([configdir targconfigname]);        %get paths{1}
    nlxchs{isess}=ncschannels;
    sesspaths{isess}=fullfile(paths{1},'matlab','xsess',filesep);
    load([sesspaths{isess} 'xinfo.mat']);
    xinfos{isess}=xinfo;
    load([sesspaths{isess} 'xinfo_burstgrps.mat']);
    xinfos_burstgrps{isess}=xinfo;
    if isess==2
        commonnlx=intersect(nlxchs{isess-1},nlxchs{isess});
    end
    if isess>2
        commonnlx=intersect(commonnlx,nlxchs{isess});
    end
end
lfpids=(~contains(commonnlx,'eyed') & ~contains(commonnlx,'lickx') & ~contains(commonnlx,'pulse'));
plotparam.trialtypes.names={'all'  'sleft'  'sright'  'aftersuccess'  'afterfail'  'phase1'  'phase4'};
plotparam.savepath=sesspaths{isess};
plotparam.lfpchs=commonnlx(lfpids);        
plotparam.dasites=dasites;
plotparam.sessnums=sessnums;
%plot mean defined lags for big vs small vs targ
%plot variance lag-defined cross-covariance wfs big vs. small vs targ
plotxsess(xinfos,plotparam);

plotparam.trialtypes.names={'burstgrp','noburstgrp'};
plotxsess(xinfos_burstgrps,plotparam,'label','bursts_');
plotxsessb(xinfos_burstgrps,plotparam);

%plot within session type (ie only big)