function plotmultiple(sessnums, varargin)
%plot xinfo multiple sessions/days
%get dir with config files
argnum=1;
sesstypes={'bigreward','smallreward','targetbreak'};
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
xinfosneg={};       %da decreases trials
nlxchs={};
commonnlx={};           %get lfp chs shared across sessions
trialgrps=[];
binfos={};
%get data/concatenate multiple sessions
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
    load([sesspaths{isess} 'binfo.mat']);
    xinfos{isess}=xinfo;
    binfos{isess}=binfo;
    if isess==2
        commonnlx=intersect(nlxchs{isess-1},nlxchs{isess});
    end
    if isess>2
        commonnlx=intersect(commonnlx,nlxchs{isess});
    end
    for itype=1:length(sesstypes)
        currpath=fullfile(paths{1},'matlab',[sesstypes{itype} '_pro'],'analyzed',filesep);
        load([currpath 'trialtypes']);
        trialgrps=setfield(trialgrps,{isess},'trialinfo',{itype},'trialtypes',trialtypes);
    end
end

plotparam.trialgrps=trialgrps;
lfpids=(~contains(commonnlx,'eyed') & ~contains(commonnlx,'lickx') & ~contains(commonnlx,'eyex') & ~contains(commonnlx,'pulse'));
plotparam.savepath=[sesspaths{isess},'mult',filesep];
if ~isdir(plotparam.savepath)
mkdir(plotparam.savepath);
end
plotparam.lfpchs=commonnlx(lfpids);        
plotparam.dasites=dasites;
plotparam.sessnums=sessnums;

params_mag={'damax','damin',...
    'lfpmax','lfpmin','lfppostmax','zlagcoef'};
params_xcov={'maxprelagts',...
    'minprelagts','maxpostlagts','minpostlagts','maxprecoef',...
    'minprecoef','maxpostcoef','minpostcoef'};
params_timing_damax={'delt_lfpmax_damax','delt_lfpfall_damax'...
    'delt_lfpmin_damax','delt_lfprise_damax',...
    'delt_lfppostmax_damax'};
params_timing_darise={'delt_lfpmax_darise','delt_lfpfall_darise'...
    'delt_lfpmin_darise','delt_lfprise_darise',...
    'delt_lfppostmax_darise'};   
params_timing_dafall={'delt_lfpmax_dafall','delt_lfpfall_dafall'...
    'delt_lfpmin_dafall','delt_lfprise_dafall',...
    'delt_lfppostmax_dafall'};
 params_timing_damin={'delt_lfpmax_damin','delt_lfpfall_damin'...
    'delt_lfpmin_damin','delt_lfprise_damin',...
    'delt_lfppostmax_damin'};

ygrp{1}=params_timing_damin;
ygrp{2}=params_timing_darise;
ygrp{3}=params_timing_damax;
ygrp{4}=params_mag;
ygrp{5}=params_xcov;

ynam{1}='timing_damin';
ynam{2}='timing_darise';
ynam{3}='timing_damax';
ynam{4}='mag';
ynam{5}='xcov';

yvar.ygrp=ygrp;
yvar.ynam=ynam;
plotxmultsessfoc(xinfos,plotparam,{'all'},yvar,'daall');
%{
%temp comment 1/21/2019
plotxmultsess(xinfos,plotparam,'daall');
plotxmultsess(xinfos,plotparam);
plotxmultsess(xinfos,plotparam,'daneg');        %plot da neg

plotxbehmultsess(xinfos,binfos,plotparam,'daall');

plotxbehmultsess(xinfos,binfos,plotparam);
plotxbehmultsess(xinfos,binfos,plotparam,'daneg');

%}

%plot multiple da neg
%plot corrs of behavior to xinfo variables
end