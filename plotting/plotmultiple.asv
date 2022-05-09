function plotmultiple(sessnums, varargin)
%plot xinfo multiple sessions/days
%get dir with config files
argnum=1;
sesstypes={'bigreward','smallreward','targetbreak'};
sesstypes={'bigreward','smallreward'};
sessnums2=[];
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'dasites'
            argnum=argnum+1;
            dasites=varargin{argnum};
        case 'cleo'
            argnum=argnum+1;
            sessnums2=varargin{argnum};
    end
    argnum=argnum+1;
end
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
%default on putamen pc in lab --> CHANGED 2021 for Pitt
configdir='C:\Users\putamen\Documents\MATLAB\fscv\analysis\config\';

d=dir(configdir);
filenames={d.name};
homed=fullfile('Y:', 'data_MIT', 'analysis', 'patra', filesep);
homed2=fullfile('Y:', 'data_MIT', 'analysis', 'cleo', filesep);

sesspaths={};
xinfosneg={};       %da decreases trials
dachs={};
nlxchs={};
commonnlx={};           %get lfp chs shared across sessions
trialgrps=[];
corrs={};
corrsbeh={};
datm={};
datmcomb={};
betatm={};
betatmcomb={};
behtm={};
behtmcomb={};
xinfotemp=[];
xinfos=[];
xbinfos=[];
respgrouped={};
resplrgrouped={};
binfos=[];
trlists=[];
sessnumsall=[sessnums sessnums2];
subj='patra';
%get data/concatenate multiple sessions
for isess=1:length(sessnumsall)
    sessnum=sessnumsall(isess);
    if isess>length(sessnums)
        subj='cleo';
        homed=homed2;
    end
    sessid=num2str(sessnum);
    disp(['loading sess # ' sessid]);
    targfiles=strfind(filenames,['chronic' sessid 'chconfigsimple.m']);
    processfiles=find(~cellfun(@isempty,targfiles));
    if strcmp(subj,'cleo')
        targfiles=strfind(filenames,['cleo_chronic' sessid 'chconfigsimple.m']);    
        processfiles=find(~cellfun(@isempty,targfiles));
    end
    if isempty(processfiles) && strcmp(subj,'cleo')
        %get chconfig file instead of chconfigsimple (for cleo)
        targfiles=strfind(filenames,['cleo_chronic' sessid '.m']);    
        processfiles=find(~cellfun(@isempty,targfiles));
    end    
    targconfigname=filenames{processfiles};
    run([configdir targconfigname]);        %get paths{1}
    pathname=fullfile(homed, ['chronic' sessid], filesep);
    nlxchs{isess}=ncschannels;
    %load trlist
    load([pathname 'trlist.mat']);
    trlists(isess).list=trlist;
    trlists(isess).sessid=sessid;
    
    %load xinfo/binfo for each session having all pertinent averaged traces
    %for signals, and defined metrics
    sesspaths{isess}=fullfile(pathname,'all','xsess',filesep);
    load([sesspaths{isess} 'binfo.mat']);
    for ib=1:length(binfo)
    binfo=setfield(binfo,{ib},'sessionid',sessid);
    end
    binfos=[binfos binfo];
    
    load([sesspaths{isess} 'xinfo.mat']);
    load([sesspaths{isess} 'xbinfo.mat']);
    dachs{isess}=unique({xinfo(1:end).siteda});
    [dapair,lfppair]=getsitepairs(dachs{isess},subj);    %get corresponding target lfp pairs
    alllfps=[xinfo.sitelfp];
    xinfotemp=xinfo;
    targrows=~contains({xinfo.event},'outcome');
    xinfotemp=xinfo(targrows);
    xinfotemp=rmfield(xinfotemp,'xcovda');
    if ~isempty(alllfps)
    targrows=contains({xinfo.sitelfp},lfppair);
    %xinfotemp=xinfo(targrows);          %get data for targeted lfp pairs
    else
        disp('no lfps');
    end
    if isfield(xinfotemp,'daneg')
    xinfotemp=rmfield(xinfotemp,'dareb');
    xinfotemp=rmfield(xinfotemp,'daneg');
    xinfotemp=rmfield(xinfotemp,'dapos');
    end
    xinfotemp=rmfield(xinfotemp,'sessionperiod');
    %put unique site identifiers in xinfo/xbinfo based on sessnum/orig da site
    sites=getsites(sessnum,dachs{isess},subj);
    %remove data corresponding to unidentified sites (ie. sites with x in
    %xcel)
    badsites=dachs{isess}(~ismember(dachs{isess},{sites.probeid}));     
    targrowsgood=~contains({xinfotemp.siteda},badsites);
    xinfotemp=xinfotemp(targrowsgood);
    xsites=find(contains({sites.site},'x'));
    badsites2={sites(xsites).probeid};
    targrowsgood=~contains({xinfotemp.siteda},badsites2);
    xinfotemp=xinfotemp(targrowsgood);
    for ix=1:length(xinfotemp)
        xinfotemp=setfield(xinfotemp,{ix},'sessionid',sessid);
        curda=getfield(xinfotemp,{ix},'siteda');
        cursite=sites(contains({sites.probeid},curda)).site;
        if ~isempty(cursite)
            xinfotemp=setfield(xinfotemp,{ix},'dasite',cursite);
        end            
    end
    
    
    %do same for xbinfo  
    xbinfotemp=[];
    if ~isempty(xbinfo)
        %remove empty cells first (chronic 25 cleo)
        empties=cellfun('isempty',{xbinfo.daall});
        xbinfo(empties)=[];
    targrowsb=~contains({xbinfo.sitelfp},{'eyex'});
    xbinfotemp=xbinfo(targrowsb);
    if isfield(xbinfo,'daneg')
    xbinfotemp=rmfield(xbinfotemp,'dareb');
    xbinfotemp=rmfield(xbinfotemp,'daneg');
    xbinfotemp=rmfield(xbinfotemp,'dapos');
    end
    xbinfotemp=rmfield(xbinfotemp,'sessionperiod');
        xbinfotemp=rmfield(xbinfotemp,'freq');

        %remove data corresponding to unidentified sites (ie. sites with x in
    %xcel)
    targrowsgood=~contains({xbinfotemp.siteda},badsites);
    xbinfotemp=xbinfotemp(targrowsgood);
    for ix=1:length(xbinfotemp)
        xbinfotemp=setfield(xbinfotemp,{ix},'sessionid',sessid);
        curda=getfield(xbinfotemp,{ix},'siteda');
        cursite=sites(strcmp({sites.probeid},curda)).site;
        xbinfotemp=setfield(xbinfotemp,{ix},'dasite',cursite);
    end
    end

    if isess==1
        xinfos=xinfotemp;
        xbinfos=xbinfotemp;
    else
        xinfos=[xinfos xinfotemp];
        xbinfos=[xbinfos xbinfotemp];
    end
       
    %load corr matrices for big/small types for each session &
    %selective task-modulated signals
    sesspathscor{isess}=fullfile(pathname,'all','corr',filesep);
    if isdir(sesspathscor{isess})
    load([sesspathscor{isess} 'corr-bigsmall-sametype.mat']);
    corrs{isess}=corrmat;
    datmcomb{isess}=datmall;
    datm{isess}=datms;
    betatmcomb{isess}=betatmall;    
    %only keep beta values for selected lfp pairs above
    betatemp={};
    if strcmp(subj,'cleo')
    for iit=1:length(betatms)
        countbeta=1;
        for iix=1:length(betatms{iit})
            betasite=betatms{iit}{iix}.site;
            targeted=find(strcmp(lfppair,betasite));
            if ~isempty(targeted)
                %site match
                betatemp{iit}{countbeta}=betatms{iit}{iix};
                countbeta=countbeta+1;
            end
        end
    end
    betatms=betatemp;
    end
    betatm{isess}=betatms;
    if exist('rewexp')>0
    respgrouped{isess}=rewexp;
    resplrgrouped{isess}=sideexp;
    end
    %{
    %don't use any ways 5/2019
    load([sesspathscor{isess} 'corr-bigsmall-beh-sametype.mat']);
    behtmcomb{isess}=behtmall;
    corrsbeh{isess}=corrmat;
    behtm{isess}=behtms;
    %}
    end
    
    if isess==2
        commonnlx=intersect(nlxchs{isess-1},nlxchs{isess});
    end
    if isess>2
        commonnlx=intersect(commonnlx,nlxchs{isess});
    end
    for itype=1:length(sesstypes)
        currpath=fullfile(pathname,[sesstypes{itype}],filesep);
        load([currpath 'trialtypes']);
        trialtypesall=trialtypes;
        if exist('trialtypeslfp')>0
            trialtypesall=[trialtypes trialtypeslfp];
        end
        %trialgrps=setfield(trialgrps,{isess},'trialinfo',{itype},'trialtypes',trialtypes);
        if isfield(trialtypes,'site')
            %new trialtypes, different across channels, merge
            % all channels since anyways selective in actual data
            %variables eg. xinfo & datm...
            temptrialtypes=trialtypesall;
            trialtypes=[];
            for iss=1:length(temptrialtypes)
                for inn=1:length(temptrialtypes(iss).names)
                    if iss==1
                        trialtypes.names{inn}=temptrialtypes(iss).names{inn};
                        trialtypes.nums{inn}=temptrialtypes(iss).nums{inn};
                    else
                        trialtypes.nums{inn}=sort(unique([trialtypes.nums{inn} temptrialtypes(iss).nums{inn}]));
                    end
                end
            end
        end
        trialgrps(isess).trialinfo(itype)=trialtypes;
    end
    trialgrps(isess).sessid=sessid;
%}
end
    
count=1;
clfp=1;
sitesda={};
siteslfp={};
dasites=unique({xinfos(1:end).siteda});
for is=1:length(xinfos)
    if str2num(xinfos(is).sessionid)<30
        xinfos(is).subj='cleo';
    else
        xinfos(is).subj='patra';
    end
end
patids=find(contains({xinfos.subj},'patra'));
cleoids=find(contains({xinfos.subj},'cleo'));
if ~isempty(patids) && ~isempty(cleoids)
    homed=fullfile('Y:', 'data_MIT', 'analysis', 'cleopatra', filesep);
    if ~isdir(homed)
        mkdir(homed)
    end
end
siteslfp=unique({xinfos(patids).sitelfp});
emptylfp=strcmp(siteslfp,'');
siteslfp=siteslfp(~emptylfp);
   plotparam.sessnums=sessnumsall;
   plotparam.sessnumsp=sessnums;
   plotparam.sessnumsc=sessnums2;
plotparam.dasites=dasites;
plotparam.trialgrps=trialgrps;
lfpids=(~contains(commonnlx,'eyed') & ~contains(commonnlx,'lickx') & ~contains(commonnlx,'eyex') & ~contains(commonnlx,'pulse'));
lfpids2=(~contains(siteslfp,'eyed') & ~contains(siteslfp,'lickx') & ~contains(siteslfp,'eyex') & ~contains(siteslfp,'pulse') & ~contains(siteslfp,'p2') & ~contains(siteslfp,'p6') & ~contains(siteslfp,'pm3'));

plotparam.savepath=[homed,'mult',filesep];
if ~isdir(plotparam.savepath)
mkdir(plotparam.savepath);
end
plotparam.lfpchs=siteslfp(lfpids2);      
plotparam.lfpchsall=siteslfp(lfpids2);
%plotparam.dasites=dasites;

%assign lfp site pairs to all da chs in patra % REMOVE IRRELEAANT lfpchs
%from plotparam before this
sitelists=assignpairs(sessnums,xinfos,plotparam);

trorgsimp=maketrorg;
trorg=maketrorg2;
trtable=gettrtable(trlists,trorg);
trtablesimp=gettrtable(trlists,trorgsimp);
trorgsimp2=maketrorg4;
trtablesimp2=gettrtable(trlists,trorgsimp2);
sessnums2=[67 68 69 75 79 80 83 84 91 92 94 95 96 100 102 109 113 114 127 179 181];

%5/20/20
%get tabdata for all windows, and cormatrix again..for specific task conditions
starts = 0.2:0.6:3.6;
ends=0.8:.6:3.8;
%starts2 = repmat(starts', length(ends));
count=1;
xmats={};
for it=1:length(starts)
    for iit=1:length(ends)
        if ends(iit)-starts(it)>0
            xmats{count}=[starts(it) ends(iit)];
            count=count+1;
        end
    end
end
[tabdata_all,negdata,mdata]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xmats,'sitelists',sitelists,'cormatrix');              
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xmats,'getcondition',{'left'},'sitelists',sitelists,'cordata',tabdata_all,'cormatrix');  
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xmats,'getcondition',{'right'},'sitelists',sitelists,'cordata',tabdata_all,'cormatrix');  
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xmats,'getcondition',{'big'},'sitelists',sitelists,'cordata',tabdata_all,'cormatrix');  
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xmats,'getcondition',{'small'},'sitelists',sitelists,'cordata',tabdata_all,'cormatrix');  
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xmats,'getcondition',{'big','left'},'sitelists',sitelists,'cordata',tabdata_all,'cormatrix');  
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xmats,'getcondition',{'big','right'},'sitelists',sitelists,'cordata',tabdata_all,'cormatrix');  
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xmats,'getcondition',{'small','left'},'sitelists',sitelists,'cordata',tabdata_all,'cormatrix');  
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xmats,'getcondition',{'small','right'},'sitelists',sitelists,'cordata',tabdata_all,'cormatrix');  

%5/20/20 do Tcov (ie. Fig. 2B) for specific task conditions
outdata=gettcov(sessnums,xinfos,datm,plotparam,[-3 4],{'big','small','all'});
outdata=gettcov(sessnums,xinfos,datm,plotparam,[-3 4],{'small','all'});
outdata=gettcov(sessnums,xinfos,datm,plotparam,[-3 4],{'big','all'});
outdata=gettcov(sessnums,xinfos,datm,plotparam,[-3 4],{'big','small','left'});
outdata=gettcov(sessnums,xinfos,datm,plotparam,[-3 4],{'big','small','right'});
outdata=gettcov(sessnums,xinfos,datm,plotparam,[-3 4],{'big','left'});
outdata=gettcov(sessnums,xinfos,datm,plotparam,[-3 4],{'big','right'});

outdata3c=gettcov([113],xinfos,datm,plotparam,[-3 4],{'big','small','all'},'lfps',{'pl2-pl3'},'dasites',{'p5'}); %fig 3C
outdata3b=gettcov([83],xinfos,datm,plotparam,[-3 4],{'big','small','all'},'lfps',{'cl4-cl6'},'dasites',{'cl5'}); %fig 3b
outdata5b=gettcov([83],xinfos,datm,plotparam,[-3 4],{'big','small','all'},'lfps',{'cl3-cl6'},'dasites',{'cl5'}); %fig 5b
outdata5c=gettcov([94],xinfos,datm,plotparam,[-3 4],{'big','small','all'},'lfps',{'pl1-p5'},'dasites',{'pl3'}); %fig 5c

[tabdata,sitesigs]=corsmultiple([67 68 69 75 83],plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','getcondition',{'left'},'sitelists',sitelists,'win',{[0.2 0.8]},'hist','constrainpairs');  
save('dabetacordata_selected_earlybeta_contraonly','tabdata','sitesigs')

%5/25/20 REDO above, instead of "getconditon" which is for histogram after
%tab data already stored and transferred. use "condition"....
[tabdata,sitesigs]=corsmultiple([67 68 69 75 83],plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','condition','left','sitelists',sitelists,'win',{[0.2 0.8]},'hist','constrainpairs');  
save('dabetacordata_selected_earlybeta_contraonly_rev2','tabdata','sitesigs')

%05/19/20, correlations DA vs beta for fixed task conditions, had already
%done this below, but now safe all p values, r values make histogram
%pos/neg, MAKE SURE TO LOAD CORRECT CONDITION DATA
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','condition','left','sitelists',sitelists,'cordata',tabdata,'hist');  
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','condition','right','sitelists',sitelists,'cordata',tabdata,'hist');  
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','condtypes',{'big','big','big'},'sitelists',sitelists,'cordata',tabdata,'hist');  
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','condtypes',{'small','small','small'},'sitelists',sitelists,'cordata',tabdata,'hist');  


%get tabdata cors for windows showing positive correlations as found in
%cormatrix_multsess folder
xts={  [0.2 1.4],[0.2 2], [0.2 2.6],[0.2 3.8], [0.8 1.4], [0.8 2.6], [0.8 3.2], [1.4 2.6],[1.4 3.8], [2 2.6],[2 3.2], [2 3.8],  [2.6 3.2],[2.6 3.8], [3.2 3.8]};
[tabdata_all,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'sitelists',sitelists);  
save('dabetacordata_targ_multiwin_posbias','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','getcondition',{'left'},'sitelists',sitelists,'cordata',tabdata_all,'hist');  
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','getcondition',{'right'},'sitelists',sitelists,'cordata',tabdata_all,'hist');  
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','getcondition',{'big'},'sitelists',sitelists,'cordata',tabdata_all,'hist');  
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','getcondition',{'small'},'sitelists',sitelists,'cordata',tabdata_all,'hist');  
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','getcondition',{'big','left'},'sitelists',sitelists,'cordata',tabdata_all,'hist');  
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','getcondition',{'big','right'},'sitelists',sitelists,'cordata',tabdata_all,'hist');  




[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'condtypes',{'big','big','big'},'sitelists',sitelists);  
save('dabetacordata_targ_multiwin_bigonly','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'condtypes',{'small','small','small'},'sitelists',sitelists);  
save('dabetacordata_targ_multiwin_smallonly','tabdata','sitesigs')

%find frequency distribution of beta
fftdata=getbetaspecmulti(plotparam,xinfos,binfos);
save('beta_spec_data_all','fftdata')
fftdata2=getbetaspecmulti(plotparam,xinfos,binfos);

%04/29/20, final history dependence of arousal/behavior
plotbehsummary(xbinfos,binfos,plotparam,'norm','mean','ttypes',{{'big','fail'},{'big','succ'}},'event','targ','beh',{'lick','pulse','eye','lrt','rrt'},'wins',{[0 2],[0 4],[0.2 .8],[0 0],[0 0]},'patra','norm','simple','plotlines');
plotbehsummary(xbinfos,binfos,plotparam,'norm','mean','ttypes',{{'big','ersm'},{'big','erbig'},{'big','fail'},{'big','succ'}},'event','targ','beh',{'lick','pulse','eye','lrt','rrt'},'wins',{[0 2],[0 4],[0.2 .8],[0 0],[0 0]},'patra','norm','simple','plotlines','normall');
plotbehsummary(xbinfos,binfos,plotparam,'norm','mean','ttypes',{{'big','ersm'},{'big','erbig'},{'big','fail'},{'big','succ'}},'event','targ','beh',{'lick','pulse','eye','lrt','rrt'},'wins',{[0 2],[0 4],[0.2 .8],[0 0],[0 0]},'patra','norm','simple','plotlines');
plotbehsummary(xbinfos,binfos,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'event','targ','beh',{'lick','pulse','eye','lrt','rrt'},'wins',{[0 2],[0 4],[0.2 .8],[0 0],[0 0]},'patra','norm','simple','plotlines');
plotbehsummary(xbinfos,binfos,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'event','targ','beh',{'lick','pulse','eye','lrt','rrt'},'wins',{[0 2],[0 4],[0.2 .8],[0 0],[0 0]},'patra','norm','simple','plotlines','normall');

%04/30 multiple regression with p values for each variable, targeted
%towards sites/sessions with sig neg cor
[tabdata,negdata,mdata]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',{[0.2 0.8]},'sitelists',sitelists,'multiregress');              

%04/30 2d matrix correlation
starts = 0.2:0.6:3.6;
ends=0.8:.6:3.8;
%starts2 = repmat(starts', length(ends));
count=1;
xmats={};
for it=1:length(starts)
    for iit=1:length(ends)
        if ends(iit)-starts(it)>0
            xmats{count}=[starts(it) ends(iit)];
            count=count+1;
        end
    end
end
%first run to get tabdata, saved 04/30 'corrdata_alldavsbetawins'
[tabdata,negdata,mdata]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xmats,'sitelists',sitelists,'cormatrix');              
%already saved data
[tabdata,negdata,mdata]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xmats,'sitelists',sitelists,'cormatrix','cordata',tabdata);              
%do multiple regression again on those with significant negative
%correlations according to best window for CN (.2-3.2) and putamen (.2-2.6)
%determined by above cormatrix
%05/01 DG: shoot down null hypothesis that da vs beta corr can be fully
%explained by bhe variables by showing that correlations exist between
%residuals of DA vs beh VS Beta vs beh
[tabdata,negdata,mdata]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',{[0.2 3.2], [.2 2.6]},'sitelists',sitelists,'multiregress','cordata',tabdata,'testcor');              
allpvals=[mdata.pvals];
%find p vals of beta coefficients of beta vs da that remain siginificant
%after multiple regression, 17/33 total sites, 51.5%, 5 significant lic (3 of which still
%davsbetacorr), 4 sig rrstd (2 still davsbet), 3 trt (1 still davsbet), 3 type (3 still davsbet), 2 side (1 still davsbet), 8 trialnum (3 still davsbet)
length(find(allpvals(1,:)<=0.05))
length(find(allpvals(2,:)<=0.05))   %lick
length(find(allpvals(3,:)<=0.05))   %rrstd
length(find(allpvals(4,:)<=0.05))   %trt
length(find(allpvals(5,:)<=0.05))   %type
length(find(allpvals(6,:)<=0.05))   %side
length(find(allpvals(7,:)<=0.05))   %trialnum

%05/02 use cormatrix results as relevant time windows to look at da vs beta
%corr as function of different task conditions, including 0.2-0.8 best for
%task modulation of individual signals..(but not amongst top 8 time windows for da vs beta)
%use top 4 windows for each region
xts={[0.2-0.8],[0.2 2.6],[0.8 2], [0.2 2], [0.2 3.2], [0.2 1.4],[0.8 1.4],       [0.2 3.8], [0.2 2], [1.4 3.8]};
xts={[0.2 0.8],[0.2 2.6],[0.8 2], [0.2 2], [0.2 3.2], [0.2 1.4],       [0.2 3.8], [0.2 2], [1.4 3.8]};
xts={[0.2 0.8],[0.2 1.4],[0.2 2],[0.8 2],[0.2 2.6],[0.2 3.2],[0.2 3.8],[1.4 3.8]};

%JUST DO CORMATRIX FOR ALL
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'sitelists',sitelists);  
save('dabetacordata_targ_multiwin','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'condition','left','sitelists',sitelists);  
save('dabetacordata_targ_multiwin_leftonly','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'condition','right','sitelists',sitelists);  
save('dabetacordata_targ_multiwin_rightonly','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'condtypes',{'big','big','big'},'sitelists',sitelists);  
save('dabetacordata_targ_multiwin_bigonly','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'condtypes',{'small','small','small'},'sitelists',sitelists);  
save('dabetacordata_targ_multiwin_smallonly','tabdata','sitesigs')


[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'sitelists',sitelists,'constrainpairs');  
save('dabetacordata_targ_multiwin_paired','tabdata','sitesigs')

maxfft{1}=plotspectrosel(92,'cl3-cl6',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});


maxfft={};
maxfft{1}=plotspectrosel(92,'cl3-cl6',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{2}=plotspectrosel(92,'cl4-cl6',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{3}=plotspectrosel(92,'cl1-cl4',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{4}=plotspectrosel(92,'cl1-cl3',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{5}=plotspectrosel(92,'pl1-p5',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
%p1-p5 too flat in 92
%maxfft{6}=plotspectrosel(92,'p1-p5',[],'targ',plotparam,'win',[-3 4],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big','small'});
maxfft{6}=plotspectrosel(68,'p1-p3',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{7}=plotspectrosel(68,'pl2-p1',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{8}=plotspectrosel(68,'p1-pl3',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{9}=plotspectrosel(68,'cl4-cl5',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{10}=plotspectrosel(113,'p1-p3',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{11}=plotspectrosel(113,'pl2-p1',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{12}=plotspectrosel(113,'p1-pl3',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{13}=plotspectrosel(113,'pl2-pl3',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{14}=plotspectrosel(113,'cl1-cl4',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{15}=plotspectrosel(179,'pl2-p3',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{16}=plotspectrosel(179,'pl2-p1',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{17}=plotspectrosel(179,'pl1-pl2',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{18}=plotspectrosel(179,'pl1-p1',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{19}=plotspectrosel(179,'p1-p3',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
%maxfft{20}=plotspectrosel(179,'cl1-cl4',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
%maxfft{21}=plotspectrosel(179,'cl3-cl4',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{22}=plotspectrosel(94,'cl3-cl4',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{23}=plotspectrosel(94,'cl4-cl6',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{24}=plotspectrosel(94,'cl3-cl6',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
%maxfft{25}=plotspectrosel(94,'pl1-p5',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{26}=plotspectrosel(100,'cl4-cl5',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{27}=plotspectrosel(100,'cl1-cl5',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{28}=plotspectrosel(100,'cl1-cl4',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
%maxfft{29}=plotspectrosel(100,'p1-p3',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
maxfft{30}=plotspectrosel(68,'cl1-cl4',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'});
plotspectrosel(68,'cl1-cl4',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'},'nomaxfft');
plotspectrosel(94,'cl3-cl6',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'},'nomaxfft');
plotspectrosel(179,'pl2-p1',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'},'nomaxfft');
plotspectrosel(179,'pl1-pl2',[],'targ',plotparam,'win',[-3 4],'freq',[10 50],'markers',{'fix','fixeye','targeye','outcome'},'xinfos',xinfos,'binfos',binfos,'types',{'big'},'nomaxfft');

%pupil diameter vs DA/beta, using same time windwos as lfp
sessnums2=[67 68 69 75 79 80 83 84 91 92 94 95 96 100 102 109 113 114 127 179 181];
xts{1}=[0.2 0.8];
xts{2}=[0.2 1];
xts{3}=[0.2 2];
xts{4}=[1 2];
xts{5}=[2 3.8];
xts{6}=[0.2 3.8];
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','eye','ewin',[0.2 0.8],'sitelists',sitelists);              %beta vs eye
save('da_eye_cordata_targ_ewin_200-800ms','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','eye','ewin',[0.2 2],'sitelists',sitelists);              %beta vs eye
save('da_eye_cordata_targ_ewin_200-2000ms','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'lfponly','beh','eye','ewin',xts,'sitelists',sitelists);            %beta vs eye
save('beta_eye_cordata_targ_multiwin','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','eye','ewin',[0.2 0.8],'sitelists',sitelists,'condtypes',{'small','small'});               %beta vs eye
save('da_eye_cordata_targ_ewin_200-800ms_small','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'lfponly','beh','eye','ewin',xts,'sitelists',sitelists,'condtypes',{'small','small'});           %beta vs eye
save('beta_eye_cordata_targ_multiwin_small','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','eye','ewin',[0.2 0.8],'sitelists',sitelists,'condtypes',{'big','big'});                %beta vs eye
save('da_eye_cordata_targ_ewin_200-800ms_big','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'lfponly','beh','eye','ewin',xts,'sitelists',sitelists,'condtypes',{'big','big'});              %beta vs eye
save('beta_eye_cordata_targ_multiwin_big','tabdata','sitesigs')

%multiple regression
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',{[0.2 0.8]},'sitelists',sitelists,'multiregress');              

%eye velocity 
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','eyev','ewin',[0 1],'sitelists',sitelists);              %beta vs eye
save('da_eyev_cordata_targ_ewin_0-1s','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','eyev','ewin',[0 1],'sitelists',sitelists,'condition','left');              %beta vs eye
save('da_eyev_cordata_targ_ewin_0-1s_left','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','eyev','ewin',[0 1],'sitelists',sitelists,'condition','right');              %beta vs eye
save('da_eyev_cordata_targ_ewin_0-1s_right','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','eyev','ewin',[0 1],'sitelists',sitelists,'condtypes',{'big','big'});              %beta vs eye
save('da_eyev_cordata_targ_ewin_0-1s_big','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','eyev','ewin',[0 1],'sitelists',sitelists,'condtypes',{'small','small'});                %beta vs eye
save('da_eyev_cordata_targ_ewin_0-1s_small','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'lfponly','beh','eyev','ewin',[0 1],'sitelists',sitelists);            %beta vs eye
save('beta_eyev_cordata_targ_multiwin','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'lfponly','beh','eyev','ewin',[0 1],'sitelists',sitelists,'condition','left');            %beta vs eye
save('beta_eyev_cordata_targ_left_multiwin','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'lfponly','beh','eyev','ewin',[0 1],'sitelists',sitelists,'condition','right');            %beta vs eye
save('beta_eyev_cordata_targ_right_multiwin','tabdata','sitesigs')

%eye distance
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','eyedist','ewin',[0 1],'sitelists',sitelists);              %beta vs eye
save('da_eyedist_cordata_targ_ewin_0-1s','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','eyedist','ewin',[0 1],'sitelists',sitelists,'condition','left');              %beta vs eye
save('da_eyedist_cordata_targ_ewin_0-1s_left','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','eyedist','ewin',[0 1],'sitelists',sitelists,'condition','right');              %beta vs eye
save('da_eyedist_cordata_targ_ewin_0-1s_right','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','eyedist','ewin',[0 1],'sitelists',sitelists,'condtypes',{'big','big'});              %beta vs eye
save('da_eyedist_cordata_targ_ewin_0-1s_big','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','eyedist','ewin',[0 1],'sitelists',sitelists,'condtypes',{'small','small'});                %beta vs eye
save('da_eyedist_cordata_targ_ewin_0-1s_small','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'lfponly','beh','eyedist','ewin',[0 1],'sitelists',sitelists);            %beta vs eye
save('beta_eyedist_cordata_targ_multiwin','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'lfponly','beh','eyedist','ewin',[0 1],'sitelists',sitelists,'condition','left');            %beta vs eye
save('beta_eyedist_cordata_targ_left_multiwin','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'lfponly','beh','eyedist','ewin',[0 1],'sitelists',sitelists,'condition','right');            %beta vs eye
save('beta_eyedist_cordata_targ_right_multiwin','tabdata','sitesigs')

%licking
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','lick','lwin',[0 2],'sitelists',sitelists);              %beta vs eye
save('da_lick_cordata_targ_lwin_0-2s','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','lick','lwin',[0 2],'sitelists',sitelists,'condtypes',{'big','big'});              %beta vs eye
save('da_lick_cordata_targ_lwin_0-2s_big','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','lick','lwin',[0 2],'sitelists',sitelists,'condtypes',{'small','small'});                %beta vs eye
save('da_lick_cordata_targ_lwin_0-2s_small','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','lick','lwin',[0 2],'sitelists',sitelists,'condition','left');              %beta vs eye
save('da_lick_cordata_targ_lwin_0-2s_left','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','lick','lwin',[0 2],'sitelists',sitelists,'condition','right');                %beta vs eye
save('da_lick_cordata_targ_lwin_0-2s_right','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'lfponly','beh','lick','lwin',[0 2],'sitelists',sitelists);            %beta vs eye
save('beta_lick_cordata_targ_multiwin','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'lfponly','beh','lick','lwin',[0 2],'sitelists',sitelists,'condtypes',{'big','big'});               %beta vs eye
save('beta_lick_cordata_targ_multiwin_big','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'lfponly','beh','lick','lwin',[0 2],'sitelists',sitelists,'condtypes',{'small','small'});             %beta vs eye
save('beta_lick_cordata_targ_multiwin_small','tabdata','sitesigs')

[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'lfponly','beh','lick','lwin',[0 2],'sitelists',sitelists,'condition','left');             %beta vs eye
save('beta_lick_cordata_targ_multiwin_left','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'lfponly','beh','lick','lwin',[0 2],'sitelists',sitelists,'condition','right');             %beta vs eye
save('beta_lick_cordata_targ_multiwin_right','tabdata','sitesigs')
%HRV, rrstd
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','rrstd','hrwin',[1 4],'sitelists',sitelists);              %beta vs eye
save('da_rrstd_cordata_targ_hrwin_1-4s','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','rrstd','hrwin',[1 4],'sitelists',sitelists,'condtypes',{'big','big'});              %beta vs eye
save('da_rrstd_cordata_targ_hrwin_1-4s_big','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','rrstd','hrwin',[1 4],'sitelists',sitelists,'condtypes',{'small','small'});                %beta vs eye
save('da_rrstd_cordata_targ_multiwin_small_hrwin_1-4s','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'lfponly','beh','rrstd','hrwin',[1 4],'sitelists',sitelists);            %beta vs eye
save('beta_rrstd_cordata_targ_multiwin_hrwin_1-4s','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'lfponly','beh','rrstd','hrwin',[1 4],'sitelists',sitelists,'condtypes',{'big','big'});              %beta vs eye
save('beta_rrstd_cordata_targ_multiwin_big_hrwin_1-4s','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'lfponly','beh','rrstd','hrwin',[1 4],'sitelists',sitelists,'condtypes',{'small','small'});              %beta vs eye
save('beta_rrstd_cordata_targ_multiwin_small_hrwin_1-4s','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','rrstd','hrwin',[0 4],'sitelists',sitelists);              %beta vs eye
save('da_rrstd_cordata_targ_hrwin_0-4s','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','rrstd','hrwin',[0 4],'sitelists',sitelists,'condtypes',{'big','big'});              %beta vs eye
save('da_rrstd_cordata_targ_hrwin_0-4s_big','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','daonly','beh','rrstd','hrwin',[0 4],'sitelists',sitelists,'condtypes',{'small','small'});                %beta vs eye
save('da_rrstd_cordata_targ_hrwin_0-4s_small','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'lfponly','beh','rrstd','hrwin',[0 4],'sitelists',sitelists);            %beta vs eye
save('beta_rrstd_cordata_targ_multiwin','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'lfponly','beh','rrstd','hrwin',[0 4],'sitelists',sitelists,'condtypes',{'big','big'});              %beta vs eye
save('beta_rrstd_cordata_targ_multiwin_big','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'lfponly','beh','rrstd','hrwin',[0 4],'sitelists',sitelists,'condtypes',{'small','small'});              %beta vs eye
save('beta_rrstd_cordata_targ_multiwin_small','tabdata','sitesigs')

sessnums2=[67 68 69 75 79 80 83 84 91 92 94 95 96 100 102 109 113 114 127 179 181];
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','metriclfp','lfpmints','xvarflag','sitelists',sitelists);              %ERD latencies vs DA peaks
save('dabetacordata_lfpmints','tabdata','sitesigs')

%how do time windows used for beta averaging affect da vs beta correlation
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'sitelists',sitelists);  
save('dabetacordata_targ_multiwin','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'condition','left','sitelists',sitelists);  
save('dabetacordata_targ_multiwin_leftonly','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'condition','right','sitelists',sitelists);  
save('dabetacordata_targ_multiwin_rightonly','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'condtypes',{'big','big','big'},'sitelists',sitelists);  
save('dabetacordata_targ_multiwin_bigonly','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'condtypes',{'small','small','small'},'sitelists',sitelists);  
save('dabetacordata_targ_multiwin_smallonly','tabdata','sitesigs')


%DG 2d PLOT 04/23/20
starts = 0:0.2:3.6;
ends=0.2:.2:3.8;
%starts2 = repmat(starts', length(ends));
count=1;
xmats={};
for it=1:length(starts)
    for iit=1:length(ends)
        if ends(iit)-starts(it)>0
            xmats{count}=[starts(it) ends(iit)];
            count=count+1;
        end
    end
end
xmatrix1=plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'small','left'},{'small','right'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xmats,'patra','gain','simple','xmatrix');
xmatrix1_con=plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'small','left'},{'small','right'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xmats,'patra','constrainpairs',sitelists,'gain','simple','xmatrix');
xmatrix2_con=plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'small','left'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xmats,'patra','constrainpairs',sitelists,'gain','simple','xmatrix');
xmatrix2=plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'small','left'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xmats,'patra','gain','simple','xmatrix');
xmatrix3_con=plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xmats,'patra','constrainpairs',sitelists,'gain','simple','xmatrix');
xmatrix3=plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xmats,'patra','gain','simple','xmatrix');
xmatrix4_con=plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','right'},{'small','right'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xmats,'patra','constrainpairs',sitelists,'gain','simple','xmatrix');
xmatrix4=plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','right'},{'small','right'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xmats,'patra','gain','simple','xmatrix');

%how do time windows used for beta averaging affect task influences

plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'small','left'},{'small','right'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'small','left'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','right'},{'small','right'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');

%how do time windows used for beta averaging affect da vs beta correlations

%CORS ALL SESSIONS & TABULATE
%with .2 offset
xts={};
xts{1}=[0.2 .8];
xts{2}=[0.2 3.8];
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'small','left'},{'small','right'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'small','left'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','right'},{'small','right'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');

plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'small','left'},{'small','right'}},'lfp','event','targeye','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'small','left'}},'lfp','event','targeye','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'lfp','event','targeye','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','right'},{'small','right'}},'lfp','event','targeye','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');

[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'constrainpairs');
save('dabetacordata_targ_earlywinvs_wholewinpadded_paired','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts);
save('dabetacordata_targ_earlywinvs_wholewinpadded','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targeye','win',xts,'constrainpairs');
save('dabetacordata_targeye_earlywinvs_wholewinpadded_paired','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targeye','win',xts);
save('dabetacordata_targeye_earlywinvs_wholewinpadded','tabdata','sitesigs')

ii=0:2:2;
wint=2;
xts={};
for it=1:length(ii)
xts{it}=[ii(it) ii(it)+wint];
end
xts{4}=[0 4];
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'small','left'},{'small','right'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'small','left'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','right'},{'small','right'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');

plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'small','left'},{'small','right'}},'lfp','event','targeye','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'small','left'}},'lfp','event','targeye','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'lfp','event','targeye','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','right'},{'small','right'}},'lfp','event','targeye','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');


[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'constrainpairs');
save('dabetacordata_targ_2swin_paired','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts);
save('dabetacordata_targ_2swin','tabdata','sitesigs')

[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targeye','win',xts,'constrainpairs');
save('dabetacordata_targeye_2swin_paired','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targeye','win',xts);
save('dabetacordata_targeye_2swin','tabdata','sitesigs')

ii=0:1:3;
wint=1;
xts={};
for it=1:length(ii)
xts{it}=[ii(it) ii(it)+wint];
end
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'small','left'},{'small','right'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'small','left'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','right'},{'small','right'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');

plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'small','left'},{'small','right'}},'lfp','event','targeye','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'small','left'}},'lfp','event','targeye','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'lfp','event','targeye','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','right'},{'small','right'}},'lfp','event','targeye','metric','xx','xinfo',xinfos,'xwin',xts,'patra','constrainpairs',sitelists,'gain','simple');

plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'small','left'},{'small','right'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xts,'patra','gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'small','left'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xts,'patra','gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xts,'patra','gain','simple');
plotdasummaryxwins(betatm,plotparam,'norm','mean','ttypes',{{'big','right'},{'small','right'}},'lfp','event','targ','metric','xx','xinfo',xinfos,'xwin',xts,'patra','gain','simple');



[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts);
save('dabetacordata_targ_1swin','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targeye','win',xts);
save('dabetacordata_targeye_1swin','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',xts,'constrainpairs');
save('dabetacordata_targ_1swin_paired','tabdata','sitesigs')
[tabdata,sitesigs]=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targeye','win',xts,'constrainpairs');
save('dabetacordata_targeye_1swin_paired','tabdata','sitesigs')




tabdata=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targeye','win',[0 1]);
tabdata_offset=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targeye','win',[0.2 1]);
tabdata_2swin=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targeye','win',[0 2]);
tabdata_targ_1swin_offsetdual=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',[0.2 .8]);
tabdata_targ_lfptm_targimwin=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','uselfptm');
tabdata_targ_1swin_offset=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',[0.2 1]);
tabdata_targ_1swin=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',[0 1]);
tabdata_targ_2swin=corsmultiple(sessnums2,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,'event','targ','win',[0 2]);






plotmultxtraces(83,xinfos,plotparam,'win',[0 4],'event','targ','ttypes',{'big','small','all'},'datm',datm,'sort','targpeak','lfps',{'pl1-p5'},'smoothlfp',[5 8]);
plotmultxtraces(83,xinfos,plotparam,'win',[0 4],'event','targ','ttypes',{'big','small','all'},'datm',datm,'sort','targpeak','lfps',{'pl1-p5'},'plotimmean','pad',.2,'ext',.8);


targetlfp={92,'cl3-cl6',{'big','small'},'targ','percentile',{'lfp',30}};
trialnumslfp=getmulttrials(plotparam,xinfos,binfos,targetlfp);

targetda={92,'cl5',{'big','small'},'targ','percentile',{'da',30}};
trialnumsda=getmulttrials(plotparam,xinfos,binfos,targetda);

trialnums=getintersecttrials(trialnumsboth,{{'top25da','bot25lfp'},{'bot25da','top25lfp'}});
trialnums=getintersecttrials(trialnumsboth,{{'top25da','bot25lfp'},{'bot25da','top25lfp'},{'top25da','top25lfp'},{'bot25da','bot25lfp'}});

targses=92;
targsite='pl2';
targlfp='pl1-p5';
pertile=25;

targses=92; %SIGNIFICANT DA/BETA ANTI- CORRELATION
targsite='cl5';
targlfp='eye';
pertile=25;
targetbeh={targses,targlfp,{'big'},'targeye','percentile',{'beh',pertile},'win',[0 4],'fixedcond',{'right'}};
trialnumsbeh=getmulttrials(plotparam,xbinfos,binfos,targetbeh);
plottracesel(targses,xbinfos,'eye',plotparam,trialnumsbeh,'win',[-4 7],'se','beh','condition','right','event','targeye');
plottracesel(92,xinfos,'cl5',plotparam,trialnumsbeh,'win',[-3 4],'se','condition','right');
plottracesel(92,xinfos,'cl3-cl6',plotparam,trialnumsbeh,'lfp','win',[-3 4],'se','condition','right');

targetbeh={targses,targlfp,{'big'},'targeye','percentile',{'beh',pertile},'win',[0 4],'fixedcond',{'left'}};
trialnumsbeh=getmulttrials(plotparam,xbinfos,binfos,targetbeh);
plottracesel(targses,xbinfos,'eye',plotparam,trialnumsbeh,'win',[-4 7],'se','beh','condition','left','event','targeye');
plottracesel(targses,xinfos,'cl5',plotparam,trialnumsbeh,'win',[-3 4],'se','condition','left');
plottracesel(targses,xinfos,'cl3-cl6',plotparam,trialnumsbeh,'lfp','win',[-3 4],'se','condition','left');

targses=80;
targsite='cl6';
[~,targbeta]=getsitepairs({targsite}); targbeta=targbeta{1};
targlfp='eye';
pertile=25;
targetbeh={targses,targlfp,{'big'},'targeye','percentile',{'beh',pertile},'win',[0 4],'fixedcond',{'right'}};
trialnumsbeh=getmulttrials(plotparam,xbinfos,binfos,targetbeh);
plottracesel(targses,xbinfos,'eye',plotparam,trialnumsbeh,'win',[-4 7],'se','beh','condition','right','event','targeye');
plottracesel(targses,xinfos,'cl5',plotparam,trialnumsbeh,'win',[-3 4],'se','condition','right');
plottracesel(targses,xinfos,'cl3-cl6',plotparam,trialnumsbeh,'lfp','win',[-3 4],'se','condition','right');

targetbeh={targses,targlfp,{'big'},'targeye','percentile',{'beh',pertile},'win',[0 4],'fixedcond',{'left'}};
trialnumsbeh=getmulttrials(plotparam,xbinfos,binfos,targetbeh);
plottracesel(targses,xbinfos,'eye',plotparam,trialnumsbeh,'win',[-4 7],'se','beh','condition','left','event','targeye');
plottracesel(targses,xinfos,targsite,plotparam,trialnumsbeh,'win',[-3 4],'se','condition','left','event','targeye');
plottracesel(targses,xinfos,targbeta,plotparam,trialnumsbeh,'lfp','win',[-3 4],'se','condition','left','event','targeye');

targses=67;     %SIGNIFICANT DA CORRELATION
targsite='cl5';
[~,targbeta]=getsitepairs({targsite}); targbeta=targbeta{1};
targlfp='lick';
pertile=15;
targetbeh={targses,targlfp,{'big'},'targ','percentile',{'beh',pertile},'win',[0 1]};
trialnumsbeh=getmulttrials(plotparam,xbinfos,binfos,targetbeh);
plottracesel(targses,xbinfos,'lick',plotparam,trialnumsbeh,'win',[-4 7],'se','beh','event','targ');
plottracesel(targses,xinfos,targsite,plotparam,trialnumsbeh,'win',[-3 4],'se','event','targ');
plottracesel(targses,xinfos,targbeta,plotparam,trialnumsbeh,'lfp','win',[-3 4],'se','event','targ');

targses=83; %SIGNIFICANT DA/BETA ANTI- CORRELATION, also discriminates lick, rt not always discriminating lick
targsite='cl5';
[~,targbeta]=getsitepairs({targsite}); targbeta=targbeta{1};
targlfp='trt';
pertile=30;
targetbeh={targses,targlfp,{'big'},'targ','percentile',{'beh',pertile},'fixedcond',{'left'}};
trialnumsbeh=getmulttrials(plotparam,xbinfos,binfos,targetbeh);
plottracesel(targses,xbinfos,'lick',plotparam,trialnumsbeh,'win',[-4 7],'se','beh','event','targ');
plottracesel(targses,xinfos,targsite,plotparam,trialnumsbeh,'win',[-3 4],'se','event','targ');
plottracesel(targses,xinfos,targbeta,plotparam,trialnumsbeh,'lfp','win',[-3 4],'se','event','targeye');

%combined behavior groups
targses=67; 
targsite='cl5';
[~,targbeta]=getsitepairs({targsite}); targbeta=targbeta{1};
targlfp='trt';
pertile1=40;
pertile2=40;

targetbeh={targses,targlfp,{'big'},'targ','percentile',{'beh',pertile1},'fixedcond',{'left'}};
trialnumsrt=getmulttrials(plotparam,xbinfos,binfos,targetbeh);
targetlick={targses,'lick',{'big'},'targ','percentile',{'beh',pertile2}};
trialnumslick=getmulttrials(plotparam,xbinfos,binfos,targetlick);
trialnumsboth=[trialnumsrt trialnumslick];
trialnums=getintersecttrials(trialnumsboth,{{['top' num2str(pertile1) 'trt'],['top' num2str(pertile2) 'lick']},...
    {['bot' num2str(pertile1) 'trt'],['bot' num2str(pertile2) 'lick']}});
plottracesel(targses,xbinfos,'lick',plotparam,trialnums,'win',[-4 7],'se','beh','event','targ');
plottracesel(targses,xinfos,targsite,plotparam,trialnums,'win',[-3 4],'se','event','targ');
plottracesel(targses,xinfos,targbeta,plotparam,trialnums,'lfp','win',[-3 4],'se','event','targeye');



% DA / BETA SIGNAL GROUPS
targetlfp={targses,targlfp,{'big','small'},'targ','percentile',{'lfp',pertile}};
targetlfp={targses,targlfp,{'big'},'targ','percentile',{'lfp',pertile},'targimwin'};
%targetlfp={targses,targlfp,{'big'},'targ','percentile',{'lfp',pertile}};
targetda={targses,targsite,{'big','small'},'targ','percentile',{'da',pertile}};
targetda={targses,targsite,{'big'},'targ','percentile',{'da',pertile}};

trialnumslfp=getmulttrials(plotparam,xinfos,binfos,targetlfp);
trialnumsda=getmulttrials(plotparam,xinfos,binfos,targetda);
trialnumsboth=[trialnumsda trialnumslfp];
trialnums=getintersecttrials(trialnumsboth,{{['top' num2str(pertile) 'da'],['bot' num2str(pertile) 'lfp']},...
    {['bot' num2str(pertile) 'da'],['top' num2str(pertile) 'lfp']},...
    {['top' num2str(pertile) 'da'],['top' num2str(pertile) 'lfp']},...
    {['bot' num2str(pertile) 'da'],['bot' num2str(pertile) 'lfp']}});
plottracesel(targses,xbinfos,'eye',plotparam,trialnums,'win',[-4 7],'se','beh','condition','right','event','targeye');
plottracesel(targses,xbinfos,'eye',plotparam,trialnums,'win',[-4 7],'se','beh','condition','left','event','targeye');



[trlist,sitemap]=gettrialseries(92,trlists,xinfos,xbinfos,binfos,datm,betatm,plotparam);

%trialnums=getintersecttrials(trialnumsboth,{{'top30da','bot30lfp'},{'bot30da','top30lfp'},{'top30da','top30lfp'},{'bot30da','bot30lfp'}});
plottracesel(targses,xbinfos,'pulse',plotparam,trialnums,'win',[0 12],'se','beh');
plottracesel(targses,xbinfos,'pulse',plotparam,trialnums,'rrstd','hist','win',[0 4],'se','beh');
plottracesel(targses,xbinfos,'pulse',plotparam,trialnums,'rmssd','hist','win',[0 4],'se','beh');

plottracesel(targses,binfos,'trt',plotparam,trialnums,'hist','condition','right');
plottracesel(targses,binfos,'trt',plotparam,trialnums,'hist','condition','left');
plottracesel(targses,xbinfos,'lick',plotparam,trialnums,'win',[0 1],'se','beh','hist','event','targ');

plottracesel(targses,xbinfos,'lick',plotparam,trialnums,'win',[-4 12],'se','beh');
plottracesel(targses,xbinfos,'lick',plotparam,trialnums,'win',[-2 12],'se','beh','event','fix');


plotspectrosel(92,'cl3-cl6',trialnums(1),'targ','win',[-3 7],'markers',{'fix','fixeye','targeye','outcome'});
plotspectrosel(targses,targlfp,trialnumsda(1),'targ',plotparam,'win',[-3 7],'markers',{'fix','fixeye','targeye','outcome'});
plotspectrosel(targses,targlfp,trialnumsda(2),'targ',plotparam,'win',[-3 7],'markers',{'fix','fixeye','targeye','outcome'});
plotspectrosel(targses,targlfp,trialnumsda(1),'targ',plotparam,'win',[-3 7],'markers',{'fix','fixeye','targeye','outcome'},'condition','left');
plotspectrosel(targses,targlfp,trialnumsda(2),'targ',plotparam,'win',[-3 7],'markers',{'fix','fixeye','targeye','outcome'},'condition','left');

plotspectrosel(targses,targlfp,trialnumsda(1),'targ',plotparam,'win',[-3 7],'markers',{'fix','fixeye','targeye','outcome'},'condition','right');
plotspectrosel(targses,targlfp,trialnumsda(2),'targ',plotparam,'win',[-3 7],'markers',{'fix','fixeye','targeye','outcome'},'condition','right');

plottracesel(92,xinfos,'cl5',plotparam,trialnums,'win',[0 5],'se','condition','right');
plottracesel(92,xinfos,'cl3-cl6',plotparam,trialnums,'lfp','win',[0 5],'label',target,'se','condition','right');
plottracesel(92,xinfos,'cl3-cl6',plotparam,trialnums,'lfp','win',[0 5],'label',target,'se','condition','left');

plottracesel(92,xinfos,'cl5',plotparam,trialnums,'win',[0 5],'se');
plottracesel(92,xinfos,'cl3-cl6',plotparam,trialnums,'lfp','win',[0 5],'se');

plottracesel(targses,xbinfos,'pulse',plotparam,trialnums,'win',[0 12],'se','beh');
plottracesel(92,xbinfos,'pulse',plotparam,trialnumsda,'win',[0 5],'se','beh');
plottracesel(92,xbinfos,'pulse',plotparam,trialnumslfp,'win',[0 5],'se','beh');


plottracesel(92,xbinfos,'eye',plotparam,trialnums,'win',[-4 7],'se','beh','condition','left','event','targeye');
plottracesel(92,xbinfos,'eyev',plotparam,trialnums,'win',[0 6],'se','beh');
plottracesel(92,xbinfos,'eyev',plotparam,trialnums,'win',[0 6],'se','beh','condition','left');

plottrialseries(92,trlists,xinfos,xbinfos,binfos,plotparam,'rt','datm','cl5',datm,'lfptm','cl3-cl6',betatm,'hr')

plotcortrials(68,trlists,datm,betatm,binfos,plotparam,'datm','cl6','lfptm','cl1-cl5','rt','hf',[0 0 -1],'hfsametype','types',{'all','all','small'},'metric','targpeak');
plotcortrials(92,trlists,datm,betatm,binfos,plotparam,'datm','cl5','lfptm','cl3-cl6','rt','hf',[0 0 0])
plotcortrials(100,trlists,datm,betatm,binfos,plotparam,'datm','cl6','lfptm','cl1-cl5','rt','types',{'all','all','small'},'metric','targpeak','metriclfp','rewshortwin');
plotcortrials(94,trlists,datm,betatm,binfos,plotparam,'datm','cl5','lfptm','cl3-cl6','eye','condition','right','types',{'all','all','all'},'metric','targpeak','metriclfp','targimwin');
plotcortrials(92,trlists,datm,betatm,binfos,plotparam,'datm','cl5','lfptm','cl3-cl6','hr','types',{'all','all','big'},'metric','targpeak','metriclfp','targimwin');

plotmultxtraces(67,xinfos,plotparam,'win',[0 4],'event','targeye','ttypes',{'big','all'},'sort','damax','bplot',xbinfos,'pulse','plotimmean','ext',4,'pad',1);
plotmultxtraces(75,xinfos,plotparam,'win',[0 4],'event','targ','ttypes',{'big','all'},'sort','damax','bplot',xbinfos,'lick','plotmax','ext',4);
plotmultxtraces(69,xinfos,plotparam,'win',[0 4],'event','targ','ttypes',{'small','all'},'sort','damax');
plotmultxtraces(92,xbinfos,plotparam,'win',[0 4],'event','targ','beh','lick','ttypes',{'big','small','all'},'sort','damax','plotimmean','ext',3);
plotmultxtraces(92,xinfos,plotparam,'win',[0 4],'event','targ','ttypes',{'big','small','all'},'sort','damax','plotimmean','ext',1);

plotmultxtraces(92,xinfos,plotparam,'win',[0 4],'event','targ','ttypes',{'big','small','all'},'datm',datm,'sort','targwin');
plotmultxtraces(92,xinfos,plotparam,'win',[0 4],'event','targ','ttypes',{'big','small','all'},'datm',datm,'sort','targpeak','lfps',{'cl4-cl6'});
plotmultxtraces(92,xinfos,plotparam,'win',[0 4],'event','targ','ttypes',{'big','small','all'},'datm',datm,'sort','targpeak','lfps',{'cl4-cl6'},'plotmax');
plotmultxtraces(83,xinfos,plotparam,'win',[0 4],'event','targ','ttypes',{'big','small','all'},'datm',datm,'sort','targpeak','smoothlfp',5);

%da traces sorted by datm peak with licking side by side
plotmultxtraces(83,xinfos,plotparam,'win',[0 4],'event','targ','ttypes',{'big','all'},'datm',datm,'sort','targpeak','bplot',xbinfos,'lick','plotimmean','ext',4);
plotmultxtraces(83,xinfos,plotparam,'win',[0 4],'event','targ','ttypes',{'big','all'},'datm',datm,'sort','targpeak','bplot',xbinfos,'lick','plotimmean','ext',3);

%da traces sorted by datm peak with hrv side by side
plotmultxtraces(75,xinfos,plotparam,'win',[0 4],'event','targ','ttypes',{'big','all'},'datm',datm,'sort','targpeak','bplot',xbinfos,'pulse','rrstd','winb',[1 3.8]);
plotcortrials(83,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','pl2','rrstd','lfptm','pl1-p5','types',{'small','small','small'},'event','targ','metric','targpeak','metriclfp','targwin','hrwin',[1 4]);
plotmultxtraces(83,xinfos,plotparam,'win',[0 4],'event','targ','ttypes',{'small','all'},'datm',datm,'sort','targpeak','bplot',xbinfos,'pulse','rrstd','winb',[1 4]);

%correlation to hr
plotcortrials(114,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','cl3','lfptm','cl1-cl5','hr','types',{'all','all','big'},'event','targ','metric','targpeak','metriclfp','targimwin','hrwin',[.1 4]);
plotcortrials(114,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','cl3','lfptm','cl1-cl5','hr','types',{'all','all','big'},'event','targ','metric','targpeak','metriclfp','targimwin','hrwin',[.1 4]);
plotcortrials(91,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','p5','rrstd','lfptm','cl1-cl5','types',{'big','big','big'},'event','targ','metric','targpeak','metriclfp','targimwin','hrwin',[1 4]);
plotcortrials(92,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','pl2','rrstd','lfptm','pl1-p5','types',{'small','small','small'},'event','targ','metric','targpeak','metriclfp','targwin','hrwin',[0 4]);

%good cors to hr
plotcortrials(75,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','pl2','rrstd','lfptm','pl1-p5','types',{'small','small','small'},'event','targ','metric','targpeak','metriclfp','targimwin','hrwin',[1 4]);
plotcortrials(69,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','cl6','rrstd','lfptm','cl1-cl5','types',{'small','small','small'},'event','targ','metric','targpeak','metriclfp','targimwin','hrwin',[0 4]);
plotcortrials(69,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','pl1','rrstd','lfptm','p1-pl3','types',{'small','small','small'},'event','targ','metric','targpeak','metriclfp','targimwin','hrwin',[0 4]);
plotcortrials(95,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','cl5','rrstd','lfptm','cl3-cl6','types',{'small','small','small'},'event','targ','metric','targpeak','metriclfp','targimwin','hrwin',[0 6]);
plotcortrials(95,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','pl2','rrstd','lfptm','pl1-p5','types',{'small','small','small'},'event','targ','metric','targpeak','metriclfp','targimwin','hrwin',[0 4]);
plotcortrials(114,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','pl1','rrstd','lfptm','p1-pl3','types',{'small','small','small'},'event','targ','metric','targpeak','metriclfp','targimwin','hrwin',[1 4]);

plotcortrials(67,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','p5','rrstd','lfptm','pl1-p5','types',{'small','small','small'},'event','targ','metric','targpeak','metriclfp','targimwin','hrwin',[0 4]);



%good example with very strong correlation pl2 and pl1-p5
plotcortrials(92,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','pl2','rrstd','lfptm','pl1-p5','types',{'small','small','small'},'event','targ','metric','targpeak','metriclfp','targimwin','hrwin',[1 4]);
%both da and beta correlated to hr
plotcortrials(83,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','pl2','rrstd','lfptm','pl1-p5','types',{'small','small','small'},'event','targ','metric','targpeak','metriclfp','targimwin','hrwin',[1 4]);
%good cors to hr adjusted win

plotcortrials(114,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','p5','rrstd','lfptm','p1-pl3','types',{'small','small','small'},'event','targ','metric','targpeak','metriclfp','targimwin','hrwin',[3 6]);
plotcortrials(67,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','pl2','rrstd','lfptm','pl1-p5','types',{'small','small','small'},'event','targ','metric','targpeak','metriclfp','targimwin','hrwin',[2 4]);

%correlation with licking
plotcortrials(114,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','cl3','lfptm','cl1-cl5','lick','types',{'all','all','big'},'event','targ','metric','targpeak','metriclfp','targimwin','lwin',[.1 3]);
plotcortrials(83,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','cl5','lick','trialseq','types',{'big','big','big'},'event','targ','metric','targpeak','metriclfp','targimwin','lwin',[0 3]);
plotcortrials(83,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'lfptm','cl3-cl6','lick','datm','cl5','types',{'big','big','big'},'event','targ','metric','targpeak','metriclfp','targimwin','lwin',[0 3]);
plotmultxtraces(83,xinfos,plotparam,'win',[0 4],'event','targ','ttypes',{'big','all'},'datm',datm,'sort','targpeak','bplot',xbinfos,'lick','plotimmean','ext',3,'zbeh');
plotmultxtraces(83,xinfos,plotparam,'win',[0 4],'event','targ','ttypes',{'big','all'},'datm',datm,'sort','targpeak','bplot',xbinfos,'lick','plotimmean','ext',4,'zbeh');
plotcortrials(83,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'lfptm','cl3-cl6','lick','datm','cl5','types',{'big','big','big'},'event','targ','metric','targpeak','metriclfp','targimwin','lwin',[0 4]);
plotcortrials(83,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','cl5','lick','trialseq','types',{'big','big','big'},'event','targ','metric','targpeak','metriclfp','targimwin','lwin',[0 4]);

%correlation pupil
plotcortrials2(94,trlists,datm,betatm,binfos,plotparam,'eye','datm','cl5','lfpx','cl3-cl6',xinfos,'win',[.2 .8],'offset',[0 0],'types',{'all','all','all'},'behx',xbinfos,'event','targ','metric','targpeak','ewin',[.2 .8]);

plotcortrials(114,trlists,datm,betatm,binfos,plotparam,'datm','cl3','lfptm','cl1-cl5','eye','types',{'all','all','big'},'event','targ','metric','targpeak','metriclfp','targimwin');
plotcortrials(114,trlists,datm,betatm,binfos,plotparam,'datm','cl3','eye','trialseq','types',{'big','big','big'},'event','targ','metric','targpeak','metriclfp','targimwin');
plotcortrials(114,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','cl3','lfptm','cl1-cl5','eye','types',{'all','all','big'},'event','targ','metric','targpeak','metriclfp','targimwin','ewin',[2 4]);
plotcortrials(114,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','cl3','lfptm','cl1-cl5','eye','types',{'all','all','all'},'event','targ','metric','targpeak','metriclfp','targimwin','ewin',[.1 .8]);
plotcortrials(114,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','cl3','eye','trialseq','types',{'all','all','all'},'event','targ','metric','targpeak','metriclfp','targimwin','ewin',[.1 .8]);

%correlation against trial progression
plotcortrials(80,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','cl6','rr','trialseq','types',{'big','big','big'},'event','targ','metric','targpeak','metriclfp','targimwin','hrwin',[1 4]);

%correlation all sites
plotcormultiple(plotparam.sessnums,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'types',{'small'},'event','targ','metric','targpeak','metriclfp','targimwin','hrwin',[1 4]);
plotcormultiple(plotparam.sessnums,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'types',{'big'},'event','targ','metric','targpeak','metriclfp','targimwin','hrwin',[1 4]);
plotcormultiple(plotparam.sessnums,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'types',{'big','small'},'event','targ','metric','targpeak','metriclfp','targimwin','xvars',{'da','lfp','trialnum'},'yvars',{'da','lfp','trialnum','hr','rrstd','lick','eye','trt'},'hrwin',[1 4]);
plotcormultiple(plotparam.sessnums,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'types',{'big','small'},'event','targ','metric','targpeak','metriclfp','targpeak','xvars',{'da','lfp','trialnum'},'yvars',{'da','lfp','trialnum','hr','rrstd','lick','eye','trt'},'hrwin',[1 4]);
plotcormultiple(plotparam.sessnums,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'types',{'big','small'},'event','targ','metric','targpeak','metriclfp','targpeak','xvars',{'da','lfp','trialnum'},'yvars',{'da','lfp','trialnum','hr','rrstd','lick','eye','trt'},'hrwin',[1 4]);

%cor da vs beta
plotcortrials(92,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','cl5','lfptm','cl4-cl6','types',{'all','all','all'},'trialseq','event','targeye','metric','targwin','metriclfp','targpeak');
plotcortrials(83,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','cl5','lfptm','cl3-cl6','types',{'all','all','all'},'trialseq','event','targ','metric','targpeak','metriclfp','targwin');
plotcortrials(67,trlists,datm,betatm,binfos,plotparam,'behx',xbinfos,'datm','cl5','lfptm','cl3-cl6','types',{'all','all','all'},'trialseq','event','targ','metric','targpeak','metriclfp','targwin');


%run loop for plotcortrials2 to get info for each session/site pairs
        %all valid da sites
        %{
        targevent='targeye';
    targdasites=plotparam.dasites;
    dasites=getsites(sessnums,targdasites,'patra');
    uniquedasites=unique({dasites(1:end).site});
    cndasites=uniquesites(contains(uniquedasites,'c'));
    pdasites=uniquesites(contains(uniquedasites,'p'));
    sites=getlfpsites(sessnums,plotparam.lfpchs,'patra');
    uniquesites=unique({sites(1:end).site});
    cnsites=uniquesites(contains(uniquesites,'c'));
    psites=uniquesites(contains(uniquesites,'p'));
    savepath=fullfile(plotparam.savepath,'tables_cor_');
    %constrain lfp pair for each da site recorded loop through each session
    %to get optimal lfp pair for each da site, already called assignpairs
    %in 'plotmultiple.m' and suppleid here
    pairedsites=unique({sitelists(1:end).site}); %constrainedpairs
    cnpairedsites=uniquesites(contains(pairedsites,'c'));
    ppairedsites=uniquesites(contains(pairedsites,'p'));
    tdata_targeye_nooffset={};
    tdata_targeye_offset={};
    tdata_targeye_2swin={};
    tdata_targ_offset_1swin={};
    tdata_targ_offsetdual_1swin={};
    tdata2={};
    count=1;
    for ise=1:length(sessnums)
        disp(['sess ' num2str(sessnums(ise))]);
        curses=sessnums(ise);
        cursessids=find([sites.sessnum]==sessnums(ise));
        curlfpchsites={sites(cursessids).probeid};
        curlfpsites={sites(cursessids).site};
        targxinfoids=find(strcmp({xinfos.sessionid},num2str(curses)));
        lfpsinxinfo=unique({xinfos(targxinfoids).sitelfp}); %lfp chs stored in xinfo
        curlfpchsites=intersect(curlfpchsites,lfpsinxinfo); %only lfp chs actually stored in xinfo
        lfpidsinxinfo=find(ismember(curlfpchsites,lfpsinxinfo));
        curlfpsites=curlfpsites(lfpidsinxinfo);
        curdasessids=find([dasites.sessnum]==sessnums(ise));
        curdachsites={dasites(curdasessids).probeid};
        curdasites={dasites(curdasessids).site};
        dasinxinfo=unique({xinfos(targxinfoids).siteda}); %da chs stored in xinfo
        daidsinxinfo=find(ismember(curdachsites,dasinxinfo));
        curdasites=curdasites(daidsinxinfo);
        curdachsites=intersect(curdachsites,dasinxinfo); %only da chs actually stored in xinfo
        for id=1:length(curdachsites)
            curda=curdachsites{id};
            disp(['da ' curda]);
            curlfps={};
            curlfps=curlfpchsites(contains(curlfpchsites,curda(1)));        %region specific lfps, ie if da ch in CN, lfp chs in CN           
            for il=1:length(curlfps) 
                disp(['lfp ' curlfps{il}]);
                tempdata=plotcortrials2(curses,trlists,datm,betatm,binfos,plotparam,'datm',curda,'lfpx',curlfps{il},xinfos,'win',[0 1],'behx',xbinfos,'types',{'all','all','all'},'trialseq','event',targevent,'metric','targpeak','noplot');
                tempdata=setfield(tempdata,{1},['siteda'],curdasites{id});
                tempdata=setfield(tempdata,{1},['sitelfp'],curlfpsites{il});
                tempdata2=plotcortrials2(curses,trlists,datm,betatm,binfos,plotparam,'datm',curda,'lfptm',curlfps{il},'xinfo',xinfos,'win',[0 1],'behx',xbinfos,'types',{'all','all','all'},'trialseq','event',targevent,'metric','targpeak','metriclfp','targimwin','noplot');
                tempdata2=setfield(tempdata2,{1},['siteda'],curdasites{id});
                tempdata2=setfield(tempdata2,{1},['sitelfp'],curlfpsites{il});
                if ~isempty(tdata)
                    tdata=[tdata tempdata];
                    tdata2=[tdata2 tempdata2];
                else
                    tdata=tempdata;
                    tdata2=tempdata2;
                end
                count=count+1;
            end
        end
        close all;
    end
        %}

%plot average time course trace based on trial history
targses=94;
targsite='cl5';
targsitelfp='cl3-cl6';
targsite='cl3';
targsitelfp='cl1-cl5';
%cleo
plottracesel(16,xinfos,'cl5b',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'all'});
plottracesel(16,xinfos,'cl6100',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'all'});
plottracesel(16,xinfos,'p6',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'all'});
plottracesel(16,xinfos,'p6',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'condtypes',{'big'},{'left','right'});
plottracesel(18,xinfos,'p6',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'all'});
plottracesel(19,xinfos,'p10',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'condtypes',{'big'},{'left','right'});
plottracesel(23,xinfos,'p3',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'condtypes',{'big'},{'left','right'});
plottracesel(23,xinfos,'p10',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'condtypes',{'big'},{'left','right'});
plottracesel(24,xinfos,'p3',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'condtypes',{'big'},{'left','right'});
plottracesel(24,xinfos,'p3',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'all'});
plottracesel(23,xinfos,'p10',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'all'});
plottracesel(19,xinfos,'p10',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'all'});
plottracesel(19,xinfos,'p3',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'all'});

plottracesel(24,xinfos,'p3',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'condtypes',{'big'},{'left','right'});
plottracesel(18,xinfos,'cl6',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'condtypes',{'big'},{'suc','fail'});
plottracesel(18,xinfos,'c5',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'condtypes',{'big'},{'suc','fail'});
plottracesel(24,xinfos,'p6',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'condtypes',{'big'},{'suc','fail'});
plottracesel(25,xinfos,'cl7-p6',plotparam,[],'win',[-3 4],'lfp','se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'all'});
plottracesel(16,xinfos,'p3-p4',plotparam,[],'win',[-3 4],'lfp','se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'all'});
plottracesel(16,xinfos,'cl3-cl1',plotparam,[],'win',[-3 4],'lfp','se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'all'});
plottracesel(14,xinfos,'p4-p2',plotparam,[],'win',[-3 4],'lfp','se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'all'});
plottracesel(14,xinfos,'p2-p0',plotparam,[],'win',[-3 4],'lfp','se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'all'});
plottracesel(14,xinfos,'cl4-cl3',plotparam,[],'win',[-3 4],'lfp','se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'all'});


plottracesel(109,xinfos,'pl2',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'condtypes',{'big'},{'left','right'});
plottracesel(109,xinfos,'pl3',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'condtypes',{'big'},{'left','right'});
plottracesel(92,xinfos,'pl3',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'condtypes',{'big'},{'left','right'});
plottracesel(83,xinfos,'pl2',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'condtypes',{'big'},{'left','right'});
plottracesel(113,xinfos,'pl1',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'condtypes',{'big'},{'left','right'});

plottracesel(83,xinfos,'pl1-p5',plotparam,[],'win',[-3 4],'ibi','se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'all'});
plottracesel(92,xinfos,'cl3-cl6',plotparam,[],'win',[-3 4],'ibi','se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'all'});
plottracesel(92,xinfos,'cl3-cl6',plotparam,[],'win',[-3 4],'lfp','se','event','targ','binfo',binfos,'condtypes',{'big'},{'phase1','phase4'});
plottracesel(92,xinfos,'cl5',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'condtypes',{'big'},{'phase1','phase4'});

plottracesel(92,xinfos,'pl1-p5',plotparam,[],'win',[-3 4],'ibi','se','event','targ','binfo',binfos,'condtypes',{'big'},{'left','right'});
plottracesel(83,xinfos,'pl1-p5',plotparam,[],'win',[-3 4],'lfp','se','event','targ','binfo',binfos,'condtypes',{'big'},{'left','right'});

plottracesel(94,xinfos,'cl5',plotparam,[],'win',[-3 4],'se','binfo',binfos,'trtable',trtable,{'big'},{{'-2fix','-1fix','-1small','-1big'}});

plottracesel(113,xinfos,'cl3',plotparam,[],'win',[-3 4],'se','binfo',binfos,'trtable',trtable,{'small'},{{'-3fix','-2fix','-1fix','-1big','-1small'}});
plottracesel(113,xinfos,'cl1-cl4',plotparam,[],'win',[-3 4],'se','lfp','binfo',binfos,'trtable',trtable,{'small'},{{'-3fix','-2fix','-1fix','-1big','-1small'}});

plottracesel(92,xinfos,'cl3-cl6',plotparam,[],'win',[-3 4],'lfp','se','event','targ','binfo',binfos,'trtable',trtable,{'big'},{'+1fix','+1reward'});
plottracesel(58,xinfos,'cl6',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'trtable',trtable,{'big'},{'-1fix','-1small','-1big'});

plottracesel(94,xinfos,'cl5',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'condtypes',{'big'},{'phase1','phase4'});
plottracesel(92,xinfos,'cl5',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'all'});

plottracesel(102,xinfos,'cl3',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'phase1'});
plottracesel(102,xinfos,'cl3',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'all'});
plottracesel(113,xinfos,'cl3',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'all'});
plottracesel(114,xinfos,'cl3',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'all'});
plottracesel(91,xinfos,'cl3',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'all'});
plottracesel(80,xinfos,'cl3',plotparam,[],'win',[-3 4],'se','event','targ','binfo',binfos,'trtypes',{'big','small'},{'all'});


plottracesel(94,xbinfos,'pulse',plotparam,[],'hist','beh','event','targ','hr','trtypes',{'big','small'},{'all'});
plottracesel(68,xbinfos,'pulse',plotparam,[],'hist','beh','event','targ','hr','trtypes',{'big','small'},{'all'});
plottracesel(68,xbinfos,'pulse',plotparam,[],'hist','beh','event','targ','rr','trtypes',{'big','small'},{'all'});

plottracesel(92,xbinfos,'eye',plotparam,[],'win',[-3 4],'beh','se','event','targ','binfo',binfos,'condition','left','condtypes',{'big'},{'phase1','phase4'});

plottracesel(67,xbinfos,'pulse',plotparam,[],'win',[-3 6],'se','beh','event','targ','pulse','trtypes',{'big','small'},{'all'});
plottracesel(83,xbinfos,'pulse',plotparam,[],'win',[-3 6],'se','beh','event','targ','pulse','trtypes',{'big','small'},{'all'});
plottracesel(92,xbinfos,'pulse',plotparam,[],'win',[-3 6],'se','beh','event','targ','pulse','trtypes',{'big','small'},{'all'});
plottracesel(91,xbinfos,'pulse',plotparam,[],'win',[-3 6],'se','beh','event','targ','pulse','trtypes',{'big','small'},{'all'});

plottracesel(92,binfos,'trt',plotparam,[],'hist','beh','condition','left','trtypes',{'big','small'},{'all'});
plottracesel(92,binfos,'trt',plotparam,[],'hist','beh','condition','right','trtypes',{'big','small'},{'all'});
plottracesel(92,binfos,'rt',plotparam,[],'hist','beh','trtypes',{'big','small'},{'all'});

plottracesel(100,xinfos,'cl1-cl5',plotparam,[],'win',[-3 4],'lfp','se','event','targ','trtable',trtable,{'big'},{'-1fix','-1small','-1big'});


%plotxmultsessfoc(xinfos,plotparam,{'all'},yvar,'daall');
%plotxmultsessfoc(xinfos,plotparam,{'all'},yvar);
%plotxmultsess(xinfos,plotparam);
%plotxbehmultsess(xinfos,xbinfos,plotparam,'evts',{'fix','targ','targeye'});
%plotxbehmultsess(xinfos,xbinfos,plotparam,'evts',{'fix','targ','targeye'},'sestypes',{'reward'});
%plotxbehmultsess(xinfos,xbinfos,plotparam,'daall','evts',{'targ','targeye'},'sestypes',{'reward'});
%plotsummarybeh(xinfos,
plottracesummary({'83'},xinfos,'cl5',plotparam,'condition','left','win',[-1 3]);
plottracesummary({'83'},xinfos,'cl3-cl6',plotparam,'condition','all','lfp','win',[-1 3]);
plottracesummary({'92'},xinfos,'cl5',plotparam,'condition','left','win',[-1 3]);
plottracesummary({'92'},xinfos,'cl3-cl6',plotparam,'condition','all','lfp','win',[-1 3]);
%plottracesummary({'83'},xinfos,'cl5',plotparam,'condition','all','win',[-1 3],'plotz');
plottracesummary({'100'},xinfos,'cl6',plotparam,'condition','left','win',[-1 3]);
plottracesummary({'94'},xinfos,'cl5',plotparam,'condition','left','win',[-1 3]);
plottracesummary({'95'},xinfos,'cl5',plotparam,'condition','all','win',[-1 3]);
plottracesummary({'75'},xinfos,'cl5',plotparam,'condition','all','win',[-1 3]);
plottracesummary({'79'},xinfos,'cl3-cl6',plotparam,'condition','all','lfp','win',[-1 3]);
plottracesummary({'67'},xinfos,'cl5',plotparam,'condition','all','win',[0 3.75],'se');
plottracesummary({'67'},xinfos,'cl3-cl6',plotparam,'condition','all','win',[0 3.75],'se','lfp');
plottracesummary({'67'},xinfos,'cl3-cl6',plotparam,'condition','all','win',[0 3.75],'se','lfp','dapos');
plottracesummary({'67'},xinfos,'cl3-cl6',plotparam,'condition','all','win',[0 3.75],'se','lfp','daneg');

plottracesummary({'92'},xinfos,'pl2',plotparam,'condition','left','win',[-1 4]);
plottracesummary({'68'},xinfos,'cl1-cl5',plotparam,'condition','all','win',[-1 3],'lfp');
plottracesummary({'68'},xinfos,'cl6',plotparam,'condition','all','win',[0 3.75],'se');
plottracesummary({'127'},xinfos,'cl6',plotparam,'condition','left','win',[0 3.75],'se');
plottracesummary({'83'},xinfos,'cl5',plotparam,'condition','all','win',[0 3.75],'se');
plottracesummary({'83'},xinfos,'cl3-cl6',plotparam,'condition','all','win',[0 3.75],'se','lfp');
plottracesummary({'83'},xinfos,'cl3-cl6',plotparam,'condition','all','win',[0 3.75],'se','lfp','dapos','targda','cl5');
plottracesummary({'83'},xinfos,'cl3-cl6',plotparam,'condition','all','win',[0 3.75],'se','lfp','daneg','targda','cl5');
plottracesummary({'92'},xinfos,'cl3-cl6',plotparam,'condition','all','win',[0 3.75],'se','lfp','dapos','targda','cl5');
plottracesummary({'92'},xinfos,'cl3-cl6',plotparam,'condition','all','win',[0 3.75],'se','lfp','daneg','targda','cl5');
plottracesummary({'100'},xinfos,'cl6',plotparam,'condition','all','win',[0 3.75],'se');
plottracesummary({'100'},xinfos,'cl1-cl5',plotparam,'condition','all','win',[0 3.75],'se','lfp','dapos','targda','cl6');
plottracesummary({'100'},xinfos,'cl1-cl5',plotparam,'condition','all','win',[0 3.75],'se','lfp','daneg','targda','cl6');
plottracesummary({'67'},xinfos,'cl5',plotparam,'condition','all','win',[0 3.75],'se');
plottracesummary({'67'},xinfos,'cl3-cl6',plotparam,'condition','all','win',[0 3.75],'se','lfp','dapos','targda','cl5');
plottracesummary({'67'},xinfos,'cl3-cl6',plotparam,'condition','all','win',[0 3.75],'se','lfp','daneg','targda','cl5');


%plot avg behaviors
plottracesummary({'113'},xbinfos,'eye',plotparam,'condition','right','win',[0 3.75],'se','bplot','event','targeye');

plottracesummary({'83'},xbinfos,'lick',plotparam,'condition','all','win',[0 6],'se','bplot','event','targ');

plottracesummary({'83'},xbinfos,'eye',plotparam,'condition','right','win',[0 3.75],'se','bplot','event','targeye','dapos','targda','cl5');
plottracesummary({'83'},xbinfos,'eye',plotparam,'condition','right','win',[0 3.75],'se','bplot','event','targeye','daneg','targda','cl5');
plottracesummary({'114'},xbinfos,'pulse',plotparam,'condition','all','win',[0 10],'se','bplot');

plottracesummary({'67'},xbinfos,'pulse',plotparam,'condition','all','win',[-3 6],'se','bplot','event','targ');
plottracesummary({'83'},xbinfos,'pulse',plotparam,'condition','all','win',[-3 6],'se','bplot','event','targ');

plotmultxtraces(94,xinfos,plotparam,'win',[0 4],'event','targ','ttypes',{'big','small','all'},'sort','damax','bplot',xbinfos,'lick');

%plotdasummary(datm,plotparam,'norm','mean');
%plotdasummary(datm,plotparam,'mean','scattypes','grpreg');

%sensitivity
davalsrew=plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'simple','gain');
davalsmot=plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','small','left'},{'big','small','right'}},'simple','gain');
groupvals=[ davalsmot davalsrew];
plotnormscatter(groupvals,plotparam);
lfpvalsrew=plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'lfp','simple','gain','metric','targimwin');
lfpvalsmot=plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','small','left'},{'big','small','right'}},'lfp','simple','gain','metric','targwin');
groupvalslfp=[ lfpvalsmot lfpvalsrew];
plotnormscatter(groupvalslfp,plotparam);
lfpvalsrew=plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'lfp','gain');

plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','left'},{'small','left'},{'big','right'},{'small','right'}},'plotlines','minmaxnorm');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','ersm'},{'big','erbig'},{'big','fail'},{'big','success'}},'plotlines','event','targ','patra','simple');
groupvals=plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','left'},{'small','left'},{'big','right'},{'small','right'}},'plotlines','event','targ','patra','simple');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'},{'small','left'},{'small','right'}},'plotlines','event','targ','patra','simple');

%whole win
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','left'},{'small','left'},{'big','right'},{'small','right'}},'plotlines','event','targ','patra','metric','targwin','simple');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'small','left'},{'big','right'},{'small','right'}},'plotlines','lfp','event','targ','metric','xxwholewin','xinfo',xinfos,'xwin',[0 4],'patra','constrainpairs',sitelists,'simple');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'},{'small','left'},{'small','right'}},'plotlines','event','targ','metric','targwin','patra','simple');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'},{'small','left'},{'small','right'}},'plotlines','lfp','event','targ','metric','xxwholewin','xinfo',xinfos,'xwin',[0 4],'patra','constrainpairs',sitelists,'simple');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','ersm'},{'big','erbig'},{'big','fail'},{'big','success'}},'plotlines','event','targ','patra','metric','targwin','simple');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','ersm'},{'big','erbig'},{'big','fail'},{'big','success'}},'plotlines','lfp','event','targ','metric','xxwholewin','xinfo',xinfos,'xwin',[0 4],'patra','constrainpairs',sitelists,'simple');

plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'},{'small','left'},{'small','right'}},'plotlines','simple','lfp','event','targ','metric','xxwholewin','xinfo',xinfos,'xwin',[0 4]);
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'},{'small','left'},{'small','right'}},'plotlines','simple','lfp','event','targ','metric','xxlatewin','xinfo',xinfos,'xwin',[1 4]);

%constrained lfp pairs
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','ersm'},{'big','erbig'},{'big','fail'},{'big','success'}},'plotlines','lfp','event','targ','metric','xxearlywin','xinfo',xinfos,'xwin',[0 1],'patra','constrainpairs',sitelists,'simple');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'},{'small','left'},{'small','right'}},'plotlines','simple','lfp','event','targ','metric','xxearlywin','xinfo',xinfos,'xwin',[0 1],'patra','constrainpairs',sitelists);
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'},{'small','left'},{'small','right'}},'plotlines','simple','lfp','event','targ','metric','xxlatewin','xinfo',xinfos,'xwin',[1 4],'patra','constrainpairs',sitelists);
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'small','left'},{'big','right'},{'small','right'}},'plotlines','simple','lfp','event','targ','metric','xxearlywin','xinfo',xinfos,'xwin',[0 1],'patra','constrainpairs',sitelists);

plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'small','left'},{'big','right'},{'small','right'}},'plotlines','lfp','event','targ','metric','xxearlywin','xinfo',xinfos,'xwin',[0 1],'patra','constrainpairs',sitelists);

plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','ersm'},{'big','erbig'},{'big','fail'},{'big','success'}},'plotlines','lfp','event','targ','metric','xxearlywin','xinfo',xinfos,'xwin',[0 1],'patra','simple');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'},{'small','left'},{'small','right'}},'plotlines','simple','lfp','event','targ','metric','xxearlywin','xinfo',xinfos,'xwin',[0 1],'patra');
groupvalslfp=plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'},{'small','left'},{'small','right'}},'plotlines','simple','lfp','event','targ','metric','xxlatewin','xinfo',xinfos,'xwin',[1 4],'patra');
groupvalslfp=plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'small','left'},{'big','right'},{'small','right'}},'plotlines','simple','lfp','event','targ','metric','xxearlywin','xinfo',xinfos,'xwin',[0 1],'patra');


%simple summaires
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'simple','plotlines');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','left'},{'small','left'}},'simple','plotlines');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'simple','plotlines');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'simple','plotlines','metric','targpeakabs');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'metric','targwin');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'lfp','metric','targwin','simple','plotlines');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'lfp','metric','targimwin','simple','plotlines');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'lfp','metric','targpeak','simple','plotlines');

plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'small','left'}},'lfp','metric','targwin','simple','plotlines');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'lfp','metric','targwin','simple','plotlines');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'lfp','metric','targpeak','simple','plotlines');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'lfp','metric','targimwin','simple','plotlines');

plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','success'},{'big','fail'}},'simple','plotlines');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','erbig'},{'big','ersm'}},'simple','plotlines');

plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','success'},{'big','fail'}},'lfp','metric','targimwin','simple','plotlines');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','success'},{'big','fail'}},'lfp','metric','targpeak','simple','plotlines');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','success'},{'big','fail'}},'lfp','metric','targwin','simple','plotlines');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','erbig'},{'big','ersm'}},'lfp','metric','targimwin','simple','plotlines');

plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','-3fix'},{'big','-2fix'},{'big','-1fix'},{'big','-1small'},{'big','-1big'}},'simple');
plotdasummary2(betatm,plotparam,'norm','mean','trtable',trtable,{{'big','-3fix'},{'big','-2fix'},{'big','-1fix'},{'big','-1small'},{'big','-1big'}},'lfp','metric','targimwin','simple');

%Da/beta mag based on reward size oncoming
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}});
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'cleo');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'cleo');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}});
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'metric','targwin');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'scattypes','grpreg');

plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'cleo','metric','targpeakabs');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'lfp','cleo','metric','targwin');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'plotlines','event','targ','metric','targwin','cleo');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'plotlines','event','targ','cleo','simple');


plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','left'},{'small','left'}});
plotdasummary2(datm,plotparam,'mean','ttypes',{{'big','all'},{'small','all'}},'scattypes','grpreg');
plotdasummary2(datm,plotparam,'mean','ttypes',{{'big','left'},{'small','left'}},'scattypes','grpreg');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'metric','targwin');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'lfp');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'small','left'}},'lfp');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'lfp','scattypes','grpreg');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'lfp','scattypes','grpreg','metric','targwin');

%Da/beta mag based on direction
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','small','left'},{'big','small','right'}});
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'metric','targwin');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}});
plotdasummary2(datm,plotparam,'mean','ttypes',{{'big','left'},{'big','right'}},'scattypes','grpreg');

plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'lfp','metric','targwin');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'lfp','metric','targwin','scattypes','grpreg');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'lfp');

%da/beta mag based on history performance failure/succ
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','success'},{'big','fail'}});
plotdasummary2(datm,plotparam,'mean','ttypes',{{'big','success'},{'big','fail'}},'scattypes','grpreg');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','success'},{'big','fail'}},'metric','targwin');
plotdasummary2(datm,plotparam,'mean','ttypes',{{'big','success'},{'big','fail'}},'scattypes','grpreg','metric','targwin');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','success'},{'big','fail'}},'lfp');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','success'},{'big','fail'}},'lfp','scattypes','grpreg');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','success'},{'big','fail'}},'lfp','metric','targwin');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','success'},{'big','fail'}},'lfp','scattypes','grpreg','metric','targwin');

plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left','success'},{'big','left','fail'}},'lfp','metric','targwin');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left','success'},{'big','left','fail'}},'lfp');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','right','success'},{'big','right','fail'}},'lfp');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','right','success'},{'big','right','fail'}},'lfp','metric','targwin');

%da/beta mag based on history performance failure/succ for Fix win
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','success'},{'big','fail'}},'metric','fixpeak','event','fix');
plotdasummary2(datm,plotparam,'mean','ttypes',{{'big','success'},{'big','fail'}},'scattypes','grpreg','metric','fixpeak','event','fix');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','success'},{'big','fail'}},'lfp','metric','fixpeak','event','fix');
plotdasummary2(betatm,plotparam,'mean','ttypes',{{'big','success'},{'big','fail'}},'lfp','scattypes','grpreg','metric','fixpeak','event','fix');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','success'},{'big','fail'}},'lfp','metric','fixwin','event','fix');
plotdasummary2(betatm,plotparam,'mean','ttypes',{{'big','success'},{'big','fail'}},'lfp','scattypes','grpreg','metric','fixwin','event','fix');

plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left','success'},{'big','left','fail'}},'lfp','metric','targwin');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','left','success'},{'big','left','fail'}},'lfp');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','right','success'},{'big','right','fail'}},'lfp');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','right','success'},{'big','right','fail'}},'lfp','metric','targwin');

%Da/beta mag based on history of reward size
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','ersm'},{'big','erbig'}});
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','ersm'},{'big','erbig'}},'metric','targwin');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','ersm'},{'big','erbig'}},'lfp');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','ersm'},{'big','erbig'}},'lfp','metric','targwin');
plotdasummary2(datm,plotparam,'mean','ttypes',{{'big','ersm'},{'big','erbig'}},'scattypes','grpreg');
plotdasummary2(betatm,plotparam,'mean','ttypes',{{'big','ersm'},{'big','erbig'}},'lfp','scattypes','grpreg');
plotdasummary2(betatm,plotparam,'mean','ttypes',{{'big','ersm'},{'big','erbig'}},'lfp','metric','targwin','scattypes','grpreg');

%Da/beta mag based on history of reward sizev for Fix win
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','ersm'},{'big','erbig'}},'metric','fixpeak','event','fix');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','ersm'},{'big','erbig'}},'metric','fixwin','event','fix');
plotdasummary2(datm,plotparam,'mean','ttypes',{{'big','ersm'},{'big','erbig'}},'metric','fixpeak','event','fix','scattypes','grpreg');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','ersm'},{'big','erbig'}},'lfp','metric','fixpeak','event','fix');
plotdasummary2(betatm,plotparam,'mean','ttypes',{{'big','ersm'},{'big','erbig'}},'lfp','metric','fixpeak','event','fix','scattypes','grpreg');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','ersm'},{'big','erbig'}},'lfp','metric','fixwin','event','fix');
plotdasummary2(betatm,plotparam,'mean','ttypes',{{'big','ersm'},{'big','erbig'}},'lfp','metric','fixwin','event','fix','scattypes','grpreg');


%Da/beta mag based on session phase for target period and fix period
%(control is fix period)
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','se1'},{'big','se4'}});
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'small','se1'},{'small','se4'}});
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','se1'},{'big','se4'}},'metric','fixpeak','event','fix');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'small','se1'},{'small','se4'}},'metric','fixpeak','event','fix');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'small','se1'},{'small','se4'}},'metric','rewpeak','event','rew');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','se1'},{'big','se4'}},'metric','rewpeak','event','rew');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','se1'},{'big','se4'}},'lfp');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'small','se1'},{'small','se4'}},'lfp');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'big','se1'},{'big','se4'}},'lfp','metric','targwin');
plotdasummary2(betatm,plotparam,'norm','mean','ttypes',{{'small','se1'},{'small','se4'}},'lfp','metric','targwin');

%timings based on reward size oncoming
plotdasummary2(datm,plotparam,'mean','ttypes',{{'big','all'},{'small','all'}},'metric','damaxts','xinfo',xinfos,'datype','postrials');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'metric','damaxts','xinfo',xinfos,'datype','postrials');
plotdasummary2(datm,plotparam,'mean','ttypes',{{'big','all'},{'small','all'}},'metric','lfpmints','xinfo',xinfos,'datype','postrials');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'metric','lfpmints','xinfo',xinfos,'datype','postrials');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'metric','delt_lfpmin_damax','xinfo',xinfos,'datype','postrials');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','success'},{'big','fail'}},'metric','lfpmints','xinfo',xinfos,'datype','postrials');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','success'},{'big','fail'}},'metric','damaxts','xinfo',xinfos,'datype','postrials');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','ersm'},{'big','erbig'}},'metric','damaxts','xinfo',xinfos,'datype','postrials');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','success'},{'big','fail'}},'metric','delt_lfpmin_damax','xinfo',xinfos,'datype','postrials');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','ersm'},{'big','erbig'}},'metric','delt_lfpmin_damax','xinfo',xinfos,'datype','postrials');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','left'},{'big','right'}},'metric','delt_lfpmin_damax','xinfo',xinfos,'datype','postrials');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','ase1'},{'big','ase4'}},'metric','damaxts','xinfo',xinfos,'datype','postrials');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','ase1'},{'big','ase4'}},'metric','lfpmints','xinfo',xinfos,'datype','postrials');

%xcov parameters based on reward size oncoming
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'metric','minprecoef','xinfo',xinfos,'datype','postrials');
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'metric','damaxts','xinfo',xinfos,'datype','postrials');
plotdasummary2(datm,plotparam,'mean','ttypes',{{'big','all'},{'small','all'}},'metric','mincoef','xcov',xinfos);
plotdasummary2(datm,plotparam,'norm','mean','ttypes',{{'big','all'},{'small','all'}},'metric','mincoefatlag','xcov',xinfos);

%da based on history trial
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','-3break'},{'big','-2break'},{'big','-1break'}});
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','-3fix'},{'big','-2fix'},{'big','-1fix'},{'big','-1big'},{'big','-1small'}});

plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'small','-3fix'},{'small','-2fix'},{'small','-1fix'},{'small','-3target'},{'small','-2target'},{'small','-1target'},{'small','-1big'},{'small','-1small'}});
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','-3fix'},{'big','-2fix'},{'big','-1fix'},{'big','-3target'},{'big','-2target'},{'big','-1target'},{'big','-1big'},{'big','-1small'}});
plotdasummary2(betatm,plotparam,'lfp','norm','mean','trtable',trtable,{{'big','-3fix'},{'big','-2fix'},{'big','-1fix'},{'big','-3target'},{'big','-2target'},{'big','-1target'},{'big','-1big'},{'big','-1small'}});
plotdasummary2(betatm,plotparam,'lfp','norm','mean','trtable',trtable,{{'big','-3fix'},{'big','-2fix'},{'big','-1fix'},{'big','-3target'},{'big','-2target'},{'big','-1target'},{'big','-1big'},{'big','-1small'}},'metric','targwin');

plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','-3fix'},{'big','-2fix'},{'big','-1fix'},{'big','-3target'},{'big','-2target'},{'big','-1target'},{'big','-1big'},{'big','-1small'}},'event','fix','metric','fixpeak');

plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'small','-3break'},{'small','-2break'},{'small','-1break'}});

plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','-3small'},{'big','-2small'},{'big','-1small'}});
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','-3small'},{'big','-2small'},{'big','-1small'},{'big','-1big'},{'big','-2big'},{'big','-3big'}});
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','-3small'},{'big','-2small'},{'big','-1small'},{'big','-1big'},{'big','-2big'},{'big','-3big'},{'big','-3break'},{'big','-2break'},{'big','-1break'}});

plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','small','-3break'},{'big','small','-2break'},{'big','small','-1break'}});
plotdasummary2(betatm,plotparam,'lfp','norm','mean','trtable',trtable,{{'big','small','-3break'},{'big','small','-2break'},{'big','small','-1break'}});
plotdasummary2(betatm,plotparam,'lfp','norm','mean','trtable',trtable,{{'big','-3break'},{'big','-2break'},{'big','-1break'}});
plotdasummary2(betatm,plotparam,'lfp','norm','mean','trtable',trtable,{{'big','-3break'},{'big','-2break'},{'big','-1break'}},'metric','targwin');
plotdasummary2(betatm,plotparam,'lfp','norm','mean','trtable',trtable,{{'big','-3small'},{'big','-2small'},{'big','-1small'},{'big','-1big'},{'big','-2big'},{'big','-3big'},{'big','-3break'},{'big','-2break'},{'big','-1break'}});

%da/beta future trials
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','+1reward'},{'big','+1break'},{'big','+2break'},{'big','+3break'}});
plotdasummary2(betatm,plotparam,'lfp','norm','mean','trtable',trtable,{{'big','+1reward'},{'big','+1break'},{'big','+2break'},{'big','+3break'}});
plotdasummary2(betatm,plotparam,'lfp','norm','mean','trtable',trtable,{{'big','+1reward'},{'big','+1break'},{'big','+2break'},{'big','+3break'}},'metric','targwin')
plotdasummary2(betatm,plotparam,'lfp','norm','mean','trtable',trtable,{{'big','+3reward'},{'big','+2reward'},{'big','+1reward'},{'big','+1break'},{'big','+2break'},{'big','+3break'}},'metric','targwin')
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','+3reward'},{'big','+2reward'},{'big','+1reward'},{'big','+1target'},{'big','+2target'},{'big','+3target'},{'big','+1fix'},{'big','+2fix'},{'big','+3fix'}})
plotdasummary2(betatm,plotparam,'lfp','norm','mean','trtable',trtable,{{'big','+3reward'},{'big','+2reward'},{'big','+1reward'},{'big','+1target'},{'big','+2target'},{'big','+3target'},{'big','+1fix'},{'big','+2fix'},{'big','+3fix'}},'metric','targwin')

%xinfo history trials
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','-3break'},{'big','-2break'},{'big','-1break'},{'big','-1big','-1small'}},'metric','damaxts','xinfo',xinfos);
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','-3small'},{'big','-2small'},{'big','-1small'},{'big','-1big'},{'big','-2big'},{'big','-3big'}},'metric','damaxts','xinfo',xinfos);
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'small','-3big'},{'small','-2big'},{'small','-1big'},{'small','-1small'},{'small','-2small'},{'small','-3small'}},'metric','damaxts','xinfo',xinfos);

plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','-3break'},{'big','-2break'},{'big','-1break'},{'big','-1big','-1small'}},'metric','lfpmints','xinfo',xinfos);
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','-3break'},{'big','-2break'},{'big','-1break'}},'metric','delt_lfpmin_damax','xinfo',xinfos);

%xinfo future tirals
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','+2reward'},{'big','+1reward'},{'big','+1break'},{'big','+2break'},{'big','+3break'}},'metric','damaxts','xinfo',xinfos);
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','+2reward'},{'big','+1reward'},{'big','+1break'},{'big','+2break'},{'big','+3break'}},'metric','lfpmints','xinfo',xinfos);
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','+3reward'},{'big','+2reward'},{'big','+1reward'},{'big','+1break'},{'big','+2break'},{'big','+3break'}},'metric','delt_lfpmin_damax','xinfo',xinfos);
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','small','+2reward'},{'big','small','+1reward'},{'big','small','+1break'},{'big','small','+2break'}},'metric','delt_lfpmin_damax','xinfo',xinfos);
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','small','+1reward'},{'big','small','+1break'}},'metric','delt_lfppostmax_damax','xinfo',xinfos);
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','small','+1reward'},{'big','small','+1break'}},'metric','delt_lfpmin_damax','xinfo',xinfos);
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','small','+1reward'},{'big','small','+1break'}},'metric','delt_lfpmin_damax','xinfo',xinfos,'grpreg');
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','+3reward'},{'big','+2reward'},{'big','+1reward'},{'big','+1target'},{'big','+2target'},{'big','+3target'},{'big','+1fix'},{'big','+2fix'},{'big','+3fix'}},'metric','delt_lfpmin_damax','xinfo',xinfos);
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','+3reward'},{'big','+2reward'},{'big','+1reward'},{'big','+1target'},{'big','+2target'},{'big','+3target'},{'big','+1fix'},{'big','+2fix'},{'big','+3fix'}},'metric','damaxts','xinfo',xinfos);
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','+3reward'},{'big','+2reward'},{'big','+1reward'},{'big','+1target'},{'big','+2target'},{'big','+3target'},{'big','+1fix'},{'big','+2fix'},{'big','+3fix'}},'metric','lfpmints','xinfo',xinfos);
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','-3big'},{'big','-2big'},{'big','-1big'},{'big','-1target'},{'big','-2target'},{'big','-3target'},{'big','-1fix'},{'big','-2fix'},{'big','-3fix'}},'metric','damaxts','xinfo',xinfos);
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','-3big'},{'big','-2big'},{'big','-1big'},{'big','-1target'},{'big','-2target'},{'big','-3target'},{'big','-1fix'},{'big','-2fix'},{'big','-3fix'}},'metric','delt_lfpmin_damax','xinfo',xinfos);
plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','left','-3big'},{'big','left','-2big'},{'big','left','-1big'},{'big','left','-1target'},{'big','left','-2target'},{'big','left','-3target'},{'big','left','-1fix'},{'big','left','-2fix'},{'big','left','-3fix'}},'metric','delt_lfpmin_damax','xinfo',xinfos);


plotdasummary2(datm,plotparam,'norm','mean','trtable',trtable,{{'big','small','+1reward'},{'big','small','+1break'}},'metric','delt_lfpmin_damax','xinfo',xinfos,'grpreg');

plotbehcorsummary(xinfos,xbinfos,plotparam,'daall','evts','targ','ttypes',{'big','small','all'});
plotbehcorsummary(xinfos,xbinfos,plotparam,'daall','evts','targeye','ttypes',{'big','small','all'});
plotbehcorsummary(xinfos,xbinfos,plotparam,'daall','evts','targeye','ttypes',{'big','small','left'});
plotbehcorsummary(xinfos,xbinfos,plotparam,'daall','evts','targeye','ttypes',{'big','small','right'});
plotbehcorsummary(xinfos,xbinfos,plotparam,'daall','evts','fix','ttypes',{'big','small','all'});
plotbehcorsummary(xinfos,xbinfos,plotparam,'daall','evts','fixeye','ttypes',{'big','small','all'});


plotdasummarylr(datm,plotparam,'norm','mean');
plotdasummarylr(datm,plotparam,'mean','scattypes','grpreg');

plotsummarycounts(respgrouped,plotparam,'mean');
plotsummarycounts(respgrouped,plotparam,'scattypes','grpreg','cond',{{'big','non'},{'sma','non'}});
plotsummarycounts(respgrouped,plotparam,'scattypes','grpreg','cond',{{'big','non'},{'big','big'}});

plotsummarycounts(respgrouped,plotparam,'scattypes','grpreg','cond',{{'big','big'},{'sma','non'}});
plotsummarycounts(respgrouped,plotparam,'scattypes','grpreg','cond',{{'big','non'},{'sma','sma'}});

plotxsummary(xinfos,plotparam);
plotxsummary(xinfos,plotparam,'scattypes','grpreg','sesstypes',{'big','small'},'conditions',{'all','all'});
plotxsummary(xinfos,plotparam,'scattypes','grpreg','sesstypes',{'big','small'},'conditions',{'left','left'});
plotxsummary(xinfos,plotparam,'scattypes','grpreg','sesstypes',{'big','big'},'conditions',{'left','right'});
plotxsummary(xinfos,plotparam,'scattypes','grpreg','sesstypes',{'small','small'},'conditions',{'left','right'});

plotxsummary(xinfos,plotparam,'norm');

plotbehsummary(datm,plotparam);
%{
 plotxtracesavgov(xinfos{plotses},trialinfo,plotparam,'evtype','targ','sesstypes',...
        {'bigreward','smallreward'},'trialtypes',{'sleft'},'daall','win',[-1 4]);
%}
%{
%temp comment 1/21/2019
plotxmultsess(xinfos,plotparam,'daall');
plotxmultsess(xinfos,plotparam,'daneg');        %plot da neg

plotxbehmultsess(xinfos,binfos,plotparam,'daall');

plotxbehmultsess(xinfos,binfos,plotparam);
plotxbehmultsess(xinfos,binfos,plotparam,'daneg');

%}

%plot multiple da neg
%plot corrs of behavior to xinfo variables
end