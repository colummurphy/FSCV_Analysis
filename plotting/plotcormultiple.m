function outdata=plotcormultiple(sessnums,trlists,datms,lfptms,binfos,plotparam,varargin)
%multiple sessions
sessiontypeids={'bigreward','smallreward','targetbreak','fixbreak'};
plottypes={'bigreward','smallreward'};
condtype='all';
trlistsnew={};

figpos=[50,50,1600,900];
hf=figure('visible','off');     %figure for each channel
if ispc
hf=figure('visible','on');     %figure for each channel
end
set(hf,'position',figpos,'color',[1 1 1]);
set(0,'CurrentFigure',hf);    %set figure handle to current figure
clf(hf,'reset');
set(hf,'color',[1 1 1]);
 
sessmark={'sq','o','^','+','v','*','p','h','x','d','<','>','o','^','v','<','>','+','*','p','h','x','d','<','>'};
sessmarktxt={' sq','\o','\Delta','+','\nabla','*','p','h','\times','d','<','>','.','\Delta','\nabla','<','>','+','*','p','h','x','d','<','>'};
marksize=[200 200 150 200 200 200 200 200 200 200 200 200 50 50 50 50 50 50 50 50 50 50 50 200 200 200 200 200 200 200 200 200 200 200 200 200];
cmapp=[1 .5 0;1 .55 0; 1 .6 0; 1 .65 0; 1 .7 0; 1 .75 0; 1 .8 0; 1 .85 0; 1 .9 0];
cmapc=[1 .5 0;1 .55 0; 1 .6 0; 1 .65 0; 1 .7 0; 1 .75 0; 1 .8 0; 1 .85 0; 1 .9 0];
cmap2=cool;
cmap1=hot;
markc=cmap1(1:3:end,:);
markp=cmap2(1:2:end,:);

xoff=100;
yoff=100;
mar=25;
axsize=[400,250];
fontsize=15;
triallim=[];
bigtype=[];
smalltype=[];
targtype=[];
fixtype=[];
plotvar={};
argnum=1;
siteda='';
sitelfp='';
bflag=0;
plotrt=0;
plothr=0;
ploteye=0;
linfos={};
xinfos={};
xbinfos={};
hfsametype=0;
metric='targpeak';
metriclfp='targwin';
event='targ';
win = [0 4];
winb= [.1 4];
rate=10;            %sample rate for xinfo data
%metriclfp='targpeakabs';
plotnum=0;
plothist=[];
plottypes=[];
wine=[];
winhr=[.1 4];
winl=[];
xvars={};
yvars={};
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'dax'
            %values from xinfo
            argnum=argnum+1;
            siteda=varargin{argnum};
            argnum=argnum+1;
            xinfos=varargin{argnum};
            plotnum=plotnum+1;
            plotvars{plotnum}='da';
        case 'datm'
            %values from datm, supplied instead of xinfos
            argnum=argnum+1;
            siteda=varargin{argnum};
            plotnum=plotnum+1;
            plotvars{plotnum}='da';
        case 'lfpx'
            argnum=argnum+1;
            sitelfp=varargin{argnum};            
            argnum=argnum+1;
            xinfos=varargin{argnum};
            plotnum=plotnum+1;
            plotvars{plotnum}='lfp';
        case 'lfptm'
            %values from betatm
            argnum=argnum+1;
            sitelfp=varargin{argnum};
            plotnum=plotnum+1;
            plotvars{plotnum}='lfp';
        case 'behx'
            %values from xbinfos
            argnum=argnum+1;
            xbinfos=varargin{argnum};
        case 'metric'
            %for datm values
            argnum=argnum+1;
            metric=varargin{argnum};
        case 'metriclfp'
            %for betatm values
            argnum=argnum+1;
            metriclfp=varargin{argnum};
        case 'trialseq'
            plotnum=plotnum+1;
            plotvars{plotnum}='trialnum';
        case 'rt'
            plotrt=1;
            plotnum=plotnum+1;
            plotvars{plotnum}='trt';
        case 'frt'
            plotrt=1;
            plotnum=plotnum+1;
            plotvars{plotnum}='frt';
        case 'hr'
            plothr=1;
            plotnum=plotnum+1;
            plotvars{plotnum}='hr';
                    case 'rr'
            plothr=1;
            plotnum=plotnum+1;
            plotvars{plotnum}='rr';
                    case 'rrstd'
            plothr=1;
            plotnum=plotnum+1;
            plotvars{plotnum}='rrstd';
                    case 'rmssd'
            plothr=1;
            plotnum=plotnum+1;
            plotvars{plotnum}='rmssd';
        case 'lick'
            plotnum=plotnum+1;
            plotvars{plotnum}='lick';
        case 'eye'
            bflag=1;
            ploteye=1;
            plotnum=plotnum+1;
            plotvars{plotnum}='eye';
        case 'event'
            argnum=argnum+1;
            event=varargin{argnum};    
        case 'hf'
            argnum=argnum+1;
            plothist=varargin{argnum};
        case 'hfsametype'
            hfsametype=1;       %constrict types of futre/past trials
        case 'types'
            argnum=argnum+1;
            plottypes=varargin{argnum};
        case 'condition'
            argnum=argnum+1;
            condtype=varargin{argnum};
        case 'bwin'
            argnum=argnum+1;
            winb=varargin{argnum};
        case 'win'
            argnum=argnum+1;
            win=varargin{argnum};
        case 'lwin'
            argnum=argnum+1;
            winl=varargin{argnum};
        case 'hrwin'
            argnum=argnum+1;
            winhr=varargin{argnum};
        case 'ewin'
            argnum=argnum+1;
            wine=varargin{argnum};
        case 'trlists'
            argnum=argnum+1;        %already obtained trlistsnew
            trlistsnew=varargin{argnum};
        case 'xvars'
            argnum=argnum+1;
            xvars=varargin{argnum};
        case 'yvars'
            argnum=argnum+1;
            yvars=varargin{argnum};
    end
    argnum=argnum+1;
end     
plotvarsx={'lfp','da','trialnum'};
plotvarsy={'trialnum','hr','rrstd','rmssd','lick','eye','trt','frt'};
if ~isempty(xvars)
    plotvarsx=xvars;
end
if ~isempty(yvars)
    plotvarsy=yvars;
end
sessiontypes=plottypes;
savepath=fullfile(plotparam.savepath);
labels=['corall_' event '_da_' metric '_lfp_' metriclfp];
sesslab=sessiontypes{1};
if length(sessiontypes)>1
for is=2:length(sessiontypes)
    sesslab=[sesslab '_' sessiontypes{is}];
end
end
labels=[labels '_' sesslab];
if ~isempty(winhr)
    labels=[labels '_winhr_' num2str(winhr(1)*rate) '_' num2str(winhr(2)*rate)];
end
if ~isempty(wine)
    labels=[labels '_wine_' num2str(wine(1)*rate) '_' num2str(wine(2)*rate)];
end
if ~isempty(winl)
    labels=[labels '_winl_' num2str(winl(1)*rate) '_' num2str(winl(2)*rate)];
end
if isempty(wine)
    wine=winb;
end
if isempty(winhr)
    winhr=winb;
end
if isempty(winl)
    winl=winb;
end
%{
for is=1:length(sessnums)
    text(axa,is*75-200,axpos(4)+50,[num2str(sessnums(is)) ' ' sessmarktxt{is}],'interpreter','tex','units','pixels','fontweight','bold','color',[0 0 0]);
    ranj(is)=rand(1,1)*.8-.5;
end
%}
axa={};
axpos={};
for ip=1:2
        axa{ip}=subplot(1,2,ip);   hold(axa{ip},'on');
        set(axa{ip},'units','pixels');
        axis('square')         
        set(axa{ip},'ytick',1:length(plotvarsy));
        set(axa{ip},'yticklabel',plotvarsy,'fontsize',fontsize);
        set(axa{ip},'ylim',[0 length(plotvarsy)+.5]);
        set(axa{ip},'TickLabelInterpreter','none');
        set(axa{ip},'xtick',1:length(plotvarsx));
        set(axa{ip},'xticklabel',plotvarsx,'fontsize',fontsize);
        set(axa{ip},'xlim',[0 length(plotvarsx)+.5]);
        %set(axa{ip},'xTickLabelRotation',90)
        axpos{ip}=get(axa{ip},'position');
        for iline=1:length(plotvarsy)       
        aa=plot(axa{ip},[0 length(plotvarsx)], ...
            [iline iline],'linewidth',1,'color',[0 0 0]);
        aa.Color(4)=.2;
        end    
                for iline=1:length(plotvarsx)
        aa=plot(axa{ip},[iline iline], ...
            [0 length(plotvarsy)],'linewidth',1,'color',[0 0 0]);
        aa.Color(4)=.2;
                end   
        set(axa{ip},'position',[axpos{ip}(1) axpos{ip}(2)-150 axpos{ip}(3) axpos{ip}(4)]);
end   
  axsesslegend=axes;   hold(axsesslegend,'on');
%  set(axsesslegend,'units','pixels');
legpos=get(axsesslegend,'position');
set(axsesslegend,'position',[.05,.9,.82,.09])
set(axsesslegend,'box','off','visible','off')
text(axsesslegend,-1,.5,'sess # :','interpreter','tex','fontweight','bold','color',[0 0 0],'fontsize',fontsize);
for is=1:length(sessnums)
   % text(axa{1},is*75,axpos(4)+150,[num2str(sessnums(is)) ' ' sessmarktxt{is}],'interpreter','tex','units','pixels','fontweight','bold','color',[0 0 0]);
    scatter(axsesslegend,is,.5,marksize(is)-5,sessmark{is},'markeredgecolor',[0 0 0]);
    text(axsesslegend,is-.5,.5,num2str(sessnums(is)),'interpreter','tex','fontweight','bold','color',[0 0 0],'fontsize',fontsize-2);
    ranj(is)=rand(1,1)*.8-.5;
end
ylim(axsesslegend,[0 1]);
titletext=[sesslab ' | ' event ' | da '  metric ' | lfp ' metriclfp ];  
text(axsesslegend,-50,legpos(2)+legpos(4)+70,titletext,'units','pixels','fontweight','bold','fontsize',fontsize,'interpreter','none');

%set up plots
if isempty(plothist)
plothist=0;       %1 means 1 trial in future value
end
if isempty(plottypes)
plottypes={'all'};
end

count=1;
if isempty(trlistsnew)
    
for ise=1:length(sessnums)
%populate trlists for all sessnums
sessnum=sessnums(ise);
sessid=num2str(sessnum);
trlistid=find(strcmp({trlists.sessid},sessid));
trlist=trlists(trlistid).list;
numtrials=length(trlist);
trialgrps=plotparam.trialgrps;
sessnums=plotparam.sessnums;
targdasites=plotparam.dasites;
sites=getsites(sessnums,targdasites);
uniquesites=unique({sites(1:end).site});
[dapair,lfppair]=getsitepairs(targdasites);
targlsites=plotparam.lfpchs;
lsites=getlfpsites(sessnums,targlsites);
targses=find(strcmp({trialgrps.sessid},sessid));
trialinfo=trialgrps(targses).trialinfo;

targdas=find([sites.sessnum]==sessnum);
for ida=1:length(targdas)
    siteda=sites(targdas(ida)).site;
    sitedalab=sites(targdas(ida)).probeid;
    sitedach=sites(targdas(ida)).ch;
    [dapair,lfppair]=getsitepairs(sitedalab);
    sitelfplab=lfppair(1);
        lsites=getlfpsites(sessnum,sitelfplab);
        if isempty(lsites)
    lsites=getlfpsites(sessnum,lfppair);
        end
    sitelfplab=lsites.probeid;
    sitelfp=lsites(1).site;
    sitelfpch=lsites(1).ch;
    
%initialize trlist
for ix=1:length(trlist)
    trlist(ix).da=nan;
    trlist(ix).lfp=nan;   
    trlist(ix).hr=nan;   
    trlist(ix).eye=nan;   
    trlist(ix).lick=nan;   
    trlist(ix).rr=nan;   
    trlist(ix).rrstd=nan;   
    trlist(ix).rmssd=nan;   
end
%first get good trials list for trial type big/small/targ
maplist={};     %listing of rows mapped by sel trials; ie. seltrials{1}(1) --> maplist{1}(1) = targ row for specific type & trial
seltrials={};
dtm={};
lfptm=[];
xtargda=[];
xtarglfp=[];
goodtrials=[];
labx='';
laby='';
for itype=1:length(sessiontypes)
    itypeid=find(contains(sessiontypeids,sessiontypes(itype)));
    if strcmp(condtype,'all')
    goodtrialgrp=find(contains(trialinfo(itypeid).names,'all')==1);
    incswgrp=find(strcmp(trialinfo(itypeid).names,'switch')==1);
    goodtrials=sort(unique([trialinfo(itypeid).nums{goodtrialgrp} trialinfo(itypeid).nums{incswgrp}]));  
    else
    goodtrialgrp=find(contains(trialinfo(itypeid).names,condtype)==1);
    goodtrials=sort(unique([trialinfo(itypeid).nums{goodtrialgrp}]));          
    end
    seltrials{itype}=goodtrials;
    tlabelsgood=goodtrials+99;              %convert to trial labels ie 100,101,etc.
    typelistids=find(strcmp({trlist.type},sessiontypes{itype}));
    maplist(itype).trlist=find(ismember([trlist.id],tlabelsgood) & ...
        contains({trlist.type},sessiontypes{itype}));
    maplist(itype).trialnums=goodtrials;
    maplist(itype).type=sessiontypes{itype}; 
   
%get signal values 
    dtm{itype}=datms{targses}{itypeid}{sitedach};           %get all da targeted values for current trial type big/small
    %if ~isequal(dtm{itype}.trialnums,maplist(itype).trialnums)
        %fewer trials in data variable than listed good trials, reset to match
        %measured good values
        maplist(itype).trlist=find(ismember([trlist.id],dtm{itype}.trialnums+99) & ...
            contains({trlist.type},sessiontypes{itype}));
        maplist(itype).trialnums=dtm{itype}.trialnums;
   % end
    targtrials=find(ismember(dtm{itype}.trialnums,maplist(itype).trialnums));
    %user original trials values in dtm which are good
    %targtrials=1:length(dtm{itype}.trialnums);
    datatemp=getfield(dtm{itype},metric);
    dadata=datatemp(targtrials);
    %puttrials=find(ismember([trlist.id],dtm{itype}.trialnums+99)& ...
    %        contains({trlist.type},sessiontypes{itype}));
    for it=1:length(dadata)
        trlist(maplist(itype).trlist(it)).da=dadata(it);
    end

%targeted beta values from betatm
    lfptm=lfptms{targses}{itypeid};
    ltm{itype}=lfptm{sitelfpch};           %get all da targeted values for current trial type big/small
    %if ~isequal(ltm{itype}.trialnums,maplist(itype).trialnums)
        %fewer trials in data variable than listed good trials, reset to match
        %measured good values
        maplist(itype).trlist=find(ismember([trlist.id],ltm{itype}.trialnums+99) & ...
            contains({trlist.type},sessiontypes{itype}));
        maplist(itype).trialnums=ltm{itype}.trialnums;
   % end
    targtrials=find(ismember(ltm{itype}.trialnums,maplist(itype).trialnums));
    ltemp=getfield(ltm{itype},metriclfp);
    ldata=ltemp(targtrials);
    for it=1:length(ldata)
        trlist(maplist(itype).trlist(it)).lfp=ldata(it);
    end   
    
%get hr  from binfos
hr=[];
rr=[];
rrstd=[];
rmssd=[];
hrtrials=[];
if isempty(xbinfos)
    btarg=find(contains({binfos.sessiontype},sessiontypes{itype}) &...
        strcmp({binfos.sessionid},sessid) &...
        strcmp({binfos.evt},event));            %event matters since where calculated
    hr=binfos(btarg).pulse;
    hrtrials=binfos(btarg).seltrials;
else 
    %get hr from xbinfos according to b window defined
     targrows=find(contains({xbinfos.sessiontype},sessiontypes{itype}) &...
         strcmp({xbinfos.sitelfp},'pulse') &...
         strcmp({xbinfos.siteda},sitedalab) &...
        strcmp({xbinfos.sessionid},sessid) &...
        strcmp({xbinfos.event},event));            %event matters since where calculated
    if ~isempty(targrows)
    temp=xbinfos(targrows(1)).daall.lfptracesaln;
    hrtrials=xbinfos(targrows(1)).daall.trials;
        alnidx=xbinfos(targrows(1)).daall.mididx;
            hhwin=alnidx+winhr(1)*rate:alnidx+winhr(2)*rate;
    hrwin=temp(:,hhwin); 
    hr=nanmean(hrwin,2);
   rrtemp=1./temp.*60; 
   rrwin=rrtemp(:,hhwin);
   rr=nanmean(rrwin,2);
    rrstd=nanstd(rrwin,[],2);
    rmssd=rms(diff(rrtemp(:,hhwin),1,2),2);    
    end    
end
for it=1:length(hrtrials)
    tid=find(ismember(maplist(itype).trialnums,hrtrials(it)));
    if ~isempty(tid)
        trlist(maplist(itype).trlist(tid)).hr=hr(it);
        if ~isempty(xbinfos)
        trlist(maplist(itype).trlist(tid)).rr=rr(it);
        trlist(maplist(itype).trlist(tid)).rrstd=rrstd(it);
        trlist(maplist(itype).trlist(tid)).rmssd=rmssd(it);
        end
    end
end
%get lick  from binfos
lick=[];
licktrials=[];
if isempty(xbinfos)
    btarg=find(contains({binfos.sessiontype},sessiontypes{itype}) &...
        strcmp({binfos.sessionid},sessid) &...
        strcmp({binfos.evt},event));            %event matters since where calculated
    lick=binfos(btarg).lickpost;
    licktrials=binfos(btarg).seltrials;  
else
    %get  from xbinfos according to b window defined
     targrows=find(contains({xbinfos.sessiontype},sessiontypes{itype}) &...
         strcmp({xbinfos.sitelfp},'lick') &...
         strcmp({xbinfos.siteda},sitedalab) &...
        strcmp({xbinfos.sessionid},sessid) &...
        strcmp({xbinfos.event},event));            %event matters since where calculated
        if ~isempty(targrows)

    temp=xbinfos(targrows(1)).daall.lfptracesaln;
    licktrials=xbinfos(targrows(1)).daall.trials;
            alnidx=xbinfos(targrows(1)).daall.mididx;

    licks=temp(:,alnidx+winl(1)*rate:alnidx+winl(2)*rate); 
    lick=nanmean(licks,2);
        end
end
for it=1:length(hrtrials)
    tid=find(ismember(maplist(itype).trialnums,licktrials(it)));
    if ~isempty(tid)
        trlist(maplist(itype).trlist(tid)).lick=lick(it);
    end
end 
    
%eye diameter, only right side valid
eye=[];
etrials=[];
btarg=find(contains({binfos.sessiontype},sessiontypes{itype}) &...
    strcmp({binfos.sessionid},sessid) &...
    strcmp({binfos.evt},event)); 
if isempty(xbinfos)    
           %event matters since where calculated
    eye=-binfos(btarg).reyed;
    etrials=binfos(btarg).seltrialsr;    
else
    %get  from xbinfos according to b window defined
     targrows=find(contains({xbinfos.sessiontype},sessiontypes{itype}) &...
         strcmp({xbinfos.sitelfp},'eye') &...
         strcmp({xbinfos.siteda},sitedalab) &...
        strcmp({xbinfos.sessionid},sessid) &...
        strcmp({xbinfos.event},event));     
        if ~isempty(targrows)
    etrialsall=xbinfos(targrows(1)).daall.trials;
    etrialsr=binfos(btarg).seltrialsr;      %right side only
    etrialsel=find(ismember(etrialsall,etrialsr));
    etrials=intersect(etrialsall,etrialsr);
    temp=xbinfos(targrows(1)).daall.lfptracesaln(etrialsel,:);
            alnidx=xbinfos(targrows(1)).daall.mididx;
    eyewin=-temp(:,alnidx+wine(1)*rate:alnidx+wine(2)*rate); 
    eye=nanmean(eyewin,2);
        end
end
for it=1:length(etrials)
    tid=find(ismember(maplist(itype).trialnums,etrials(it)));
    if ~isempty(tid)
        trlist(maplist(itype).trlist(tid)).eye=eye(it);
    end
end


end
%finished scrolling session types for given session num & given site
%save to origianl trlists
trlistsnew(count).list=trlist;
trlistsnew(count).sessid=sessid;
trlistsnew(count).dasite=siteda;
trlistsnew(count).lfpsite=sitelfp;
trlistsnew(count).dalab=sitedalab;
trlistsnew(count).lfplab=sitelfplab;
count=count+1;
end
%finish scrolling all sites in sess
end
%finish compililng all sessions
end

sitecolor=[1 0 0];
sitecolorn=[0 .5 1];
boxpos=0:1/length(trlistsnew):1;
axleg={};
  axleg{1}=axes;   hold(axleg{1},'on');
  axleg{2}=axes;   hold(axleg{2},'on');
%  set(axsesslegend,'units','pixels');
legpos=get(axleg{1},'position');
set(axleg{1},'position',[.05,.7,.45,.2])
set(axleg{2},'position',[.52,.7,.45,.2])

set(axleg{1},'box','off','visible','off')
set(axleg{2},'box','off','visible','off')

ylim(axleg{1},[0 1]);
ylim(axleg{2},[0 1]);

countc=0;
countp=0;
for ip=1:length(trlistsnew)
sessid=trlistsnew(ip).sessid;
siteda=trlistsnew(ip).dasite;
sitelfp=trlistsnew(ip).lfpsite;
ise=find(ismember(sessnums,str2num(sessid)));
m=sessmark{ise};
iip=2;
if contains({trlistsnew(ip).dasite},'c')
    iip=1;
    countc=countc+1;
    count=countc;
else
    countp=countp+1;
    count=countp;
end
scatter(axleg{iip},count-.5,.85,marksize(ise),sessmark{ise},'markeredgecolor',[0 0 0]);
    text(axleg{iip},count-.75,.6,siteda,'interpreter','tex','fontweight','bold','color',[0 0 0],'fontsize',fontsize-5,'rotation',90);
    text(axleg{iip},count-.75,0.15,sitelfp,'interpreter','tex','fontweight','bold','color',[0 0 0],'fontsize',fontsize-5,'rotation',90);
end
    text(axleg{1},-2.5,.65,'da site','interpreter','tex','fontweight','bold','color',[0 0 0],'fontsize',fontsize-2);
    text(axleg{1},-2.5,0.25,'lfp site','interpreter','tex','fontweight','bold','color',[0 0 0],'fontsize',fontsize-2);
  

%plotting
for ip=1:length(trlistsnew)
%populate trlists for all sites all sess
sessid=trlistsnew(ip).sessid;
ise=find(ismember(sessnums,str2num(sessid)));
trlist=trlistsnew(ip).list;
numtrials=length(trlist);

bigtype=double(strcmp({trlist.type},'bigreward'));
bigtype(bigtype==0)=nan;
smalltype=double(strcmp({trlist.type},'smallreward'));
smalltype(smalltype==0)=nan;
targtype=double(strcmp({trlist.type},'targetbreak'));
targtype(targtype==0)=nan;
fixtype=double(strcmp({trlist.type},'fixbreak'));
fixtype(fixtype==0)=nan;
sidetype=[trlist.side];

for ix = 1:length(plotvarsx)
curx=plotvarsx{ix};
for iy=1:length(plotvarsy)
cury=plotvarsy{iy};
datax=[];
datay=[];
if strcmp(cury,'trialnum')
    datay=1:length(trlist);
else
    datay=[trlist.(cury)];
end  
if strcmp(curx,'trialnum')
    datax=1:length(trlist);
else
    datax=[trlist.(curx)];
end 
%datax=[trlist.(curx)];       
if strcmp(cury,'eye')
    normval=abs(nanmedian(datay));
    datay=datay./normval;
end
%datax=circshift(datax,[0 plothist(iip)]);       %shift -1 for future trial, +1 for past trial
datay=circshift(datay,[0 plothist]); %shift beh signal -1 for future trial, +1 for past trial
databigy=bigtype;
databigx=bigtype;
datasmally=smalltype;
datasmallx=smalltype;
if hfsametype
 %constrict to same type for future/past trial
 databigy=circshift(bigtype,[0 plothist]);
 datasmally=circshift(smalltype,[0 plothist]);
end  
xoutthres=nanmean(datax)+4*nanstd(datax);
youtthres=nanmean(datay)+4*nanstd(datay);
remx=[];
remy=[];
if ~strcmp(cury,'eye') && ~strcmp(cury,'lfp') && ~strcmp(cury,'da')
remy=find(abs(datay)>youtthres );
end
if strcmp(cury,'frt')
    remy2=find(datay<.05);
    remy=unique([remy remy2]);
end
if ~isempty(remy)
datay(remy)=nan;
end
nonnan=find(~isnan(datax) & ~isnan(datay));
if hfsametype
    nonnan=find(~isnan(datax) & ~isnan(datay) & (~isnan(datasmallx) | ~isnan(databigx)) & (~isnan(datasmally) | ~isnan(databigy)));
end
datax=datax(nonnan);
datay=datay(nonnan);
[r,psig]=corr(datax',datay');
%plot if sig (p<=0.05)
if psig<=0.05
m=sessmark{ise};
c=sitecolor;
fsize=40*log10(1/psig);
if r<0
c=sitecolorn;
end
iip=2;
if contains(trlistsnew(ip).dasite,'c')
    iip=1;
end
a=scatter(axa{iip},ix-boxpos(ip),iy-boxpos(ip),...
fsize,m,'markeredgecolor',c,'MarkerEdgeAlpha',.3,'linewidth',1.5,'markerfacecolor','none');
end

end
      
end
end

savename=[savepath '_' labels];
outdata=trlistsnew;
trlists=trlistsnew;
save([savepath 'trlists_' labels],'trlists','winhr','wine','winl','winb','-v7.3');
%savename=[savepath '_' event];
savefig(hf,savename);
saveas(hf,savename,'jpg')
%print(hf,savename,'-painters','-depsc');
end


