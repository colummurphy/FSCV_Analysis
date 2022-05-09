function [tabdata,sitesigs,mdata]=corsmultiple(sessnums,plotparam,xinfos,binfos,xbinfos,trlists,datm,betatm,varargin)
%CORS ALL SESSIONS & TABULATE
%run loop for plotcortrials2 to get info for each session/site pairs
    %all valid da sites
argnum=1;
flaglfptm=0;
targevent='targeye';
metricda='targpeak';
metriclfp='targimwin';
xwin{1}=[0 1];
offset=[0 0];
figpos=[50,50,1500,600];
axa={};
axsiz=[350 300];
off=75;
mar=100;
defaultmarksize=30;
condtypes={'all','all','all'};
fontsize=14;
rate=10;
 alphac=1;
lwid=1;
ranj=[];
constrainpairs=0;
condition='all';
xvarflag=0;
behtype='eye';
beh=0;
relfp=0;
daonly=0;
lfponly=0;
ewin=[0 0];
hrwin=[0 0];
lwin=[0 0];
sitelists={};
mregress=0;
countm=1;
mdata=[];
cormatrix=0;
cordata={};
testcor=0;
mdata={};
sitesigs=[];
plothist=0;
getcondition='';
xmat={};

while argnum<=length(varargin)
    switch varargin{argnum}
        case 'hist'
            plothist=1; %make histogram s of r / p values, pos/neg r's
        case 'uselfptm'
            %values from lfptm
            flaglfptm=1;
        case 'event'
            argnum=argnum+1;
            targevent=varargin{argnum};
        case 'win'
            argnum=argnum+1;
            relfp=1;
            xwin=varargin{argnum};
        case 'constrainpairs'
            constrainpairs=1;
        case 'condtypes'
            argnum=argnum+1;
            condtypes=varargin{argnum};
        case 'condition'
            argnum=argnum+1;
            condition=varargin{argnum};
        case 'getcondition'
            %must also provide tabdata variable with tabdata(1).big ,
            %tabdata(1).small etc trial indeces for task conditions
            argnum=argnum+1;
            getcondition=varargin{argnum};
        case 'xvarflag'
            %use xinfo variable as defined by metriclfp
            xvarflag=1;
            xwin{1}=[0 0 ];
        case 'metriclfp'
            argnum=argnum+1;
            metriclfp=varargin{argnum};
        case 'beh'
            argnum=argnum+1;
            beh=1;
            behtype=varargin{argnum};
        case 'daonly'
            daonly=1;
        case 'lfponly'
            lfponly=1;
        case 'ewin'
            %window for eye behavior
            argnum=argnum+1;
            ewin=varargin{argnum};
        case 'lwin'
            %window for eye behavior
            argnum=argnum+1;
            lwin=varargin{argnum};
        case 'hrwin'
            %window for eye behavior
            argnum=argnum+1;
            hrwin=varargin{argnum};
        case 'sitelists'
            argnum=argnum+1;
            sitelists=varargin{argnum}; %already generated sitelists, otherwise takes long to retrieve below
        case 'multiregress'
            mregress=1;
            condtypes=repmat({'all'},1,10);
        case 'cormatrix'
            %cor matrix da/beta with beta time windows
            cormatrix=1;
                        condtypes=repmat({'all'},1,10);
        case 'cordata'
            argnum=argnum+1;
            cordata=varargin{argnum};       %already got tabdata on last run provided here
        case 'testcor'
            %do multiple regression on those with significant neg cor using wins for CN/Put determined by cormatrix CN (.2-3.2) putamen (.2-2.6)
            %05/01 DG: shoot down null hypothesis that da vs beta corr can be fully
            %explained by bhe variables by showing that correlations exist between
            %residuals of DA vs beh VS Beta vs beh
            testcor=1;
            mregress=1;
            condtypes=repmat({'all'},1,10);

    end
    argnum=argnum+1;
end
if isempty(sitelists)
sitelists=assignpairs(sessnums,xinfos,plotparam);
end
targdasites=plotparam.dasites;
dasites=getsites(sessnums,targdasites,'patra');
uniquedasites=unique({dasites(1:end).site});
cndasites=uniquedasites(contains(uniquedasites,'c'));
pdasites=uniquedasites(contains(uniquedasites,'p'));
sites=getlfpsites(sessnums,plotparam.lfpchs,'patra');
uniquesites=unique({sites(1:end).site});
cnsites=uniquesites(contains(uniquesites,'c'));
psites=uniquesites(contains(uniquesites,'p'));
if constrainpairs
     sites=sitelists;   
end
savepath=fullfile(plotparam.savepath,'cors_multsess',filesep);

if ~isdir(savepath)
mkdir(savepath);
end
savepath=[savepath 'cor_'];

if mregress
    savepath=fullfile(plotparam.savepath,'regress_multsess',filesep,'mreg_');
    if ~isdir(fullfile(plotparam.savepath,'regress_multsess',filesep))
        mkdir(fullfile(plotparam.savepath,'regress_multsess',filesep));
    end
end
if cormatrix
    savepath=fullfile(plotparam.savepath,'cormatrix_multsess',filesep,'cormat_');
    if ~isdir(fullfile(plotparam.savepath,'cormatrix_multsess',filesep))
        mkdir(fullfile(plotparam.savepath,'cormatrix_multsess',filesep));
    end
end
if testcor
       savepath=fullfile(plotparam.savepath,'mregress_multsess',filesep,'cormat_');
    if ~isdir(fullfile(plotparam.savepath,'mregress_multsess',filesep))
        mkdir(fullfile(plotparam.savepath,'mregress_multsess',filesep));
    end
end
%constrain lfp pair for each da site recorded loop through each session
%to get optimal lfp pair for each da site, already called assignpairs
%in 'plotmultiple.m' and suppleid here
pairedsites=unique({sitelists(1:end).site}); %constrainedpairs
tabdata={};
mdata={};
negdata={};
count=1;
if isempty(cordata)
for ix=1:length(xwin)
    disp(['xwin ' num2str(xwin{ix}(1)) '-' num2str(xwin{ix}(2)) 's']);
    tdata={};
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
        curdachsites{1}='c';
        curdachsites{2}='p';
        curdasites{1}='c';
        curdasites{2}='p';
    if ~lfponly
        curdasessids=find([dasites.sessnum]==sessnums(ise));
        curdachsites={dasites(curdasessids).probeid};
        curdasites={dasites(curdasessids).site};
        dasinxinfo=unique({xinfos(targxinfoids).siteda}); %da chs stored in xinfo
        daidsinxinfo=find(ismember(curdachsites,dasinxinfo));
        curdasites=curdasites(daidsinxinfo);
        curdachsites=intersect(curdachsites,dasinxinfo); %only da chs actually stored in xinfo
    end
    for id=1:length(curdachsites)
        curda=curdachsites{id};
        disp(['da ' curda]);
        curlfps={};
        pairedlfps=find(contains(curlfpchsites,curda(1)));  
        curlfps=curlfpchsites(pairedlfps);        %region specific lfps, ie if da ch in CN, lfp chs in CN   
        pairedlfpsites=curlfpsites(pairedlfps);
        
        %CONSTRAINED PAIRS FOR EACH DA
        if constrainpairs
            daid=find(strcmp({sites(cursessids).daid},curda));
            curlfps={};
            curlfps{1}=sites(cursessids(daid)).probeid;
            pairedlfpsites={};
            pairedlfpsites{1}=sites(cursessids(daid)).site;
        end
        if daonly
            %otherwise redundant DA sites in tabdata for all non-analyzed
            %lfp pairs
            curlfps=curlfps(1);
        end
        for il=1:length(curlfps) 
          %  disp(['lfp ' curlfps{il}]);
          curewin=[]; curlwin=[]; curhrwin=[];
          if contains(behtype,'eye')
            curewin=ewin;
          if iscell(ewin)
              %multiple wins (e.g. eye to follow lfp win)
              curewin=ewin{ix};
          end
          end
          if contains(behtype,'lick')
              curlwin=lwin;
          end
          if (contains(behtype,'hr') || contains(behtype,'rrstd'))
              curhrwin=hrwin;
          end
            tempdata={};
            if ~mregress
            if ~beh
                if ~flaglfptm && ~xvarflag 
                    tempdata=plotcortrials2(curses,trlists,datm,betatm,binfos,plotparam,'datm',curda,'lfpx',curlfps{il},xinfos,'win',xwin{ix},'offset',offset,'behx',xbinfos,'types',condtypes,'condition',condition,'event',targevent,'metric',metricda,'noplot');
                    %tempdata=plotcortrials2(curses,trlists,datm,betatm,binfos,plotparam,'datm',curda,'lfpx',curlfps{il},xinfos,'win',xwin{ix},'offset',offset,'behx',xbinfos,'types',condtypes,'trialseq','event',targevent,'metric',metricda,'noplot','nodata');
                elseif flaglfptm 
                   tempdata=plotcortrials2(curses,trlists,datm,betatm,binfos,plotparam,'datm',curda,'lfptm',curlfps{il},'xinfo',xinfos,'behx',xbinfos,'types',condtypes,'trialseq','event',targevent,'condition',condition,'metric',metricda,'metriclfp',metriclfp,'noplot','nodata');
                elseif ~flaglfptm && xvarflag 
                    tempdata=plotcortrials2(curses,trlists,datm,betatm,binfos,plotparam,'datm',curda,'lfpx',curlfps{il},xinfos,'metriclfp','lfpmints','xvar','behx',xbinfos,'types',condtypes,'condition',condition,'event',targevent,'metric',metricda);
                end
            end
            if beh
                if daonly && contains(behtype,'eye')
                tempdata=plotcortrials2(curses,trlists,datm,betatm,binfos,plotparam,'datm',curda,behtype,'behx',xbinfos,'types',condtypes,'condition',condition,'event',targevent,'metric',metricda,'ewin',curewin);                
                elseif lfponly && contains(behtype,'eye')
                tempdata=plotcortrials2(curses,trlists,datm,betatm,binfos,plotparam,'lfpx',curlfps{il},xinfos,behtype,'win',xwin{ix},'offset',offset,'behx',xbinfos,'types',condtypes,'condition',condition,'ewin',curewin);                
                end
                if daonly && contains(behtype,'lick')
                tempdata=plotcortrials2(curses,trlists,datm,betatm,binfos,plotparam,'datm',curda,behtype,'behx',xbinfos,'types',condtypes,'condition',condition,'event',targevent,'metric',metricda,'lwin',lwin);                
                elseif lfponly && contains(behtype,'lick')
                tempdata=plotcortrials2(curses,trlists,datm,betatm,binfos,plotparam,'lfpx',curlfps{il},xinfos,behtype,'win',xwin{ix},'offset',offset,'behx',xbinfos,'types',condtypes,'condition',condition,'lwin',lwin);                
                end
                if daonly && (contains(behtype,'hr') || contains(behtype,'rrstd'))
                tempdata=plotcortrials2(curses,trlists,datm,betatm,binfos,plotparam,'datm',curda,behtype,'behx',xbinfos,'types',condtypes,'condition',condition,'event',targevent,'metric',metricda,'hrwin',hrwin);                
                elseif lfponly && (contains(behtype,'hr') || contains(behtype,'rrstd'))
                tempdata=plotcortrials2(curses,trlists,datm,betatm,binfos,plotparam,'lfpx',curlfps{il},xinfos,behtype,'win',xwin{ix},'offset',offset,'behx',xbinfos,'types',condtypes,'condition',condition,'hrwin',hrwin);                
                end
                if daonly && contains(behtype,'rt') 
                tempdata=plotcortrials2(curses,trlists,datm,betatm,binfos,plotparam,'datm',curda,behtype,'types',condtypes,'condition',condition,'event',targevent,'metric',metricda);                
                elseif lfponly && contains(behtype,'rt') 
                tempdata=plotcortrials2(curses,trlists,datm,betatm,binfos,plotparam,'lfpx',curlfps{il},xinfos,behtype,'win',xwin{ix},'offset',offset,'types',condtypes,'condition',condition);                
                end
            end
            else
                %multiple regression analysis rather than single pair
                %correlations
               % tempdata=multiregress(curses,trlists,datm,betatm,binfos,plotparam,'datm',curda,'lfpx',curlfps{il},xinfos,'lick','rrstd','rt','behx',xbinfos,'event',targevent,'metric',metricda,'regressvars',{'type','side','trialnum'});                     
                [tempdata, mdl]=multiregress(curses,trlists,datm,betatm,binfos,plotparam,'datm',curda,'lfpx',curlfps{il},xinfos,'win',xwin{ix},'lick','rrstd','rt','behx',xbinfos,'event',targevent,'metric',metricda,'regressvars',{'type','side','trialnum'},'negcortarg');                     
               if ~isempty(mdl)
                   mdata{countm}=mdl;
                   countm=countm+1;
                        tempdata=setfield(tempdata,{1},['siteda'],curdasites{id});
                        tempdata=setfield(tempdata,{1},['sitelfp'],pairedlfpsites{il});
                   if isempty(negdata)
                       negdata=tempdata;
                   else
                    negdata=[negdata tempdata];
                   end
               end
            end
            if ~isempty(tempdata)
            count=count+1;
            tempdata=setfield(tempdata,{1},['siteda'],curdasites{id});
            tempdata=setfield(tempdata,{1},['sitelfp'],pairedlfpsites{il});
            if ~isempty(tdata) && ~isempty(tempdata)
                tdata=[tdata tempdata];
            elseif ~isempty(tempdata)
                tdata=tempdata;
            end
            end
        end
    end
    close all;
end
tabdata{ix}=tdata;
end
else
    tabdata=cordata;
end

if plothist
    %make histograms of r values p values grouped by pos/neg r's
figposxmat=[50,50,1900,900];
figsess=figure('visible','off');     %figure for each channel
if ispc
figsess=figure('visible','on');     %figure for each channel
end
set(0,'CurrentFigure',figsess);    %set figure handle to current figure
axa={};
axsiz=[350 300];
%set up plot
set(figsess,'position',figposxmat,'color',[1 1 1]);
        offx=100;
        axpos={};
        offy=100;
        count=1;
    for ip=1:2
        for ix=1:4
            axa{ip}{ix}=subplot(2,4,count);   hold(axa{ip}{ix},'on');  
            set(axa{ip}{ix},'units','pixels');
           % axpos{ip}{ix}=get(axa{ip}{ix},'position');
           % set(axa{ip}{ix},'xlim',[1 length(xids)]);
          %  set(axa{ip}{ix},'ylim',[1 length(yids)]);
          
            set(axa{ip}{ix},'xlim',[-.025 1.025]);
            set(axa{ip}{ix},'position',[(ix-1)*axsiz(1)+offx*ix,figposxmat(4)-axsiz(2)*ip-offy*ip,axsiz(1),axsiz(2)]);
           % xlabel(axa{ip}{ix},'start time (s)');
          %  ylabel(axa{ip}{ix},'end time (s)');
            count=count+1;
        end
    end
         set(axa{1}{1},'units','pixels')
    axpos=get(axa{1}{1},'position');

    xwinstab={};
    tid=[];
   targr='cp';
   nids={};
   binedges=0:.05:1;
   bine=1:50;
   logbins=logspace(-10,0,20);
   logbins=sort([logbins .05]);
   
            ticks=0:.2:1;
            ticksnew=sort([ticks .05]);
                      nids={};
          pids={};
          prs=[];
          pps=[];
          nrs=[];
          nps=[];
    for ix=1:length(tabdata)
        sitesigs={};
        xwinstab{ix}=tabdata{ix}.metriclfp;
         for ip=1:2
           % if isequal(round(xwinstab{ix}(1)*10),round(xwin{ip}(1)*10)) && isequal(round(xwinstab{ix}(2)*10),round(xwin{ip}(2)*10))
                %cn targwindow
            %    tid(ip)=ix;
          %  end
          pids{ip}=[];
          nids{ip}=[];
          if isempty(getcondition)
              %all trials
            nids{ip}=find([tabdata{ix}.r]<0 & contains({tabdata{ix}.siteda},targr(ip)));
            pids{ip}=find([tabdata{ix}.r]>0 & contains({tabdata{ix}.siteda},targr(ip)));
            [prs,bb]=histc([tabdata{ix}(pids{ip}).r],binedges);
            [pps,bb]=histc([tabdata{ix}(pids{ip}).p],binedges);
            [nrs,bb]=histc(-[tabdata{ix}(nids{ip}).r],binedges);
            [nps,bb]=histc([tabdata{ix}(nids{ip}).p],binedges);
          else
              %cond specific
              for iix=1:length(tabdata{ix})
                  %all sessions
                  if contains({tabdata{ix}(iix).siteda},targr(ip))
                      %for curren region
                                        %all sessions
                  targids=[];
                  condition=[];
                      for ic=1:length(getcondition)
                          targcond=getcondition{ic};
                          if contains(getcondition{ic},'right') || contains(getcondition{ic},'left') 
                              targcond='side';
                          end                          
                        condids=getfield(tabdata{ix},{iix},targcond);
                        if contains(getcondition{ic},'left')
                            condids=~condids;
                        end
                        if ic==1
                            targids=condids;
                            condition=getcondition{ic};
                        else
                            targids=[targids & condids];
                            condition=[condition '_' getcondition{ic}];
                        end
                      end
                   
                    xval=[tabdata{ix}(iix).datax(targids)];
                    yval=[tabdata{ix}(iix).datay(targids)];
                    nonnan=find(~isnan(xval) & ~isnan(yval) & ~isoutlier(xval) & ~isoutlier(yval)); %based on median outlier removal rather than mean 04/25/2020
                     datax=xval(nonnan);
                     datay=yval(nonnan);
                    [r,psig]=corr(datax',datay');
                    tabdata{ix}(iix).r=r;          %update values in tabdata according to new condition
                    tabdata{ix}(iix).p=psig;
                    tabdata{ix}(iix).type=getcondition;
                    if r>0 
                        pids{ip}=[pids{ip} iix];
                    elseif r<0 
                        nids{ip}=[nids{ip} iix];
                    end
                  end
              end

             [prs,bb]=histc([tabdata{ix}(pids{ip}).r],binedges);
            [pps,bb]=histc([tabdata{ix}(pids{ip}).p],logbins);
            [nrs,bb]=histc(-[tabdata{ix}(nids{ip}).r],binedges);
            [nps,bb]=histc([tabdata{ix}(nids{ip}).p],logbins);
            sigidspos=find([tabdata{ix}(pids{ip}).p]<=0.05 & [tabdata{ix}(pids{ip}).r]>0);
            sitesigs(ip).sitespos={tabdata{ix}(pids{ip}(sigidspos)).sitelfp};
                sitesigs(ip).sitesposDA={tabdata{ix}(pids{ip}(sigidspos)).siteda};
        sitesigs(ip).sessposDA=[tabdata{ix}(pids{ip}(sigidspos)).sess];
               sitesigs(ip).p=[tabdata{ix}(pids{ip}(sigidspos)).p];
               sitesigs(ip).r=[tabdata{ix}(pids{ip}(sigidspos)).r];
 negidspos=find([tabdata{ix}(nids{ip}).p]<=0.05 & [tabdata{ix}(nids{ip}).r]<0);
            sitesigs(ip).sitesneg={tabdata{ix}(nids{ip}(negidspos)).sitelfp};
                sitesigs(ip).sitesnegDA={tabdata{ix}(nids{ip}(negidspos)).siteda};
        sitesigs(ip).sessnegDA=[tabdata{ix}(nids{ip}(negidspos)).sess];
               sitesigs(ip).pn=[tabdata{ix}(nids{ip}(negidspos)).p];
               sitesigs(ip).rn=[tabdata{ix}(nids{ip}(negidspos)).r];
          end
            cla(axa{ip}{1});
            cla(axa{ip}{2});
            cla(axa{ip}{3});
            cla(axa{ip}{4});
            
              titletext=[targevent ' | ' condition ' | ' condtypes{1}];  
            plotsavename=['hist_' targevent '_' condition '_' condtypes{1}];        
            %text(axa{1}{1},0,axpos(4)+75,titletext,'units','pixels','fontweight','bold','fontsize',fontsize,'interpreter','none');
            
                axpos=get(axa{ip}{1},'position');
                wintitle=[titletext ' | ' num2str(xwinstab{ix}(1)) '-' num2str(xwinstab{ix}(2)) ' s' ];
                region='caudate nucleus';
                if ip==2
                    region='putamen';
                    text(axa{1}{1},0,axpos(4)+75,wintitle,'units','pixels','fontweight','bold','fontsize',fontsize,'interpreter','none');
                end
            text(axa{ip}{1},0,axpos(4)+35,region,'units','pixels','fontweight','bold','fontsize',12,'interpreter','none');
            histplot= bar(axa{ip}{1},binedges+(binedges(2)-binedges(1))/2,prs,'FaceColor',[0 0 0],'FaceAlpha',.2); %,'BarWidth',1,'FaceAlpha',.7,'linestyle','none');
            title(axa{ip}{1},['positive correlations']);
            xlabel(axa{ip}{1},'R');
            ylabel(axa{ip}{1},'counts');
           % histplot2= bar(axa{ip}{2},binedges+(binedges(2)-binedges(1))/2,pps,'FaceColor',[0 0 0],'FaceAlpha',.2);
            histplot2=semilogx(axa{ip}{2},log10(logbins),pps,'linestyle','none','marker','o','markersize',6,'markeredgecolor',[0 0 0],'markerfacecolor',[.5 .5 .5]);
            hold(axa{ip}{2},'on'); plot(axa{ip}{2},[log10(.05) log10(.05)],[get(axa{ip}{2},'ylim')],'color',[0 0 0],'linestyle','--');
            set(axa{ip}{2}, 'XTickMode', 'auto', 'XTickLabelMode', 'auto','xlimmode','auto')
            title(axa{ip}{2},['positive correlations']);
            xlabel(axa{ip}{2},'P');
            ylabel(axa{ip}{2},'counts');
                      %  set(axa{ip}{2},'xtick',ticksnew)
            histplot= bar(axa{ip}{3},binedges+(binedges(2)-binedges(1))/2,nrs,'FaceColor',[0 0 0],'FaceAlpha',.2);
            title(axa{ip}{3},['negative correlations']);
            xlabel(axa{ip}{3},'R');
            ylabel(axa{ip}{3},'counts');
           % histplot2= bar(axa{ip}{4},binedges+(binedges(2)-binedges(1))/2,nps,'FaceColor',[0 0 0],'FaceAlpha',.2);
            histplot2=semilogx(axa{ip}{4},log10(logbins),nps,'linestyle','none','marker','o','markersize',6,'markeredgecolor',[0 0 0],'markerfacecolor',[.5 .50 .50]);
            hold(axa{ip}{4},'on'); plot(axa{ip}{4},[log10(.05) log10(.05)],[get(axa{ip}{4},'ylim')],'color',[0 0 0],'linestyle','--');
            set(axa{ip}{4}, 'XTickMode', 'auto', 'XTickLabelMode', 'auto','xlimmode','auto')
            title(axa{ip}{4},['negative correlations']);
            xlabel(axa{ip}{4},'P');
            ylabel(axa{ip}{4},'counts');   
                                %    set(axa{ip}{2},'xtick',ticksnew)
         end 
         plotsaveout=[savepath plotsavename '_' num2str(xwinstab{ix}(1)*10) '-' num2str(xwinstab{ix}(2)*10) 'ds'];
         save(plotsaveout,'tabdata','sitesigs')
         savefig(figsess,plotsaveout);
saveas(figsess,plotsaveout,'jpg')
print(figsess,plotsaveout,'-painters','-depsc');
%save(plotsaveout,'tabdata','sitesigs')
    end  
    
end

if testcor
%already supplied tabdata & testing multi for sig neg's
%find sit sites at targeted windows defined for cn and put xwins{ip}
    xwinstab={};
    tid=[];
    for ix=1:length(tabdata)
        xwinstab{ix}=tabdata{ix}.metriclfp;
        for ip=1:2
        if isequal(round(xwinstab{ix}(1)*10),round(xwin{ip}(1)*10)) && isequal(round(xwinstab{ix}(2)*10),round(xwin{ip}(2)*10))
            %cn targwindow
            tid(ip)=ix;
        end
        end
    end
    %get sigdata for cn and p
    %redo multiregress and test cor of beh residuals on negatively
    %correlated sites from tabdata alaready generated/supplied
    %TRY subset_regression_helen1 (Supplied by Ken Amemori 05/01/20)
    sids={};
    targr='cp';
    mdata=[];
    count=1;
    for ip=1:2
        sids{ip}=find([tabdata{tid(ip)}.r]<0 & [tabdata{tid(ip)}.p]<=0.05 & contains({tabdata{tid(ip)}.siteda},targr(ip)));
        for id=1:length(sids{ip})
            curlfpch=tabdata{tid(ip)}(sids{ip}(id)).sitechlfp;
            curses=tabdata{tid(ip)}(sids{ip}(id)).sess;
            curda=tabdata{tid(ip)}(sids{ip}(id)).sitechda;
            targevent=tabdata{tid(ip)}(sids{ip}(id)).event;
            [tempdata, mdl]=multiregress(curses,trlists,datm,betatm,binfos,plotparam,'datm',curda,'lfpx',curlfpch,xinfos,'win',xwin{ip},'lick','rrstd','rt','behx',xbinfos,'event',targevent,'metric',metricda,'regressvars',{'type','side','trialnum'},'multicor');    
            if ~isempty(tempdata)
                count=count+1;
                tempdata=setfield(tempdata,{1},['siteda'],tabdata{tid(ip)}(sids{ip}(id)).siteda);
                tempdata=setfield(tempdata,{1},['sitelfp'],tabdata{tid(ip)}(sids{ip}(id)).sitelfp);
                if ~isempty(mdata) && ~isempty(tempdata)
                    mdata=[mdata tempdata];
                elseif ~isempty(tempdata)
                    mdata=tempdata;
                end
            end

        end            
    end    

else



%%%%%%
if mregress
    sitesigs=negdata;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PLOTTING bar graphs of correlations CN vs PUT
if ~mregress
%plotting setup for bar graphs of da vs beta or da/beta vs beh correlations as a function of time window
figsess=figure('visible','off');     %figure for each channel
if ispc
figsess=figure('visible','on');     %figure for each channel
end
set(figsess,'position',figpos,'color',[1 1 1]);
set(0,'CurrentFigure',figsess);    %set figure handle to current figure
clf(figsess,'reset');
set(figsess,'color',[1 1 1]);
axpos={};
%cn vs put
for ip=1:2
    axa{ip}=subplot(1,2,ip);   hold(axa{ip},'on');
    set(axa{ip},'units','pixels');
    axpos{ip}=get(axa{ip},'position');
   % set(axa{ip},'xtick',1:length(condlabels),'xticklabel',condlabels,'ticklabelinterpreter','none','xticklabelrotation',45);
    xlab1{1}=metriclfp;
    set(axa{ip},'xtick',1,'xticklabel',xlab1,'ticklabelinterpreter','none');
    if ~flaglfptm && ~xvarflag
    for ix=1:length(xwin)
        xlab1{ix}=[num2str(xwin{ix}(1)+offset(1)) '-' num2str(xwin{ix}(2)+offset(2)) 's']; 
    end
    set(axa{ip},'xtick',1:length(xwin),'xticklabel',xlab1,'ticklabelinterpreter','none');
    else
                set(axa{ip},'xtick',1,'xticklabel',metriclfp,'ticklabelinterpreter','none');
    end
    set(axa{ip},'xlim',[0 length(xwin)+1]);
    if length(xwin)>5
        set(axa{ip},'xTickLabelRotation',90)
    end
    set(axa{ip},'position',[(ip-1)*axsiz(1)+ip*(mar/2)+off,125,axsiz(1),axsiz(2)])
      %  ax=text(axa{ip},10,axsiz(2),plotvar{ip},'units','pixels','fontweight','bold');
end

yaxlabel=['R (p<0.05), DA ' metricda ' vs Beta '  ];   %yaxis label
if beh && daonly
    yaxlabel=['R (p<0.05), DA ' metricda ' vs ' behtype  ];   %yaxis label
end
if beh && lfponly
    yaxlabel=['R (p<0.05), beta ' metriclfp ' vs ' behtype  ];   %yaxis label
    if relfp
        yaxlabel=['R (p<0.05), beta win vs ' behtype  ];   %yaxis label
    end
end
ylabel(axa{1},yaxlabel,'interpreter','none');
set(axa{1},'units','pixels')
axpos=get(axa{1},'position');
titletext=[targevent ' | ' condtypes{:} ' | ' condition ];  
text(axa{1},-100,axpos(4)+160,titletext,'units','pixels','fontweight','bold','fontsize',fontsize,'interpreter','none');

%scatter Rs
rsig={};
siteids{1}=find(contains({tabdata{1}.siteda},'c'));
siteids{2}=find(contains({tabdata{1}.siteda},'p'));
sitesigs={};
for ix=1:length(xwin)        
        sigids{1}=find([tabdata{ix}.p]<0.05 & contains({tabdata{ix}.siteda},'c'));
        sigids{2}=find([tabdata{ix}.p]<0.05 & contains({tabdata{ix}.siteda},'p'));
        rsig{1}{ix}=[tabdata{ix}(sigids{1}).r];
        rsig{2}{ix}=[tabdata{ix}(sigids{2}).r];
        for ip=1:2 
        randomnumbers=rand(length(sigids{ip}),1);
        ranj=randomnumbers*.5-.3;
         scatter(axa{ip},ix+ranj,...
                    rsig{ip}{ix},defaultmarksize,'o','markeredgecolor',[0 0 0],...
                        'MarkerEdgeAlpha',alphac,'linewidth',lwid,'markerfacecolor',[0 0 0],'markerfacealpha',.25);
    end            
end

%bar averages R
for ix=1:length(xwin)
    for ip=1:2
         rdatas=[tabdata{ix}.r];
      pdatas=[tabdata{ix}.p];
      sitedata={tabdata{ix}.siteda};
      siteids{1}=(contains(sitedata,'c'));
      siteids{2}=(contains(sitedata,'p'));
       rdata{ip}=rdatas(siteids{ip});
        pdata{ip}=pdatas(siteids{ip});
                sigids{ip}=(pdata{ip}<=0.05);
       % sigids{ip}=find([tabdata{ix}.p]<=0.05 & siteids{ip});
        %plot bar graph
        ylimax=get(axa{ip},'ylim');
       % ylimax=[-max(abs(ylimax)) max(abs(ylimax))];
        ax2 = axes('Position',get(axa{1},'Position'),'XAxisLocation','top',...
        'YAxisLocation','right','Color','none','XColor','k','YColor','k');
        hold(ax2, 'all');  %   <--------------------------------
        %barsimp=bar(axa{1},(1:length(meansda)),meansda,.5,'FaceAlpha',0.05,'FaceColor',[0 0 0],'EdgeColor',[0 0 0]);
        %        set(axa{1},'ylim',ylimax)
        meansdapos(ip,ix)=nanmean(rdata{ip}(sigids{ip} & rdata{ip}>0));       %groups with positive values, nan if no positive R's
                meansdaneg(ip,ix)=nanmean(rdata{ip}(sigids{ip} & rdata{ip}<0));     %groups with neg values, nan if no negative R's
                stdsdapos(ip,ix)=nanstd(rdata{ip}(sigids{ip} & rdata{ip}>0)); 
                stdsdaneg(ip,ix)=nanstd(rdata{ip}(sigids{ip} & rdata{ip}<0)); 
                cisdapos(ip,ix)=stdsdapos(ip,ix)./sqrt(length(rdata{ip}(sigids{ip} & rdata{ip}>0))).*1.96;
                cisdaneg(ip,ix)=stdsdaneg(ip,ix)./sqrt(length(rdata{ip}(sigids{ip} & rdata{ip}<0))).*1.96;               
                totalsigpos=length(rdata{ip}(sigids{ip} & rdata{ip}>0));
                totalsigneg=length(rdata{ip}(sigids{ip} & rdata{ip}<0));    
                fractpossig=totalsigpos./length(find(siteids{ip}));
                fractnegsig=totalsigneg./length(find(siteids{ip}));
                
       % meanrpos=nanmean(rsig{ip}{ix}(rsig{ip}{ix}>0));       
       % meanrneg=nanmean(rsig{ip}{ix}(rsig{ip}{ix}<0));       
        barsimp=bar(axa{ip},ix,meansdapos(ip,ix),.5,'FaceAlpha',0.05,'FaceColor',[0 0 0],'EdgeColor',[0 0 0]);
        barsimp=bar(axa{ip},ix,meansdaneg(ip,ix),.5,'FaceAlpha',0.05,'FaceColor',[0 0 0],'EdgeColor',[0 0 0]);
        ylimmax=max(abs(ylimax));
        set(axa{ip},'ylim',[-ylimmax ylimmax])
        
      %  totalsigpos=length(find([rsig{ip}{ix}]>0));
      %  totalsigneg=length(find([rsig{ip}{ix}]<0));
        if length(xwin)<5 || ix==length(xwin)
            text(axa{ip},ix,ylimmax,['+' num2str(totalsigpos) '/' num2str(length(find(siteids{ip}))) sprintf('\n') '-' num2str(totalsigneg) '/' num2str(length(find(siteids{ip})))]);   
        else
            text(axa{ip},ix,ylimmax,['+' num2str(totalsigpos) sprintf('\n') '-' num2str(totalsigneg)]);   
        end
        sigidspos=find(siteids{ip} & pdatas<=0.05 & rdatas>0);
        sigidsneg=find(siteids{ip} & pdatas<=0.05 & rdatas<0);
       % intersect(sigids{ip},find([tabdata{ix}.r]>0));
        %sigidsneg=intersect(sigids{ip},find([tabdata{ix}.r]<0));
        
        sitesigs(ix).xwin=xwin{ix};
        if ~isempty(sigidspos)
        sitesigs(ix).sitespos={tabdata{ix}(sigidspos).sitelfp};
                sitesigs(ix).sitesposDA={tabdata{ix}(sigidspos).siteda};
        sitesigs(ix).sessposDA=[tabdata{ix}(sigidspos).sess];
                if isfield(tabdata{ix}(sigidspos),'trials')
        sitesigs(ix).trialsposDA=[tabdata{ix}(sigidspos).trials];
                else
                    sitesigs(ix).trialsposDA=nan;
                end
        end
        if ~isempty(sigidsneg)
        sitesigs(ix).sitesneg={tabdata{ix}(sigidsneg).sitelfp};
        sitesigs(ix).sitesnegDA={tabdata{ix}(sigidsneg).siteda};
        sitesigs(ix).sessnegDA=[tabdata{ix}(sigidsneg).sess];    
        if isfield(tabdata{ix}(sigidsneg),'trials')
        sitesigs(ix).trialsnegDA=[tabdata{ix}(sigidsneg).trials];   
        else
            sitesigs(ix).trialsnegDA=nan;
        end
        end
    end  
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2D MATRIX PLOT of CORs function of xwins 
if cormatrix
        %2d matrix da vs beta as a function of beta avg time windows
          %make 2d matrix from start end time points
          figposxmat=[50,50,1900,900];
figsess=figure('visible','off');     %figure for each channel
if ispc
figsess=figure('visible','on');     %figure for each channel
end
set(0,'CurrentFigure',figsess);    %set figure handle to current figure
axa={};
axsiz=[350 300];

          %set up plot
              set(figsess,'position',figposxmat,'color',[1 1 1]);
             %make 2d matrix from start end time points
         maxts=max([xwin{:}]);
         mints=min([xwin{:}]);
         step=xwin{1}(2)-xwin{1}(1);    %assume step size equal for all cells
         xids=mints:step:maxts-step;
         yids=mints+step:step:maxts;
        count=1;
        offx=100;
        axpos={};
        offy=100;
        titletext=[targevent ];  
        
        for ip=1:2
            for ix=1:4
                axa{ip}{ix}=subplot(2,4,count);   hold(axa{ip}{ix},'on');  
                set(axa{ip}{ix},'units','pixels');
               % axpos{ip}{ix}=get(axa{ip}{ix},'position');
                set(axa{ip}{ix},'xtick',1:1:length(xids),'xticklabel',xids(1:1:end),'ticklabelinterpreter','none');
                set(axa{ip}{ix},'xlim',[1 length(xids)]);
                set(axa{ip}{ix},'xTickLabelRotation',90)
                set(axa{ip}{ix},'ytick',1:1:length(yids),'yticklabel',yids(1:1:end),'ticklabelinterpreter','none');
                set(axa{ip}{ix},'ylim',[1 length(yids)]);
                set(axa{ip}{ix},'position',[(ix-1)*axsiz(1)+offx*ix,figposxmat(4)-axsiz(2)*ip-offy*ip,axsiz(1),axsiz(2)]);
                xlabel(axa{ip}{ix},'start time (s)');
                ylabel(axa{ip}{ix},'end time (s)');
                count=count+1;
            end
        end
         set(axa{1}{1},'units','pixels')
    axpos=get(axa{1}{1},'position');
for ip=1:2
               xmat(ip).pcor=nan(length(xids),length(yids));        %correlation coefficients positive, for sig correlations
               xmat(ip).ncor=nan(length(xids),length(yids));
               xmat(ip).pratiosig=nan(length(xids),length(yids));   %fraction sites positive correlation
               xmat(ip).nratiosig=nan(length(xids),length(yids));
               xmat(ip).cnsites=nan(length(xids),length(yids));
               xmat(ip).putsites=nan(length(xids),length(yids));
end
          for ix=1:length(xwin)
             % rdatas=[tabdata{ix}.r];
            %  pdatas=[tabdata{ix}.p];
             % sitedata={tabdata{ix}.siteda};
             % siteids{1}=find(contains(sitedata,'c'));
            %  siteids{2}=find(contains(sitedata,'p'));
            
            if ~isempty(getcondition)
                %cond specific, recompute R's and P's for given data set,
                %update values in tabdata
                condition=getcondition;
              for iix=1:length(tabdata{ix})
                  %all sessions
                  targids=[];
                  for ic=1:length(getcondition)
                      targcond=getcondition{ic};
                      if contains(getcondition{ic},'right') || contains(getcondition{ic},'left') 
                          targcond='side';
                      end                          
                    condids=getfield(tabdata{ix},{iix},targcond);
                    if contains(getcondition{ic},'left')
                        condids=~condids;
                    end
                    if ic==1
                        targids=condids;
                        condition=getcondition{ic};
                    else
                        targids=[targids & condids];
                        condition=[condition '_' getcondition{ic}];
                    end
                  end
                    xval=[tabdata{ix}(iix).datax(targids)];
                    yval=[tabdata{ix}(iix).datay(targids)];
                    nonnan=find(~isnan(xval) & ~isnan(yval) & ~isoutlier(xval) & ~isoutlier(yval)); %based on median outlier removal rather than mean 04/25/2020
                     datax=xval(nonnan);
                     datay=yval(nonnan);
                    [r,psig]=corr(datax',datay');
                    tabdata{ix}(iix).r=r;          %update values in tabdata according to new condition
                    tabdata{ix}(iix).p=psig;
                    tabdata{ix}(iix).type=getcondition;
              end
            end
              
              rdatas=[tabdata{ix}.r];
              pdatas=[tabdata{ix}.p];
              sitedata={tabdata{ix}.siteda};
              siteids{1}=(contains(sitedata,'c'));
              siteids{2}=(contains(sitedata,'p'));
              
              
              for ip=1:2
               rdata{ip}=rdatas(siteids{ip});
              pdata{ip}=pdatas(siteids{ip});
              sigids{ip}=(pdata{ip}<=0.05);
                        meansdapos(ip,ix)=nanmean(rdata{ip}(sigids{ip} & rdata{ip}>0));       %groups with positive values, nan if no positive R's
                meansdaneg(ip,ix)=nanmean(rdata{ip}(sigids{ip} & rdata{ip}<0));     %groups with neg values, nan if no negative R's
                stdsdapos(ip,ix)=nanstd(rdata{ip}(sigids{ip} & rdata{ip}>0)); 
                stdsdaneg(ip,ix)=nanstd(rdata{ip}(sigids{ip} & rdata{ip}<0)); 
                cisdapos(ip,ix)=stdsdapos(ip,ix)./sqrt(length(rdata{ip}(sigids{ip} & rdata{ip}>0))).*1.96;
                cisdaneg(ip,ix)=stdsdaneg(ip,ix)./sqrt(length(rdata{ip}(sigids{ip} & rdata{ip}<0))).*1.96;               
                totalsigpos=length(rdata{ip}(sigids{ip} & rdata{ip}>0));
                totalsigneg=length(rdata{ip}(sigids{ip} & rdata{ip}<0));    
                fractpossig=totalsigpos./length(find(siteids{ip}));
                fractnegsig=totalsigneg./length(find(siteids{ip}));
                
                xmat(ip).pcor(find(ismember(round(xids.*10),round(xwin{ix}(1).*10))),find(ismember(round(yids.*10),round(xwin{ix}(2).*10))))=meansdapos(ip,ix);
                xmat(ip).ncor(find(ismember(round(xids.*10),round(xwin{ix}(1).*10))),find(ismember(round(yids.*10),round(xwin{ix}(2).*10))))=meansdaneg(ip,ix);
                xmat(ip).pratiosig(find(ismember(round(xids.*10),round(xwin{ix}(1).*10))),find(ismember(round(yids.*10),round(xwin{ix}(2).*10))))=totalsigpos;
                xmat(ip).nratiosig(find(ismember(round(xids.*10),round(xwin{ix}(1).*10))),find(ismember(round(yids.*10),round(xwin{ix}(2).*10))))=totalsigneg;   
                xmat(ip).cnsites(find(ismember(round(xids.*10),round(xwin{ix}(1).*10))),find(ismember(round(yids.*10),round(xwin{ix}(2).*10))))=length(find(siteids{1}));
                xmat(ip).putsites(find(ismember(round(xids.*10),round(xwin{ix}(1).*10))),find(ismember(round(yids.*10),round(xwin{ix}(2).*10))))=length(find(siteids{2}));
              end
          end
          for ip=1:2
          plotnames=fieldnames(xmat(ip));
          for iplot=1:4
              plotdata=getfield(xmat(ip),plotnames{iplot});
              plotdata=plotdata';
              cla(axa{ip}{iplot});
               if ip==1
                textCN=text(axa{1}{iplot},10,axsiz(2)+25,'CN','units','pixels','fontweight','bold');
                                                textCNsites=text(axa{1}{iplot},10,axsiz(2)+40,[num2str(xmat(ip).cnsites(1)) ' sites'],'units','pixels','fontweight','bold');
                else
                textPut=text(axa{2}{iplot},10,axsiz(2)+25,'Putamen','units','pixels','fontweight','bold');
                                textPutsites=text(axa{2}{iplot},10,axsiz(2)+40,[num2str(xmat(ip).putsites(1)) ' sites'],'units','pixels','fontweight','bold');
                end
              matplot{iplot}=image(abs(plotdata),'parent',axa{ip}{iplot},'cdatamapping','scaled');
                  hold(axa{ip}{iplot},'on')
                title(axa{ip}{iplot},plotnames{iplot});
                if contains(plotnames{iplot},'ncor')
                    title(axa{ip}{iplot},'|ncor|');
                end
                 %   colormap(hax,colortrial)
                %set artifact (on nan values) mask ontop of image using alpha data properties
                artTime=isnan(plotdata);   %find artifact points (nan periods)
                artTime=abs(artTime-1);         %make alpha data mask by inverting 1/0's
                artTime2=artTime;
                maskGray=artTime2==0;             %find Zero indices representing artifact mask
                maskGray=maskGray*.15;            %make gray rather than white default by making non-zero
                artTime=artTime+maskGray;
                set(matplot{iplot}, 'AlphaData', artTime);
                origpos=getpixelposition(axa{ip}{iplot});  
                %if iplot==length(plotnames)
                    h1=colorbar(axa{ip}{iplot},'northoutside');
                    cpos = getpixelposition(h1);
                   % ylabelbar=ylabel(h1,clabel); 
                   % ypos = getpixelposition(ylabelbar);
                    %cpos(4) = 15; 
                    set(h1,'Units','Pixels','Position', [cpos(1)+origpos(3)*.75 origpos(2)+origpos(4)+10 cpos(3)*.25 cpos(4)*.5],'fontsize',10);
               % end
          end
          end
          if ~isempty(getcondition)
            titletext=[targevent ' | ' condition];
          end
        
text(axa{1}{1},0,axpos(4)+75,titletext,'units','pixels','fontweight','bold','fontsize',fontsize,'interpreter','none');

end

end
if ~plothist
savename=[savepath targevent '_' condition '_' condtypes{1} ];
if sum(abs(xwin{1}))>0
     if length(xwin)<=5
for ix=1:length(xwin)
    savename=[savename '_' num2str((xwin{ix}(1)+offset(1))*rate) '-' num2str((xwin{ix}(2)+offset(2))*rate) ];
end
     else
for ix=1:5
    savename=[savename '_' num2str((xwin{ix}(1)+offset(1))*rate) '-' num2str((xwin{ix}(2)+offset(2))*rate) ];
end
     end
end
if xvarflag
    savename=[savename '_' metriclfp ];
end

if beh && daonly
        savename=[savename '_' metricda '_' behtype];
        if strcmp(behtype,'eye')
            savename=[savename '_ewin_' num2str(ewin(1)*rate) '-' num2str(ewin(2)*rate)];
        end
                if strcmp(behtype,'lick')
            savename=[savename '_lwin_' num2str(lwin(1)*rate) '-' num2str(lwin(2)*rate)];
                end
         if (contains(behtype,'hr') || contains(behtype,'rrstd'))
            savename=[savename '_hrwin_' num2str(hrwin(1)*rate) '-' num2str(hrwin(2)*rate)];
        end     
end
if beh && lfponly
    if ~relfp
        savename=[savename '_' metriclfp '_' behtype];
    else
        savename=[savename '_winlfp_' behtype];
    end
                    if strcmp(behtype,'lick')
            savename=[savename '_lwin_' num2str(lwin(1)*rate) '-' num2str(lwin(2)*rate)];
                    end
        if (contains(behtype,'hr') || contains(behtype,'rrstd'))
            savename=[savename '_hrwin_' num2str(hrwin(1)*rate) '-' num2str(hrwin(2)*rate)];
        end            
end
if relfp
    if length(xwin)<=5
for ix=1:length(xwin)
    savename=[savename '_' num2str((xwin{ix}(1)+offset(1))*rate) '-' num2str((xwin{ix}(2)+offset(2))*rate) ];
end
    end
end
if constrainpairs
    savename=[savename '_pairedsites'];
end

if ~mregress
savefig(figsess,savename);
saveas(figsess,savename,'jpg')
print(figsess,savename,'-painters','-depsc');
save(savename,'tabdata','sitesigs')
if cormatrix
    save(savename,'tabdata','sitesigs','xmat')
end
else
    save(savename,'tabdata','negdata','mdata')
end
end

