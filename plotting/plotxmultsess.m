function plotdata=plotxmultsess(xinfos,plotparam,varargin)
%2/17/2019 pre, now updated xinfos so not arrayed by sess num, all in one
%chart
%plot bar scatter plots of properties of cross-covariance lags/waveforms
%12/31/2018 updated for timing characteristics from new xclust 
%same as plotx, but plot big/small/targ on same plots
%for plotting multiple sessions
%1/3/2018, udpates as in plotx/plotxallsess
argnum=1;
datype='dapos';
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'daneg'
            datype='daneg';
        case 'dareb'
            datype='dareb';
        case 'daall'
            datype='daall';
    end
    argnum=argnum+1;
end
fontsize=9;
savepath=plotparam.savepath;
sessnums=plotparam.sessnums;
trialgrps=plotparam.trialgrps;
trialinfo=trialgrps(1).trialinfo;
lfpchs=plotparam.lfpchs;
%get lfp ch groups p & c
pgroup=find(contains(lfpchs,'p')==1);
cgroup=find(contains(lfpchs,'c')==1);
targdasites=plotparam.dasites;
figpos=[50,50,1800,900];
figsess=figure('position',figpos,'color',[1 1 1]);
set(0,'CurrentFigure',figsess);    %set figure handle to current figure
axa={};

xinfo=xinfos;
eventrows={xinfo(1:end).event};
eventtypes=unique(eventrows);

sessionrows={xinfo(1:end).sessiontype};
sessiontypes=unique({xinfo(1:end).sessiontype});
params={'damax','damaxts','damin','damints','dafallts','darisets',...
    'lfpmax','lfpmaxts','lfpmin','lfpmints','lfprisets','lfpfallts',...
    'lfppostmax','lfppostmaxts'};
params_mag={'damax','damin',...
    'lfpmax','lfpmin','lfppostmax','zlagcoef'};
params_xcov={'maxprelagts',...
    'minprelagts','maxpostlagts','minpostlagts','maxprecoef',...
    'minprecoef','maxpostcoef','minpostcoef'};
params_timing_damax={'delt_lfpmax_damax','delt_lfpfall_damax'...
    'delt_lfpmin_damax','delt_lfprise_damax',...
    'delt_lfppostmax_damax','lfppostmaxts','damaxts','lfpmints'}; 
params_timing_darise={'delt_lfpmax_darise','delt_lfpfall_darise'...
    'delt_lfpmin_darise','delt_lfprise_darise',...
    'delt_lfppostmax_darise','lfppostmaxts','darisets','lfprisets'};   
params_timing_dafall={'delt_lfpmax_dafall','delt_lfpfall_dafall'...
    'delt_lfpmin_dafall','delt_lfprise_dafall',...
    'delt_lfppostmax_dafall','lfppostmaxts','dafallts','lfpfallts'};
 params_timing_damin={'delt_lfpmax_damin','delt_lfpfall_damin'...
    'delt_lfpmin_damin','delt_lfprise_damin',...
    'delt_lfppostmax_damin','lfppostmaxts','damints','lfpmints'};

params=[params_mag params_xcov params_timing_damax...
   params_timing_darise params_timing_dafall params_timing_damin];   

ygrp{1}=params_mag;
ygrp{2}=params_xcov;
ygrp{3}=params_timing_damax;

ynam{1}='mag';
ynam{2}='xcov';
ynam{3}='timing_damax';



axpos={};
axsiz=[750,320];
axoff=100;
mark=[0 .1 .9; 0 .5 .3; 0 .9 .1];
mark=[.7 0 0; 0.6 0 .7; 0 .7 .1];
mark=[0 0.1 .7; 0.8 0 .5; 0 .7 0;.7 0 0; 0.4 0 .7];

markr=[1 .5 .0; .3 .3 0; .9 .1 0];
typeoff=[-0.3 -0.15 0];
sessmark={'.','+','o','^','s','*','v','p','h','x'};
sessmarktxt={'.','+','o','\delta','sq','*','\nabla','p','h','x'};
for ida=1:length(targdasites)
    %each da separate subplot  
    lfpsites={};
    daregion=contains(targdasites(ida),'c');     %1=='c'
    if daregion==1
        lfpsites=lfpchs(cgroup);
    else
        lfpsites=lfpchs(pgroup);
    end
for ilfp=1:length(lfpsites)    
for ievent=1:length(eventtypes)
for yid=1:length(ygrp)    
    %timing/magnitudes/xcov parameters to plot in diff figs    
    clf(figsess,'reset');
    set(figsess,'color',[1 1 1]);
    trialtypes=trialinfo(1); %get trial groups for session
    axpos={};
    for ip=1:4
        axa{ip}=subplot(2,2,ip);   hold(axa{ip},'on');
        set(axa{ip},'units','pixels');
        axpos{ip}=get(axa{ip},'position');
        if mod(ip,2)
            %left plot            
            set(axa{ip},'position',[axoff axpos{ip}(2) axsiz(1) axsiz(2)]);
        else
            set(axa{ip},'position',[axsiz(1)+axoff*2 axpos{ip}(2) axsiz(1) axsiz(2)]);
        end            
        set(axa{ip},'xtick',1:length(trialtypes.names),'xticklabel',trialtypes.names);
        set(axa{ip},'xlim',[0 length(trialtypes.names)+1]);
        set(axa{ip},'xTickLabelRotation',90)
    end
    
    labels={};
    titletext=[eventtypes{ievent} ' | ' targdasites{ida} ' x ' lfpsites{ilfp}];  
    set(axa{1},'units','pixels')
    axpos=get(axa{1},'position');
    text(axa{1},-30,axpos(4)+30,titletext,'units','pixels','fontweight','bold');
    %set legend
    labelx={};
    for xleg=1:length(sessnums)
        labelx{xleg}=[num2str(sessnums(xleg))];  
        axpos=get(axa{1},'position');
        text(axa{1},100+xleg*50,axpos(4)+30,labelx{xleg},'units','pixels','fontweight','bold','color',mark(xleg,:));
        text(axa{1},100+xleg*50+18,axpos(4)+30,sessmarktxt{xleg},'interpreter','tex','units','pixels','fontweight','bold','color',[0 0 0]);
       % axpos=get(axc{1},'position');
      %  text(axc{1},100+xleg*50,axpos(4)+30,label{xleg},'units','pixels','fontweight','bold','color',mark(xleg,:));

    end
    curgrp=ygrp{yid};
    count=1;
    plotdata=[];    
    sesslab=[];
for xx=1:length(sessnums)
sesslab=[sesslab '_' num2str(sessnums(xx))];
trialinfo=trialgrps(xx).trialinfo;
for itype=1:length(sessiontypes)
    %each da-lfp pair separate figure
    %& each window event type & session type same figure   
    targrow=find((contains({xinfo.siteda},targdasites(ida)) & ...
        contains({xinfo.sitelfp},lfpsites(ilfp)) & ...
        strcmp({xinfo.event},eventtypes(ievent))) & ... 
        contains({xinfo.sessionid},num2str(sessnums(xx))) &...
        contains({xinfo.sessiontype},sessiontypes(itype))==1);
    if ~isempty(targrow)
    %curdata=xinfo(targrow).dapos;
    curdata=getfield(xinfo(targrow),datype);    %'dapos' or 'daneg' types
    %trial ids for this data group for all/phase1/phase4..
    trialtypes=trialinfo(itype); %get trial groups for session
    grptids={};     %trial nums selected to index for this group of data
    ranj={};    %ran jitter, depends on trial nums each group
    for itrial=1:length(trialtypes.names)
        grptids(itrial).name=trialtypes.names{itrial};
        grptids(itrial).nums=...
            find(ismember(curdata.trials,trialtypes.nums{itrial})==1);
        %remove switch trials for all trial conditions except switch itself
        %1/3/2018
        if ~strcmp(trialtypes.names{itrial},'switch')
            %take out switch trials
            swtrials=trialtypes.nums{strcmp(trialtypes.names,'switch')};
            allbutsw=curdata.trials(ismember(curdata.trials,trialtypes.nums{itrial})...
                & ~ismember(curdata.trials,swtrials));
            grptids(itrial).nums=...
                find(ismember(curdata.trials,allbutsw)==1);
        end
        ranj{itrial}=rand(1,length(grptids(itrial).nums))*.1;
    end
    
    cntp=1;
for did=1:length(curgrp)
    %scroll parameters for group (ie. all timing variables)
    pid=ceil(cntp/2);        %plot axes, pairs
    grpdata=getfield(curdata,curgrp{did});
    xdata=[];
    xdata={grptids(:).name};
    gids={grptids(:).nums};
    ydata=[];
    xoff=0;
    if mod(cntp,2)
        %odd
        yyaxis(axa{pid},'left');
        xoff=-0.25+typeoff(itype);
        cmark=mark(itype,:);
    else
        %even
        yyaxis(axa{pid},'right');
        xoff=0.25+typeoff(itype);
        cmark=markr(itype,:);
    end
    set(axa{pid},'ycolor',cmark);
    for ix=1:length(xdata)
        plotdata(count).sessnum=sessnums(xx);
        plotdata(count).sesstype=sessiontypes{itype};
        plotdata(count).x=grptids(ix).name;
        plotdata(count).yvar=curgrp{did};
        plotdata(count).yavg=nan;
        plotdata(count).yci=nan;
        plotdata(count).ydata=nan;
        if ~isempty(grpdata(gids{ix}))
            ydata{ix}=grpdata(gids{ix});
            xplot=repmat(ix,1,length(ydata{ix}))+xoff+ranj{ix};
            %scatter(axa{pid},xplot,...
            %    ydata{ix},15,'.','markeredgecolor',cmark,...
            %    'MarkerEdgeAlpha',.3,'linewidth',.5);
            ciy=nanstd(ydata{ix})./sqrt(length(ydata{ix}(~isnan(ydata{ix}))))*1.96;
            avgy=nanmean(ydata{ix});
            %ciline=[avgy-ciy  avgy+ciy];
            %plot(axa{pid},[ix ix]+xoff-.05,ciline,'-','linewidth',1,'color',cmark);
            scatter(axa{pid},ix+xoff,...
                avgy,25,sessmark{xx},'markeredgecolor',cmark,...
                'MarkerEdgeAlpha',.5,'linewidth',.5);

            plotdata(count).yavg=avgy;
            plotdata(count).yci=ciy;
            plotdata(count).ydata=ydata{ix};
        end
        count=count+1;
    end
    if itype==1 && xx==1
        hy=ylabel(axa{pid},curgrp{did},'interpreter','none');
        if ~mod(cntp,2)
            %right axis        
            set(hy,'rotation',270,'units','pixels');
            hypos=get(hy,'position');
            set(hy,'position',[hypos(1)+15,hypos(2),hypos(3)]);  
        end
    end
    if itype==length(sessiontypes) && xx==length(sessnums)       
        targys=find(contains({plotdata.yvar},curgrp{did})==1 ...
         & contains({plotdata.x},'all')==1);
        getys=[plotdata(targys).ydata];    %all data for curren setup
        ystd=nanstd(getys);
        ymean=nanmean(getys);       
        if ~isnan(ystd) && ~isempty(ystd)
            set(axa{pid},'ylim',[ymean-ystd ymean+ystd]);  
        end
    end
    cntp=cntp+1;
    set(findall(axa{pid},'-property','FontSize'),'FontSize',fontsize)
    set(axa{pid},'xlim',[0.15 length(xdata)+0.75])
end
end
end
end
%plot lines between sessions big/small/targ to better visualize differences
%yvars=unique({plotdata(:).yvar}); %problem not sorted in orderplotted already
%xvars=unique({plotdata(:).x});
%{
yvars=ygrp{yid};
xvars=xdata;
cntp=1;
for iy=1:length(yvars)
for ix=1:length(xvars)
for isess=1:length(sessnums)
    pid=ceil(cntp/2);        %plot axes, pairs
    if mod(cntp,2)
        %odd
        yyaxis(axa{pid},'left');
        xoff=-0.25+typeoff;
        cmark=mark(isess,:);
    else
        %even
        yyaxis(axa{pid},'right');
        xoff=0.25+typeoff;
        cmark=mark(isess,:);
    end
    curids=find(contains({plotdata.yvar},yvars{iy})==1 ...
         & ismember({plotdata.x},xvars{ix})==1 ...
         & ismember([plotdata.sessnum],sessnums(isess))==1);
    ploty=[plotdata(curids).yavg];
    lp=plot(axa{pid},repmat(ix,1,length(ploty))+xoff,ploty,'-','linewidth',.5,'color',cmark);
    lp.Color(4)=0.3;
end
end
    cntp=cntp+1;
end
%}
    
save([savepath ynam{yid} '_' datype '_' eventtypes{ievent} ...
    '_' targdasites{ida} 'x' lfpsites{ilfp}],'plotdata');
savefig(figsess,[savepath ynam{yid} '_' datype '_' eventtypes{ievent} '_' targdasites{ida} 'x' lfpsites{ilfp}]);
saveas(figsess,[savepath ynam{yid} '_' datype '_' eventtypes{ievent} '_' targdasites{ida} 'x' lfpsites{ilfp}],'tif')
delete(findall(figsess,'type','text')) 
        
end
end
end
end

end