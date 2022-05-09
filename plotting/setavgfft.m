function evffts=setavgfft(data,fs,plotparam,hax,seltrials,varargin)
evffts={};
fontsize=14;
plotnum=0;
numtrials=size(data,1);
if isempty(seltrials)
    seltrials=1:numtrials;
end
argnum=1;
markers={};
iscsc=0;        %not csc data, ie fscv time stamps
cscname=[];
pulse=0;
lfp2=0;
logscale=0;
eye=0;
blink=0;
noplot=0;
sitename={};
freq=[];
cscale=[];
subbl=0;
pinkbl=0;
alignrelids=[]; %in form of [ts width; ts2 width2; etc.]
alignmarks={};      %alignment labels, origins to use
eids=[];
rmdc=0;     %flags for removing dc
zflag=0;     %flags for zflagalization'
%eec=0;       %# events with markers supplied % SOMEHOW THIS VARIABLE
%CAUSES mcr_is_external_function_available NO MATTER WHAT NAME
eventtypes={};      %strings of event labels
events={};          %event markers ts's
logplot=0;
aligntype=[];
alignevent={};
basets=0;
fax=[];
alnevt='targeye';

while argnum <= length(varargin)
    switch varargin{argnum}     
        case 'align'
            argnum=argnum+1;
            alignevent=varargin{argnum};
        case 'events'
            %event marker ids in fscv sample rates ids
            eec=eec+1;
            argnum=argnum+1;
            eventtypes{eec}=varargin{argnum};    %event label
            argnum=argnum+1;
            events{eec}=varargin{argnum};        %fscv ts's
        case 'cscale'
            argnum=argnum+1;
            plotparam.cscale=varargin{argnum};
        case 'plotnum'
            argnum=argnum+1;
            plotnum=varargin{argnum};
        case 'sitename'
            argnum=argnum+1;        %fscv ch site name
            cscname=varargin{argnum};
        case 'alignrelavgwin'
            %get avg spectrum specific time point relative alignidx
            argnum=argnum+1;
            alignmarks=varargin{argnum};    %alignment labels for origins
            %can be 'outcome' eg) alignidx, 'fix', 'target','fixeye'
            argnum=argnum+1;
            alignrelids=varargin{argnum}; %user provides time(s) relative alignidx and avg width (s)
        case 'freq'
            argnum=argnum+1;        %fscv ch site name
            freq=varargin{argnum};     %flag to plot lfp2            
        case 'log'
            logscale=1;
            logplot=1;
        case 'noplot'
            noplot=1;   %don't plot
        case 'rmdc'
            rmdc=1;
        case 'norm'
            zflag=1;
        case 'subbl'
            subbl=1;
        case 'pinkbl'
            pinkbl=1;
            subbl=1;
        case 'fax'
            argnum=argnum+1;        %user provides ax input for fft's
            fax=varargin{argnum};
        case 'subwin'
            argnum=argnum+1;    %user provides baseline sub ts
            basets=varargin{argnum};

           
    end
    argnum = argnum + 1;
end

fslabels=plotparam.samplespersec;           %sample rate for markers/win
win=plotparam.win;      %in fscv samples
rewardidx=plotparam.alignidx;
idsfixeye=plotparam.samplesfixeye{1};
idsfix=plotparam.samplesfix{1};
idstarg=plotparam.samplestarg{1};
idstargeye=plotparam.samplestargeye{1};

alignidx=rewardidx.*fs./plotparam.samplespersec; %convert to nlx sample rate
ts1=win(1)./plotparam.samplespersec;        %get absolute time
ts2=win(end)./plotparam.samplespersec;
tslfp=0:1/fs:(length(data(1,:))-1)/fs;
winids=find(tslfp>=ts1 & tslfp<=ts2);
idsfixeye=round(plotparam.tfixeye.*fs);
idsfix=round(plotparam.tfix.*fs);
idstarg=round(plotparam.ttarg.*fs);
idstargeye=round(plotparam.ttargeye.*fs);
    
shiftdata=[];
if isfield(plotparam,'shift')
    %shift data to targ or fix marker
    shiftdata=plotparam.shift;
end
for ii=1:length(events)
    %convert marker id's from fscv samples to nlx sample rate
    markers{ii}=events{ii}.*fs./plotparam.samplespersec;
end

idsrew=repmat(alignidx,1,size(data,1));
alnevt=plotparam.alnevt;
idsaln=[];
switch alnevt
    case 'fix'
        idsaln=idsfix;
    case 'fixeye'
        idsaln=idsfixeye;
    case 'targ'
        idsaln=idstarg;
    case 'targeye'
        idsaln=idstargeye;
    case 'outcome'
        %do not shift
        idsaln=idsrew;
end
info.alnevt=alnevt;
info.idsaln=idsaln;
alnidx=alignidx;
%shift signals according to alnevt idsaln
badtrls=[];
idsfixeye2=[];
idsfix2=[];
idstarg2=[];
idstargeye2=[];
idsrew2=[];

if ~isempty(idsaln) && ~strcmp(alnevt,'outcome')
    meanaln=nanmean(idsaln);
    stdaln=nanstd(idsaln);
    badtrls=(idsaln<meanaln-stdaln*20 | idsaln>meanaln+stdaln*20 | isnan(idsaln));
    data(badtrls,:)=nan;          %set signal to nan for bad trials
   % alnidx=nanmean(idsaln(~badtrls));       %new alignment idx for all trials
    daaln=[];
    for itrial=1:size(data,1)
        %numshift to align to meanalnwin
        if ~isnan(idsaln(itrial))
        shiftpos=round((alnidx-idsaln(itrial)));
        %shift da signal to indicated alignment idx
        daaln(itrial,:)=circshift(data(itrial,:),shiftpos,2);
        %shift all markers also based on alignment idx specified
        else
            daaln(itrial,:)=nan(1,size(data,2));
        end
        idsfixeye2(itrial)=idsfixeye(itrial)+shiftpos;
        idsfix2(itrial)=idsfix(itrial)+shiftpos;
        idstarg2(itrial)=idstarg(itrial)+shiftpos;
        idstargeye2(itrial)=idstargeye(itrial)+shiftpos;
        idsrew2(itrial)=idsrew(itrial)+shiftpos;
        if isnan(nanmean(data(itrial,:)))
            %if all signal is nan, designate bad trial so don't process
            badtrls(itrial)=1;
        end
    end
    %replace original markers with newly shifted events
    idsfixeye=idsfixeye2;
    idsfix=idsfix2;
    idstarg=idstarg2;
    idstargeye=idstargeye2;
    idsrew=idsrew2;
    %replace original data with shifted data
    data=daaln;
end


if noplot==0
    cla(hax);
fftwin=[1 .25];     %dan gibson sliding win default
if isfield(plotparam,'fftwin')
    fftwin=plotparam.fftwin;
end
interval=plotparam.interval;        %interval xtick marks def 2.5s
cminshow=plotparam.cminshow;
nantrials=find(isnan(data(:,1))==1);
%remove nantrials from seltrials
seltrials(find(ismember(seltrials,nantrials)==1))=[];   %remove bad trials
winidsextended=winids(1)-ceil((fftwin(1)-fftwin(2))*fs):winids(end)+ceil((fftwin(1)-fftwin(2))*fs);       %extended window to account for window chopped off in boundaries
winidsextended=winidsextended(winidsextended>0 & winidsextended<size(data,2));
winids=winidsextended;          %04/09/2019
datawin=data(seltrials,winids);     %data within window/selected trials
alldatasel=data(seltrials,:);     %all unwindowed data/selected trials
    
if isempty(freq)
    freq=[2.5 fs/2];
end
%[S,f]=chronux_mtspectrumc(datawin.*1e6,[1.8 1],0,fs,freq,0,0);  %convert to microvolts first
%[S,f]=chronux_mtspectrumc(datawin.*1e6,[3 5],0,fs,[freq(1)-5 freq(2)],0,0);  %convert to microvolts first
tapers=[1.8 1];     %dan gibson default
if isfield(plotparam,'ffttapers')
    tapers=plotparam.ffttapers;
end
dataforfft=datawin'.*1e6; %sel trials, windowed, conver to microvolts
alldataforfft=alldatasel'.*1e6; %selected trials, unwindowed, conver to microvolts

[S,winmid,f]=lfp_mtspecgram2(dataforfft,fftwin,tapers,0,fs,[freq(1) freq(2)],0,1,rmdc,zflag);
sfftwin=[3 0.25];   %wider window for more frequency resolution to get min/max pass bands
sfftwin=fftwin;
[Sh,winmidh,fh]=lfp_mtspecgram2(dataforfft,sfftwin,tapers,0,fs,[freq(1) freq(2)],0,1,rmdc,zflag);

fb=[];
bfft=[];
if subbl==1
    %get baseline spectrum, ie average across time & then across trials
   % for itrial=1:length(seltrials)
    %    tid=seltrials(itrial);
    %    tdata=dataforfft(tid,:);
           % baselinefft=alldataforfft(1:fftwin(1)*fs,:); %baseline data same size as samples        
        if basets==0
            baseid=10;      %arbitrary 10s into file recording
                      else
                   baseid=basets;
        end
                           baselineids=round(baseid*fs-fftwin(1)*fs/2):round(baseid*fs+fftwin(1)*fs/2);
            baselineids=baselineids(1:fftwin(1)*fs);
            baselinefft=alldataforfft(baselineids,:);
            baselineids=round(baseid*fs-sfftwin(1)*fs/2):round(baseid*fs+sfftwin(1)*fs/2);
            baselineids=baselineids(1:sfftwin(1)*fs);
            baselineffts=alldataforfft(baselineids,:);
      %  end
         [bfft,fb]=chronux_mtspectrumc(baselinefft,tapers,0,fs,[freq(1) freq(2)],0,1);  
                  [sbfft,sfb]=chronux_mtspectrumc(baselineffts,tapers,0,fs,[freq(1) freq(2)],0,1);  

   % [bfft,fb]=chronux_mtspectrumc(dataforfft,tapers,0,fs,[freq(1) freq(2)],0,1);  
   % bfft=resample(bfft,size(S,2),size(bfft,1));     %convert to resolution of S
   % fb=resample(fb,size(S,2),size(fb,2));
    if pinkbl==1
        bfft(bfft<0)=0;     %set negative #'s to 0 incase created during resample
        bl.sum=bfft;
        bl.f=fb;
        bl.N=1;
        bls.sum=sbfft;
        bls.f=sfb;
        bls.N=1;
        [pinkblf, a, b] = dg_fitPinkSpec(bl);
        bfft=pinkblf.sum;
        [pinkblfs, a, b] = dg_fitPinkSpec(bls);
        sbfft=pinkblfs.sum;
    end
    bldata=repmat(bfft', [size(S,1) 1]);
    S = S - bldata;
     sbldata=repmat(sbfft', [size(Sh,1) 1]);
    Sh = Sh - sbldata;   
end
plotf=find(f>freq(1));
%hI = image(hax,tslfp(winids(winmid)), f(plotf), S(plotf,:),'cdatamapping','scaled');
specdata=S(:,plotf)';
if logplot
    specdata=10*log10(specdata);
end
%get histogram of lfp powers/frequencies to find passbands at specific time
%period
targf=find(fh>12);       %above dc's
aa=var(Sh(:,targf),1);       %get variance across time to find highest variable f
maxf=fh(targf(find(aa==max(aa))));
fid=find(fh==maxf);
if isempty(fid)
    fid=round(median(find(fh>20 & fh<30))); %default beta range
    warning('peak center freq not found, default used');
end
maxv=max(Sh(:,fid));     %max power at this center freq
maxtsi=find(Sh(:,fid)==maxv);
%S2=S(maxtsi,:)-min(S(maxtsi,:));    %get spectrum at max point and shift above zero for log
%logS=10.*log10(S2);
maxS=Sh(maxtsi,:);
stdS=std(maxS);
thresS=find(maxS<=(max(maxS)-1.5*stdS));
lowf=fh(thresS(thresS<fid));
if ~isempty(lowf)
    lowf=lowf(end);     %low band
    if lowf<10
        lowf=fh(find(fh>=10 & fh<15));       %default freq low if <10
        lowf=lowf(1);
    end
else
    lowf=fh(find(fh>=10 & fh<15));       %default freq low if <10
    lowf=lowf(1);
end
highf=fh(thresS(thresS>fid));
if ~isempty(highf)
    highf=highf(1);
end
maxts=winmid(maxtsi)/fs+winids(1)/fs-alignidx/fs;       %ts relative to outcomeidx
evffts.maxts=maxts;     %store
evffts.maxf=maxf;
evffts.maxfl=lowf;
evffts.maxfh=highf;
evffts.maxS=maxS;
if ~isempty(fax)
    %plot on user supplied handle fft for each window
    plot(fax,fh,maxS)
    plot(fax,lowf, maxS(find(fh==lowf)),'ro')
    plot(fax, highf, maxS(find(fh==highf)),'ro')
    %pass band defined as 25% decrease below peak or 1.25 db
    legend(fax,[evffts.type 'max' 'fl' 'fh']);
    ylim([-stdS max(maxS)+stdS]);
    xlabel(fax,'freq (Hz)')
    title(fax,[plotparam.trialtype ' trials | lfp \times 10^6 (V^2) | ' cscname ...
        ' | tapers ' num2str(tapers(1)) ' / ' num2str(tapers(2)) ])
    ylabel(fax,'\muV^2'); 
    set(findall(fax,'-property','FontSize'),'FontSize',fontsize)
    hold(fax,'off');
end

[~,mintid]=min(abs(round(tslfp(winids(winmid))*fslabels)-plotparam.win(1)));
[~,maxtid]=min(abs(round(tslfp(winids(winmid))*fslabels)-plotparam.win(end)));
tdata=tslfp(winids(winmid(mintid:maxtid)));
pdata=specdata(:,mintid:maxtid);
hI = image(hax,tdata, f(plotf), pdata,'cdatamapping','scaled');

set(hax,'YDir','normal')        %flip y-axis values

%imagetrials=image(hax, da(seltrials,win(1):win(end)),'cdatamapping','scaled');
hold(hax,'on')
origpos=getpixelposition(hax);      %get original position 
%set color map to specified 
%colormap(hax,colortrial)

%plot event markers
targf=find(plotf>6);
means=nanmean(pdata(targf,:),1);  %average whole spectrum across frequencies
stds=nanstd(pdata(targf,:),0,1);     %std's across frequencies
means=nanmean(means);
stds=nanmean(stds);

if ~isempty(cscale) 
    caxis(hax,cscale); 
else
    
    cscale(1)=means-3*stds;
    cscale(1)=0;
    cscale(2)=means+4*stds;       
    caxis(hax,cscale); 
    
end
%plot event markers
colorm=plotparam.markercolors;
marks=tslfp(round([nanmean(idsfix(seltrials)) nanmean(idsfix(seltrials))]));
line(hax,marks,[f(1) f(end)],'Color',colorm(1,:),'LineWidth',1,'LineStyle','--')
marks=tslfp(round([nanmean(idsfixeye(seltrials)) nanmean(idsfixeye(seltrials))]));
line(hax,marks,[f(1) f(end)],'Color',colorm(2,:),'LineWidth',1,'LineStyle','--')
if ~strcmp(plotparam.trialtype,'fixbreak')
    if ~strcmp(plotparam.alnevt,'targ')
        marks=tslfp(round([nanmean(idstarg(seltrials)) nanmean(idstarg(seltrials))]));
        line(hax,marks,[f(1) f(end)],'Color',colorm(3,:),'LineWidth',1,'LineStyle','--')
    end
    if ~strcmp(plotparam.alnevt,'targeye')
        marks=tslfp(round([nanmean(idstargeye(seltrials)) nanmean(idstargeye(seltrials))]));
        line(hax,marks,[f(1) f(end)],'Color',colorm(4,:),'LineWidth',1,'LineStyle','--')
    end
end
if ~strcmp(plotparam.alnevt,'outcome')
   %not aligned to outcome, mark
    marks=tslfp(round([nanmean(idsrew(seltrials)) nanmean(idsrew(seltrials))]));
    line(hax,marks,[f(1) f(end)],'Color',colorm(5,:),'LineWidth',1,'LineStyle','--')
end

%organize plot
xticks=round(tslfp(winids(winmid(mintid)))):interval:round(tslfp(winids(winmid(maxtid))));
mintime=xticks(1)-round(alnidx/fs);
%mintime=round(tslfp((win(1)-alnidx)))/plotparam.samplespersec;
maxtime=xticks(end)-round(alnidx/fs);
%maxtime=round((win(end)-alnidx))/plotparam.samplespersec;
xticklabels=mintime:interval:maxtime;
xticklabels=round(xticklabels.*10)./10;
xticklabels=num2str(xticklabels');

set(hax,'xlim',[tslfp(winids(winmid(mintid))) tslfp(winids(winmid(maxtid)))])
set(hax,'ylim',[f(plotf(1)) f(plotf(end))])
set(hax,'XTick',xticks)
set(hax,'xticklabel',xticklabels)
set(hax,'tickdir','out','box','off')
xlabel(hax,'time (s)')
title(hax,[plotparam.trialtype ' | ' plotparam.triallabel ' | lfp \times 10^6 (V^2) | ' cscname ...
    ' | tapers ' num2str(tapers(1)) ' / ' num2str(tapers(2)) ' | win/step (s) ' ...
    num2str(fftwin(1)) ' / ' num2str(fftwin(2))])
ylabelleft=ylabel(hax,'freq'); 
h1=colorbar(hax,'eastoutside');
cpos = getpixelposition(h1);
ylabelbar=ylabel(h1,'V^2' ,'Rotation',270,'fontsize',fontsize); ypos = getpixelposition(ylabelbar);
cpos(3) = 15; 
set(h1,'Units','Pixels','Position', [cpos(1)+75 origpos(2) cpos(3) origpos(4)]);
%set(hax,'Units','Pixels','Position', [origpos(1) origpos(2) origpos(3) origpos(4)]);
set(ylabelbar,'Units','Pixels','Position', [ypos(1)+65 ypos(2)+cpos(4)/2 ]...
    ,'fontsize',fontsize);

 set(findall(hax,'-property','FontSize'),'FontSize',fontsize)
hold(hax,'off')

end


end
%}

