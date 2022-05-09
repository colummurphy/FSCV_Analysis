function [evffts, maxfft]=getavgspec(data,fs,ts,alnts,hax,varargin)
evffts={};
fs=1000;
fontsize=14;
maxfft={};
numtrials=size(data,1);
argnum=1;
markers={};
cscale=[];
subbl=0;
pinkbl=0;
rmdc=0;     %flags for removing dc
zflag=0;     %flags for zflagalization'
eventtypes={};      %strings of event labels
logplot=0;
basets=0;
fax=[];
freq=[];
fftwin=[1 .25];     %dan gibson sliding win default
win=[];
tapers=[1.8 1];     %dan gibson default
markers={};
data2=[];
alnts2=[];
markers2={};
subfft=0;
targts=0;
getmax=0;
while argnum <= length(varargin)
    switch varargin{argnum}          
        case 'cscale'
            argnum=argnum+1;
            cscale=varargin{argnum};         
        case 'freq'
            argnum=argnum+1;        %fscv ch site name
            freq=varargin{argnum};     %flag to plot lfp2      
        case 'fftwin'
            argnum=argnum+1;
            fftwin=varargin{argnum};
        case 'tapers'
            argnum=argnum+1;
            tapers=varargin{argnum};
        case 'log'
            logplot=1;        
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
        case 'win'
            argnum=argnum+1;
            win=varargin{argnum};
        case 'marks'
            argnum=argnum+1;
            markers=varargin{argnum};
        case 'sub'
            %subtract another spectra to show difference among conditions
            %supply signals x trials, and alnts for these data
            argnum=argnum+1;
            data2=varargin{argnum};
            argnum=argnum+1;
            alnts2=varargin{argnum};
            subfft=1;
        case 'marks2'
            %for sub data
            argnum=argnum+1;
            markers2=varargin{argnum};
        case 'targts'
            %for single point spectra
            argnum=argnum+1;
            targts=varargin{argnum};
        case 'getmaxfft'
            %get fft at point of max variance across time for broad beta
            %power, store in basefft variable
            getmax=1;
    end
    argnum = argnum + 1;
end
outcomets=ts(median(1:length(ts)));
idsaln=round(alnts.*fs);
shiftids=round(outcomets.*fs-alnts.*fs);     %numshift to align to meanalnwin    
daaln=nan(size(data));
daaln2=nan(size(data2));
markts={};
markts2={};
for itrial=1:size(data,1)
    %shift signal to alignment idx
    daaln(itrial,:)=circshift(data(itrial,:),shiftids(itrial),2);
    if ~isempty(markers)
        for ix=1:length(markers)
            markts{ix}(itrial)=markers{ix}(itrial)+outcomets-alnts(itrial);
        end
    end
end
if subfft
    shiftids2=round(outcomets.*fs-alnts2.*fs);     %numshift to align to meanalnwin    
    for itrial=1:size(data2,1)
    %shift signal to alignment idx
    daaln2(itrial,:)=circshift(data2(itrial,:),shiftids2(itrial),2);
    if ~isempty(markers2)
        for ix=1:length(markers2)
            markts2{ix}(itrial)=markers2{ix}(itrial)+outcomets-alnts2(itrial);
        end
    end
    end
end
data=daaln;
data2=daaln2;
cla(hax);
interval=1;
relts=ts-outcomets;
winids=find(relts>=win(1) & relts<=win(2));
winidsextended=winids(1)-ceil((fftwin(1)-fftwin(2))*fs):winids(end)+ceil((fftwin(1)-fftwin(2))*fs);       %extended window to account for window chopped off in boundaries
winidsextended=winidsextended(winidsextended>0 & winidsextended<size(data,2));
winids=winidsextended;          %04/09/2019
datawin=data(:,winids);     %data within window/selected trials
alldatasel=data;     %all unwindowed data/selected trials
if isempty(freq)
    freq=[2.5 fs/10];
end
%[S,f]=chronux_mtspectrumc(datawin.*1e6,[1.8 1],0,fs,freq,0,0);  %convert to microvolts first
%[S,f]=chronux_mtspectrumc(datawin.*1e6,[3 5],0,fs,[freq(1)-5 freq(2)],0,0);  %convert to microvolts first
dataforfft=datawin'.*1e6; %sel trials, windowed, conver to microvolts
alldataforfft=alldatasel'.*1e6; %selected trials, unwindowed, conver to microvolts

[S,winmid,f]=lfp_mtspecgram2(dataforfft,fftwin,tapers,0,fs,[freq(1) freq(2)],0,1,rmdc,zflag);
if subfft
    dataforfft=data2(:,winids)'.*1e6; %sel trials, windowed, conver to microvolts
    [S2,winmid2,f2]=lfp_mtspecgram2(dataforfft,fftwin,tapers,0,fs,[freq(1) freq(2)],0,1,rmdc,zflag);
    S=S-S2;
end
fb=[];
bfft=[];
if subbl==1
    %get baseline spectrum, ie average across time & then across trials
    if basets==0
        baseid=win(1)+round(outcomets)-3;      %arbitrary window
        else
        baseid=basets+round(outcomets);
    end
    baselineids=round(baseid(1)*fs-fftwin(1)*fs/2):round(baseid(1)*fs+fftwin(1)*fs/2);
    baselineids=baselineids(baselineids>0 & baselineids<=size(data,2));
    baselineffts=alldataforfft(baselineids,:);
    [bfft,fb]=chronux_mtspectrumc(baselineffts,tapers,0,fs,[freq(1) freq(2)],0,1);  
    if pinkbl==1
        bfft(bfft<0)=0;     %set negative #'s to 0 incase created during resample
        bl.sum=bfft;
        bl.f=fb;
        bl.N=1;      
        [pinkblf, a, b] = dg_fitPinkSpec(bl);
        bfft=pinkblf.sum;      
    end
    bldata=repmat(bfft', [size(S,1) 1]);
    S = S - bldata; 
end
if getmax
    %get spectrum at time point of max broad beta across averaged trial
    %spectrogram
    offset=.5;
    windowids=find(((winids(winmid)-winids(winmid(1)))./fs)>offset & ((winids(winmid)-winids(winmid(1)))./fs)<(winids(winmid(end))-winids(winmid(1)))./fs-offset);
    [maxts maxtsid]=max(nanmean(S(windowids,find(f>=13 & f<=33)),2));  %time point of max broad beta
    maxtsid=maxtsid+windowids(1)-1;
    baseid=winids(winmid(maxtsid));    
    baselineids=round(baseid(1)-fftwin(1)*fs/2):round(baseid(1)+fftwin(1)*fs/2);
    baselineids=baselineids(baselineids>0 & baselineids<=size(data,2));
    baselineffts=alldataforfft(baselineids,:);
    [bfft,fb]=chronux_mtspectrumc(baselineffts,tapers,0,fs,[freq(1) freq(2)],0,1);  
    bfft(bfft<0)=0;     %set negative #'s to 0 incase created during resample
    bl.sum=bfft;
    bl.f=fb;
    bl.N=1;      
    [pinkblf, a, b] = dg_fitPinkSpec(bl);
    pfft=pinkblf.sum;  
    normfft = bfft - pfft; 
    maxfft.freq=fb';
    maxfft.bfft=bfft;
    maxfft.pfft=pfft;
    maxfft.normfft=normfft;
    maxfft.maxts=baseid./fs-round(outcomets);       %ts relative to alignment event
    %find peak and width of beta spectrum
    fsearchids=find(fb' >= 13 & fb' <=33);
    [peakspec,fpeakid]=max(normfft(fsearchids));
    fpeakid=fpeakid+fsearchids(1)-1;
    fpeak=fb(fpeakid);
    
    %differentiate spectra to find "hump" defined by decreasing slopes
    %around peak
    diffspec=diff(normfft);
    diffspec=[diffspec(1); diff(normfft)]; 
    winboundarymin=find(diffspec(1:fpeakid)>0)-1;       %rising slop frequencies
    winboundarymin=winboundarymin(winboundarymin>0);
    winboundarymax=find(diffspec(fpeakid:end)<0)+fpeakid-1;       %falling slope
   % betaedges=find(normfft<0 & fb' >= 10 & fb' <=40);
    discontinuity=find(diff(winboundarymax)>1);     %find id until next rising slope
    if ~isempty(discontinuity)
    winboundarymax=winboundarymax(1:discontinuity(1));
    end
    discontinuity=find(diff(winboundarymin)>1);  
        if ~isempty(discontinuity)
    winboundarymin=winboundarymin(discontinuity(end)+1:end);
        end
    [minspecrising,minfid]=min(normfft(winboundarymin));
    minfid=minfid+winboundarymin(1)-1;
    fmin=fb(minfid);
    [minspecfalling,maxfid]=min(normfft(winboundarymax));
    maxfid=maxfid+winboundarymax(1)-1;
    fmax=fb(maxfid);
    
   % fmin=fb(betaedges(discontinuity(1)));
   % fmax=fb(betaedges(discontinuity(end)+1));
    maxfft.fmin=fmin;
    maxfft.fmax=fmax;
    maxfft.fpeak=fpeak;
    
    bldata=repmat(pfft', [size(S,1) 1]);
    S = S - bldata; 
end
plotf=find(f>freq(1));
specdata=S(:,plotf)';
if logplot
    specdata=10*log10(specdata);
end

%[~,mintid]=min(abs(round(ts(winids(winmid))*fs)-win(1)*fs));
%[~,maxtid]=min(abs(round(ts(winids(winmid))*fs)-win(end)*fs));
%tdata=tslfp(winids(winmid(mintid:maxtid)));
tswin=relts(winids(winmid));
pdata=specdata;

hI = image(hax,tswin, f(plotf), pdata,'cdatamapping','scaled');

set(hax,'YDir','normal')        %flip y-axis values

%imagetrials=image(hax, da(seltrials,win(1):win(end)),'cdatamapping','scaled');
hold(hax,'on')
origpos=getpixelposition(hax);      %get original position 
%set color map to specified 
%colormap(hax,colortrial)

%plot event markers
if ~isempty(markts)
    for ix=1:length(markts)
        plotmark=repmat(nanmean(markts{ix}),1,2)-outcomets;
        af=line(hax,plotmark,[f(1) f(end)],'Color',[1 .4 1],'LineWidth',1,'LineStyle','--');
        af.Color(4)=0.75;
    end
end
if getmax
    %plot lines showing spectra width
     plotmark=[maxfft.fmin maxfft.fmin];
     af=line(hax,[tswin(1) tswin(end)],plotmark,'Color',[0 0 0],'LineWidth',1,'LineStyle','--');
      af.Color(4)=.5;
     af=line(hax,[tswin(1) tswin(end)],[maxfft.fmax maxfft.fmax],'Color',[0 0 0],'LineWidth',1,'LineStyle','--');
    af.Color(4)=.5;
   % af=line(hax,[tswin(1) tswin(end)],[maxfft.fpeak maxfft.fpeak],'Color',[0 1 0],'LineWidth',1,'LineStyle','--');
  %  af.Color(4)=.5;
end

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
    if subfft
        cscale(1)=means-4*stds;
        cscale(2)=means+4*stds;
    end
    caxis(hax,cscale);     
end
evffts.cscale=cscale;
%organize plot
xticks=win(1):interval:win(end);
xticklabels=round(xticks.*10)./10;
xticklabels=num2str(xticklabels');

rbmap=redblue(128);
colormap(rbmap);
set(hax,'xlim',[win(1) win(end)])
set(hax,'ylim',[f(plotf(1)) f(plotf(end))])
set(hax,'XTick',xticks)
set(hax,'xticklabel',xticklabels)
set(hax,'tickdir','out','box','off')
xlabel(hax,'time (s)')
ylabelleft=ylabel(hax,'freq'); 
set(findall(hax,'-property','FontSize'),'FontSize',fontsize)
h1=colorbar(hax,'eastoutside');
cpos = getpixelposition(h1);
ylabelbar=ylabel(h1,'\muV^2' ,'Rotation',270,'fontsize',fontsize); ypos = getpixelposition(ylabelbar);
cpos(3) = 15; 
set(hax,'units','pixels');
posax=get(hax,'position');
set(h1,'Units','Pixels','Position', [cpos(1)+65 posax(2) cpos(3) posax(4)]);
%set(hax,'Units','Pixels','Position', [origpos(1) origpos(2) origpos(3) origpos(4)]);
set(ylabelbar,'Units','Pixels','Position', [ypos(1)+65 ypos(2)+cpos(4)/2 ]...
    ,'fontsize',fontsize);

hold(hax,'off')

end
