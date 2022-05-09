function processed = spectralSub(samples,samplerate,varargin)
%spectral subtraction to remove FSCV harmonic noise
%07/01/2018 - make chunk smaller and increase overlap so more like slideing
%window to preven tartifacts between pasting windows
%03/17/2018 - basic zeroing of peaks to remove harmonic artifacts
%03/21/2018 - convert fft signal to abs/phase components & interpolate
%around artifact harmonics depending on width of peaks
%using sliding windo w 07/02/2018
argnum=1;
cleo=0;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'cleo'
            cleo=1;

    end
    argnum=argnum+1;
end
overlap=15;             % percent overlap of signal in time, increased form 5 to 10 07/01/2018
overlap=50;             % critical param percent overlap of signal in time, increased form 5 to 10 07/01/2018
if cleo==1
    overlap=90;             % cleo
end
chunkwidth=1e5;             %num samples per chunk before 07/01/2018
padwidth=2;
gap=3;      %gap between signals for linear interpolation to prevent glitches where pasted together
gap=0;
%change again 07/01/2018, seems like if too short, cannot fully
%characterize noisy spectrum to remove, maybe increase overlap
chunkwidth=ceil(samplerate)*10;             %num samples per chunk = samples per second
chunkwidth=ceil(samplerate)*5;             %num samples per chunk = samples per second
if cleo==1
chunkwidth=ceil(samplerate)*3;             %cleo
end
%added to calculate smoothing width for fillgap based on length of samples
avgSpectralWidth=ceil(chunkwidth/samplerate*1.25);             %number of samples +/- peak to average across for interpolation
if avgSpectralWidth<=2
    avgSpectralWidth=3;
end
padwidth=floor(avgSpectralWidth/2);
padwidth=0;
%avgSpectralWidth=5;             %number of samples +/- peak to average across for interpolation
chunkwidthoverlap=chunkwidth*overlap/100;
chunkstep=chunkwidth-chunkwidthoverlap;     %step size for sliding window
harmonics=10:10:1000;            %harmonics of FSCV signal
%harmonics=10:10:8000;            %harmonics of FSCV signal, almost no
%difference with higher freq limit.

harmonicWidth=.1;           %needed to discretize rounded frequences estimated after FFT to find peak freq's
%harmonics=[harmonics -harmonics];   %along negative freq axis as well
harmonicRange=[];
peakWidthLimit=5/samplerate*chunkwidth;     %if more than 5 Hz +/- width do not use (sample rate/chunkwidth is freq resolutiOn)

%Identify peaks and 0.5 Hz around them
for ii=1:length(harmonics)
    if ii==1
        harmonicRange=harmonics(ii)-harmonicWidth:0.1:harmonics(ii)+harmonicWidth;
    else
        harmonicRange=[harmonicRange harmonics(ii)-harmonicWidth:0.1:harmonics(ii)+harmonicWidth];
    end
end
%spectral subtraction on all selected channels
%chunks=1:chunkwidth:length(samples);       %perform subtraction in chunks
chunks=1:chunkstep:length(samples);       %perform subtraction in chunks
freq=0;
processed=[];
midchunkovup=floor(chunkwidthoverlap/2);
midchunkovdown=ceil(chunkwidthoverlap/2);

if midchunkovup+midchunkovdown~=chunkwidthoverlap
    error
end
winsb=[chunks+chunkwidth-1-midchunkovup];
winsa=[chunks+midchunkovdown];

for chunkID=1:length(chunks)
    samplesChunk=[];
    if chunkID==1
        win=1:chunkwidth;
        win=win(win>=1 & win<=length(samples));
        samplesChunk=samples(win);
    end
    if chunkID>1
        win=chunks(chunkID):chunks(chunkID)+chunkwidth-1;
        win=win(win>=1 & win<=length(samples));
        samplesChunk=samples(win);
    end    
    window = tukeywin(length(samplesChunk),.1);      %tukey window in time domain to avoid spectral leakage at edges, 
    %.1 r parameter means 10% signals (5% each side) attenuated, but rest is 1
    windowedSamples=window'.*samplesChunk;
    samplesFFT = reshape(fft(windowedSamples), [], 1);  %get fft
    samplesFFTshift=fftshift(samplesFFT);       %shift spectra to normal frequency scale (ie. mirroring 0)
        df=samplerate./length(samplesFFT);
        freq = (-samplerate/2:df:samplerate/2-df)' + mod(length(samplesFFT),2)*df/2;
   freq0ID=fix(length(samplesFFT)/2)+1;
    denoiseShift=samplesFFTshift;
    denoiseShiftMag=abs(denoiseShift);
    denoiseShiftPhase=angle(denoiseShift);

    freqRounded=round(freq*10)/10;      %get rounded to nearest .1  place
    harmonicFind=ismember(freqRounded,harmonics);   %find integer peaks only
    harmonicFind=ismember(freqRounded,harmonicRange);   %find frequency
    %IDs of noisy harmonics (Default width set above)
    harmonicIDs=find(harmonicFind>0); 
    %harmonicWidthSamps=harmonicIDs(1)-freq0ID-lowerBound(end);
    %for each predicted harmonic noise peak, find individual widths (based
    %on 5% attenuation of mag(peak) and interpolate around this
    
    %smooth fft to find peak widths
   % smoothfft=smoothts(denoiseShiftMag,'g',avgSpectralWidth);       
    smoothfft=smooth(denoiseShiftMag,avgSpectralWidth);
    %/add 07/01/2018 to find peak widths because sometimes glitch in peaks 
    %makes the width shorter than it should be for removal
    for ii=1:length(harmonicIDs)
        %lowerBound=find(denoiseShiftMag(freq0ID:harmonicIDs(ii))<.05*denoiseShiftMag(harmonicIDs(ii))); %determine width of noisy peaks, find border
        %upperBound=find(denoiseShiftMag(harmonicIDs(ii):end)<.05*denoiseShiftMag(harmonicIDs(ii)));
        lowerBound=find(smoothfft(freq0ID:harmonicIDs(ii))<.75*max(smoothfft(harmonicIDs(ii)-avgSpectralWidth:harmonicIDs(ii)+avgSpectralWidth))); %determine width of noisy peaks, find border
        upperBound=find(smoothfft(harmonicIDs(ii):end)<.75*max(smoothfft(harmonicIDs(ii)-avgSpectralWidth:harmonicIDs(ii)+avgSpectralWidth)));
        %if do not found bounds skip.
        if isempty(lowerBound) || isempty(upperBound)
            continue
        end
        harmonicWidthDown=harmonicIDs(ii)-freq0ID-lowerBound(end)+padwidth;
        harmonicWidthUp=upperBound(1)+padwidth;        
      %  harmonicWidthDownP=harmonicIDs(ii)-freq0ID-lowerBoundP(end)+2;
       % harmonicWidthUpP=upperBoundP(1)+2;        

        %if bounds are too wide skip. maybe peak is already too low
        if harmonicWidthDown>peakWidthLimit || harmonicWidthUp>peakWidthLimit
            continue
        end

        bounds=[harmonicIDs(ii)-harmonicWidthDown harmonicIDs(ii)+harmonicWidthUp];
        %remove nan since no longer user fillgaps 07/03/2018
       %  denoiseShiftMag(bounds(1):bounds(2))=nan;      %05/05/2018, replace linear interpolation
      % denoiseShiftPhase(bounds(1):bounds(2))=nan;
        %denoiseShiftMag(bounds(1):bounds(2))=smoothts(denoiseShiftMag,'g',avgSpectralWidth);
        xx2=bounds(2):bounds(2)+avgSpectralWidth;
        xx=bounds(1)-avgSpectralWidth:bounds(1);
        xq=xx(1):xx2(end);
        xxx=[xx xx2];
        xxx=unique(xxx);
        interpolatedpoints=interp1(xxx, denoiseShiftMag(xxx), xq,'pchip');
        interpolatedphase=interp1(xxx, denoiseShiftPhase(xxx), xq,'pchip');
        denoiseShiftMag(xq)=interpolatedpoints';
        denoiseShiftPhase(xq)=interpolatedphase';
        %denoiseShiftMag(xq)=interpolatedpointsp';
    end
    
    %fill gaps takes too long 07/03/208, changed to interpolate above
 % denoiseShiftMag=fillgaps(denoiseShiftMag,avgSpectralWidth);     %interpolate around nanned harmonic peaks in a smooth manner 05/05/2018
  %  denoiseShiftPhase=fillgaps(denoiseShiftPhase,avgSpectralWidth);     %interpolate around nanned harmonic peaks in a smooth manner 05/05/2018

    
    denoiseShift=denoiseShiftMag.*exp(1j.*denoiseShiftPhase);    %retrive complex format 
    denoiseShiftMirror=flipud(denoiseShift(freq0ID:end));  %mirror of pos spectrum
    conjDenoise=conj(denoiseShiftMirror);
    %denoiseShift(1:freq0ID)=conj(denoiseShiftMirror);  %complex conjugate to get neg spectrum
    if mod(length(samplesChunk),2)==1
        %if odd
        denoiseShift(1:freq0ID-1)=conjDenoise(1:end-1);
    else
        %if even
        denoiseShiftMirror=flipud(denoiseShift(freq0ID-1:end));
        conjDenoise=conj(denoiseShiftMirror);
        denoiseShift(2:freq0ID)=conjDenoise(1:end-1);
    end
    %denoiseShift(harmonicIDs)=0;        %set spectra to 0 at harmonic frequencies of artifact
    %denoiseShift(harmonicIDs)=denoiseShift(harmonicIDs)*.01;
    denoiseFFT=ifftshift(denoiseShift); %shift denoised signal back to normal frequency scale
    processedChunk=ifft(denoiseFFT);     %convert spectrum back to timescale
    %check imag components should be very small or something ws
    if sum(imag(processedChunk))>1e-10
        %if greater than .1 nanovolt
        warning(['imaginary components from inverse FFT too high at chunk ' num2str(chunkID) '/' num2str(length(chunks))]);
    end
    if chunkID==1
        processed=real(processedChunk);
    elseif chunks(chunkID)+chunkwidth-1<=length(samples)
        %current chunk period below sample length, nothing cut off
       processed1=processed(1:end-midchunkovup-gap);
        processed2=real(processedChunk(midchunkovdown+gap+1:end));
       % pastinginterval(1:gap*2,1)=nan;
       % processed=[processed1; pastinginterval; processed2];
       % processed=fillgaps(processed,avgSpectralWidth);     %interpolate around nanned harmonic peaks in a smooth manner 05/05/2018
        %fill gaps function too time consuming, live with glitches between
        %chunks
        processed=[processed1; processed2];
    elseif chunks(chunkID)+chunkwidth-1>length(samples)
        %over sample length, stop loop
        processed1=processed(1:end-midchunkovup-gap);
        processed2=real(processedChunk(midchunkovdown+gap+1:end));
                processed=[processed1; processed2];

       % pastinginterval(1:gap*2,1)=nan;
       % processed=[processed1; pastinginterval; processed2];
        %  processed=fillgaps(processed,avgSpectralWidth);  
          break
       % processed=[processed(1:end-chunkwidthoverlap); processedChunk(chunkwidthoverlap+1:end)];
    end
end
processed=real(processed);      %only save real components, 
if length(processed)~=length(samples)
    warning('not equal lengths')
    %may not equal each other if last chunk had even # of samples
    processed=[processed; 0];        %add zero at end
end
%check spec sub fft
%{
fig1=figure;
set(fig1,'color',[1 1 1]);
hold on;
plot(freq,abs(samplesFFTshift))
xlabel('frequency (Hz)')
ylabel('magnitude')
xlim([0 100])
title('LFP fft');
%legend({'original signal' ,'FSCV harmonics subtracted'});
ylim([0 5]); ylim([0 1]);

xlim([-1 100]);
plot(freq,abs(denoiseShift))
hold off;

fig2=figure;
set(fig2,'color',[1 1 1], 'units','pixels','position',[100 100 800 400]);
hold on;
ts=0:length(samples)-1;
ts=ts./samplerate;
plot(ts,samples)
xlabel('time (s)')
ylabel('amplitude (V)')
title('LFP (time-domain)');
ylim([-2e-4 2e-4]);
xlim([35 36]);
plot(ts,processed)
hold off;
%}


%Imaginary componets should be small (should add some check) because of sample rate jitter
end