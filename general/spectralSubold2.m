function processed = spectralSub(samples,samplerate)
%spectral subtraction to remove FSCV harmonic noise
%07/01/2018 - make chunk smaller and increase overlap so more like slideing
%window to preven tartifacts between pasting windows
%03/17/2018 - basic zeroing of peaks to remove harmonic artifacts
%03/21/2018 - convert fft signal to abs/phase components & interpolate
%around artifact harmonics depending on width of peaks
overlap=15;             % percent overlap of signal in time, increased form 5 to 10 07/01/2018
chunkwidth=1e5;             %num samples per chunk before 07/01/2018
padwidth=2;
%change again 07/01/2018, seems like if too short, cannot fully
%characterize noisy spectrum to remove, maybe increase overlap
chunkwidth=ceil(samplerate)*10;             %num samples per chunk = samples per second
chunkwidth=ceil(samplerate)*5;             %num samples per chunk = samples per second

%added to calculate smoothing width for fillgap based on length of samples
avgSpectralWidth=ceil(chunkwidth/samplerate*1.25);             %number of samples +/- peak to average across for interpolation
if avgSpectralWidth<=2
    avgSpectralWidth=3;
end
padwidth=round(avgSpectralWidth/2);
%avgSpectralWidth=5;             %number of samples +/- peak to average across for interpolation
chunkwidthoverlap=chunkwidth*overlap/100;
%harmonics=10:10:1000;            %harmonics of FSCV signal
harmonics=10:10:8000;            %harmonics of FSCV signal

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
chunks=1:chunkwidth:length(samples);       %perform subtraction in chunks
processed=[];
processedChunks=[];
freq=[];
skiplast=0;
for chunkID=1:length(chunks)
    samplesChunk=[];
    if chunkID==1
        samplesChunk=samples(1:chunks(chunkID+1)+chunkwidthoverlap);
    end
    if chunkID>1 && chunkID<length(chunks)
        if chunks(chunkID+1)+chunkwidthoverlap>length(samples)
                    samplesChunk=samples(chunks(chunkID)-chunkwidthoverlap:end);
                    skiplast=1;
        else
        samplesChunk=samples(chunks(chunkID)-chunkwidthoverlap:chunks(chunkID+1)+chunkwidthoverlap);
        end
    end
    if chunkID==length(chunks)
        if skiplast==1
            %already stored data processed to end because chunk overlap
            %overlapped to end of recording, quit loop & finish process
            continue;
        end
        samplesChunk=samples(chunks(chunkID(end))-chunkwidthoverlap:end);
    end
    window = tukeywin(length(samplesChunk),.1);      %tukey window in time domain to avoid spectral leakage at edges, 
    %.1 r parameter means 10% signals (5% each side) attenuated, but rest is 1
    windowedSamples=window'.*samplesChunk;
    samplesFFT = reshape(fft(windowedSamples), [], 1);  %get fft
    samplesFFTshift=fftshift(samplesFFT);       %shift spectra to normal frequency scale (ie. mirroring 0)
    if mod(length(samplesChunk),2)==0
        %if even
        %get frequency scale
        %usually first number in FFT is not equal to end of spectra so take
        %out
      % warning('even');
        samplesFFT=samplesFFT(2:end);
        samplesFFTshift=samplesFFTshift(2:end);
        freq=samplerate*(-(fix(length(samplesFFT)/2)):(length(samplesFFT)/2))/length(samplesFFT);
    else
        %if odd
        %get frequency scale
        freq=samplerate*(-(fix(length(samplesFFT)/2)):(length(samplesFFT)/2))/length(samplesFFT);   
    end
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
    smoothfft=smoothts(denoiseShiftMag,'g',avgSpectralWidth);       
    %/add 07/01/2018 to find peak widths because sometimes glitch in peaks 
    %makes the width shorter than it should be for removal
    for ii=1:length(harmonicIDs)
        %lowerBound=find(denoiseShiftMag(freq0ID:harmonicIDs(ii))<.05*denoiseShiftMag(harmonicIDs(ii))); %determine width of noisy peaks, find border
        %upperBound=find(denoiseShiftMag(harmonicIDs(ii):end)<.05*denoiseShiftMag(harmonicIDs(ii)));
        lowerBound=find(smoothfft(freq0ID:harmonicIDs(ii))<.05*smoothfft(harmonicIDs(ii))); %determine width of noisy peaks, find border
        upperBound=find(smoothfft(harmonicIDs(ii):end)<.05*smoothfft(harmonicIDs(ii)));
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
     %   boundsP=[harmonicIDs(ii)-harmonicWidthDownP harmonicIDs(ii)+harmonicWidthUpP];
       % avgBefore=mean(denoiseShiftMag(bounds(1)-avgSpectralWidth:bounds(1)),'omitnan');    %mag's
       % avgAfter=mean(denoiseShiftMag(bounds(2):bounds(2)+avgSpectralWidth),'omitnan');
       % avgBeforeP=mean(denoiseShiftPhase(boundsP(1)-avgSpectralWidth:boundsP(1)),'omitnan');   %phases
       % avgAfterP=mean(denoiseShiftPhase(boundsP(2):boundsP(2)+avgSpectralWidth),'omitnan');
      %  avgBeforeP=mean(denoiseShiftPhase(bounds(1)-avgSpectralWidth:bounds(1)),'omitnan');    %mag's
      %  avgAfterP=mean(denoiseShiftPhase(bounds(2):bounds(2)+avgSpectralWidth),'omitnan');

       % denoiseShiftMag(bounds(1):bounds(2))=linspace(avgBefore,avgAfter,bounds(2)-bounds(1)+1);
       % denoiseShiftPhase(bounds(1):bounds(2))=linspace(avgBeforeP,avgAfterP,bounds(2)-bounds(1)+1);
         denoiseShiftMag(bounds(1):bounds(2))=nan;      %05/05/2018, replace linear interpolation
       denoiseShiftPhase(bounds(1):bounds(2))=nan;

    end
    
  denoiseShiftMag=fillgaps(denoiseShiftMag,avgSpectralWidth);     %interpolate around nanned harmonic peaks in a smooth manner 05/05/2018
    denoiseShiftPhase=fillgaps(denoiseShiftPhase,avgSpectralWidth);     %interpolate around nanned harmonic peaks in a smooth manner 05/05/2018

    
    denoiseShift=denoiseShiftMag.*exp(1j.*denoiseShiftPhase);    %retrive complex format 
    denoiseShiftMirror=flipud(denoiseShift(freq0ID:end));  %mirror of pos spectrum
    conjDenoise=conj(denoiseShiftMirror);
    %denoiseShift(1:freq0ID)=conj(denoiseShiftMirror);  %complex conjugate to get neg spectrum
    denoiseShift(1:freq0ID-1)=conjDenoise(1:end-1);
    %denoiseShift(harmonicIDs)=0;        %set spectra to 0 at harmonic frequencies of artifact
    %denoiseShift(harmonicIDs)=denoiseShift(harmonicIDs)*.01;
    denoiseFFT=ifftshift(denoiseShift); %shift denoised signal back to normal frequency scale
    processedChunk=ifft(denoiseFFT);     %convert spectrum back to timescale
    %check imag components should be very small or something wrong with
    %ifft
    if sum(imag(processedChunk))>1e-10
        %if greater than .1 nanovolt
        warning(['imaginary components from inverse FFT too high at chunk ' num2str(chunkID) '/' num2str(length(chunks))]);
    end
    if chunkID==1
        processed=processedChunk;
    else
        processed=[processed(1:end-chunkwidthoverlap); processedChunk(chunkwidthoverlap+2:end)];
    end
end
processed=real(processed);      %only save real components, 
if length(processed)~=length(samples)
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