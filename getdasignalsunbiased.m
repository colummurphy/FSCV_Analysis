function getdasignals(timewin)
%input argument is fscv timestamp id's  to process
%updates processed.glitchids & processed.detected

global processed parameters plotParam

t1=timewin(1);
t2=timewin(2);
selch=plotParam.selch;
allchdata={};
processed.glitchids={};
processed.detected={};

for ii=1:length(selch)
    allchdata{selch(ii)}=processed.Iread(selch(ii)).data(:,t1:t2);
end

%find periods of recording where large & long-lasting movement
%artifacts & shared across channels
badids=findsigartifacts(allchdata,parameters);  %timestamps
badids=round(badids.*parameters.samplerate)+1;  %conver to samples

%chunk file
for ii=1:length(selch)
    disp(['detecting ch' num2str(selch(ii))])
    %find ids where HF glitches, to skip these data points in further
    %processing, 
    glitchids=findglitches(...
        processed.Iread(selch(ii)).rawdata(:,t1:t2),...
        parameters.glitchThres);
    %combine hf glitches w/ movement artifact periods & store
    processed.glitchids{selch(ii)}=sort([glitchids'; badids])';

    %find local "transient' increases in DA & their parameters
    processed.detected{selch(ii)}=detectdatransientschunks(...
        processed.Iread(selch(ii)).data(:,t1:t2),parameters,...
        processed.glitchids{selch(ii)},plotParam,'timestamp',...
        processed.LFPread.LFPts(1)+(t1-1)/parameters.samplerate); 

end
        
end