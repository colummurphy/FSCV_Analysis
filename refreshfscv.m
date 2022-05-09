function [Isub, BG]=refreshfscv(Iread,plotParam)
Isub=[];
BG=[];
selch=plotParam.selch;
BGstartpoint=plotParam.xsel;
BGavg=plotParam.BGavg;
t_end=plotParam.t_end;
orig_t_range=t_end;
%find boundaries, same for all channels
sizeMeas=size(Iread);
%set time scale
t_range=sizeMeas(2);        %total # timepoints
orig_t_range=t_range; %display(orig_t_range)
if t_end>t_range
    t_end=t_range;
end
%plot current subtracted from avg of ##BGavg points from BGstartpoint
if BGstartpoint<=0
    BGstartpoint=BGavg+1;
end
if (BGstartpoint-BGavg)<=0
    BGstartpoint=BGstartpoint+BGavg;
end
BGend=(round(BGstartpoint)+BGavg);
if BGend>t_range
    BGend=t_range;
end
for ii=1:length(selch)
    %initialize original plot based on loaded data
    BGvectors=Iread(:,(round(BGstartpoint)-BGavg):BGend); 
    BG=mean(BGvectors,2);         %mean ## sample points BG for subtraction
    BGrefmatrix=repmat(BG,1,orig_t_range);   %tile BGavg vector for entire meas matrix span
    Isub=Iread-BGrefmatrix;     %subtracted current data
end
%processed.Iread=Iread;

end