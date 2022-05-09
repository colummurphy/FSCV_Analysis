function info = calcdata(Iread, Isub, cvs,selch)
global plotParam parameters
%calc noise SNR etc for display
numch=length(selch);
info=[];
xsel=plotParam.xsel;
avgbg=plotParam.BGavg;


for ii=1:numch
    subData=Isub(selch(ii)).data;
    measData=Iread(selch(ii)).data;
    if ~isempty(subData)
        if xsel>size(subData,2)
            xsel=size(subData,2);
        end
        Imax=max(subData(:,xsel));      %select max current along single time point as selected
        V_max_ind=find(flipud(subData(:,xsel))>=max(subData(:,xsel)));     
        index_DA=V_max_ind;        %index for current along oxidation potential
        ind1=xsel-avgbg;
        if ind1<1
            ind1=1;
        end
        ind2=xsel+avgbg;
        if ind2>size(subData,2)
            ind2=size(subData,2);
        end
        noise_DA=subData(:,(ind1:ind2));                    %current noise around BG cursor (subtracted out current that is) ie across various potentials
        noise_DA=abs(noise_DA(index_DA,:));             %current noise around BG cursor (time) at oxidation potential
        noise_DA=rms(noise_DA);
        noise=noise_DA;         
        snr=Imax/noise_DA(1);
        abs_I=max(measData(:,xsel));
        ImaxVIndex=find(subData(:,xsel)==Imax);
        ImaxVIndex=ImaxVIndex(1);
        ind1=xsel-5;
        ind2=xsel+5;
        if ind1<1
            ind1=1;
        end
        if ind2>size(subData,2)
            ind2=size(subData,2);
        end
        Ipeak=max(subData(ImaxVIndex,ind1:ind2));
        Eox=subData(ImaxVIndex,xsel);
        ImaxTIndex=find(subData(ImaxVIndex,:)==Imax);
        ImaxTIndex1=find(ImaxTIndex>=xsel,1,'first'); ImaxTIndex2=find(ImaxTIndex<=xsel,1,'first');
        ImaxTIndex=ImaxTIndex(min([ImaxTIndex1 ImaxTIndex2]));
        IhalfTIndex=find(subData(ImaxVIndex,ImaxTIndex:end)<=Ipeak/2);
        if length(IhalfTIndex)<=1
            IhalfTIndex=[1 1];
        end
        IhalfTIndex=ImaxTIndex+IhalfTIndex(1);
        thalfSamples=IhalfTIndex-ImaxTIndex;
        reverseIhalfTIndex=find(subData(ImaxVIndex,1:ImaxTIndex)<=Ipeak);
        if length(reverseIhalfTIndex)<=1
            reverseIhalfTIndex=[1 1];
        end
        reverseIhalfTIndex=max([reverseIhalfTIndex(end) 1]);
        reversethalfSamples=reverseIhalfTIndex-ImaxTIndex;
        tRiseSamples=reversethalfSamples;
        D=cvs{selch(ii)};          %[ianodal icathodal]
        Dproj=parameters.Vc'*D;        %projections of unknown data set onto relevant PC's of training set matrix Ats (Lxn * nxm= Lxm)
        Cu=parameters.F*Dproj;         %concentrartions predicted of each L component (ie DA/pH/background)
        E=D-(parameters.Vc*Dproj);     %residuals
        Q=diag(E'*E)';      %nA^2
        KDA=parameters.K(:,1);
        KPH=parameters.K(:,2);
        KBG=parameters.K(:,3);
        KM=zeros(length(KBG),1);
        try
            %only if movement component exists in template
            KM=parameters.K(:,4);
        catch ME
        end
        R=corr2(KDA(:,1),D);
        RPH=corr2(KPH(:,1),D);
        RBG=corr2(KBG(:,1),D);
        RM=corr2(KM(:,1),D);

        info(selch(ii)).VmaxID=V_max_ind;
    info(selch(ii)).Eox=Eox;
    info(selch(ii)).noise=noise;
    info(selch(ii)).Imax=Imax;
    info(selch(ii)).Ipeak=Ipeak;
    info(selch(ii)).snr=snr;
    info(selch(ii)).Iabs=abs_I;
    info(selch(ii)).IabsBG=mean(abs(measData(:,xsel)));      %average background current
    info(selch(ii)).thalfSamples=thalfSamples;
    info(selch(ii)).tRiseSamples=tRiseSamples;
    info(selch(ii)).Dproj=Dproj;
    info(selch(ii)).E=E;
    info(selch(ii)).Q=Q;
    info(selch(ii)).R=R;
    info(selch(ii)).RPH=RPH;
    info(selch(ii)).RBG=RBG;
    info(selch(ii)).RM=RM;

    else
        info(selch(ii)).VmaxID=1;
        info(selch(ii)).Eox=0;
        info(selch(ii)).noise=0;
        info(selch(ii)).Imax=0;
        info(selch(ii)).Ipeak=0;
        info(selch(ii)).snr=0;
        info(selch(ii)).Iabs=0;
        info(selch(ii)).IabsBG=0;
        info(selch(ii)).thalfSamples=0;
        info(selch(ii)).tRiseSamples=0;
        info(selch(ii)).Dproj=0;
        info(selch(ii)).E=0;
        info(selch(ii)).Q=0;
        info(selch(ii)).R=0;
        info(selch(ii)).RPH=0;
        info(selch(ii)).RBG=0;
        info(selch(ii)).RM=0;
    end
end
                
end