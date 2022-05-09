function pct = getpct(rawdata,subdata,bg, parameters,varargin)
%04/08/2018
%similar to getpctdirect, but includes subfunction to remove hf glitches
smoothwidth=4;        %# of samples +/- to interpolate with nan around glitches/artifacts
hfwidth=2;      %# samplse +/- hf glitch (not as longa s movement artifacts so usu smaller
hfwidth=3;      %# 07/27/2018 changed to be same as in autotrialdir 
%was set to 5 pre-04/18/2018, then set to 2 04/18/2018 to look more closely
%at da changes around reward time
BGavgi=bg;
saturationthres=parameters.saturationLevel;     %add 07/04/2018
saturationthres=[];
glitchThresT=parameters.glitchThresT;
glitchThres=parameters.glitchThres;
Rthres=parameters.RThres;
RthresOut=parameters.RThresOut;
QaDAPH=parameters.Qa;
CV_mat_DA=parameters.DAcvs;
ats=parameters.CV;
cts=parameters.Cts;
[aa,oxid]=max(parameters.K(:,1));   %DA oxidation current

mstdlim=2; %1/2019
numpcs=3;       %da,ph,bg
if isfield(parameters,'numpcs')
    %how many components to use, just ph & da, or inlcude bg/move?
   numpcs=parameters.numpcs; 
end
includem=parameters.includepcam;        %include movement pca component or just use for artiact removal;
argnum=1;
rembgph=0;
rthresp=0.95;    %rthresout for ph
rthresb=0.95;  %rthresout for bg
rthresb=0.85;  %rthresout for bg

satwidth=[-5 30];    %atleast 2 seconds needed for saturated current to restablilize
bgref=0;
KDA=[];
KBG=[];
KM=[];
KPH=[];
Qa=0;
while argnum <= length(varargin)
    switch varargin{argnum}
        case 'removebgph'
            %make correlated bg/ph signals also glitches like movement
            rembgph=1;
       case 'nanwidth'
            %width of nanned signal padding
            argnum=argnum+1;
            smoothwidth=varargin{argnum};
       case 'glitchwidth'
            %width of nanned signal padding
            argnum=argnum+1;
            if ~isempty(varargin{argnum})
            hfwidth=varargin{argnum};
            end
        case 'satwidth'
            %width of large saturation current -/+ samples
            argnum=argnum+1;
            satwidth=varargin{argnum};
       case 'rthresp'
            %ph component correlation coefficient threshold exclusion
            argnum=argnum+1;
            rthresp=varargin{argnum};
        case 'rthresb'
            %bg component correlation coefficient threshold exclusion
            argnum=argnum+1;
            rthresb=varargin{argnum};
    end
    argnum = argnum + 1;
end

isub=subdata;
rawexists=0;
if ~isempty(rawdata)
    bgref=rawdata(:,6:15);
    bgref=mean(bgref,2);         %mean ## sample points BG for subtraction
    refmatrix=repmat(bgref,1,size(rawdata,2));   %tile BGavg vector for entire meas matrix span

    rawdata=rawdata-refmatrix;      %raw data is subtracted data
    rawexists=1;
end

if isempty(saturationthres)
    saturationthres=mean(abs(bgref)).*.015;
end

ImaxBG=max(BGavgi);
IDAproj=[];
pHplot=[];
IDsaboveRout=[];        %idx of bad signals (ie. correlated to movement template or glitch)
  Dproj=[];
  F=[];
KscaledM=[];        %movement templates repeated to find correlated data (ie. artifacts)
Mplot=[];           %also calculate M current component if we were to include all 4 signal categories
Vc=[];
%if do not want to use movement component, but exists in loaded pca
%templates, remove from templates cts and ats variables
%still calculated KscaledM for artifact removal
%{
if includem==0 && numpcs==4
    [Vc, F, ~, ~]=generatePCR(ats,cts);
    Dproj=Vc'*isub;
    Cu=F*Dproj;         %current contributions predicted of each component 
    K=pinv(F*transpose(Vc));        %na/uM, ideal CV's for each analyte
    Mplot=[];
    KM=zeros(length(Vc),1);
    KscaledM=KM*Dproj(3,:);
    try
        Mplot=Cu(4,:);
        KM=K(:,4);
        KscaledM=KM*Dproj(4,:);     %ideal CV template to correlate to
    catch
    end
    KscaledM=KM*Dproj(4,:);     %ideal CV template to correlate to
    
    %now reduce templates provided to without movement templates from 
    %predict concentrations of main non-overlapping components
    midx=find(cts(4,:)==0);
    midx=midx(end);
    cts=cts(1:3,1:midx);
    ats=ats(:,1:midx);
elseif includem==1 
    [Vc, F, ~, ~]=generatePCR(ats,cts);
    Dproj=Vc'*isub;
    Cu=F*Dproj;         %current contributions predicted of each component 
    K=pinv(F*transpose(Vc));        %na/uM, ideal CV's for each analyte
    Mplot=[];
    KM=zeros(length(Vc),1);
    KscaledM=KM*Dproj(3,:);    
end
%}

%get current contributions of each component
if numpcs==4
    %{
if size(cts,1)>3
        %reduce impact of bg (cts(4,:)), DOES NOT MAKE DIFFERENCE
    midx=find(cts(4,:)~=0);
    midx=midx(1);
    cts4=cts(4,midx:end);
    cts4=cts4.*4;
    cts(4,midx:end)=cts4;
end
    %}
  [Vc, F, Aproj, Qa]=generatePCR(ats,cts);
Dproj=Vc'*isub;
Cu=F*Dproj;         %concentrartions predicted of each L component (ie DA/pH/background)
E=isub-(Vc*Dproj);     %residuals
Q=diag(E'*E)';      %nA^2
IDAproj=Cu(1,:);            %current da projection
pHplot=Cu(2,:);
Mplot=Cu(3,:);          %1/2019, changed from BGplot
BGplot=Cu(4,:);
K=pinv(F*transpose(Vc));        %na/uM, ideal CV's for each analyte
KDA=K(:,1);
KPH=K(:,2);
%KBG=K(:,3);
KM=K(:,3);          %3RD COMPONENT NOW MOVEMENT, NOT BG 01/2019
KscaledM=KM*Dproj(3,:);     %ideal CV template to correlate to
KBG=K(:,4);
  
elseif numpcs==3
%REMOVE LAST COMPONENT< BG, 1 /2019 after getting BG component for
%reference
if size(cts,1)==4
[Vc, F, Aproj, Qa]=generatePCR(ats,cts);
Dproj=Vc'*isub;
Cu=F*Dproj;         %concentrartions predicted of each L component (ie DA/pH/background)
E=isub-(Vc*Dproj);     %residuals
Q=diag(E'*E)';      %nA^2
BGplot=Cu(4,:);          %1/2019, changed from BGplot
K=pinv(F*transpose(Vc));        %na/uM, ideal CV's for each analyte
KBG=K(:,4);
end
if size(cts,1)>3
    midx=find(cts(4,:)==0);
    midx=midx(end);
    cts=cts(1:3,1:midx);
    ats=ats(:,1:midx);
end
[Vc, F, Aproj, Qa]=generatePCR(ats,cts);
Dproj=Vc'*isub;
Cu=F*Dproj;         %concentrartions predicted of each L component (ie DA/pH/background)
E=isub-(Vc*Dproj);     %residuals
Q=diag(E'*E)';      %nA^2
IDAproj=Cu(1,:);            %current da projection
pHplot=Cu(2,:);
Mplot=Cu(3,:);          %1/2019, changed from BGplot
K=pinv(F*transpose(Vc));        %na/uM, ideal CV's for each analyte
KDA=K(:,1);
KPH=K(:,2);
KM=K(:,3);          %3RD COMPONENT NOW MOVEMENT, NOT BG 01/2019, except for rodent
KscaledM=KM*Dproj(3,:);     %ideal CV template to correlate to
%    KBG=nan(length(KM),1);
%BGplot=nan(1,size(Cu,2));
elseif numpcs==2
    %reduce to 2 components only
    if size(cts,1)>2
        midx=find(cts(3,:)~=0);
        midx=midx(1)-1;
        cts=cts(1:2,1:midx);
        ats=ats(:,1:midx);
    end
       
    [Vc, F, Aproj, Qa]=generatePCR(ats,cts);
    Dproj=Vc'*isub;
    Cu=F*Dproj;         %concentrartions predicted of each L component (ie DA/pH/background)
    E=isub-(Vc*Dproj);     %residuals
    Q=diag(E'*E)';      %nA^2
    IDAproj=Cu(1,:);            %current da projection
    pHplot=Cu(2,:);
    BGplot=nan(1,size(Cu,2));     %make bg zeros
    K=pinv(F*transpose(Vc));        %na/uM, ideal CV's for each analyte
    KDA=K(:,1);
    KPH=K(:,2);
    KBG=nan(size(K,1),1);
end
VcPH=Vc(:,2);
DprojPH=VcPH'*isub;
EPH=isub-(VcPH*DprojPH);       %residuals without pH only 
VcM=[];
DprojM=[];
EM=[];
VcPHM=[];
DprojPHM=[];
EPHM=[];
if size(Vc,2)>2
VcM=Vc(:,3);
DprojM=VcM'*isub;
EM=isub-(VcM*DprojM);       %residuals without M only 
VcPHM=Vc(:,2:3);
DprojPHM=VcPHM'*isub;
EPHM=isub-(VcPHM*DprojPHM);       %residuals without pH & M
end
pct.EPH=EPH;
pct.EM=EM;
pct.EPHM=EPHM;
pct.F=F;
pct.Dproj=Dproj;
pct.Vc=Vc;
%find cv's that have hf glitch within them
if rawexists==1
    %get differential along cvs (I vs V that is) 
    maxamps=max(abs(rawdata),[],1);
    diffrawcvs=diff(rawdata,1,1);   
    maxdiff=max(abs(diffrawcvs),[],1)-min(max(abs(diffrawcvs),[],1));   %subtraction 07/22/2018
    tbb=find(maxdiff>glitchThres);
    idxglitch=tbb;
    idxglitch=unique(idxglitch);  %remove duplicates
    %get raw signal along time for v-selected
    diffIrawt=diff(rawdata,2,1); %take out hf glitches along time by looking at unfiltered data
    maxdifft=max(abs(diffIrawt),[],1)-mean(max(abs(diffIrawt),[],1));   %add 07/20/2018
    %get differential of this signal to detect HF glitches
    rawglitch=find(maxdifft>glitchThres);
    if ~isempty(idxglitch)
        idxglitch=unique([idxglitch rawglitch]);
    else
        idxglitch=rawglitch;
    end
    
    %for each detected sample add +/- smoothwidth samples around
    for ii=1:length(idxglitch)
        idsaroundglitch=idxglitch(ii)-hfwidth:idxglitch(ii)+hfwidth;
        %very large current changes (ie large movements)
        %want to take out signal more widely, more of forward current
        if abs(mean(rawdata(:,idxglitch(ii))))>saturationthres
            idsaroundglitch=idxglitch(ii)+satwidth(1):idxglitch(ii)+satwidth(2);           
        end
        if ii==1
            IDsaboveRout=idsaroundglitch;
        else
            IDsaboveRout=[IDsaboveRout idsaroundglitch];
        end
    end
    IDsaboveRout=sort(unique(IDsaboveRout));    
    IDsaboveRout=IDsaboveRout(IDsaboveRout>0 & ...
        IDsaboveRout<size(rawdata,2)); %only keep ids in original window

end

dsdata=[];
%smooth data & downsample data to find vectors correlated
%to movement artifact template
   % dsdata=smoothts(isub,'g',smoothwidth);  Matlab 2013 and earlier
if exist('smoothdata')
    %if is function
dsdata=smoothdata(isub,2,'gaussian',smoothwidth);  %maybe not needed, anyways process entire dataset below
else
    %use <2013 version of smooth data, need to check more carefully
    dsdata=smoothts(isub,'g',smoothwidth); 
end

%get correlation coefficients between data all time and movementtemplate
rr=[];
rm=[];
if ~isempty(KscaledM)
rr=abs(corr(dsdata,KscaledM(:,1)));
rm=find(rr>RthresOut);
end
if isempty(IDsaboveRout)
    IDsaboveRout=rm';
else
    if ~isempty(rm)
    IDsaboveRout=[IDsaboveRout rm'];        %idx of bad signals (ie. correlated to movement template or glitch)
    end
end

if rembgph==1
    rb=abs(corr(isub,KBG));
    routb=find(rb>rthresb);
    rp=abs(corr(isub,KPH));
    routp=find(rp>rthresp);
    IDsaboveRout=[IDsaboveRout routb' routp'];
end

% ITplot=ITplot./(0.0341*ImaxBG)*1000;    %normalize to conc
ITplot=IDAproj./(0.0282*ImaxBG*1.24)*1000;    %normalize to conc
DAproj=ITplot;

%DEFINE LIMIT HERE BASED on first 50% of lower Q in data 01/2019
sortedQ=sort(Q);
thresQ=nanmean(sortedQ(1:round(length(sortedQ)/2)));
QaDAPH=thresQ*7;        %2/9/2019, previously used default QaDAPH, now changed
if isfield(parameters,'traditionalQthres')
    if parameters.traditionalQthres==1
    QaDAPH=Qa*.25;
    end
end
%QaDAPH=2.5;                 %2/21/2019
IDsaboveQ=find(Q>(QaDAPH));     %Qa forDA & PH
Mthres=nanstd([pHplot IDAproj BGplot])*mstdlim;     %1/20/2019 Thres for value of M component (not just corr looking for);
IDsaboveMthres=find(Mplot>Mthres);
IDsaboveRout=[IDsaboveRout find(abs(diff(Mplot))>(glitchThresT))];     %if above Rthres then must be artifact
ITplot(IDsaboveQ)=nan;     ITplot(IDsaboveRout)=nan; 
pct.pHplot=pHplot;         %current's before nan
pHplot(IDsaboveQ)=nan;      pHplot(IDsaboveRout)=nan;
ITplot(IDsaboveMthres)=nan;
pct.Mplot=Mplot;          %current's before nan
pct.BGplot=BGplot;        %current's before nan
pct.DAproj=DAproj;          %current's before nan
pct.DAiso=ITplot;           %current's after  nan
pct.rawIT=rawdata;
pct.Qplot=Q;
pct.Iox=nanmean(isub(oxid-2:oxid+2,:),1);


end
