function [parameters,csc_map,event_codes]=CM_getparams(subjectid,recordid,ncschannels,varargin)
%get plotParam and parameters variables as well as maps for subject &
%recording setting
csc_map={};
event_codes={};

%default fscv params
anodal_lim=1.3;
cathodal_lim=-0.4;
vlim=anodal_lim-cathodal_lim;
Vox=0.6;
sampling_freq=10;   %10 hz sampling freq 
offset=0;
argnum=1;
sessionnum=[];
while argnum<=length(varargin)
    switch varargin{argnum}
        case 'sessnum'
            argnum=argnum+1;
            sessionnum=varargin{argnum};
    end
    argnum=argnum+1;
end
if isempty(sessionnum)
    %delete variable, so does not get fed into patra_map_biplar to det
        %which csc_map to use
        clear sessionnum;
end
%get csc_map, PCR file dir, event_codes for subject
switch subjectid
    case 'patra'
        patra_map
    case 'patrabipolar'
        CM_patra_map_bipolar
    case 'patrabichunky'
        patra_map_bipolar_chunky
    case 'cleo'
        cleo_map_bipolar
    case 'cleo16'
        cleo_map_bipolar16
    case 'cleo17'
        cleo_map_bipolar17
    case 'cleo25'
        cleo_map_bipolar25
    case 'cleo14'
        cleo_map_bipolar14   
    case 'cleo18'
        cleo_map_bipolar18
    case 'cleo19b'
        cleo_map_bipolar19b
    case 'cleo19'
        cleo_map_bipolar19
    case 'cleo20'
        cleo_map_bipolar20
    case 'cleo22'
        cleo_map_bipolar22
    case 'cleo24'
        cleo_map_bipolar24
    case 'cleo28'
        cleo_map_bipolar28    
    case 'cleo15'
        cleo_map_bipolar15
    case 'cleo23'
        cleo_map_bipolar23    
    case 'cleo22b'
        cleo_map_bipolar22b
    case 'cleo19b'
        cleo_map_bipolar19b
    case 'cfmea'
        rodent_map
    otherwise
        %no ncs channels, maybe da recording only
end

if exist('recordid')
    if ~isempty(recordid)
%get settings for recording session, for artifact removal
%load parameters variables
switch recordid
    case 'default'
        settings_default
    case 'noisycleo'
        settings_cleo_noisy
    case 'cleo'
        settings_cleo_noisy
        if exist(['settings_' subjectid])>0
            run(['settings_' subjectid]);
        end
    case 'cleo14'
        settings_cleo14
    case 'rodent'
        settings_rodent
        parameters.differentialChannel=diffChannel;
        parameters.badChannels =badChannels;        
        parameters.resort=resort;
        parameters.IVconv=1e9/4.99e6; 
end
end
else
    settings_default
end

if exist('ncschannels')
    if ~isempty(ncschannels)
    %targeted ncs to display
    %get cscid's in the order they are given in ncschannels
    for ii=1:length(ncschannels)
        cscid=find(ismember(csc_map,ncschannels{ii}))-1;
        parameters.NCSchannels(ii)=str2num(csc_map{cscid});
    end
    else
        parameters.NCSchannels=[];
    end
end
    
parameters.templatesrcdir=PCR_srcdir;   %from subjectid loaded map file
parameters.PCR_template=PCR_template; 

%default fscv paramters based on templates
load([parameters.templatesrcdir parameters.PCR_template])      %loads Ats_mat & Cts_mat
[index1,~]=find(Cts_mat~=0); index_da=find(index1==1); index_ph=find(index1==2);
Cts_mat_DAPH=Cts_mat(1:2,1:index_ph(end)); CV_mat_DA=Ats_mat(:,index_da); CV_mat_pH=Ats_mat(:,index_ph);

parameters.DAcvs=CV_mat_DA;
parameters.Cts=Cts_mat;
parameters.CV=Ats_mat;

if strcmp(recordid,'rodent')
    %scale pca templates
    [index1,index2]=find(Cts_mat~=0); index_da=find(index1==1); index_ph=find(index1==2);
Cts_mat_DA=Cts_mat(1,index_da); Cts_mat_pH=Cts_mat(2,index_ph); Cts_mat_DAPH=Cts_mat(1:2,1:index_ph(end));

resampled=resampleCV(Ats_mat,length(Ats_mat),parameters.samplesperscan);       %resample templates to match recorded length
Ats_mat=resampled;
CV_mat_DA=Ats_mat(:,index_da);
CV_mat_pH=Ats_mat(:,index_ph);
%calculate PCR parameters from templates
[Vc F Aproj Qa]=generatePCR(Ats_mat,Cts_mat);
[VcDAPH FDAPH AprojDAPH QaDAPH]=generatePCR([CV_mat_DA CV_mat_pH],Cts_mat_DAPH);
QaDAPH=Qa*5;            %just for CFMEA
K=pinv(FDAPH*transpose(VcDAPH));        %na/uM, ideal CV's for each analyte
parameters.DAcvs=CV_mat_DA;
parameters.Cts=Cts_mat_DAPH;
parameters.CV=[CV_mat_DA CV_mat_pH];


end

%calculate PCR parameters from templates
[Vc, F, Aproj, Qa]=generatePCR(Ats_mat,Cts_mat);
[VcDAPH, FDAPH, AprojDAPH, QaDAPH]=generatePCR([CV_mat_DA CV_mat_pH],Cts_mat_DAPH);
K=pinv(F*transpose(Vc));        %na/uM, ideal CV's for each analyte

parameters.Vc=Vc;  parameters.F=F;    parameters.Aproj=Aproj;    
parameters.Qa=Qa;  parameters.QaDAPH=QaDAPH;  parameters.K=K;

%default info about fscv data based on templates
sizeData=size(Ats_mat); 
samplesperscan=sizeData(1); 
stepsizeCV=((vlim-vlim/samplesperscan)*2/samplesperscan);
Vrange=(cathodal_lim:stepsizeCV:(anodal_lim-vlim*2/samplesperscan))+offset;
Vrange_cathodal=(anodal_lim:-stepsizeCV:(cathodal_lim+vlim*2/samplesperscan))+offset;
if (length(Vrange)+length(Vrange_cathodal))<samplesperscan
    Vrange=[Vrange 1.3+offset];
    Vrange_cathodal=(anodal_lim:-stepsizeCV:(cathodal_lim+vlim*2/samplesperscan))+offset;
end
Voxid=find(Vrange<=Vox+0.01); Voxid=Voxid(end);
Voxid=length(Vrange_cathodal)+length(Vrange)-Voxid;

parameters.Voxid=Voxid;
parameters.Vrange=Vrange;
parameters.Vrange_cathodal=Vrange_cathodal;



parameters.LFPterm='.ncs';
parameters.includeMovementThres=1;          %include this threshold to remove singals
parameters.samplesperscan=samplesperscan;
parameters.samplerate=sampling_freq;
end