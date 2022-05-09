%Patra map pinout CSC channels to electrode names
PCR_srcdir='A:\mit\injectrode\experiments\fscv\templates\pcr_templates\';
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
if ~ispc
    %chunky dir
    PCR_srcdir=fullfile('home','schwerdt','matlab','analysis','analysis','pcr_templates');
    PCR_srcdir=[PCR_srcdir filesep];
end
PCR_template='patra_clean_titanium_da_ph_bg_movep'; %generated cleaner with average of DA in vitro's and non-noisy BG/movement from patra 04/01/2018 for Patra titanium reference causing redox shifts
PCR_template='patra_clean_titanium_da_ph_bg_movep2'; %4/18/18 made cts for movement correlate with da-like ox peak instead of switch
PCR_template='patra_clean_3pcs_nobg_frommovep2_orig4182018'; %1/28/2019 no bg, to include all 3
PCR_template='patra_clean_pcs_da_ph_m_bg'; %1/29/2019 all 4, order is now da, ph, m, bg
PCR_template='patra_clean_mscaled_da_ph_m_bg.mat';     %2/9/2019, scaled m Cts coeff to match second hump rather than domiannt m hump
PCR_template='pcr_invitrophbgda_comb_forCFMEA_amp.mat';     %11/2017 rodent map