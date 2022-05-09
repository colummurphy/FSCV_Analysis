%Patra map pinout CSC channels to electrode names
PCR_srcdir='A:\mit\injectrode\experiments\fscv\templates\pcr_templates\';
pctype=computer;
ispc=strcmpi(pctype,'pcwin64');
if ~ispc
    PCR_srcdir=fullfile('home','schwerdt','matlab','analysis','analysis','pcr_templates');
    PCR_srcdir=[PCR_srcdir filesep];
end
%PCR_template='patra_clean_titanium_da_ph_bg_movep'; %generated cleaner with average of DA in vitro's and non-noisy BG/movement from patra 04/01/2018 for Patra titanium reference causing redox shifts
%PCR_template='PCR_invitrophbgda_movement_comb';     %includes in vivo mvmtn from chronic 23, 27, 18 for Cleo analysis
PCR_template='cleo_mscaled_da_ph_m_bg';             %3/1/2019, from  PCR_invitrophbgda_movement_comb, used 2nd peak as M scale, now in order of da, ph, m, bg, as with patra

csc_map={
    '3' 'cl6'...
    '4' 'c5'...
    '304'  'cl6-c5'...
    '5' 'c4' ...   
    '504'  'c4-c5'...
     '7'     'c3' ...
     '8'    'c2'...
     '9'    'c1'...
     '805'  'c2-c4'...
     '907'  'c1-c3'...
     '10'   'cl5'...
     '3010' 'cl6-cl5'...
     '11'   'cl4'...
     '11010'    'cl4-cl5'...
     '12'   'cl3'...
     '12011'    'cl3-cl4'...
     '14'   'cl2'...
     '14012'    'cl2-cl3'...
     '803'     'c2-cl6'...
     '15'   'p6'...     
     '16'   'p5'...
     '32'   'p4'...
     '31'   'p3'...
     '15031'    'p6-p3'...
     '30'   'p2'...
     '28'   'm1'...
     '27'   'm3'...
     '26'   'm4'...
     '25'   'm2'...
   '33'    'eyed'   ...
    '34'    'eyex'  ...
    '35'    'eyey'  ...    
    };

event_codes={
    '4'     'display_fix' ...
    '5'     'start_fix' ...
    '6'     'break_fix' ...
    '10'    'display_target' ...
    '11'    'start_target'  ...
    '12'    'break_target'  ...
    '14'    'error' ...
    '29'    'left_condition'    ...
    '30'    'right_condition'   ...
    '45'    'reward_big'    ...
    '46'    'reward_small'  ...
    };
%07/27/2018 - found out left-condition/right_condition not proper fix spot
%on, need to use "display_fix" instead

%error_trial means fixation break, target fixation break, or something else
%ie # error > fix or target breaks
%but error < fix + target breaks + no enter trials (37)
%# correct trials (ie. id = 13) = = # 45 + #46
%however if we calculate total  trials based on #3 start trial cue) - successful
%trials(#13), get # error trials (14)
%If use #4 (fixspot on) as indicating total trials then this is equal to
%#37 (no enter) + #13 (correct) + #14 (error)
%29 & 30 presented at same as 4

