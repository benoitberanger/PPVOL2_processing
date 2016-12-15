close all
clear all
fclose('all')
clc

%%

chemin={[pwd filesep 'PRISMA_PPVOL2_processing']}
% suj = get_subdir_regex(chemin)
%ou
suj = get_subdir_regex(chemin,'Sujet'); %to get all subdir that start with 2
%to see the content
char(suj)


%functional and anatomic subdir
par.dfonc_reg='(\d|Meridian)$';
par.dfonc_reg_oposit_phase = 'blip$'; % [AZER]
par.danat_reg='t1_mpr';

%for the preprocessing : Volume selecytion
par.anat_file_reg  = '^s.*nii'; %le nom generique du volume pour l'anat
par.file_reg  = '^f.*nii'; %le nom generique du volume pour les fonctionel

% par.TR = 2.1; %for slice timing
par.run=1;par.display=0;


%%


%anat segment
anat = get_subdir_regex(suj,par.danat_reg)
fanat = get_subdir_regex_files(anat,par.anat_file_reg,1)

par.GM   = [1 0 1 0]; % Unmodulated / modulated / native_space dartel / import
par.WM   = [1 0 1 0];
j = job_do_segment(fanat,par)

%apply normalize on anat
fy = get_subdir_regex_files(anat,'^y',1)
fanat = get_subdir_regex_files(anat,'^ms',1)
j=job_apply_normalize(fy,fanat,par)


%anat brain extract

ff=get_subdir_regex_files(anat,'^c[123]',3);
fo=addsufixtofilenames(anat,'/mask_brain');
do_fsl_add(ff,fo)
fm=get_subdir_regex_files(anat,'^mask_b',1); fanat=get_subdir_regex_files(anat,'^s.*nii',1);
fo = addprefixtofilenames(fanat,'brain_');
do_fsl_mult(concat_cell(fm,fanat),fo);



%% get subdir

dfonc = get_subdir_regex_multi(suj,par.dfonc_reg)
dfonc_op = get_subdir_regex_multi(suj,par.dfonc_reg_oposit_phase)
dfoncall = get_subdir_regex_multi(suj,{par.dfonc_reg,par.dfonc_reg_oposit_phase })
anat = get_subdir_regex_one(suj,par.danat_reg) %should be no warning



%%

%realign and reslice
fprintf('realign and reslice : dfonc \n')
par.file_reg = '^f.*nii'; par.type = 'estimate_and_reslice';
j = job_realign(dfonc,par)

%realign and reslice opposite phase
fprintf('realign and reslice opposite phase: dfonc_op \n')
par.file_reg = '^f.*nii'; par.type = 'estimate_and_reslice';
j = job_realign(dfonc_op,par)

%topup and unwarp
fprintf('topup and unwarp \n')
par.file_reg = {'^rf.*nii'}; par.sge=0;
do_topup_unwarp_4D(dfoncall,par)

%coregister mean fonc on brain_anat
fprintf('coregister mean fonc on brain_anat \n')
fanat = get_subdir_regex_files(anat,'^brain_',1)

par.type = 'estimate';
for nbs=1:length(suj)
    fmean(nbs) = get_subdir_regex_files(dfonc{nbs}(1),'^utmeanf');
end

fo = get_subdir_regex_files(dfonc,'^utrf.*nii',1)
j=job_coregister(fmean,fanat,fo,par)

%apply normalize
fy = get_subdir_regex_files(anat,'^y',1)
j=job_apply_normalize(fy,fo,par)

%smooth the data
ffonc = get_subdir_regex_files(dfonc,'^wutrf')
par.smooth = [6 6 6];
j=job_smooth(ffonc,par);


%% first level

onsetdir = get_subdir_regex(suj,'onset')

sta=r_mkdir(suj,'stat')
st_localizer =r_mkdir(sta,'Localizer')
st_localizer_1run =r_mkdir(sta,'Localizer_1run')
st_meridian =r_mkdir(sta,'Meridian')

fons_localizer_1run = get_subdir_regex_files(onsetdir,'^loca',1);
fons_localizer = {repmat(fons_localizer_1run{1},[2 1])};
fons_meridian = get_subdir_regex_files(onsetdir,'^meridian',1);

par.file_reg = '^sw.*nii';
dfoncdir_localizer = get_subdir_regex_multi(suj,'Localizer_\d$')
dfoncdir_localizer_1run = get_subdir_regex_multi(suj,'Localizer_2$')
dfoncdir_meridian = get_subdir_regex_multi(suj,'Meridian$')

par.TR=1.000;
par.delete_previous=1;


%% Specify model

par.run = 1
par.display = 0

j_specify_localizer = job_first_level12(dfoncdir_localizer,st_localizer,fons_localizer,par)
j_specify_localizer_1run = job_first_level12(dfoncdir_localizer_1run,st_localizer_1run,fons_localizer_1run,par)
j_specify_meridian = job_first_level12(dfoncdir_meridian,st_meridian,fons_meridian,par)


%% Estimate model

fspm_localizer = get_subdir_regex_files(st_localizer,'SPM',1)
j_estimate_localizer = job_first_level12_estimate(fspm_localizer,par)

fspm_localizer_1run = get_subdir_regex_files(st_localizer_1run,'SPM',1)
j_estimate_localizer_1run = job_first_level12_estimate(fspm_localizer_1run,par)

fspm_meridian = get_subdir_regex_files(st_meridian,'SPM',1)
j_estimate_meridian = job_first_level12_estimate(fspm_meridian,par)


%% Define contrast

contrast_localizer.names = {
    'stimleft'
    'blank'
    'stimright'
    'stimleft  - blank'
    'stimright - blank'
    'stimleft  - stimright'
    'stimright - stimleft'
    'stimleft  + stimright - 2*blanck'
    }';
contrast_localizer.values = {
    [ 1 0 0]
    [ 0 1 0]
    [ 0 0 1]
    [1 -1 0]
    [0 -1 1]
    [1 0 -1]
    [-1 0 1]
    [1 -2 1]
    }';
contrast_localizer.types = cat(1,repmat({'T'},[1 length(contrast_localizer.names)]));

contrast_meridian.names = {
    'horz'
    'blank'
    'vert'
    'horz - blank'
    'vert - blank'
    'horz - vert'
    'vert - horz'
    'horz + vert - 2*blanck'
    }';
contrast_meridian.values = {
    [ 1 0 0]
    [ 0 1 0]
    [ 0 0 1]
    [1 -1 0]
    [0 -1 1]
    [1 0 -1]
    [-1 0 1]
    [1 -2 1]
    }';
contrast_meridian.types = cat(1,repmat({'T'},[1 length(contrast_meridian.names)]));

par.delete_previous=1;


%% Estimate contrast

par.run = 1;
par.display = 0;

j_contrast_localizer_rep = job_first_level12_contrast_rep(fspm_localizer,contrast_localizer,par)
j_contrast_localizer_1_run = job_first_level12_contrast(fspm_localizer_1run,contrast_localizer,par)
j_contrast_meridian = job_first_level12_contrast(fspm_meridian,contrast_meridian,par)


