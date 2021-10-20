function expt = rivindiff_experiment

if ~ispc
    %     expt.session_dir = '~/google_drive/Projects/Experiments/Rivalry_IndivDiff/Experiment/';
    expt.session_dir = '~/OneDrive/Projects/Experiments/Rivalry_IndivDiff/Experiment/';
else
    expt.session_dir = 'C:\OneDrive - UC Davis\Projects\Experiments\Rivalry_IndivDiff\Experiment\';
    %     expt.session_dir = 'C:\OneDrive\Projects\Experiments\Rivalry_IndivDiff\Experiment\';
end
expt.eeg_datadir = 'Data/eeglab';
expt.eeg_epochdir = 'Epochs';

subj_codes = [3577 3578 3581 3582 3583 3584 3585 3586 3587 3588 ...
    3589 3591 3592 3593 3597 3598 3600 3601 3602 3605 ...
    3607 3608 3609 3610 3613 3614 3615 3616 ...
    3618 3619 3620 3621 3623 3624 3626 3627 3628 3629 ...
    3630 3632 3633 3634 3635 3636 3639 3641 3642 3643 ...
    3644 3646 3647 3648 3650 3651 3652 3653 3654 3656 ...
    3657 3659 3661 3663 3664 3665 3666 3668 3669 3670 ...
    3672 3673 3674 3675 3676 3677 3678 3679 3680 3683 ...
    3689 3690 3692 3693 3694 3695 3697];

subj_fam_millness = [3 6 0 2 0 0 0 0 2 0 ...
    0 2 0 0 0 0 0 0 6 0 ...
    0 0 0 0 0 0 0 0 ...
    0 2 0 0 0 0 0 0 0 2 ... % 3618 3619 3620 3621 3623 3624 3626 3627 3628 3629
    0 0 0 5 0 9 0 0 0 0 ... % 3630 3632 3633 3634 3635 3636 3639 3641 3642 3643 ..
    0 0 3 0 0 0 0 0 0 2 ... % 3644 3646 3647 3648 3650 3651 3652 3653 3654 3656 ...
    0 10 0 0 0 0 0 0 0 0 ... % 3657 3659 3661 3663 3664 3665 3666 3668 3669 3670 ...
    0 0 0 0 0 2 0 0 0 2 ... % 3672 3673 3674 3675 3676 3677 3678 3679 3680 3683 ...
    0 0 0 0 0 0 0]; % 3689 3690 3692 3693 3694 3695 3697];
% 0-none, 2-depression, 3-anxiety, 4-depress+anxiety, 5-OCD, 6-ADD, 7-Epilepsy, 8-PTSD, 9-Bipolar, 10-asperger

subj_millness = [0 2 0 0 0 0 0 2 0 5 ...
    0 0 0 0 0 0 0 0 6 0 ...
    0 0 0 0 0 0 3 0 ...
    0 2 0 0 0 0 0 0 0 3 ...
    0 0 3 8 0 0 0 0 0 3 ...
    0 0 3 0 0 0 2 0 0 0 ...
    3 0 0 0 0 0 0 4 0 0 ...
    0 0 0 0 0 2 0 0 0 2 ...
    3 0 0 0 0 0 0]; % 0-none, 2-depression, 3-anxiety, 4-depress+anxiety, 5-OCD, 6-ADD, 7-Epilepsy, 8-PTSD

subj_age = [20    23    40    19    23    21    18    21    21    23,...
    23    19    19    19    21    19    18    20    23    18,...
    21    19    18    20    18    23    18    18,...
    19    18    18    18    18    23    21    19    21    18,...
    22    31    20    27    19    20    20    20    20    20,...
    22    21    20    20    19    18    18    21    21    19,...
    22    23    20    19    22    27    23    19    19    20,...
    19    19    20    19    21    19    18    21    28    19,...
    32    20    23    20    19    19    25];

subj_sex = ['M','F', 'F', 'F', 'F', 'F', 'M', 'F', 'M', 'F', ...
    'M','F','F','M','F','M','F','M','F','M',...
    'F','F','F','F','F','F','F','F',...
    'F','F','F','F','M','F','F','F','F','F',...
    'F','M','F','F','M','F','F','F','M','F',...
    'F','F','F','F','F','F','F','F','F','F',...
    'F','M','M','F','F','F','M','F','M','F',...
    'F','M','M','M','M','F','F','F','F','F',...
    'M','F','F','F','F','F','F'];

invalid_riv = [{[]; []}, {[]; []}, {[]; [3]}, {[2, 7]; []}, {[]; []}, ... % mixed percept for most of the run
    {[1:5]; []}, {[]; []}, {[]; []}, {[]; []}, {[]; []},...
    {[]; [1 2]}, {[7 8]; []}, {[]; []}, {[]; []}, {[]; []},...
    {[]; []}, {[]; []}, {[]; []}, {[]; []}, {[]; []},...
    {[]; []}, {[]; []}, {[]; []}, {[]; []}, ....
    {[10]; []}, {[]; []}, {[]; []}, {[]; []}, ...
    {[]; []}, {[]; []}, {[]; []}, {[]; []}, {[]; []}, ...
    {[]; []}, {[]; []}, {[1]; []}, {[]; []}, {[]; []}, ...
    {1:9; []}, {[]; []}, {[]; []}, {[]; []}, {[]; []}, ...
    {[]; []}, {[8 12]; []}, {[]; []}, {[]; []}, {[]; []}, ...
    {[]; []}, {[]; []}, {[]; []}, {[]; []}, {[]; []}, ...
    {[]; []}, {[]; []}, {[1 2]; []}, {[]; []}, {[]; []}, ...
    {[]; []}, {[]; []}, {[]; []}, {[]; []}, {[]; []}, ...
    {[]; []}, {[]; []}, {[]; []}, {[]; []}, {[]; []}, ...
    {[]; []}, {[]; []}, {[]; []}, {[]; []}, {[]; []}, ...
    {[]; []}, {[]; []}, {[]; []}, {[]; []}, {[]; []}, ...
    {[]; []}, {[]; []}, {[]; []}, {[]; []}, {[]; []}, ...
    {[]; []}, {[]; []}];
invalid_sfm = [{[]; []}, {[]; []}, {[]; []}, {[]; []}, {[]; []},...%3577-3583 mixed percept for most of the run
    {[3,4]; []}, {[]; []}, {[]; []}, {[]; [3]}, {[]; []},...% 3584-3588
    {[]; []}, {[]; []}, {[]; []}, {[]; []}, {[]; []},    ...% 3589-3597
    {[]; []}, {[3]; []}, {[]; []}, {[]; []}, {[]; []},   ...% 3598-3605
    {[]; []}, {[]; []}, {[]; []}, {[]; []},   ... % 3607-3000
    {[]; []}, {[]; []}, {[]; []}, {[]; []},   ... % 3001-3616
    {[]; []}, {[]; []}, {[]; []}, {[3]; []}, {[]; []},   ...% 3618-3623
    {[]; []}, {[]; []}, {[]; []}, {[1]; []}, {[]; []},   ...% 3624-3629
    {1:6; []}, {[]; []}, {[]; []}, {[]; []}, {[]; []},   ...% 3630-3635
    {[]; []}, {[]; []}, {[]; []}, {[]; []}, {[]; []},   ... % 3636-3643
    {[]; []}, {[]; []}, {[]; []}, {[]; []}, {[]; []},   ... % 3644-3650
    {[]; []}, {[]; []}, {[]; []}, {[]; []}, {[]; []},   ... % 3651-3656
    {[]; []}, {[]; []}, {[]; []}, {[]; []}, {[]; []},   ... % 3657-3664
    {[]; []}, {[]; []}, {[]; []}, {[]; []}, {[]; []},   ... % 3665-3670
    {[]; []}, {[]; []}, {[]; []}, {[]; []}, {[]; []},   ... % 3672-3676
    {[]; []}, {[]; []}, {[]; []}, {[]; []}, {[]; []},   ... % 3677-3683
    {[]; []}, {[]; []}, {[]; []}, {[]; []}, {[]; []},   ... % 3689-3697
    {[]; []}, {[1 7]; []}];

reverse_rivbutton = false(1, numel(subj_codes));
rev_subjs = [3586, 3605, 3613, 3623, 3627, 3633, 3639, ...
    3663, 3672, 3690];
rev_subjs = [3586, 3597, 3605, 3613, 3623, 3627, 3633, 3639, ...
    3663, 3672, 3690];
for rs = 1:numel(rev_subjs)
    reverse_rivbutton(subj_codes==rev_subjs(rs)) = true;
end

bad_fix_runs = {[3597 3], [3630 1], [3635 1], [3648 1], [3666 3], [3695 4], [3689 3], [3679 3], ...
    [3636 1], [3697 4]};
badfixruns = cell(1, numel(subj_codes));
for rs = 1:numel(bad_fix_runs)
    badfixruns{subj_codes==bad_fix_runs{rs}(1)} = bad_fix_runs{rs}(2:end);
end

bad_rest_runs = {[3591 6], [3602 10], [3605 6], [3615 10], [3626 10], [3627 6], ...
    [3636 6], [3634 10], [3643 10], [3665 8], [3673 8], [3697 10]};
bad_rest_runs = {};
badrestruns = cell(1, numel(subj_codes));
for rs = 1:numel(bad_rest_runs)
    badrestruns{subj_codes==bad_rest_runs{rs}(1)} = bad_rest_runs{rs}(2:end);
end


subj_codes = num2strcell(subj_codes);
subj_valid = [1 ones(1, 9),...
    ones(1, 10), ...
    ones(1, 4), ones(1, 4), ... 
    1 1 1 1 1 1 1 1 1 1, ...
    1 1 1 1 1 1 1 1 1 1, ...
    1 1 1 1 1 1 1 1 1 1, ...
    1 1 1 1 1 1 1 1 1 1, ...
    1 1 1 1 1 1 1 1 1 1, ...
    1 1 1 1 1 1 1];

subj_valid_riv = [ones(1, 10), ...
    ones(1, 10), ...
    1 1 1 1 ones(1, 3) 1, ...
    1 1 1 1 1 1 1 1 1 1, ...
    1 1 1 1 1 1 1 1 1 1, ...
    1 1 1 1 1 1 1 1 1 1, ...
    1 1 1 1 1 1 1 1 1 1, ...
    1 1 1 1 1 1 1 1 1 1, ...
    1 1 1 1 1 1 1]; % 3608, 3620, 3628, 3629 with imbalanced rivalry
subj_valid_sfm = [ones(1, 10), ...
    ones(1, 10), ...
    1 1 1 1 ones(1, 3) 1, ...
    1 1 1 1 1 1 1 1 1 1, ... % 3618, 3628, 3648- weird sfm
    1 1 1 1 1 1 1 1 1 1, ...
    1 1 1 1 1 1 1 1 1 1, ...
    1 1 1 1 1 1 1 1 1 1, ...
    1 1 1 1 1 1 1 1 1 1, ...
    1 1 1 1 1 1 1];

subj_valid_eeg = [1 ones(1, 5) 1 1 1 1, ... %3577- noisy session due to coughing, 3587- weird signal lots of low frequency stuff, 3588- just one session valid
    1 1 1 1 1 1 1 1 1 1, ...
    1 1 1 1 1 1 1 1, ... % two only behavioral sessions
    1 1 1 1 1 1 1 1 0 1, ...
    1 1 1 1 1 1 1 1 1 1, ...
    1 1 1 1 1 1 1 1 1 1, ...
    1 1 1 1 1 1 1 1 1 1, ...
    1 1 1 1 1 1 1 1 1 1, ...
    1 1 1 1 1 1 1]; % exclude subject 3577 for excessive noise. 

% subj_valid_eeg = zeros(1, numel(subj_valid_eeg)); % exclude subject 3577 for excessive noise. 
% subj_valid_eeg(9) = 1;

subj_cols = [[14 41 84] ./ 255; [101 53 113] ./ 255; [208 3 121] ./ 255; ...
    [30 194 124] ./ 255; [198 136 19] ./ 255; [48 118 233] ./ 255;...
    [217 131 190] ./ 255; [126 184 130] ./ 255; [98 212 251] ./ 255;...
    [229 188 72] ./ 255; [112 101 163] ./ 255; [164 162 227] ./ 255;...
    [201 169 176] ./ 255; [239 2 35] ./ 255; [242 141 158] ./ 255;...
    [55 16 162] ./ 255; [160 172 103] ./ 255; [178 42 125] ./ 255;...
    [209 176 123] ./ 255; [22 8 216] ./ 255; [99 62 21] ./ 255;...
    [176 33 3] ./ 255; [31 204 75] ./ 255; [12 137 183] ./ 255;...
    [25 128 168] ./ 255; [202 147 132] ./ 255;[158 21 16] ./ 255;...
    [28 20 39] ./ 255; [218 92 15] ./ 255;[94 72 152] ./ 255;...
    [160 131 127] ./ 255; [138 155 170] ./ 255;[98 49 234] ./ 255;...
    [111 82 24] ./ 255;[57 69 44] ./ 255;[106 145 120] ./ 255;...
    [232 132 226] ./ 255;[140 77 247] ./ 255;[82 95 130] ./ 255;...
    [240 4 140] ./ 255;[148 137 249] ./ 255;[27 132 149] ./ 255;...
    [230 164 134] ./ 255;[230 217 157] ./ 255;[45 137 132] ./ 255;...
    [115 44 147] ./ 255;[255 132 119] ./ 255;[194 78 141] ./ 255; ...
    [151 152 38] ./ 255;[168 132 168] ./ 255;[6 226 240] ./ 255; ...
    [218 222 185] ./ 255;[250 235 172] ./ 255;[91 135 195] ./ 255;...
    [184 253 81] ./ 255;[191 105 156] ./ 255;[154 145 162] ./ 255;...
    [14 41 84] ./ 255; [101 53 113] ./ 255; [208 3 121] ./ 255; ...
    [30 194 124] ./ 255; [198 136 19] ./ 255; [48 118 233] ./ 255;...
    [217 131 190] ./ 255; [126 184 130] ./ 255; [98 212 251] ./ 255;...
    [229 188 72] ./ 255; [112 101 163] ./ 255; [164 162 227] ./ 255;...
    [201 169 176] ./ 255; [239 2 35] ./ 255; [242 141 158] ./ 255;...
    [55 16 162] ./ 255; [160 172 103] ./ 255; [178 42 125] ./ 255;...
    [209 176 123] ./ 255; [22 8 216] ./ 255; [99 62 21] ./ 255;...
    [176 33 3] ./ 255; [31 204 75] ./ 255; [12 137 183] ./ 255;...
    [25 128 168] ./ 255; [202 147 132] ./ 255;[158 21 16] ./ 255; ...
    [25 128 168] ./ 255; [202 147 132] ./ 255;[158 21 16] ./ 255];

cd(expt.session_dir)
datadirs = dir;
datadirs = datadirs([datadirs.isdir]);
datadirs = datadirs(3:end);

nsubj = numel(subj_codes);
expt.nsubj = nsubj;

expt.flick_freqs = [14.4 18];

expt.trigs.rest = {'51', '56'};
expt.trigs.fix = {'41', '46'};
expt.trigs.riv = {'12', '17'; '13', '18'};
expt.trigs.sfm = {'22', '27'; '23', '28'};
expt.trigs.mono = {'30', '36'; '33', '39'; '32', 8; '31', 8; '35', 8; '34', 8};
expt.trigs.artif = {'rej1', 'rej2'};

expt.occip_electrodes = {'O1', 'Oz', 'O2', 'POz', 'PO3', 'PO4', 'PO7', 'PO8'};
expt.occiponly_electrodes = {'O1', 'Oz', 'O2'};
expt.pariet_electrodes = {'POz', 'PO3', 'PO4', 'PO7', 'PO8', 'P7', 'P3', ...
    'Pz', 'P4', 'P8', 'CP1', 'CP2', 'CP5', 'CP6'};
expt.parietonly_electrodes = {'P7', 'P3', 'Pz', 'P4', 'P8'};
expt.leftpariet_electrodes = {'P7', 'P3', 'CP5'};
expt.rightpariet_electrodes = {'P4', 'P8', 'CP6'};
expt.front_electrodes = {'F7', 'F3', 'Fz', 'F4', 'F8', 'FC5', 'FC1', 'FC2', 'FC6'};
expt.leftfront_electrodes = {'F7', 'F3', 'FC5', 'FC1'};
expt.rightfront_electrodes = {'F4', 'F8','FC2', 'FC6'};
expt.parietoccip_electrodes = {'O1', 'Oz', 'O2', 'POz', 'PO3', 'PO4', 'PO7', ...
    'PO8', 'Pz', 'P3', 'P4', 'P7', 'P8'};
expt.rightpar_electrodes = {'CP6', 'P8', 'P4'};
expt.leftpar_electrodes = {'CP5', 'P7', 'P3'};
expt.latpar_electrodes = {'PO7', 'PO3', 'P7', 'P3', 'P4', 'P8', 'PO8', 'PO4'};
expt.poboth_electrodes = {'PO3', 'PO4'}; % 1.7
expt.po3_electrodes = {  'POz'     'PO3'    'PO4'}; % 1.5
expt.po4_electrodes = { 'Pz', 'P3', 'P4', 'CP2', 'CP3'}; %1.6
expt.po_electrodes = { 'POz', 'PO3', 'PO4', 'PO8', 'PO7'}; %1.6

expt.fix_cluster= {'O1', 'O2', 'Oz', 'POz', 'PO3', 'PO4', 'PO7', 'Pz', 'P4', ...
    'P8', 'CP1', 'C3'};
expt.rest_cluster = {'O1', 'PO3'};
expt.riv_cluster = {'O2', 'POz', 'Pz', 'CP5', 'CP2'};
expt.onebyf_fix_cluster = {'O1', 'PO4'};
expt.onebyf_riv_cluster = {'PO4'};

expt.fix_cluster_whole = {'O1', 'O2', 'PO3', 'PO4', 'PO7', 'PO8'};
expt.rest_cluster_whole = {'CP6', 'PO8'};
expt.rest_cluster_whole = {'PO8'};
expt.riv_cluster_whole = {'O2', 'POz', 'Pz', 'P3', 'CP5', 'CP2'};


expt.onebyf_eeg_filename = 'rejevs_icacomprem_gaprem_filt_rivindiff.set';

alldirs = [datadirs.name];

for ii = 1:nsubj
    
    findirs = strfind(alldirs, ['_' subj_codes{ii}]);
    snum = 1;
    expt.data(ii).subj_code = subj_codes{ii};
    expt.data(ii).subj_age = subj_age(ii);
    expt.data(ii).subj_sex = subj_sex(ii);
    expt.data(ii).subj_valid = subj_valid(ii);
    expt.data(ii).subj_valid_riv = subj_valid_riv(ii);
    expt.data(ii).subj_valid_sfm = subj_valid_sfm(ii);
    expt.data(ii).subj_valid_eeg = subj_valid_eeg(ii);
    expt.data(ii).subj_millness = subj_millness(ii);
    expt.data(ii).subj_fam_millness = subj_fam_millness(ii);
    expt.data(ii).subj_cols = subj_cols(ii, :);
    expt.data(ii).invalid_riv = invalid_riv(:, ii)';
    expt.data(ii).invalid_sfm = invalid_sfm(:, ii)';
    
    expt.data(ii).session_name{snum} = alldirs(findirs(snum)-6:findirs(snum)+4);
    sdate = alldirs(findirs(snum)-6:findirs(snum)-1); sdate = sdate([3 4 1 2 5 6]);
    expt.data(ii).session_date{snum} = sdate;
    
    expt.data(ii).reverse_rivbuttons = reverse_rivbutton(ii);
    expt.data(ii).bad_fix_runs = badfixruns(ii);
    expt.data(ii).bad_rest_runs = badrestruns(ii);

    expt.data(ii).nsess = numel(findirs);
    if numel(findirs)==2
        snum = 2;
        expt.data(ii).session_name{snum} = alldirs(findirs(2)-6:findirs(2)+6);

        sdate = alldirs(findirs(snum)-6:findirs(snum)-1); 
        sdate = sdate([3 4 1 2 5 6]);
        expt.data(ii).session_date{snum} = sdate;
    end
    
end
end