%% Analysis script for rivalry individual differences

if ~ispc
    addpath(genpath('~/OneDrive - UC Davis/Projects/Analysis_Such'))
    addpath ~/'OneDrive - UC Davis'/Projects/NeuroMatlabToolboxes/Paths/
    addpath ~/'OneDrive - UC Davis'/Projects/Experiments/Rivalry_IndivDiff/Analysis/
    AddEEGPaths
    eeglab
else
    addpath(genpath('C:\OneDrive - UC Davis\Projects\Analysis_Such\'))
    addpath C:\Analysis\NeuroMatlabToolboxes\Paths
    addpath C:\'OneDrive - UC Davis'\Projects\Experiments\Rivalry_IndivDiff\Analysis\
    AddEEGPaths
    eeglab
end

%% initialize experiment

exp_rivind = rivindiff_experiment;

% compute rivalry and sfm durations
pre_post = 3; % 0-60 test subjects, 1-prelim subjects session 1, ...
% 2-prelim subjects session 2, 3-combine subjects, 4-combine subjects w/o
% anxiety and depression

do_gamfits = false;
exclude_mixed = true;
exp_rivind = rivindiff_durations(exp_rivind, pre_post, 1.5, 1, do_gamfits, ...
    exclude_mixed);

%% behavioral analyses
% correlation between sessions of riv and sfm durations and riv and sfm

pre_post = 3; % 0-60 test subjects, 1-prelim subjects session 1, ...
% 2-prelim subjects session 2, 3-combine subjects, 4-combine subjects w/o
% anxiety and depression
rivindiff_plotdurationcorrs(exp_rivind, pre_post);


%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Analysis on full subject set  %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% correlate rivalry and sfm rates with peak alpha frequency
plot_type = [3];
runtype = 'fix';
chanset = 15; % 0-all, 1-occipital, 2-parietal, 3-L-parietal, 4-R-parietal,
% 5-frontal, 6-L-frontal, 7-R-frontal, 9-right parietal, 10-left pariet,
% 11-occip only, 13- par only, 14- PO3, PO4, 15- po3, po4, poz, 16- po4, 17-fix
% cluster 18-rest cluster 19-riv cluster, 20-fix whole cluster, 21-rest
% whole cluster
gauss_model = 2;
pre_post = 01; % 0-60 test subjects, 1-prelim subjects session 1, ...
% 2-prelim subjects session 2, 3-combine subjects, 4-combine subjects w/o
% anxiety and depression
rem_mh_subj = false;
freqband_names = {[]};
analysis_filestart = 'csd_';
iaf_method = 3; % 2-centroid, 3-peak
fano_flag = 0; % 0- mean, 1-fano factor, 2-std, 3-cv
take_log = true;
exclude_mix = false;

for fbn = 1:numel(freqband_names)
    [rivdur_af, subjcode_af, alphafreq, sfmdur_af, subjcode_af_sfm, alphafreq_sfm, ...
        valsubj_af_riv, iaf_allsubj, mean_alpha_aroundpeak, valriv_alphamp, rivdur_all, ...
        subjcode_all]...
        = rivindiff_alphafreq(exp_rivind, runtype, plot_type, chanset, ...
        gauss_model, freqband_names{fbn}, analysis_filestart, pre_post, ...
        iaf_method, rem_mh_subj, fano_flag, take_log, exclude_mix);
end

mean_alpha_aroundpeak = mean_alpha_aroundpeak';

%% robust correlation with peak alpha frequency
do_plots = false;
corr_type = 'Pearson';

rc = robust_correlation(alphafreq, rivdur_af, do_plots);

nol = ~rc.outliers.bivariate.boxplot;
% [rr, pp] = corrcoef(alphafreq, rivdur_af);
% [rr2, pp2] = corrcoef(alphafreq(nol), rivdur_af(nol));
[rr, pp] = corr(alphafreq', rivdur_af', 'Type', corr_type, 'tail', 'both');
[rr2, pp2] = corr(alphafreq(nol)', rivdur_af(nol)', 'Type', corr_type, ...
    'tail', 'both');

figure, scatter(alphafreq, rivdur_af, 80,...
    'f', 'MarkerFaceColor','b')
hold on, scatter(alphafreq(~nol), rivdur_af(~nol), 80,...
    'f', 'MarkerFaceColor','r')
p = polyfit(alphafreq, rivdur_af, 1);
sfit = polyval(p, alphafreq);
hold on, plot(alphafreq, sfit, 'k', 'LineWidth', 3)
text(alphafreq, rivdur_af+.1, subjcode_af)
if fano_flag
    axis([7.2 12.5 0 8])
else
    axis([7.2 12.5 1 6])
end
title(sprintf('r=%g, p=%g; robust r=%g, p=%g', rr, pp, rr2, pp2))

[sort_rivdur, sortinds] = sort(rivdur_af);
sort_af = alphafreq(sortinds);
[h, p, ~, stats] = ...
    ttest2((sort_af(1:floor(end/2+1))), (sort_af(ceil(end/2+1):end)));
fprintf('alpha frequency difference between fast and slow half switchers, p=%g\n',...
    p)

%% robust correlation with alpha amplitude
do_plots = false;

rc = robust_correlation(mean_alpha_aroundpeak, rivdur_all, do_plots);
nol = ~rc.outliers.bivariate.boxplot;
[rr, pp] = corrcoef(mean_alpha_aroundpeak, rivdur_all);
[rr2, pp2] = corrcoef(mean_alpha_aroundpeak(nol), rivdur_all(nol));

figure, scatter(mean_alpha_aroundpeak, rivdur_all, 80,...
    'f', 'MarkerFaceColor','b')
hold on, scatter(mean_alpha_aroundpeak(~nol), rivdur_all(~nol), 80,...
    'f', 'MarkerFaceColor','r')
p = polyfit(mean_alpha_aroundpeak, rivdur_all, 1);
sfit = polyval(p, mean_alpha_aroundpeak);
hold on, plot(mean_alpha_aroundpeak, sfit, 'k', 'LineWidth', 3)
if strcmp(runtype, 'fix')
    axis([0 2.5 1 6])
elseif strcmp(runtype, 'rest')
    axis([0 8.5 1 6])
elseif strcmp(runtype, 'riv')
    axis([0 1.8 1 6])
end
title(sprintf('r=%g, p=%g; robust r=%g, p=%g', rr(2), pp(2), rr2(2), pp2(2)))

%% plot session 1 vs session 2 of alpha freq and rivalry dur
rivdur_af2 = rivdur_af;
subjcode_af2 = subjcode_af;
alphafreq2 = alphafreq;
mean_alpha_aroundpeak2 = mean_alpha_aroundpeak;
rivdur_all2 = rivdur_all;
subjcode_all2 = subjcode_all;
if strcmp(runtype, 'fix')
    valsubj2 = ~strcmp(subjcode_af2, '3587');
else
    valsubj2 = true(1, numel(subjcode_af2));
end

%%
% session 1 rivalry vs. session 2
figure, scatter(rivdur_all, rivdur_all2, 80,...
    'f', 'MarkerFaceColor','b')
% text(rivdur_all, rivdur_all2, subjcode_all)
mm = [min([rivdur_all rivdur_all2]) max([rivdur_all rivdur_all2])];
line(mm, mm, 'Color', 'k', 'LineWidth', 3)
[rr, ~, ~, ~, ~, ~, pp] = ICC([rivdur_all; rivdur_all2]', 'A-1');
% [rr, pp] = corrcoef(rivdur_all, rivdur_all2);
cc = sprintf('riv sess1-2 icc %.3g, p = %.3g', rr, pp);
title(cc)
xlabel('session 1'), ylabel('session 2')
axis square

% session 1 alphafreq vs. session 2 alphafreq
figure, scatter(alphafreq(valsubj2), alphafreq2(valsubj2), 80,...
    'f', 'MarkerFaceColor','b')
% text(alphafreq(valsubj2), alphafreq2(valsubj2), subjcode_af(valsubj2))
mm = [min([alphafreq(valsubj2) alphafreq2(valsubj2)]) max([alphafreq(valsubj2) alphafreq2(valsubj2)])];
line(mm, mm, 'Color', 'k', 'LineWidth', 3)
[rr, ~, ~, ~, ~, ~, pp] = ICC([alphafreq(valsubj2); alphafreq2(valsubj2)]', 'A-1');
% [rr, pp] = corrcoef(alphafreq, alphafreq2);
cc = sprintf('alphafreq sess1-2 icc %.3g, p = %.3g', rr, pp);
title(cc)
xlabel('session 1'), ylabel('session 2')
axis square

% % session 1 rivalry vs. session 2 alphafreq
% figure, scatter(rivdur_af, alphafreq2)
% [rr, pp] = corrcoef(rivdur_af, alphafreq2);
% cc = sprintf('iaf sess1-riv sess2-alphafreq %.3g, p = %.3g', rr(2), pp(2));
% title(cc)

% session 1 alpha amp vs. session 2
figure, scatter(mean_alpha_aroundpeak, mean_alpha_aroundpeak2, 80,...
    'f', 'MarkerFaceColor','b')
% text(mean_alpha_aroundpeak, mean_alpha_aroundpeak2, subjcode_all)
mm = [min([mean_alpha_aroundpeak mean_alpha_aroundpeak2]) max([mean_alpha_aroundpeak mean_alpha_aroundpeak2])];
line(mm, mm, 'Color', 'k', 'LineWidth', 3)
% [rr, pp] = corrcoef(mean_alpha_aroundpeak, mean_alpha_aroundpeak2);
[rr, ~, ~, ~, ~, ~, pp] = ICC([mean_alpha_aroundpeak; mean_alpha_aroundpeak2]', 'A-1');
cc = sprintf('sess-by-sess alpha amp ICC %.3g, p = %.3g', rr, pp);
title(cc)
xlabel('session 1'), ylabel('session 2')
axis square

%% reliability between alpha fequency of fix and rest
alpha2s = []; nval = 0; rd2s = [];
for ap = 1:numel(alphafreq)
    elnum = strcmp(subjcode_af{ap}, subjcode_af2);
    if any(elnum)
        nval = nval+1;
        alpha2s(:, nval) = [alphafreq(ap) alphafreq2(elnum)];
        rd2s(nval) = rivdur_af(ap);
    end
end
figure, scatter(alpha2s(1, :), alpha2s(2, :))
[rr, pp] = corrcoef(alpha2s(1, :), alpha2s(2, :));
mm = [min([alpha2s(1, :) alpha2s(2, :)]) max([alpha2s(1, :) alpha2s(2, :)])];
line(mm, mm)

[rr, ~, ~, ~, ~, ~, pp] = ICC(alpha2s', 'A-1');
cc = sprintf('icc %.3g, p = %.3g', rr, pp);
title(cc)

[B,BINT,R,RINT,STATS] = regress(rd2s', [ones(size(alpha2s, 2), 1) alpha2s']);

%% predict rivalry duration variance from mean and paf

exp_data = [exp_rivind.data];
riv_data = [exp_data.riv_durations];
valsubj = logical([exp_data.subj_valid_eeg]);
nsess = [exp_data.nsess];
svr = logical([exp_data.subj_valid_riv]);
subj_code = {exp_data.subj_code};

sesnum = 1;

rdm = {riv_data.mean_all}; nmrdm = numel(rdm); rivdurmed = zeros(1, nmrdm);
for nrdm = 1:nmrdm, rivdurmed(nrdm) = rdm{nrdm}(sesnum); end
rdm = {riv_data.std_all}; nmrdm = numel(rdm); rivstd = zeros(1, nmrdm);
for nrdm = 1:nmrdm, rivstd(nrdm) = rdm{nrdm}(sesnum); end

% select valid eeg
rivdurmed = rivdurmed(valsubj);
rivstd = rivstd(valsubj);
nsess = nsess(valsubj);
svr = svr(valsubj);
subj_code = subj_code(valsubj);

% select valid rivalry
rivdurmed = rivdurmed(svr);
rivstd = rivstd(svr);
nsess = nsess(svr);
subj_code = subj_code(svr);
valalpha = valsubj_af_riv(svr);

% select valid alpha
rivdurmed = rivdurmed(valalpha);
rivstd = rivstd(valalpha);
nsess = nsess(valalpha);
subj_code = subj_code(valalpha);

for nsamp = 2%[2 1]
    rd = rivdurmed(nsess==nsamp);
    rs = (rivstd(nsess==nsamp).^2)./rd;
    af = alphafreq(nsess==nsamp);
    
    [~,~,~,~,STATS] = regress(rs', [ones(size(rd, 2), 1) rd' af'])
    
    [~,~,~,~,stats] = stepwisefit([rd' af'], rs');
    
    %     [~,~,~,~,STATS] = regress(af', [ones(size(rd, 2), 1) rd' rs'])
    %
    %     [~,~,~,~,stats] = stepwisefit([rd' rs'], af');
    figure, scatter(rd, rs)
end

%% for the full dataset look at the 1/f spectra

plot_type = [2]; % 1-individual subject plots, 2-summary
runtype = 'rest';
chanset = 1; % 0-all, 1-occipital, 2-parietal, 3-L-parietal, 4-R-parietal,
% 5-frontal, 6-L-frontal, 7-R-frontal, 8-all occipital and parietal
% 9-right parietal, 10-left pariet,
% 11-occip only, 12- lateral parietal 13- par only, 14- po both, 15- po3, 16- po4
% 17-fix cluster, 19-riv cluster
% freqband_names = {'slow', 'delta', 'theta', 'alpha', 'beta', 'gamma1', 'gamma2'};

freqband_names = {'alpha'};
freqband_names = {[]};
analysis_filestart = 'csd_';
pre_post = 1; % 0-60 test subjects, 1-prelim subjects session 1, ...
% 2-prelim subjects session 2, 3-combine subjects, 4-combine subjects w/o
% anxiety and depression
fano_flag = 0; % 0- mean, 1-fano factor, 2-std, 3-cv

for fbn = 1:numel(freqband_names)
    if any(plot_type==2)
        [rivdur_eeg1f, subjcode_eeg1f, eeg1f, sfmdur_eeg1f, subjcode_eeg1f_sfm, ...
            eeg1f_sfm, nrettrans]...
            = rivindiff_eeg1byf(exp_rivind, runtype, plot_type, chanset, ...
            freqband_names{fbn}, analysis_filestart, pre_post, fano_flag);
    else
        rivindiff_eeg1byf(exp_rivind, runtype, plot_type, chanset, ...
            freqband_names{fbn}, analysis_filestart, pre_post);
    end
end
eeg1f = eeg1f';

%% robust 1/f slope corelation with mean rivalry

[rr, pp] = corrcoef(eeg1f, rivdur_eeg1f);

rc = robust_correlation(eeg1f, rivdur_eeg1f, false);
nol = ~rc.outliers.bivariate.boxplot;
[rr2, pp2] = corrcoef(eeg1f(nol), rivdur_eeg1f(nol));

figure, scatter(eeg1f, rivdur_eeg1f, 80,...
    'f', 'MarkerFaceColor','b')
hold on, scatter(eeg1f(~nol), rivdur_eeg1f(~nol), 80,...
    'f', 'MarkerFaceColor','r')
p = polyfit(eeg1f, rivdur_eeg1f, 1);
sfit = polyval(p, eeg1f);
hold on, plot(eeg1f, sfit, 'k', 'LineWidth', 3)
axis([-2 -.2 1 6])

title(sprintf('r=%g, p=%g; robust r=%g, p=%g', rr(2), pp(2), rr2(2), pp2(2)))

%% robust 1/f slope corelation with num return transitions

[rr, pp] = corrcoef(eeg1f, nrettrans);

rc = robust_correlation(eeg1f, nrettrans, false);
nol = ~rc.outliers.bivariate.boxplot;
[rr2, pp2] = corrcoef(eeg1f(nol), nrettrans(nol));

figure, scatter(eeg1f, nrettrans, 80,...
    'f', 'MarkerFaceColor','b')
hold on, scatter(eeg1f(~nol), nrettrans(~nol), 80,...
    'f', 'MarkerFaceColor','r')
p = polyfit(eeg1f, nrettrans, 1);
sfit = polyval(p, eeg1f);
hold on, plot(eeg1f, sfit, 'k', 'LineWidth', 3)

title(sprintf('r=%g, p=%g; robust r=%g, p=%g', rr(2), pp(2), rr2(2), pp2(2)))

%%
rivdur_eeg1f2 = rivdur_eeg1f;
subjcode_eeg1f2 = subjcode_eeg1f;
eeg1f2 = eeg1f;

%%
% session 1 alphafreq vs. session 2

% alpha of fixation with rest epochs
eeg1f2s = []; nval = 0;
for ap = 1:numel(eeg1f)
    elnum = strcmp(subjcode_eeg1f{ap}, subjcode_eeg1f2);
    if any(elnum)
        nval = nval+1;
        eeg1f2s(:, nval) = [eeg1f(ap) eeg1f2(elnum)];
    end
end

figure, scatter(eeg1f2s(1, :), eeg1f2s(2, :), 80,...
    'f', 'MarkerFaceColor','b')
% text(eeg1f2s(1, :), eeg1f2s(2, :), subjcode_eeg1f)
mm = [min([eeg1f eeg1f2]) max([eeg1f eeg1f2])];
line(mm, mm, 'Color', 'k', 'LineWidth', 3)
[rr, ~, ~, ~, ~, ~, pp] = ICC(eeg1f2s', 'A-1');
% [rr, pp] = corrcoef(eeg1f, eeg1f2);
cc = sprintf('sess-by-sess 1/f slope icc %.3g, p = %.3g', rr, pp);
title(cc)
axis square
% [B,BINT,R,RINT,STATS] = regress(rd2s', [ones(size(alpha2s, 2), 1) alpha2s']);
% 
% [b,bint,r,rint,stats] = stepwisefit(alpha2s', rd2s');

%% correlation of peak alpha frequency and 1/f slope

nsubjs = numel(subjcode_eeg1f);
nvalsubj = 0;
clear rivdurs varint subjcode
for ns = 1:nsubjs
    sn = find(strcmp(subjcode_eeg1f{ns}, subjcode_af));
    if ~isempty(sn)
        nvalsubj = nvalsubj + 1;
        rivdurs(nvalsubj, :) = [rivdur_eeg1f(ns), rivdur_af(sn)];
        varint(nvalsubj, :) = [eeg1f(ns), alphafreq(sn)];
        subjcode{nvalsubj} = subjcode_eeg1f{ns};
    end
end
[rr, pp] = corrcoef(varint(:, 1), varint(:, 2));
% title(sprintf('corrcoeff = %.2g, p-val = %.2g', rr(2), pp(2)))
% [b,bint,r,rint,stats] = regress(rivdurs(:, 1), [ones(nvalsubj, 1), varint]);
% [rr, pp] = corrcoef(varint(:, 2), varint(:, 1));

rc = robust_correlation(varint(:, 1), varint(:, 2));
nol = ~rc.outliers.bivariate.boxplot;
[rr2, pp2] = corrcoef(varint(nol, 1), varint(nol, 2))
figure, scatter(varint(:, 1), varint(:, 2), 'r')
text(varint(:, 1), varint(:, 2), subjcode)
% pf = polyfit(varint(:, 1), varint(:, 2), 1);
% pv = polyval(pf, [min(varint(:, 1)), max(varint(:, 1))]);
% hold on, line([min(varint(:, 1)), max(varint(:, 1))], pv, 'Color', 'k', 'LineWidth', 2)
hold on, scatter(varint(nol, 1), varint(nol, 2), 'b')
title(sprintf('r=%g p=%g; robust r=%g p=%g', rr(2), pp(2), rr2(2), pp2(2)))



%%
nsubjs = numel(subjcode_eeg1f_sfm);
nvalsubj = 0;
clear sfmdurs varint_sfm subjcode_sfm
for ns = 1:nsubjs
    sn = find(strcmp(subjcode_eeg1f_sfm{ns}, subjcode_af_sfm));
    if ~isempty(sn)
        nvalsubj = nvalsubj + 1;
        sfmdurs(nvalsubj, :) = [sfmdur_eeg1f(ns), sfmdur_af(sn)];
        varint_sfm(nvalsubj, :) = [eeg1f_sfm(ns), alphafreq_sfm(sn)];
        subjcode_sfm{nvalsubj} = subjcode_eeg1f_sfm{ns};
    end
end
figure, scatter(varint_sfm(:, 1), varint_sfm(:, 2))
pf = polyfit(varint_sfm(:, 1), varint_sfm(:, 2), 1);
pv = polyval(pf, [min(varint_sfm(:, 1)), max(varint_sfm(:, 1))]);
hold on, line([min(varint_sfm(:, 1)), max(varint_sfm(:, 1))], pv, 'Color', 'k', 'LineWidth', 2)
[rr, pp] = corrcoef(varint_sfm(:, 1), varint_sfm(:, 2));
title(sprintf('corrcoeff = %.2g, p-val = %.2g', rr(2), pp(2)))
% [b,bint,r,rint,stats] = regress(rivdurs(:, 1), [ones(nvalsubj, 1), varint]);

[b,bint,r,rint,stats] = stepwisefit(varint_sfm, sfmdurs(:, 1));

%% button press rivalry crossover full dataset

plot_type = [8]; % 1-indiv run plots, 2-indiv subject 2x2 transition plots
% 6-indiv subj combined trans plot, 7- plot each subj transition in big plot
runtype = 'riv';
chanset = 0; % 0-all, 1-occipital, 2-parietal, 3-L-parietal, 4-R-parietal, 
% §5-frontal, 6-L-frontal, 7-R-frontal, 8-pariet-occipital,
% 9-parietal-only, 10-occipital-only

% elect_flag: 1-peak, 2-percentile-90, 3-most modulated peak, 4-most
% modulated average around transition, 5-percentile most peaked
% 6-percentile most averaged around trans, 7-, 0-all specified
elect_flag = 01;

analysis_filestart = 'csd_';
pre_post = 01; % 0-60 test subjects, 1-prelim subjects session 1, ...
% 2-prelim subjects session 2, 3-combine subjects, 4-combine subjects w/o
% anxiety, and depression

rivindiff_btn_trans(exp_rivind, runtype, plot_type, chanset, elect_flag, ...
    analysis_filestart, pre_post);

% [subj_codes_btp, validriv_btp, max_mdiff_btp, peakdiff_btp, trans_subj_btp, ...
%     badsubj] = ...
%     rivindiff_btn_trans(exp_rivind, runtype, plot_type, chanset, elect_flag, ...
%     analysis_filestart, pre_post);

%% button press rivalry crossover with FFT tseries full dataset

plot_type = [8]; % 1-indiv run plots, 2-, 7- plot each subj transition in big plot
runtype = 'riv';
chanset = 7; % 0-all, 1-occipital, 2-parietal, 3-L-parietal, 4-R-parietal, 5-frontal, 6-L-frontal, 7-R-frontal, 8-pariet-occipital

% elect_flag: 0-use chans of group, 1-peak, 2-percentile-90, 3-most modulated peak, 4-most
% modulated average around transition, 5-percentile most peaked
% 6-percentile most averaged around trans, 7-
elect_flag = 0;

analysis_filestart = 'alpha_';
pre_post = 3; % 0-60 test subjects, 1-prelim subjects session 1, ...
% 2-prelim subjects session 2, 3-combine subjects, 4-combine subjects w/o
% anxiety, and depression
% dur_before_after = [.7 .7];
dur_before_after = [1 1];
do_erp = true; % this will set highpass filter to .05=20sec, otherwise 0.5=sec
if ~exist('badsubj', 'var'), badsubj = []; end

if do_erp
    [rivtrans, rivdurs, valsubj, subjcode] = rivindiff_btn_trans_rawfft(exp_rivind, runtype, plot_type, chanset, elect_flag, ...
        analysis_filestart, pre_post, [], dur_before_after, badsubj, do_erp);
    
    tt3 = linspace(-2.5, 2.5, size(rivtrans, 2));
    nsubj = size(rivtrans, 1);
    rivtransmooth = rivtrans;
    for ns = 1:nsubj
        rivtransmooth(ns, :) = smooth(rivtrans(ns, :), 100);
    end
    figure, plot(tt3, rivtransmooth)
    legend(num2strcell(rivdurs))
    mxmnrange = [-.6 0; 0 1];
    mnval = min(rivtransmooth(:, tt3>mxmnrange(1, 1)&tt3<mxmnrange(1, 2)), [], 2);
    mxval = max(rivtransmooth(:, tt3>mxmnrange(2, 1)&tt3<mxmnrange(2, 2)), [], 2);
    %         mnval = mean(rivtransmooth(:, tt3>mxmnrange(1, 1)&tt3<mxmnrange(1, 2)), 2);
    %         mxval = mean(rivtransmooth(:, tt3>mxmnrange(2, 1)&tt3<mxmnrange(2, 2)), 2);
    diffval = mxval-mnval;
    
    subjcoderiv = {exp_rivind.data.subj_code};
    nsess = [exp_rivind.data.nsess];
    
    switch pre_post
        case 1
            subjcoderiv = subjcoderiv(nsess==2);
        case 0
            subjcoderiv = subjcoderiv(nsess==1);
    end
    figure, scatter(rivdurs, diffval)
    text(rivdurs, diffval, subjcoderiv)
    %     figure, scatter(rivdurs, iaf_allsubj)
    %     text(rivdurs, iaf_allsubj, subjcoderiv)
    
    %     y = rivdurs(valsubj_af_riv);
    %     x1 = iaf_allsubj(valsubj_af_riv);
    %     x2 = diffval(valsubj_af_riv);
    %     nsubj2 = numel(y);
    %     [b,bint,r,rint,stats] = regress(y', [ones(nsubj2, 1), x1' x2]);
    
else
    % [subj_codes_btp, validriv_btp, ] =
    [subj_codes_btp, validriv_btp, peakdiff_btp] = ...
        rivindiff_btn_trans_rawfft(exp_rivind, runtype, plot_type, chanset, elect_flag, ...
        analysis_filestart, pre_post, [], dur_before_after, badsubj, do_erp);
end

%% stepwise regression for rivalry of the three DVs

nsubjs = numel(subjcode);
nvalsubj = 0;
clear rivdurs2 varint2 subjcode2
for ns = 1:nsubjs
    sn = find(strcmp(subjcode{ns}, subj_codes_btp));
    if ~isempty(sn)
        nvalsubj = nvalsubj + 1;
        rivdurs2(nvalsubj, :) = [rivdurs(ns, :), validriv_btp(sn)];
        varint2(nvalsubj, :) = [varint(ns, :), peakdiff_btp(sn)];
        subjcode2{nvalsubj} = subjcode{ns};
    end
end
figure, scatter(varint2(:, 1), varint2(:, 3))
pf = polyfit(varint2(:, 1), varint2(:, 3), 1);
pv = polyval(pf, [min(varint2(:, 1)), max(varint2(:, 1))]);
hold on, line([min(varint2(:, 1)), max(varint2(:, 1))], pv, 'Color', 'k', 'LineWidth', 2)
[rr, pp] = corrcoef(varint2(:, 1), varint2(:, 3));
title(sprintf('1/f with countermod corrcoeff = %.2g, p-val = %.2g', rr(2), pp(2)))
% [b,bint,r,rint,stats] = regress(rivdurs(:, 1), [ones(nvalsubj, 1), varint]);
figure, scatter(varint2(:, 2), varint2(:, 3))
pf = polyfit(varint2(:, 2), varint2(:, 3), 1);
pv = polyval(pf, [min(varint2(:, 2)), max(varint2(:, 2))]);
hold on, line([min(varint2(:, 2)), max(varint2(:, 2))], pv, 'Color', 'k', 'LineWidth', 2)
[rr, pp] = corrcoef(varint2(:, 2), varint2(:, 3));
title(sprintf('iaf with countermod corrcoeff = %.2g, p-val = %.2g', rr(2), pp(2)))

rdurs = rivdurs2(:, 1);
[b,~,~,~,stats] = stepwisefit(varint2, rdurs);

x1 = varint2(:, 2); x2 = varint2(:, 3);
[b,~,~,~,stats] = regress(rdurs, [ones(numel(rdurs), 1) x1 x2]);

figure, scatter3(x1, x2, rdurs,'filled')
hold on
text(double(x1), double(x2), rdurs, subjcode2)
x1fit = linspace(min(x1), max(x1), 50);
x2fit = linspace(min(x2), max(x2), 50);
[X1FIT,X2FIT] = meshgrid(x1fit,x2fit);
YFIT = b(1) + b(2)*X1FIT + b(3)*X2FIT;
mesh(X1FIT,X2FIT,YFIT)
xlabel('IAF')
ylabel('countermod amp')
zlabel('riv dur')
view(50,10)


%% stepwise regression between alpha and countermod

nsubjs = numel(subjcode_af);
nvalsubj = 0;
clear rivdurs3 varint3 subjcode3
for ns = 1:nsubjs
    sn = find(strcmp(subjcode_af{ns}, subj_codes_btp));
    if ~isempty(sn)
        nvalsubj = nvalsubj + 1;
        rivdurs3(nvalsubj, :) = [rivdur_af(ns), validriv_btp(sn)];
        varint3(nvalsubj, :) = [alphafreq(ns), peakdiff_btp(sn)];
        subjcode3{nvalsubj} = subjcode_af{ns};
    end
end
figure, scatter(varint3(:, 1), varint3(:, 2))
pf = polyfit(varint3(:, 1), varint3(:, 2), 1);
pv = polyval(pf, [min(varint3(:, 1)), max(varint3(:, 1))]);
hold on, line([min(varint3(:, 1)), max(varint3(:, 1))], pv, 'Color', 'k', 'LineWidth', 2)
[rr, pp] = corrcoef(varint3(:, 1), varint3(:, 2));
title(sprintf('iapf with countermod corrcoeff = %.2g, p-val = %.2g', rr(2), pp(2)))

rdurs = rivdurs3(:, 1);
[b,~,~,~,stats] = stepwisefit(varint3, rdurs);

out_sub = false(1, numel(rdurs)); out_sub(3) = true;
[b,~,~,~,stats] = stepwisefit(varint3(~out_sub, :), rdurs(~out_sub));


%% split each individual to slow and fast and look at alpha

plot_type = [];
runtype = 'riv';
chanset = 0; % 0-all, 1-occipital, 2-parietal, 3-L-parietal, 4-R-parietal, 5-frontal, 6-L-frontal, 7-R-frontal, 8-pariet-occipital

analysis_filestart = 'rejevs_';
pre_post = 3; % 0-60 test subjects, 1-prelim subjects session 1, ...
% 2-prelim subjects session 2, 3-combine subjects, 4-combine subjects w/o
% anxiety, and depression

rivindiff_split_alphafreq(exp_rivind, runtype, plot_type, chanset, ...
    analysis_filestart, pre_post);

%% alpha amplitude analysis

plot_type = [3];
runtype = 'rest';

analysis_filestart = 'csd_';
chanset = 0; % 0-all, 1-occipital, 2-parietal, 3-L-parietal, 4-R-parietal,
% 5-frontal, 6-L-frontal, 7-R-frontal, 9-right parietal, 10-left pariet,
% 11-occip only, 13- par only
pre_post = 0; % 0-60 test subjects, 1-prelim subjects session 1, 2-prelim subjects session 2, 3-combine subjects

rivindiff_alphaamp(exp_rivind, runtype, plot_type, chanset, ...
    analysis_filestart, pre_post)


%% age and alpha
valriv = logical([exp_rivind.data.subj_valid_riv]);
subjage = [exp_rivind.data.subj_age];
valage = subjage(valriv);
valsubj = logical([exp_rivind.data.subj_valid])&logical([exp_rivind.data.subj_valid_eeg]);
valriv2 = valriv(valsubj);
valsubj_af = valsubj_af_riv(valriv2);

figure, subplot(121), scatter(valage(valsubj_af), alphafreq)
[rr, pp] = corrcoef(valage(valsubj_af), alphafreq);
sprintf('age vs paf corr = %g, p = %g', rr(2), pp(2))

rivdur = [exp_rivind.data.riv_durations]; meandur = zeros(1, numel(rivdur));
for ma = 1:numel(rivdur), meandur(ma) = rivdur(ma).mean_all(1); end
meandur = meandur(logical(valriv));
subplot(122), scatter(valage, meandur)
[rr, pp] = corrcoef(valage, meandur);
sprintf('age vs rivdur corr = %g, p = %g', rr(2), pp(2))


%% extract epochs

run_types = {'mono'};
% run_types = {'riv', 'sfm'};
% run_types = {'rest'};
file_prefix = 'csd_';
pre_post = 3; % 0-prelim subjects, 1-test subjects

rivindiff_extractepochs(exp_rivind, run_types, file_prefix, ...
    pre_post)

%%
run_types = {'riv'};
freqband_names = {'slow', 'delta', 'theta', 'alpha', 'beta', 'gamma1', 'gamma2'};
freqband_names = {'alpha'};
nbands = numel(freqband_names);
pre_post = 0; % 0-prelim subjects, 1-test subjects
for nb = 1:nbands
    rivindiff_extractepochs(exp_rivind, run_types, [freqband_names{nb} '_'], ...
        pre_post)
end

%% Behavioural figures

exp_data = [exp_rivind.data];
riv_data = [exp_data.riv_durations];
sn = 1;
nsubj = numel(riv_data);
eye_dom = zeros(1, nsubj);
meanall_dur = eye_dom;
mean_oddeve = zeros(nsubj, 2, 3);
mean_oddeve_all = zeros(nsubj, 2);
mean_half = mean_oddeve;
mean_half_all = mean_oddeve_all;
mix_dur = meanall_dur;

for ns = 1:nsubj
    eye_dom(ns) = riv_data(ns).eye_dominance(sn);
    meanall_dur(ns) = riv_data(ns).mean_all(sn);
    mix_dur(ns) = riv_data(ns).means(sn, 3);
    mean_oddeve(ns, :, :) = squeeze(riv_data(ns).mean_oddeve(sn, :, :));
    mean_oddeve_all(ns, :)= squeeze(riv_data(ns).mean_all_oddeve(sn, :));
    mean_half(ns, :, :) = squeeze(riv_data(ns).mean_half(sn, :, :));
    mean_half_all(ns, :)= squeeze(riv_data(ns).mean_all_half(sn, :));
end
svr = logical([exp_data.subj_valid_riv]);
nsess = [exp_data.nsess];
subj_code = {exp_data.subj_code};

%%% figure 1B - histogram of mean durations across subjects + odd-even
%%% split of the same thing
[nh, bin] = hist(meanall_dur(svr), .5:.5:6.5);
figure, plot(bin, nh, 'LineWidth', 4)
title(['mean=' num2str(mean(meanall_dur(svr))),...
    'std=' num2str(std(meanall_dur(svr))/sqrt(sum(svr)))])
[nh, bin] = hist(mean_oddeve_all(svr, 1), .5:.5:6.5);
hold on, plot(bin, nh, 'LineWidth', 2)
[nh, bin] = hist(mean_oddeve_all(svr, 2), .5:.5:6.5);
hold on, plot(bin, nh, 'LineWidth', 2)
xlabel('mean duration (s)')
ylabel('frequency')

%%% measuring reliability via odd even split of the dominance durations
figure, scatter(mean_oddeve_all(svr, 1), mean_oddeve_all(svr, 2))
mm = [min([mean_oddeve_all(svr, 1); mean_oddeve_all(svr, 2)]), ...
    max([mean_oddeve_all(svr, 1); mean_oddeve_all(svr, 2)])];
line(mm, mm)
[rr, LB, UB, F, df1, df2, pp] = ICC(mean_oddeve_all, 'A-1');
cc = sprintf('even odd runs dominance icc %.3g, p = %.3g', rr, pp);
title(cc)

%%% measuring reliability via odd even split of the mix durations
figure, scatter(mean_oddeve(svr, 1, 3), mean_oddeve(svr, 2, 3))
mm = [min([mean_oddeve(svr, 1,3); mean_oddeve(svr, 2,3)]), ...
    max([mean_oddeve(svr, 1,3); mean_oddeve(svr, 2,3)])];
line(mm, mm)
mixnotnan = ~isnan(mean_oddeve(:, 1, 3)) & ~isnan(mean_oddeve(:, 2, 3)); mixnotnan = mixnotnan';
[rr, LB, UB, F, df1, df2, pp] = ICC(mean_oddeve(svr&mixnotnan, :, 3), 'A-1');
cc = sprintf('even odd runs mixed icc %.3g, p = %.3g', rr, pp);
title(cc)


%%% correlate dominance and mixed durations
mixnan = isnan(mix_dur);

meansamp = meanall_dur(svr & nsess==2 & ~mixnan); mixsamp = mix_dur(svr & nsess==2& ~mixnan);
[rr, pp] = corrcoef(meansamp, mixsamp);
rc = robust_correlation(meansamp, mixsamp, false);
nol = ~rc.outliers.bivariate.boxplot;
figure, scatter(meansamp, mixsamp, 80,...
    'f', 'MarkerFaceColor','b')
hold on, scatter(meansamp(~nol), mixsamp(~nol), 80,...
    'f', 'MarkerFaceColor','r')
p = polyfit(meansamp, mixsamp, 1);
sfit = polyval(p, meansamp);
hold on, plot(meansamp, sfit, 'k', 'LineWidth', 3)
[rr2, pp2] = corrcoef(meansamp(nol), mixsamp(nol));
title(sprintf('r=%g, p=%g; robust r=%g, p=%g', rr(2), pp(2), rr2(2), pp2(2)))
axis([1.5 6 0 16])

meansamp = meanall_dur(svr & nsess==1 & ~mixnan); mixsamp = mix_dur(svr & nsess==1 & ~mixnan);
[rr, pp] = corrcoef(meansamp, mixsamp);
rc = robust_correlation(meansamp, mixsamp, false);
nol = ~rc.outliers.bivariate.boxplot;
figure, scatter(meansamp, mixsamp, 80,...
    'f', 'MarkerFaceColor','b')
hold on, scatter(meansamp(~nol), mixsamp(~nol), 80,...
    'f', 'MarkerFaceColor','r')
p = polyfit(meansamp, mixsamp, 1);
sfit = polyval(p, meansamp);
hold on, plot(meansamp, sfit, 'k', 'LineWidth', 3)
[rr2, pp2] = corrcoef(meansamp(nol), mixsamp(nol));
title(sprintf('r=%g, p=%g; robust r=%g, p=%g', rr(2), pp(2), rr2(2), pp2(2)))
axis([1.5 6 0 16])


%%% fano-factor intersession reliability
rivdur2sess = zeros(nsubj, 2); rivstd2sess = rivdur2sess; 
rdm = {riv_data.mean_all}; nmrdm = numel(rdm); 
for nrdm = 1:nmrdm, if nsess(nrdm)==2, rivdur2sess(nrdm, :) = rdm{nrdm}; end, end
rivdur2sess = rivdur2sess(nsess==2, :);
rdm = {riv_data.std_all}; nmrdm = numel(rdm); 
for nrdm = 1:nmrdm, if nsess(nrdm)==2, rivstd2sess(nrdm, :) = rdm{nrdm}; end, end
rivstd2sess = rivstd2sess(nsess==2, :);
fano2sess = rivstd2sess.^2./rivdur2sess;
figure, scatter(fano2sess(:, 1), fano2sess(:, 2))
[rr, ~, ~, ~, ~, ~, pp] = ICC(fano2sess, 'A-1');
title(sprintf('Fano factor session1-2 reliability, ICC: r=%g, p=%g', rr, pp))


rdm = {riv_data.mean_all}; nmrdm = numel(rdm); rivdurmed = zeros(1, nmrdm);
for nrdm = 1:nmrdm, rivdurmed(nrdm) = rdm{nrdm}(1); end
rdm = {riv_data.std_all}; nmrdm = numel(rdm); rivstd = zeros(1, nmrdm);
for nrdm = 1:nmrdm, rivstd(nrdm) = rdm{nrdm}(1); end
rivdurmed = rivdurmed(svr); rivstd = rivstd(svr);
nsess2 = nsess(svr);

%%% rivalry duration mean with std
for ns = [2 1]
    rd = rivdurmed(nsess2==ns);
    rs = rivstd(nsess2==ns);

    [rr, pp] = corrcoef(rd, rs);

    rc = robust_correlation(rd, rs);
    nol = ~rc.outliers.bivariate.boxplot;

    figure, scatter(rd, rs, 80,...
        'f', 'MarkerFaceColor','b')
    hold on, scatter(rd(~nol), rs(~nol), 80,...
        'f', 'MarkerFaceColor','r')
    p = polyfit(rd, rs, 1);
    sfit = polyval(p, rd);
    hold on, plot(rd, sfit, 'k', 'LineWidth', 3)
    
    hold on, scatter(rd(nol), rs(nol), 'b')
    [rr2, pp2] = corrcoef(rd(nol), rs(nol));
    
    title(sprintf('dur mean vs. std r=%g, p=%g; robust r=%g, p=%g', rr(2), pp(2), rr2(2), pp2(2)))
end

%%% rivalry duration mean with fano factor
for ns = [2 1]
    rd = rivdurmed(nsess2==ns);
    rs = (rivstd(nsess2==ns).^2)./rivdurmed(nsess2==ns);

    [rr, pp] = corrcoef(rd, rs);

    rc = robust_correlation(rd, rs);
    nol = ~rc.outliers.bivariate.boxplot;

    figure, scatter(rd, rs, 80,...
        'f', 'MarkerFaceColor','b')
    hold on, scatter(rd(~nol), rs(~nol), 80,...
        'f', 'MarkerFaceColor','r')
    p = polyfit(rd, rs, 1);
    sfit = polyval(p, rd);
    hold on, plot(rd, sfit, 'k', 'LineWidth', 3)
    
    hold on, scatter(rd(nol), rs(nol), 'b')
    [rr2, pp2] = corrcoef(rd(nol), rs(nol));
    
    title(sprintf('fano factor r=%g, p=%g; robust r=%g, p=%g', rr(2), pp(2), rr2(2), pp2(2)))
end

%% look at the ICC at each electrode for paf between fix and rest to see why they are so different

cd(exp_rivind.session_dir)

samp = 2;

fd = load(['indiff_corr_data_s' num2str(samp) '_fix']);
rd = load(['indiff_corr_data_s' num2str(samp) '_rest']);
nchans = numel(fd.chanloc34);

fd.subjind = find(fd.valriv2&~fd.noalphasubjriv2);
fd.alliaf = NaN(numel(fd.valriv2), nchans);
fd.alliaf(fd.subjind, :) = fd.iaf5;

rd.subjind = find(rd.valriv2&~rd.noalphasubjriv2);
rd.alliaf = NaN(numel(rd.valriv2), nchans);
rd.alliaf(rd.subjind, :) = rd.iaf5;

valsubj = ~fd.noalphasubjriv2 & ~rd.noalphasubjriv2 & fd.valriv2 & rd.valriv2;
rr = zeros(1, nchans); pp = rr; tt = pp;
fd.malphaval = fd.alliaf(valsubj, :);
rd.malphaval = rd.alliaf(valsubj, :);

for nc = 1:nchans
    [rr(nc), ~, ~, ~, ~, ~, pp(nc)] = ...
        ICC([rd.alliaf(valsubj, nc) fd.alliaf(valsubj, nc)], 'A-1');
end
figure, topoplot(rr, rd.chanloc34, 'maplimits', [.52 .9]);

% [h, p] = ttest(rd.malphaval-fd.malphaval);
% figure, subplot(221), topoplot(mean(fd.malphaval), rd.chanloc34, 'maplimits', [8.6 10.2]);
% title('mean alpha freq - fix')
% subplot(222), topoplot(mean(rd.malphaval), rd.chanloc34, 'maplimits', [8.6 10.2]);
% title('mean alpha freq - rest')
% subplot(223), topoplot(mean(fd.iaf5)-mean(rd.iaf5), rd.chanloc34, 'maplimits', [-.3 .3]);
% title('mean alpha freq diff fix-rest')
% subplot(224), topoplot(h, rd.chanloc34, 'maplimits', [0 1]);
% title('p<.05 if red')

%%