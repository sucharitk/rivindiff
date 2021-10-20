function varargout = rivindiff_alphafreq(exp_rivind, runtype, plot_type, chanset, ...
    gauss_mod, freq_band, analysis_filestart, pre_post_subj, iaf_method, ...
    rem_mh_subj, fano_flag, take_log, exclude_mix)
%
% function rivindiff_alphafreq(exp_rivind, runtype, plot_type, chanset, ...
%     gauss_mod, pre_post_subj)
%
% exp_rivind: exp datastructure
% runtype:  type of run
% plot_type:
% chanset
% gauss_mod
% pre_post_subj- 0-prelim subjects, 1-test subjects
%


cd(exp_rivind.session_dir)

nsubj = exp_rivind.nsubj;
eegdatapath = exp_rivind.eeg_datadir;
epochdir = exp_rivind.eeg_epochdir;

subj_ii = 0;
sess_ii = 0;
if ~exist('runtype', 'var'), runtype = 'fix'; end

% gauss_model = ['gauss' num2str(gauss_mod)];

trigs = exp_rivind.trigs.(runtype);
ntrigs = size(trigs, 1);
freq_range = [.1 50];

alpha_freq_range = [7 14];
% alpha_freq_range = [2 7];
% alpha_freq_range = [25 59];
% alphit_freq_range = [1 5];
alphit_freq_range = [.5 5.5];

alpha_amp_calc = [-.75 .75]; % range of frequencies over which to compute alpha amplitude
% alpha_amp_calc = [-1 1]; % range of frequencies over which to compute alpha amplitude

corr_type = 'Pearson';
% corr_type = 'Spearman';

if strcmp(runtype, 'riv')
    alpha_smooth = 100;
else
    alpha_smooth = 150;
end

valsubj = false(1, nsubj);

if strcmp(runtype, 'rest')
    lastn = 10;
else
    lastn = 0;
end

switch chanset
    case 1
        chans_to_an = exp_rivind.occip_electrodes;
    case 2
        chans_to_an = exp_rivind.pariet_electrodes;
    case 3
        chans_to_an = exp_rivind.leftpariet_electrodes;
    case 4
        chans_to_an = exp_rivind.rightpariet_electrodes;
    case 5
        chans_to_an = exp_rivind.front_electrodes;
    case 6
        chans_to_an = exp_rivind.leftfront_electrodes;
    case 7
        chans_to_an = exp_rivind.rightfront_electrodes;
    case 8
        chans_to_an = exp_rivind.parietoccip_electrodes;
    case 9
        chans_to_an = exp_rivind.rightpar_electrodes;
    case 10
        chans_to_an = exp_rivind.leftpar_electrodes;
    case 11
        chans_to_an = exp_rivind.occiponly_electrodes;
    case 12
        chans_to_an = exp_rivind.latpar_electrodes;
    case 13
        chans_to_an = exp_rivind.parietonly_electrodes;
    case 14
        chans_to_an = exp_rivind.poboth_electrodes;
    case 15
        chans_to_an = exp_rivind.po3_electrodes;
    case 16
        chans_to_an = exp_rivind.po4_electrodes;
    case 17
        chans_to_an = exp_rivind.fix_cluster;
    case 18
        chans_to_an = exp_rivind.rest_cluster;
    case 19
        chans_to_an = exp_rivind.riv_cluster;
    case 20
        chans_to_an = exp_rivind.fix_cluster_whole;
    case 21
        chans_to_an = exp_rivind.rest_cluster_whole;
end

if isempty(freq_band)
    compstr = ['^' analysis_filestart];
else
    compstr = [freq_band '_' analysis_filestart];
end

sessnum = 1;
if pre_post_subj == 2
    sessnum = 2;
end

for ns = 1:nsubj
    % subject loop
    
    subj_data = exp_rivind.data(ns);
    
    if subj_data.subj_valid && subj_data.subj_valid_eeg
        if (pre_post_subj && subj_data.nsess>1) || ...
                (~pre_post_subj && subj_data.nsess==1) || ...
                pre_post_subj>=3
            % for subjects that have two sessions
            subj_ii = subj_ii + 1;
            valsubj(ns) = true;
            
            if any(plot_type==1)
                figure('Name', subj_data.subj_code)
                hold on
            end
            if any(plot_type==2) && ~mod(subj_ii-1, 16)
                f2 = figure;
                subplot(441)
            end
            
            for snum = sessnum:sessnum
                sess_ii = sess_ii + 1;
                
                sessdir = fullfile(exp_rivind.session_dir, subj_data.session_name{snum});
                if isdir(fullfile(sessdir, eegdatapath))
                    
                    cd(sessdir)
                    
                    ff = dir(epochdir);
                    ff2 = find(~isemptycell(regexp({ff.name}, compstr)), 1, 'last');
                    if ~isempty(ff2)
                        rivn=strfind(ff(ff2).name, 'rivindiff');
                        fname = ff(ff2).name(1:rivn+8);
                        
                        for ntr = 1:ntrigs
                            resfname = fullfile('Epochs', sprintf('%s_%s_%s', ...
                                fname, runtype, trigs{ntr}));
                            
                            load(resfname)
                            
                            if ~exist('NFFT', 'var')
                                NFFT = 2^16;
                                freqres = 2/fs;
                                freqs = linspace(0, fs/2, NFFT/2+1);
                                ff_range = freqs>=freq_range(1) & freqs<=freq_range(2);
                                numfel = numel(ff_range);
                                
                                if any(alpha_freq_range)
                                    if numel(alpha_freq_range)==2
                                        ff_alpha = freqs>=alpha_freq_range(1) & freqs<=alpha_freq_range(2);
                                    elseif numel(alpha_freq_range)==4
                                        ff_alpha = freqs>=(alpha_freq_range(1) & freqs<=alpha_freq_range(2))...
                                            | (freqs>=(alpha_freq_range(3) & freqs<=alpha_freq_range(4)));
                                    end
                                end
                                if chanset
                                    chanlab34 = chans_to_an;
                                else
                                    chanlab34 = {chanlocs.labels};
                                end
                                
                                nchanan = 34;
                            end
                            
                            chanlabels = {chanlocs.labels};
                            if chanset
                                chanloc_inds = get_channels_from_labels(chanlabels, chans_to_an);
                                nbchan = sum(chanloc_inds);
                                chanlabels = chanlabels(chanloc_inds);
                                nchanan = numel(chans_to_an);
                            else
                                % average over all available channels
                                nbchan = size(epochs{1}, 1);
                                chanloc_inds = true(1, nbchan);
                            end
                            
                            nepochs = numel(epochs);
                            
                            ft_eps = zeros(nepochs, nbchan, numfel);
                            for ne = 1:nepochs
                                valid_epoch = epochs{ne}(chanloc_inds, :);
                                %                             valid_epoch = epochs{ne}(chanloc_inds, artif_mask{ne});
                                ft_ep = AbsFFT_2(valid_epoch, fs, NFFT);
                                ft_eps(ne, :, :) = ft_ep(:, 1:numfel);
                            end
                            noeog = ~(strcmp(chanlabels, ...
                                {'VEOG'}) | strcmp(chanlabels, {'HEOG'}));
                            if ntr==1
                                % in case of two different types of
                                % triggers as in case of the rivalry runs
                                % save the first one and append it
                                % outside the loop
                                ft_eps2 = ft_eps;
                            end
                        end
                        if ntrigs>1
                            ft_eps = [ft_eps2; ft_eps];
                        else
                            clear ft_eps2
                        end
                        
                        switch iaf_method
                            case 1
                                
                            case 2
                                % centroid method
                                fteps = squeeze(mean(ft_eps(:, :, ff_range), 1));
                                fff = freqs(ff_range)';
                                fteps2 = squeeze(mean(ft_eps(:, :, ff_alpha), 1));
                                fff2 = freqs(ff_alpha)';
                                
                                lval = zeros(size(fteps2));
                                for nc = 1:nbchan
                                    lfit = polyfit(fff, fteps(nc, :)', 1); % fit trend line first
                                    lval(nc, :) = polyval(lfit, fff2);
                                end
                                
                                fteps = fteps2-lval; % subtract trend get alpha`
                                for nc = 1:nbchan
                                    fteps(nc, :) = smooth(fteps(nc, :), 100);
                                end
                                fteps2 = fteps-min(fteps(:));
                                fff2 = fff2';
                                
                                if nbchan<nchanan
                                    alp = indiv_alpha_freq(fff2, fteps2);
                                    chanloc_inds = get_channels_from_labels(chanlab34, chanlabels);
                                    iaf(subj_ii, chanloc_inds) = alp;
                                elseif nbchan>nchanan
                                    alp = indiv_alpha_freq(fff2, fteps2);
                                    iaf(subj_ii, :) = alp(noeog);
                                else
                                    iaf(subj_ii, :) = indiv_alpha_freq(fff2, fteps2);
                                end
                                
                                if any(plot_type==2)
                                    figure(f2)
                                    subplot(4, 4, mod(subj_ii-1, 16)+1)
                                    plot(fff2, fteps)
                                    title([subj_data.subj_code ': ' ...
                                        num2str(subj_data.riv_durations.mean_all(snum))])
                                end
                                
                                
                            case 3
                                % smooth and find the absolute peak
                                alpha_freq_range = [6.5 13.5];
                                ff_alpha = freqs>=alpha_freq_range(1) & freqs<=alpha_freq_range(2);
                                ff_alphit = freqs>=alphit_freq_range(1) & freqs<=alphit_freq_range(2);
                                
                                fteps = squeeze(mean(ft_eps(:, noeog, ff_range), 1));
                                fff = freqs(ff_range)';
                                fteps2 = squeeze(mean(ft_eps(:, noeog, ff_alpha), 1));
                                fff2 = freqs(ff_alpha)';
                                if strcmp(runtype, 'riv')
                                    interp_pts = (fff2>10.77 & fff2<10.83)|...
                                    (fff2>7.97 & fff2<8.03)|(fff2>11.97 & fff2<12.03)...
                                    |(fff2>7.17 & fff2<7.23);
                                    finterp = fff2(interp_pts);
                                    for nc = 1:nbchan
                                        fteps_interp = interp1(fff2(~interp_pts), ...
                                            fteps2(nc, ~interp_pts), finterp);
                                        fteps2(nc, interp_pts) = fteps_interp;
                                    end
                                end
                                fteps3 = squeeze(mean(ft_eps(:, noeog, ff_alphit), 1));
                                fff3 = freqs(ff_alphit)';
                                fteps6 = fteps2;
                                fff6 = fff2;
                                nbchan = sum(noeog);
                                if nbchan==1, fteps = fteps'; fteps2 = fteps2'; fteps3 = fteps3'; end
                                
                                if (strcmp(subj_data.subj_code, '3609')) ||...
                                        (strcmp(runtype, 'riv')&&(strcmp(subj_data.subj_code, '3679'))) ||...
                                        (strcmp(runtype, 'riv')&&(strcmp(subj_data.subj_code, '3653'))) ||...
                                        (strcmp(runtype, 'riv')&&(strcmp(subj_data.subj_code, '3674'))) ||...
                                        (strcmp(runtype, 'riv')&&(strcmp(subj_data.subj_code, '3690'))) ||...
                                        ~take_log
                                    fteps5 = fteps2;
                                else
                                    lval = zeros(size(fteps2));
                                    lval2 = zeros(size(fteps));
                                    for nc = 1:nbchan
                                        lfit = polyfit(log10(fff3), log10(fteps3(nc, :))', 1); % fit trend line first
                                        lval(nc, :) = polyval(lfit, log10(fff2));
                                        lval2(nc, :) = polyval(lfit, log10(fff));
                                    end
                                    fteps5 = 10.^(log10(fteps2)-lval); % subtract trend get alpha
                                    %                                     fteps4 = 10.^(log10(fteps)-lval2);
                                end
                                
                                
                                for nc = 1:nbchan
                                    fteps5(nc, :) = smooth(fteps5(nc, :), alpha_smooth);
                                end
                                fteps2 = fteps5-min(fteps5(:));
                                fff2 = fff2';
                                
                                within_alpha = [7 13];
                                fff4 = (fff2>within_alpha(1) & fff2<within_alpha(2));
                                fteps2 = fteps2(:, fff4);
                                fff2 = fff2(fff4);
                                mfteps2 = mean(fteps2, 1);
                                
                                [~, mxind] = max(mfteps2);
                                iaf(subj_ii) = fff2(mxind); % record the iaf averaged across channels for each subject
                                mx_alpha(subj_ii) = mfteps2(mxind); % record the alphamp averaged across channels for each subject
                                [malp, mxind] = max(fteps2, [], 2);
                                alp = fff2(mxind);
                                
                                smind = find(fff2(1)==fff);
                                %                                 indrange = [mxind+smind+alpha_amp_calc(1)*fs mxind+smind+fs]; % calculate alpha amp between peak apha±2hz
                                indrange = mxind+smind+fs*alpha_amp_calc; % calculate alpha amp between peak apha±2hz
                                clear alpha_amp
                                for nc = 1:nbchan
                                    alpha_amp(nc) = mean(fteps(nc, indrange(nc, 1):indrange(nc, 2)));
                                end
                                                                
                                if nbchan<nchanan
                                    chanloc_inds = get_channels_from_labels(chanlab34, chanlabels);
                                    iaf_chan(subj_ii, chanloc_inds) = alp;
                                    mx_alpha_chan(subj_ii, chanloc_inds) = malp;
                                    alpha_amp_chan(subj_ii, chanloc_inds) = alpha_amp;
                                else
                                    iaf_chan(subj_ii, :) = alp;
                                    mx_alpha_chan(subj_ii, :) = malp;
                                    alpha_amp_chan(subj_ii, :) = alpha_amp;
                                end
                                
                                if any(plot_type==5) && strcmp(subj_data.subj_code, '3598')
                                    % plot figure 2 for the paper
                                    % describing methods to calculate paf,
                                    % alpha and 1/f
                                    alpha_plotrange = [.5 20];
                                    
                                    ff_alphit = freqs>=alpha_plotrange(1) & freqs<=alpha_plotrange(2);
                                    ftep_plot = squeeze(mean(ft_eps(:, noeog, ff_alphit), 1));
                                    fff_plot = freqs(ff_alphit)';
                                    lftep_plot = log10(mean(ftep_plot, 1));
                                    lfff_plot = log10(fff_plot);
                                    figure
                                    subplot(221)
                                    plot(lfff_plot, lftep_plot)
                                    lfff3 = log10(fff3);
                                    lfit = polyfit(lfff3, log10(mean(fteps3, 1))', 1); % fit trend line first
                                    lval3 = polyval(lfit, lfff3);
                                    lval4 = polyval(lfit, lfff_plot);
                                    hold on, line([lfff3(1) lfff3(end)], ...
                                        [lval3(1) lval3(end)], 'Color', 'r', 'LineWidth', 3)
                                    ax = axis;
                                    ax(1:2) = log10(alpha_plotrange(1:2));
                                    axis(ax)
                                    
                                    subplot(222)
                                    lftep_plot = lftep_plot-lval4';
                                    plot(lfff_plot, lftep_plot)
                                    
                                    ax(3) = -.3;
                                    ax(4) = 1.1;
                                    axis(ax)
                                    subplot(223)
                                    plot(fff2, mfteps2, 'LineWidth', 2)
                                    [~, mxind] = max(mfteps2);
                                    mxf = fff2(mxind); 
                                    mmf = fff2([mxind-fs, mxind+fs]);
                                    hold on, line([mxf mxf; mmf(1) mmf(1);...
                                        mmf(2) mmf(2)]', [0 mfteps2(mxind); ...
                                        0 mfteps2(mxind-fs); 0 mfteps2(mxind+fs)]')
                                    
                                end
                                
                                if any(plot_type==2)
                                    figure(f2)
                                    subplot(4, 4, mod(subj_ii-1, 16)+1)
                                    plot(fff2, mfteps2)
                                    %                                     plot(fff6, mean(fteps6, 1))
                                    line([iaf(subj_ii) iaf(subj_ii)], [min(mfteps2) max(mfteps2)])
                                    title([subj_data.subj_code ': ' ...
                                        num2str(subj_data.riv_durations.mean_all(snum)) ': '...
                                        num2str(iaf(subj_ii))])
                                end
                        end
                        
                    else
                        sprintf('Session data not available %s', subj_data.subj_code)
                        nosess = true;
                    end
                else
                    nosess = true;
                end
            end
        end
    end
end

% alphampcutoff = .0012;
% alpha_freq_range = [7 14];
noalpha = noalphasubj(lastn, analysis_filestart, chanset, runtype);

if any(plot_type==3)
    
    valriv = [exp_rivind.data.subj_valid_riv];
    valriv2 = logical(valriv(valsubj));
    millness_riv = [exp_rivind.data.subj_millness];
    millness_riv = millness_riv(valsubj);
    millness_riv = millness_riv(valriv2);
    fam_millness_riv = [exp_rivind.data.subj_fam_millness];
    fam_millness_riv = fam_millness_riv(valsubj);
    fam_millness_riv = fam_millness_riv(valriv2);

    valsfm = [exp_rivind.data.subj_valid_sfm];
    valsfm2 = logical(valsfm(valsubj));
    %     millness_sfm = [exp_rivind.data.subj_millness];
    %     millness_sfm = millness_sfm(valsubj);
    %     millness_sfm = millness_sfm(valsfm2);
    
    % first sets of sessions
    rivdur = [exp_rivind.data(valsubj&valriv).riv_durations];
    %     if pre_post_subj>1
    
    %%% read mean durations for each subject
    if exclude_mix
        rdm = {rivdur.nomix_meanall_remoutlier}; nmrdm = numel(rdm); rivdurmed = zeros(1, nmrdm);
    else
        rdm = {rivdur.mean_all}; nmrdm = numel(rdm); rivdurmed = zeros(1, nmrdm);
    end
    for nrdm = 1:nmrdm, rivdurmed(nrdm) = rdm{nrdm}(sessnum); end
    %%% read standard deviation of durations for each subject
    rdm = {rivdur.std_all}; nmrdm = numel(rdm); rivdurstd = zeros(1, nmrdm);
    for nrdm = 1:nmrdm, rivdurstd(nrdm) = rdm{nrdm}(sessnum); end
    switch fano_flag
        case 1
            rivdurmed = (rivdurstd.^2)./rivdurmed; % fano factor
        case 2
            rivdurmed = rivdurstd; % stdev
        case 3
            rivdurmed = rivdurstd./rivdurmed; % cv
    end
    
    subjcoderiv = {exp_rivind.data(valsubj&valriv).subj_code};
    
    rdm = {rivdur.means};
    for nrdm = 1:nmrdm, mixdur(nrdm) = rdm{nrdm}(sessnum, 3); end

    sfmdur = [exp_rivind.data(valsubj&valsfm).sfm_durations];
    if pre_post_subj>1
        rdm = {sfmdur.mean_all}; nmrdm = numel(rdm); sfmdurmed = zeros(1, nmrdm);
        for nrdm = 1:nmrdm, sfmdurmed(nrdm) = rdm{nrdm}(sessnum); end
    else
        sfmdurmed = [sfmdur.mean_all];
    end
    subjcodesfm = {exp_rivind.data(valsubj&valsfm).subj_code};
    
    switch iaf_method
        case 1
            
        case 2
            iaf(~iaf) = NaN;
            
            iaf2 = nanmean(iaf, 2);
            iaf3 = iaf2(valriv2);
            
            noalphasubjriv = false(1, size(iaf3, 1));
            for nos = 1:numel(noalpha)
                noalphasubjriv(strcmp(subjcoderiv, noalpha(nos))) = true;
            end
            rivdurmed = rivdurmed(:, ~noalphasubjriv);
            iaf3 = iaf3(~noalphasubjriv);
            subjcoderiv = subjcoderiv(~noalphasubjriv);
            
            %             iaf4 = iaf2(valriv2);
            if size(rivdurmed, 1)==2, rivdurmed = rivdurmed(sessnum, :); end
            figure
            %             subplot(121)
            scatter(rivdurmed,iaf3)
            [rr, pp] = corrcoef(rivdurmed, iaf3);
            cc = sprintf('iaf riv corr %.3g, p = %.3g', rr(2), pp(2));
            title(cc)
            text(rivdurmed, iaf3+.025, subjcoderiv)
            p = polyfit(iaf3', rivdurmed, 1);
            sfit = polyval(p, iaf3');
            hold on, plot(sfit, iaf3, 'LineWidth', 3)
            
            %             iaf4 = iaf2(valsfm2);
            %             noalphasubjsfm = false(1, size(iaf4, 1));
            %             for nos = 1:numel(noalpha)
            %                 noalphasubjsfm(strcmp(subjcodesfm, noalpha(nos))) = true;
            %             end
            %             sfmdurmed = sfmdurmed(:, ~noalphasubjsfm);
            %             iaf4 = iaf4(~noalphasubjsfm);
            %             subjcodesfm = subjcodesfm(~noalphasubjsfm);
            %
            %             iaf4 = iaf2(valsfm2);
            %             if size(sfmdurmed, 1)==2, sfmdurmed = sfmdurmed(sessnum, :); end
            %             subplot(122), scatter(sfmdurmed, iaf4)
            %             [rr, pp] = corrcoef(sfmdurmed, iaf4);
            %             cc = sprintf('iaf sfm corr %.3g, p = %.3g', rr(2), pp(2));
            %             title(cc)
            %             text(sfmdurmed, iaf4+.025, subjcodesfm)
            %             p = polyfit(iaf4', sfmdurmed, 1);
            %             sfit = polyval(p, iaf4');
            %             hold on, plot(sfit, iaf4, 'LineWidth', 3)
            
            chanlabels = {chanlocs.labels};
            chanloc34 = chanlocs(~(strcmp(chanlabels, ...
                {'VEOG'}) | strcmp(chanlabels, {'HEOG'})));
            iaf5 = iaf(valriv2, :); iaf5 = iaf5(~noalphasubjriv, :);
            %             iaf6 = iaf(valsfm2, :); iaf6 = iaf6(~noalphasubjsfm, :);
            rpchanrivsfm = zeros(2, 2, nchanan);
            for nc = 1:nchanan
                % riv
                nansubj = isnan(iaf5(:, nc));
                [rr, pp] = corrcoef(iaf5(~nansubj, nc), rivdurmed(~nansubj));
                rpchanrivsfm(1, 1, nc) = rr(2); rpchanrivsfm(1, 2, nc) = pp(2);
                
                % sfm
                %                 nansubj = isnan(iaf6(:, nc));
                %                 [rr, pp] = corrcoef(iaf6(~nansubj, nc), sfmdurmed(~nansubj));
                %                 rpchanrivsfm(2, 1, nc) = rr(2); rpchanrivsfm(2, 2, nc) = pp(2);
            end
            if ~chanset
                figure
                mval = squeeze(rpchanrivsfm(1, 1, :));
                if abs(min(mval))>abs(max(mval))
                    mval = -mval;
                    mmval= max(mval);
                else
                    mmval=max(mval);
                end
                subplot(121), topoplot(mval, chanloc34, ...
                    'maplimits', [0 mmval]);
                title('riv correlation coeff'), colorbar
                pval = -log(squeeze(rpchanrivsfm(1, 2, :)));
                if max(pval)>3
                    subplot(122), topoplot(pval, chanloc34, ...
                        'maplimits', [3 max(pval)]);
                end
                title('riv -log(pval)'), colorbar
            end
            
            %             mval = squeeze(rpchanrivsfm(2, 1, :));
            %             if abs(min(mval))>abs(max(mval))
            %                 mval = -mval;
            %                 mmval= max(mval);
            %             else
            %                 mmval=max(mval);
            %             end
            %             subplot(223), topoplot(mval, chanloc34, ...
            %                 'maplimits', [0 mmval]);
            %             title('sfm correlation coeff'), colorbar
            %             pval = -log(squeeze(rpchanrivsfm(2, 2, :)));
            %             if max(pval)>3
            %                 subplot(224), topoplot(pval, chanloc34, ...
            %                     'maplimits', [3 max(pval)]);
            %             end
            %             title('sfm -log(pval)'), colorbar
            
            varargout{1} = rivdurmed;
            varargout{2} = subjcoderiv;
            varargout{3} = iaf3;
            
            varargout{4} = sfmdurmed;
            varargout{5} = subjcodesfm;
            %             varargout{6} = iaf4;
            varargout{7} = ~noalphasubjriv;
            varargout{8} = iaf;
            
        case 3
            % calculate the smoothed peak
            
            iaf2 = mean(iaf, 1);
            iaf3 = iaf2(valriv2); % select from the valid rivalry subjects
            mx_alpha2 = mean(mx_alpha, 1);
            mx_alpha3 = mx_alpha2(valriv2);
            
            subjcoderiv2 = {exp_rivind.data(valsubj).subj_code};
            
            %remove the subjects who don't have a good alpha
            noalphasubjriv2 = false(1, numel(iaf2));
            for nos = 1:numel(noalpha)
                noalphasubjriv2(strcmp(subjcoderiv2, noalpha(nos))) = true;
            end
            noalphasubjriv = noalphasubjriv2(valriv2);
            rivdurmed2 = rivdurmed;
            rivdurmed = rivdurmed(:, ~noalphasubjriv);
            mixdur = mixdur(~noalphasubjriv);
            
            if pre_post_subj==3
                nsess = [exp_rivind.data.nsess];
                nsess = nsess(valsubj);
                nsess = nsess(valriv2);
                nsess = nsess(~noalphasubjriv);
                sprintf('n-prelim = %g, n-test = %g', sum(nsess==2), sum(nsess==1))
            end
            iaf3 = iaf3(~noalphasubjriv);
            mx_alpha3 = mx_alpha3(~noalphasubjriv);
            
            subjcoderiv = subjcoderiv(~noalphasubjriv);
            
            %             iaf4 = iaf2(valriv2);
            if size(rivdurmed, 1)==2, rivdurmed = rivdurmed(sessnum, :); end
            
            %%% plot peak alpha frequency vs. rivalry
            %             figure
            %             %             subplot(121),
            %             scatter(iaf3, rivdurmed)
            %             [rr, pp] = corrcoef(iaf3, rivdurmed);
            %             cc = sprintf('iaf riv corr %.3g, p = %.3g', rr(2), pp(2));
            %             title(cc)
            %             text(iaf3, rivdurmed+.05, subjcoderiv)
            %             %             p = polyfit(iaf3, rivdurmed, 1);
            %             p = polyfit(iaf3, rivdurmed, 1);
            %             sfit = polyval(p, iaf3);
            %             hold on, plot(iaf3, sfit, 'LineWidth', 3)
            %             %             max([rivdurmed; sfit]')
            %             %             min([rivdurmed; sfit]')
            %             if rem_mh_subj
            %                 millness_riv = millness_riv(~noalphasubjriv);
            %                 fam_millness_riv = fam_millness_riv(~noalphasubjriv);
            %                 scatter(iaf3(logical(fam_millness_riv)), rivdurmed(logical(fam_millness_riv)), 'g')
            %                 scatter(iaf3(logical(millness_riv)), rivdurmed(logical(millness_riv)), 'r')
            %             end
            
            
            %             %%% plot peak alpha frequency vs. mixed duration
            %             figure
            %             scatter(iaf3(mixdur<10), mixdur(mixdur<10))
            %             [rr, pp] = corrcoef(iaf3(mixdur<10), mixdur(mixdur<10));
            %             cc = sprintf('iaf mix riv corr %.3g, p = %.3g', rr(2), pp(2));
            %             title(cc)
            %             text(iaf3, mixdur+.05, subjcoderiv)
            %             %             p = polyfit(iaf3, rivdurmed, 1);
            %             p = polyfit(iaf3, mixdur, 1);
            %             sfit = polyval(p, iaf3);
            %             hold on, plot(iaf3, sfit, 'LineWidth', 3)
            
            
            %             %%% plot peak alpha amplitude vs. rivalry
            %             figure
            %             %             subplot(121),
            %             scatter(rivdurmed, mx_alpha3)
            %             [rr, pp] = corrcoef(rivdurmed, mx_alpha3);
            %             cc = sprintf('max alpha amp riv corr %.3g, p = %.3g', rr(2), pp(2));
            %             title(cc)
            %             text(rivdurmed, mx_alpha3+.05, subjcoderiv)
            %             %             p = polyfit(iaf3, rivdurmed, 1);
            %             p = polyfit(rivdurmed, mx_alpha3, 1);
            %             sfit = polyval(p, rivdurmed);
            %             hold on, plot(rivdurmed, sfit, 'LineWidth', 3)
            %             max([rivdurmed; sfit]')
            %             min([rivdurmed; sfit]')

            %             %%% plot alpha amplitude around peak frequency vs. rivalry
            malpharoundpeak = mean(alpha_amp_chan(valriv2, :), 2);
            %             figure
            %             %             subplot(121),
            %             subjcoderiv2 = subjcoderiv2(valriv2);
            %             scatter(malpharoundpeak, rivdurmed2)
            %             [rr, pp] = corrcoef(rivdurmed2, malpharoundpeak);
            %             cc = sprintf('alpha amp around peak freq riv corr %.3g, p = %.3g', rr(2), pp(2));
            %             title(cc)
            %             text(malpharoundpeak, rivdurmed2+.05, subjcoderiv2)
            %             p = polyfit(malpharoundpeak, rivdurmed2', 1);
            %             sfit = polyval(p, malpharoundpeak);
            %             hold on, plot(malpharoundpeak, sfit, 'LineWidth', 3)
            
            
            %             max([rivdurmed2; sfit]')
            %             min([rivdurmed2; sfit]')
            
            %             figure
            %             iaf4 = iaf2(valsfm2);
            %             noalphasubjsfm = false(1, numel(iaf4));
            %             for nos = 1:numel(noalpha)
            %                 noalphasubjsfm(strcmp(subjcodesfm, noalpha(nos))) = true;
            %             end
            %             sfmdurmed = sfmdurmed(:, ~noalphasubjsfm);
            %             iaf4 = iaf4(~noalphasubjsfm);
            %             subjcodesfm = subjcodesfm(~noalphasubjsfm);
            %
            %             if size(sfmdurmed, 1)==2, sfmdurmed = sfmdurmed(sessnum, :); end
            %             subplot(122), scatter(sfmdurmed, iaf4)
            %             [rr, pp] = corrcoef(sfmdurmed, iaf4);
            %             cc = sprintf('iaf sfm corr %.3g, p = %.3g', rr(2), pp(2));
            %             title(cc)
            %             text(sfmdurmed, iaf4+.025, subjcodesfm)
            %             p = polyfit(iaf4, sfmdurmed, 1);
            %             sfit = polyval(p, iaf4);
            %             hold on, plot(sfit, iaf4, 'LineWidth', 3)
            
            chanlabels = {chanlocs.labels};
            chanloc34 = chanlocs(~(strcmp(chanlabels, ...
                {'VEOG'}) | strcmp(chanlabels, {'HEOG'})));
            iaf5 = iaf_chan(valriv2, :); iaf5 = iaf5(~noalphasubjriv, :);
            %             iaf6 = iaf_chan(valsfm2, :); iaf6 = iaf6(~noalphasubjsfm, :);
            
            if ~chanset
                rpchanrivsfm = zeros(2, 2, nchanan);
                figure, subplot(ceil(sqrt(nchanan)), ceil(sqrt(nchanan)), 1)
                for nc = 1:nchanan
                    % riv
                    nansubj = isnan(iaf5(:, nc)) | ~iaf5(:, nc);
                    %                     if corr_type
                    %                         [rr, pp] = corrcoef(iaf5(~nansubj, nc), rivdurmed(~nansubj));
                    %                         rpchanrivsfm(1, 1, nc) = rr(2); rpchanrivsfm(1, 2, nc) = pp(2);
                    %                     else
                    [rr, pp] = corr(iaf5(~nansubj, nc), ...
                        rivdurmed(~nansubj)', 'Type', corr_type, ...
                        'tail', 'both');
                    rpchanrivsfm(1, 1, nc) = rr; rpchanrivsfm(1, 2, nc) = pp;
                    %                     end
                    
                    subplot(ceil(sqrt(nchanan)), ceil(sqrt(nchanan)), nc)
                    scatter(rivdurmed(~nansubj), iaf5(~nansubj, nc))
                    text(rivdurmed(~nansubj), iaf5(~nansubj, nc), subjcoderiv(~nansubj))
                    ax=axis; ax(3:4) = [6.5 13.5]; axis(ax);
                    title([chanlabels{nc} 'rr=' num2str(rr) ' pp=' num2str(pp)])
                    
                    %                     % sfm
                    %                     nansubj = isnan(iaf6(:, nc));
                    %                     [rr, pp] = corrcoef(iaf6(~nansubj, nc), sfmdurmed(~nansubj));
                    %                     rpchanrivsfm(2, 1, nc) = rr(2); rpchanrivsfm(2, 2, nc) = pp(2);
                end
                mval = squeeze(rpchanrivsfm(1, 1, :));
                if abs(min(mval))>abs(max(mval))
                    mval = -mval;
                    mmval= max(mval);
                else
                    mmval=max(mval);
                end
                figure, topoplot(mval, chanloc34, ...
                    'maplimits', [0 mmval]);
                title('riv correlation coeff'), colorbar
                pval = -log(squeeze(rpchanrivsfm(1, 2, :)));
                if pre_post_subj==1
                    pthresh = 0.05;
                else
                    pthresh = .00147;
                end
                figure
                if max(pval)>-log10(pthresh)
                    %                     hold on, topoplot(pval, chanloc34, ...
                    %                         'maplimits', [-log(.05) max(pval)]);
                    topoplot(squeeze(rpchanrivsfm(1, 2, :))<pthresh, chanloc34, ...
                        'maplimits', [0 1]);
                end
                title('riv -log(pval)'), colorbar
                %                 save('indiff_corr_data_s3_fix', 'iaf5', 'nansubj', 'rivdurmed', 'chanloc34', 'valriv2', 'noalphasubjriv2')
                
                mac = alpha_amp_chan(valriv2, :); 
                mac = mac(~noalphasubjriv, :);
                figure('Name', 'amplitude around peak alpha')
                if strcmp(runtype, 'fix')
                    topoplot(mean(mac, 1), chanloc34, ...
                        'maplimits', [.3 1.8]);  colorbar
                elseif strcmp(runtype, 'rest')
                    topoplot(mean(mac, 1), chanloc34, ...
                        'maplimits', [1 5]);  colorbar
                else
                    topoplot(mean(mac, 1), chanloc34, ...
                        'maplimits', 'minmax');  colorbar
                end
                
                rpchanrivsfm = zeros(2, 2, nchanan);
                figure, subplot(ceil(sqrt(nchanan)), ceil(sqrt(nchanan)), 1)
                for nc = 1:nchanan
                    % riv
                    %                     if corr_type
                    %                         [rr, pp] = corrcoef(alpha_amp_chan(valriv2, nc), rivdurmed2);
                    %                         rpchanrivsfm(1, 1, nc) = rr(2); rpchanrivsfm(1, 2, nc) = pp(2);
                    %                     else
                    [rr, pp] = corr(alpha_amp_chan(valriv2, nc), ...
                        rivdurmed2', 'Type', corr_type, 'tail', 'both');
                    rpchanrivsfm(1, 1, nc) = rr; rpchanrivsfm(1, 2, nc) = pp;
                    %                     end
                    subplot(ceil(sqrt(nchanan)), ceil(sqrt(nchanan)), nc)
                    scatter(rivdurmed2, alpha_amp_chan(valriv2, nc))
                    title([chanlabels{nc} 'rr=' num2str(rr) ' pp=' num2str(pp)])
                end
                figure
                mval = squeeze(rpchanrivsfm(1, 1, :));
                if abs(min(mval))>abs(max(mval))
                    mval = -mval;
                    mmval= max(mval);
                else
                    mmval=max(mval);
                end
                subplot(121), topoplot(mval, chanloc34, ...
                    'maplimits', [0 mmval]);
                title('riv correlation coeff'), colorbar
                pval = -log(squeeze(rpchanrivsfm(1, 2, :)));
                pthresh = 3;
                if max(pval)>pthresh
                    subplot(122), topoplot(pval, chanloc34, ...
                        'maplimits', [pthresh max(pval)]);
                end
                title('riv -log(pval)'), colorbar;
                %                 mval = squeeze(rpchanrivsfm(2, 1, :));
                %                 if abs(min(mval))>abs(max(mval))
                %                     mval = -mval;
                %                     mmval= max(mval);
                %                 else
                %                     mmval=max(mval);
                %                 end
                %                 subplot(223), topoplot(mval, chanloc34, ...
                %                     'maplimits', [0 mmval]);
                %                 title('sfm correlation coeff'), colorbar
                %                 pval = -log(squeeze(rpchanrivsfm(2, 2, :)));
                %                 if max(pval)>3
                %                     subplot(224), topoplot(pval, chanloc34, ...
                %                         'maplimits', [3 max(pval)]);
                %                 end
                %                 title('sfm -log(pval)'), colorbar
            end
            
            varargout{1} = rivdurmed;
            varargout{2} = subjcoderiv;
            varargout{3} = iaf3;
            
            %             varargout{4} = sfmdurmed;
            %             varargout{5} = subjcodesfm;
            %             varargout{6} = iaf4;
            
            varargout{7} = ~noalphasubjriv2;
            varargout{8} = iaf;
            varargout{9} = malpharoundpeak; % return the alpha amplitude around peak frequency
            varargout{10} = valriv2; % return the alpha amplitude around peak frequency
            varargout{11} = rivdurmed2; % all subjects' rivalry
            varargout{12} = subjcoderiv2; % all subject codes
    end
end

end
