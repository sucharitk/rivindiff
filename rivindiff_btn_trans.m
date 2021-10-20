function varargout = rivindiff_btn_trans(exp_rivind, runtype, plot_type, chanset, ...
    elect_flag, analysis_filestart, pre_post_subj, flick_freqs, dur_before_after)
%
%
%

cd(exp_rivind.session_dir)

nsubj = exp_rivind.nsubj;
eegdatapath = exp_rivind.eeg_datadir;

subj_ii = 0;
sess_ii = 0;
if ~exist('runtype', 'var'), runtype = 'riv'; end

trigs = exp_rivind.trigs.(runtype);

% flick_freqs = [exp_rivind.flick_freqs diff(exp_rivind.flick_freqs)];
if ~exist('flick_freqs', 'var')||isempty(flick_freqs)
    flick_freqs = exp_rivind.flick_freqs;
end

fm = .1;
nff = numel(flick_freqs);

valsubj = false(1, nsubj);

indiv_plots = any(plot_type==1);
avg_chans = true;

ntrigs = 2;

%%% envelope parameters
trans_to_plot = [2.5 2.5];
trans_to_plot = [3 3];
smoothdur = [2.5];

mx_percept_btn = 3;
snr_cutoff = 3;

if ~exist('dur_before_after', 'var')||isempty(dur_before_after)
    if strcmp(runtype, 'riv')
        dur_before_after = [ 1 1]; %use for rivalry
    else
        % dur_before_after = [1 1];
        % dur_before_after = [3.5 2 2 3.5];
        dur_before_after = [2.5 2.5];
        %     dur_before_after = [1 .5 .5 1]; %use for rivalry
    end
end

switch elect_flag
    case {1, 2, 3, 5}
        peak_detect_duration = [-.6 .6];
    case {4, 6, 7}
        peak_detect_duration = [-.4 .4];
end

elsel_prctile = 85;

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
        chans_to_an = exp_rivind.po_electrodes;
    case 8
        chans_to_an = exp_rivind.parietoccip_electrodes;
    case 9
        chans_to_an = exp_rivind.parietonly_electrodes;
    case 10
        chans_to_an = exp_rivind.occiponly_electrodes;
end

compstr = ['^' analysis_filestart];

sessnum = 1;
if pre_post_subj == 2
    sessnum = 2;
end

if any(plot_type==8)
    all_trans_subj = cell(1, 2);
end
nocountermodsubj = {};
countermodbetween = [-2 2];

for ns = 1:nsubj
    % subject loop
    
    subj_data = exp_rivind.data(ns);
    
    if subj_data.subj_valid && subj_data.subj_valid_eeg
        
        if  (pre_post_subj && subj_data.nsess>1) || ...
                (~pre_post_subj && subj_data.nsess==1) || ...
                pre_post_subj>=3
            % for subjects that have two sessions
            
            subj_ii = subj_ii + 1;
            valsubj(ns) = true;
            
            if any(plot_type==6) && ~mod(subj_ii-1, 16)
                f6 = figure;
                %                 subplot(441)
            end
            
            for snum = sessnum:sessnum
                sess_ii = sess_ii + 1;
                subj_snr = true;
                
                sessdir = fullfile(exp_rivind.session_dir, subj_data.session_name{snum});
                if isfolder(fullfile(sessdir, eegdatapath))
                    
                    cd(sessdir)
                    
                    ff = dir(eegdatapath);
                    %                     ff2 = find(~isemptycell(strfind({ff.name}, 'csd_')), 1, 'last');
                    ff2 = find(~isemptycell(regexp({ff.name}, compstr)), 1, 'last');
                    if ~isempty(ff2)
                        fname = ff(ff2).name;
                        
                        if any(plot_type == 3)
                            figure('Name', subj_data.subj_code)
                        end
                        for numtrigs = 1:ntrigs
                            
                            resfname = fullfile('Epochs', sprintf('%s_%s_%s', ...
                                fname(1:end-4), runtype, trigs{numtrigs}));
                            
                            load(resfname)
                            
                            if ~exist('NFFT', 'var')
                                NFFT = 2^16;
                            end
                            
                            if chanset
                                chanloc_inds = get_channels_from_labels({chanlocs.labels}, chans_to_an);
                                nbchan = sum(chanloc_inds);
                            else
                                % average over all available channels
                                nbchan = size(epochs{1}, 1);
                                chanloc_inds = true(1, nbchan);
                            end
                            chanlifinds = find(chanloc_inds);
                            
                            nepochs = numel(epochs);
                            
                            if elect_flag
                                %%% select electrodes for analysis
                                clear freq_amps freq_phase
                                for ne = 1:nepochs
                                    freq_amps(ne, :, :, :) = FFT_at_freq(...
                                        epochs{ne}, fs, ...
                                        NFFT, flick_freqs, fm, 1);
                                end
                                freq_amps = squeeze(mean(freq_amps, 1));
                                if numel(size(freq_amps))==2
                                    % if only one frequency
                                    freq_amps2 = freq_amps;
                                    clear freq_amps
                                    freq_amps(1, :, :) = freq_amps2;
                                    [~, mxind] = max(squeeze(mean(freq_amps, 2)));
                                else
                                    [~, mxind] = max(squeeze(mean(freq_amps, 2)), [], 2);
                                end
                                
                                
                                for flf = 1:nff
                                    noiseamp = mean(freq_amps(flf, :, ...
                                        [1:10 end-10:end]), 3);
                                    sigamp = freq_amps(flf, :, mxind(flf));
                                    snr = sigamp./noiseamp;
                                    %                                 sigph = freq_phase(flf, :, mxind(flf));
                                    %                                 snr = snr.*sigph;
                                    
                                    switch elect_flag
                                        case 1
                                            [~, ind] = max(snr(chanloc_inds));
                                            electrind{numtrigs, flf} = chanlifinds(ind);
                                            electrind2{numtrigs, flf} = ind;
                                        case {2, 5, 6, 7}
                                            ind = snr(chanloc_inds)>prctile(snr, elsel_prctile);
                                            electrind{numtrigs, flf} = chanlifinds(ind);
                                            electrind2{numtrigs, flf} = find(ind);
                                        case {3, 4}
                                            resp_elect = snr(chanloc_inds)>snr_cutoff;
                                            if ~sum(resp_elect), subj_snr = false; end
                                            electrind{numtrigs, flf} = ...
                                                chanlifinds(resp_elect);
                                            electrind2{numtrigs, flf} = resp_elect;
                                    end
                                    
                                    if any(plot_type == 3)
                                        subplot(ntrigs, nff, (numtrigs-1)*nff+flf)
                                        topoplot(snr, chanlocs, 'maplimits', [0 max(snr)]);
                                        %                                     title(chanlocs(chanlifinds(ind)).labels)
                                        colorbar
                                    end
                                end
                                
                            else
                                for flf = 1:nff
                                    electrind{numtrigs, flf} = find(chanloc_inds);
                                end
                            end
                            
                            if numtrigs==1
                                epochs1 = epochs;
                                latencies1 = latencies;
                                artif_mask1 = artif_mask;
                            end

                        end
                        
                    else
                        nosess = true;
                    end
                else
                    nosess = true;
                end
            end
            
            % sort runs according to onset order
            nruntype{1} = 1:numel(epochs1); nruntype{2} = nruntype{1}(end)+(1:numel(epochs)); % put the two types of triggered runs in two separate cells in case some subjects have an extra one of one or the other
            epochs = [epochs1 epochs];
            latencies = [latencies1 latencies];
            artif_mask = [artif_mask1 artif_mask];
            [~,srt_eps] = sort(latencies(1, :));
            [~,srt_inds] = sort(srt_eps); % reverse sort and do it on the behavioral file
            
            switch runtype
                % load behavior
                case 'riv'
                    fname = sprintf('BR_Rivalry_%s_%s', subj_data.subj_code,...
                        subj_data.session_date{snum});
                    load(fullfile('Behavior', fname))
                    subj_valid = subj_data.subj_valid_riv;
                    results = results(srt_inds);
                    valid_runs = true(1, numel(results));
                    valid_runs(subj_data.invalid_riv{1}) = false;
                    subj_snr = subj_snr & subj_data.subj_valid_riv;
                    if subj_snr
                        rivdur = subj_data.riv_durations.mean_all(snum);
                    end
                case 'sfm'
                    fname = sprintf('StructureFromMotion_%s_%s', subj_data.subj_code,...
                        subj_data.session_date{snum});
                    load(fullfile('Behavior', fname))
                    subj_valid = subj_data.subj_valid_sfm;
                    results = results(srt_inds);
                    valid_runs = true(1, numel(results));
                    valid_runs(subj_data.invalid_sfm{1}) = false;
                    subj_snr = subj_snr & subj_data.subj_valid_sfm;
                    if subj_snr
                        rivdur = subj_data.sfm_durations.mean_all(snum);
                    end
            end
            if subj_valid && subj_snr
                valid_runs = valid_runs(srt_inds);
                nruntype{1}(~valid_runs(nruntype{1})) = [];
                nruntype{2}(~valid_runs(nruntype{2})) = [];
                
                plot_cols = 'bgrk';
                
                switch elect_flag
                    case {0, 1, 2}
                        
                        btn_trans{1} = btnpress_alternations(epochs(nruntype{1}), fs, ...
                            results(nruntype{1}), flick_freqs, ...
                            electrind(1, :), 1:2, artif_mask(nruntype{1}), ...
                            indiv_plots, dur_before_after, smoothdur, [], [], avg_chans, [], ...
                            trans_to_plot);
                        btn_trans{2} = btnpress_alternations(epochs(nruntype{2}), fs, ...
                            results(nruntype{2}), flick_freqs, electrind(2, :), 1:2, ...
                            artif_mask(nruntype{2}), indiv_plots, ...
                            dur_before_after, smoothdur, [], [], avg_chans, [], ...
                            trans_to_plot);
                        
                        %%% btn_trans: first cell dimension = 2 different
                        %%% types of bistable runs counterbalanced, second
                        %%% dim = n different types of transition (e.g.,
                        %%% dom1->dom2, dom2->dom1 etc., dim 3 = n
                        %%% different types of frequency envelopes
                        %%% returned, 3 dims of array = [ntrials
                        %%% nelectrodes tpoints]
                        
                        sizind = sum(trans_to_plot*fs)+1;
                        tt3 = linspace(-trans_to_plot(1), trans_to_plot(2), sizind);
                        %                         selt = tt3>peak_detect_duration(1) & ...
                        %                             tt3<peak_detect_duration(2);
                        
                        if any(plot_type==2)
                            figure('Name', sprintf('%s: %.3g', subj_data.subj_code, ...
                                rivdur))
                        end
                        
                        all_trans = cell(1, 2);
                        for tt = 1:2 % counterbalanced runs
                            bt1 = btn_trans{tt}{1};
                            if ~isempty(bt1{1})
                                nonan{1} = ~isnan(squeeze(mean(btn_trans{tt}{1}{1}(:, 1, :), 3)));
                            else
                                nonan{1} = [];
                            end
                            
                            bt2 = btn_trans{tt}{2};
                            if ~isempty(bt2{1})
                                nonan{2} = ~isnan(squeeze(mean(btn_trans{tt}{2}{1}(:, 1, :), 3)));
                            else
                                nonan{2} = [];
                            end
                            
                            for nf = 1:2
                                nn1 = nonan{1};
                                if ~isempty(nn1)
                                    co1 = btn_trans{tt}{1}{nf}(nn1, :, :);
                                    aft1 = squeeze(mean(mean(co1, 1), 2)); % A->B
                                    co1f = squeeze(mean(co1, 2));
                                    if size(co1f, 2)==1, co1f = co1f'; end % if there is only a single trial
                                    
                                    co2 = btn_trans{tt}{2}{nf}(nonan{2}, :, :);
                                    aft2 = squeeze(mean(mean(co2, 1), 2)); % A->B
                                    co2f = squeeze(mean(co2, 2));
                                    if size(co2f, 2)==1, co2f = co2f'; end
                                   
                                    if any(plot_type==2)
                                        subplot(2, 2, (2*tt-1)), hold on
                                        plot(tt3, aft1, ...
                                            plot_cols(nf))
                                        subplot(2, 2, 2*tt), hold on
                                        plot(tt3, aft2, ...
                                            plot_cols(nf))
                                    end
                                end
                            end
                        end
                        
                        clear all_trans
                        all_trans{1} = squeeze(mean([btn_trans{1}{1}{1}; ...
                            btn_trans{1}{2}{2}; btn_trans{2}{1}{2}; ...
                            btn_trans{2}{2}{1}], 2));
                        all_trans{2} = squeeze(mean([btn_trans{1}{1}{2}; ...
                            btn_trans{1}{2}{1}; btn_trans{2}{1}{1}; ...
                            btn_trans{2}{2}{2}], 2));
                        if subj_data.reverse_rivbuttons && strcmp(runtype, 'riv')
                            atemp = all_trans{1};
                            all_trans{1} = all_trans{2}; all_trans{2} = atemp;
                        end
                       
                        
                    case {3, 4}
                        
                        btn_trans{1} = btnpress_alternations(epochs(nruntype{1}), fs, ...
                            results(nruntype{1}), flick_freqs, ...
                            [], 1:2, artif_mask(nruntype{1}), indiv_plots, ...
                            dur_before_after);
                        btn_trans{2} = btnpress_alternations(epochs(nruntype{2}), fs, ...
                            results(nruntype{2}), ...
                            flick_freqs, [], 1:2, ...
                            artif_mask(nruntype{2}), indiv_plots, dur_before_after);
                        
                        sizind = max(size(btn_trans{1}{1}{1}));
                        tt3 = linspace(-trans_to_plot(1), trans_to_plot(2), sizind);
                        selt = tt3>peak_detect_duration(1) & ...
                            tt3<peak_detect_duration(2);
                        
                        if any(plot_type==4)
                            %%% plot transitions for each electrode to see
                            %%% if there are both up and down going ones
                            figure('Name', sprintf('%s: %.3g', subj_data.subj_code, ...
                                rivdur))
                            
                            for tt = 1:2
                                nonan{1} = ~isnan(squeeze(mean(btn_trans{tt}{1}{1}(:, 1, :), 3)));
                                nonan{2} = ~isnan(squeeze(mean(btn_trans{tt}{2}{1}(:, 1, :), 3)));
                                for nf = 1:2
                                    co1 = btn_trans{tt}{nf}{1}(nonan{nf}, :, :);
                                    aft1 = squeeze(mean(co1, 1)); % A->B
                                    
                                    co2 = btn_trans{tt}{nf}{2}(nonan{nf}, :, :);
                                    aft2 = squeeze(mean(co2, 1)); % A->Bsubplot(2, 2, (2*tt-1)), hold on,
                                    
                                    subplot(2, 2, 2*tt-1), hold on
                                    plot(tt3, aft1, ...
                                        plot_cols(nf))
                                    %                                     plot(tt3, mean(aft1(electrind2{tt, nf}, :), 1),...
                                    %                                         plot_cols(nf+2), 'LineWidth', 2)
                                    subplot(2, 2, 2*tt), hold on
                                    plot(tt3, aft2, ...
                                        plot_cols(nf))
                                    %                                     plot(tt3, mean(aft2(electrind2{tt, nf}, :), 1),...
                                    %                                         plot_cols(nf+2), 'LineWidth', 2)
                                end
                            end
                            
                        end
                        
                        if any(plot_type==5) && ~any(plot_type==4)
                            %%% plot the peak selected electrodes
                            figure('Name', sprintf('%s: %.3g', subj_data.subj_code, ...
                                rivdur))
                        end
                        
                        crossing_type = [1 2 2 1; 2 1 1 2];
                        if subj_data.reverse_rivbuttons
                            crossing_type = [2 1 1 2; 1 2 2 1];
                        end
                        all_trans = cell(1, 2);
                        for tt = 1:2
                            % two types of triggers with frequencies switched
                            nonan{1} = ~isnan(squeeze(mean(btn_trans{tt}{1}{1}(:, 1, :), 3)));
                            nonan{2} = ~isnan(squeeze(mean(btn_trans{tt}{2}{1}(:, 1, :), 3)));
                            
                            for nf = 1:2
                                co1 = btn_trans{tt}{nf}{1}(nonan{nf}, :, :);
                                af1 = squeeze(mean(co1, 1)); % A->B
                                zc_el1 = zero_crossing(af1, selt, ...
                                    crossing_type(tt, 2*nf-1), elect_flag-3); % get zero crossing points
                                
                                co1f = squeeze(co1(:, zc_el1, :)); if size(co1f, 2)==1, co1f = co1f'; end % if there is only a single trial
                                
                                all_trans{crossing_type(tt, 2*nf-1)} = ...
                                    [all_trans{crossing_type(tt, 2*nf-1)}; ...
                                    co1f];
                                
                                co2 = btn_trans{tt}{nf}{2}(nonan{nf}, :, :);
                                af2 = squeeze(mean(co2, 1)); % A->
                                zc_el2 = zero_crossing(af2, selt, ...
                                    crossing_type(tt, 2*nf), elect_flag-3); % get zero crossing points
                                
                                co2f = squeeze(co2(:, zc_el2, :)); if size(co2f, 2)==1, co2f = co2f'; end % if there is only a single trial
                                
                                all_trans{crossing_type(tt, 2*nf)} = ...
                                    [all_trans{crossing_type(tt, 2*nf)}; ...
                                    co2f];
                                
                                if any(plot_type==5)
                                    subplot(2, 2, (2*tt-1)), hold on
                                    plot(tt3, af1(zc_el1, :), ...
                                        plot_cols(nf), 'LineWidth', 2)
                                    subplot(2, 2, 2*tt), hold on
                                    plot(tt3, af2(zc_el2, :), ...
                                        plot_cols(nf), 'LineWidth', 2)
                                end
                            end
                        end
                        
                    case {5, 6}
                        btn_trans{1} = btnpress_alternations(epochs(nruntype{1}), fs, ...
                            results(nruntype{1}), flick_freqs, ...
                            electrind(1, :), 1:2, artif_mask(nruntype{1}), indiv_plots, ...
                            dur_before_after);
                        btn_trans{2} = btnpress_alternations(epochs(nruntype{2}), ...
                            fs, results(nruntype{2}), ...
                            flick_freqs, electrind(2, :), 1:2, ...
                            artif_mask(nruntype{2}), indiv_plots, dur_before_after);
                        
                        sizind = max(size(btn_trans{1}{1}{1}));
                        tt3 = linspace(-trans_to_plot(1), trans_to_plot(2), sizind);
                        selt = tt3>peak_detect_duration(1) & ...
                            tt3<peak_detect_duration(2);
                        
                        if any(plot_type==4)
                            %%% plot transitions for each electrode to see
                            %%% if there are both up and down going ones
                            figure('Name', sprintf('%s: %.3g', subj_data.subj_code, ...
                                rivdur))
                            
                            for tt = 1:2
                                nonan{1} = ~isnan(squeeze(mean(btn_trans{tt}{1}{1}(:, 1, :), 3)));
                                nonan{2} = ~isnan(squeeze(mean(btn_trans{tt}{2}{1}(:, 1, :), 3)));
                                for nf = 1:2
                                    co1 = btn_trans{tt}{nf}{1}(nonan{nf}, :, :);
                                    aft1 = squeeze(mean(co1, 1)); % A->B
                                    
                                    co2 = btn_trans{tt}{nf}{2}(nonan{nf}, :, :);
                                    aft2 = squeeze(mean(co2, 1)); % A->Bsubplot(2, 2, (2*tt-1)), hold on,
                                    
                                    subplot(2, 2, 2*tt-1), hold on
                                    plot(tt3, aft1, ...
                                        plot_cols(nf))
                                    %                                     plot(tt3, mean(aft1(electrind2{tt, nf}, :), 1),...
                                    %                                         plot_cols(nf+2), 'LineWidth', 2)
                                    subplot(2, 2, 2*tt), hold on
                                    plot(tt3, aft2, ...
                                        plot_cols(nf))
                                    %                                     plot(tt3, mean(aft2(electrind2{tt, nf}, :), 1),...
                                    %                                         plot_cols(nf+2), 'LineWidth', 2)
                                end
                            end
                        end
                        
                        if any(plot_type==5) && ~any(plot_type==4)
                            %%% plot the peak selected electrodes
                            figure('Name', sprintf('%s: %.3g', subj_data.subj_code, ...
                                rivdur))
                        end
                        
                        crossing_type = [1 2 2 1; 2 1 1 2];
                        if subj_data.reverse_rivbuttons
                            crossing_type = [2 1 1 2; 1 2 2 1];
                        end
                        all_trans = cell(1, 2);
                        for tt = 1:2
                            % two types of triggers with frequencies switched
                            nonan{1} = ~isnan(squeeze(mean(btn_trans{tt}{1}{1}(:, 1, :), 3)));
                            nonan{2} = ~isnan(squeeze(mean(btn_trans{tt}{2}{1}(:, 1, :), 3)));
                            
                            for nf = 1:2
                                co1 = btn_trans{tt}{nf}{1}(nonan{nf}, :, :);
                                af1 = squeeze(mean(co1, 1)); % A->B
                                zc_el1 = zero_crossing(af1, selt, ...
                                    crossing_type(tt, 2*nf-1), elect_flag-3); % get zero crossing points
                                
                                co1f = squeeze(co1(:, zc_el1, :)); if size(co1f, 2)==1, co1f = co1f'; end % if there is only a single trial
                                
                                all_trans{crossing_type(tt, 2*nf-1)} = ...
                                    [all_trans{crossing_type(tt, 2*nf-1)}; ...
                                    co1f];
                                
                                co2 = btn_trans{tt}{nf}{2}(nonan{nf}, :, :);
                                af2 = squeeze(mean(co2, 1)); % A->
                                zc_el2 = zero_crossing(af2, selt, ...
                                    crossing_type(tt, 2*nf), elect_flag-3); % get zero crossing points
                                
                                co2f = squeeze(co2(:, zc_el2, :)); if size(co2f, 2)==1, co2f = co2f'; end % if there is only a single trial
                                
                                all_trans{crossing_type(tt, 2*nf)} = ...
                                    [all_trans{crossing_type(tt, 2*nf)}; ...
                                    co2f];
                                
                                if any(plot_type==5)
                                    subplot(2, 2, (2*tt-1)), hold on
                                    plot(tt3, af1(zc_el1, :), ...
                                        plot_cols(nf), 'LineWidth', 2)
                                    subplot(2, 2, 2*tt), hold on
                                    plot(tt3, af2(zc_el2, :), ...
                                        plot_cols(nf), 'LineWidth', 2)
                                end
                            end
                        end
                  
                end
                
                if any(plot_type==8)
                    all_trans_subj{1} = [all_trans_subj{1}; ...
                        nanmean(all_trans{1}, 1)];
                    all_trans_subj{2} = [all_trans_subj{2}; ...
                        nanmean(all_trans{2}, 1)];
                end
                if any(plot_type==6)
                    figure(f6)
                    subplot(4, 4, mod(subj_ii-1, 16)+1)
                    plot(tt3, nanmean(all_trans{1}, 1), ...
                        tt3, nanmean(all_trans{2}, 1))
                    axis([-trans_to_plot(1) trans_to_plot(2) -3 3])
                    title(sprintf('%s: %.3g', ...
                        subj_data.subj_code, rivdur))
                end

            else
                fprintf('Subject did not qualify\n')
            end
            
        end
    end
end

plotcols = 'brg';
if any(plot_type==8)
    nsubj = size(all_trans_subj{1}, 1);
    figure
        subplot(121), hold on
    for pp = 1:2
        subjstd = std(all_trans_subj{pp})/sqrt(nsubj);
        plot(tt3, mean(all_trans_subj{pp}, 1), plotcols(pp))
        shadedErrorBar(tt3, mean(all_trans_subj{pp}, 1), subjstd, ...
            plotcols(pp), true)
    end
        title('average of subjects (µv)')
        subplot(122), 
    stdsubj = nanstd([all_trans_subj{1} all_trans_subj{2}], [], 2);
    stdsubj = repmat(stdsubj, [1, numel(tt3)]);
    hold on
    for pp = 1:2
        zsc = all_trans_subj{pp}./stdsubj;
        subjstd = std(zsc)/sqrt(nsubj);
        plot(tt3, mean(zsc, 1), plotcols(pp))
        shadedErrorBar(tt3, nanmean(zsc, 1), subjstd, ...
            plotcols(pp), true)
        %     subplot(122), plot(tt3, mean(, 1), ...
        %         tt3, mean(all_trans_subj{2}./subjstd, 1, 'LineWidth', 2))
    end
    title('average of subjects z-score')
    
end


end

function zc_el = zero_crossing(aft1, selt, crossing_type, peak_flag)
% choose the electrode that either goes lowest or highest around the
% transition point

if ~peak_flag
    daft1 = diff(aft1, [], 2); % differentiate the curve
    switch crossing_type
        case 1
            %  up
            zdaft1 = daft1(:, 1:end-1)>0 & daft1(:, 2:end)<0;
            ou = zdaft1.*aft1(:, 1:end-2); % scale zero crossings with peak of curve
            [~, zc_el] = max(max(ou(:, selt(1:end-2)), [], 2)); % find the electrode with the highest peak
        case 2
            %  down
            zdaft1 = daft1(:, 1:end-1)<0 & daft1(:, 2:end)>0;
            ou = zdaft1.*aft1(:, 1:end-2); % scale zero crossings with peak of curve
            [~, zc_el] = min(min(ou(:, selt(1:end-2)), [], 2)); % find the electrode with the highest peak
    end
else
    maft = mean(aft1(:, selt), 2);
    switch crossing_type
        case 1
            %  up
            [~, zc_el] = max(maft); % find the electrode with the highest peak
        case 2
            %  down
            [~, zc_el] = min(maft); % find the electrode with the highest peak
    end
    
end
end