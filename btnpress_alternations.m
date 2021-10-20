function varargout = btnpress_alternations(epochs, fs, results, freqs_to_plot, ...
    selectrodes, trans_to_eval, artif_mask, do_plots, dur_before_after, ...
    smoothdur, durplot, highpass_filt, avg_all_chans, mixdur_to_exclude, ...
    trans_to_plot, raw_eeg_trans, smooth_env, eeg_or_env, mix_btn)
%
% function varargout = btnpress_alternations(epochs, fs, results, freqs_to_plot, ...
%     selectrodes, trans_to_eval, artif_mask, do_plots, dur_before_after, ...
%     smoothdur, durplot, highpass_filt, avg_all_chans, mixdur_to_exclude)
%


if ~exist('smoothdur', 'var')||isempty(smoothdur), smoothdur = .55555; end
if ~exist('durplot', 'var')||isempty(durplot), durplot = [1 119]; end
if ~exist('trans_to_plot', 'var')||isempty(trans_to_plot), trans_to_plot = [4 4]; end
if ~exist('mixdur_to_exclude', 'var')||isempty(mixdur_to_exclude), mixdur_to_exclude = 1; end
if ~exist('raw_eeg_trans', 'var')||isempty(raw_eeg_trans), raw_eeg_trans = 0; end
if ~exist('smooth_env', 'var')||isempty(smooth_env), smooth_env = true; end
if ~exist('mix_btn', 'var'), mix_btn = 3; end
if ~exist('eeg_or_env','var')||isempty(eeg_or_env), eeg_or_env = 1; end % eeg by default

nepochs = numel(epochs);
nsamps = 36000;

numtrans = numel(trans_to_eval);
do_behav = true;
numf = numel(freqs_to_plot);
btn_trans = cell(1, numtrans);
for ntrans = 1:numtrans
    btn_trans{ntrans} = cell(1, numf);
end
if nargout>1
    trans_durations = cell(1, numtrans);
end

% trans_num = zeros(1, numtrans);

if ~exist('trans_to_eval', 'var')
    trans_to_eval = 1:2; % 1:2 - dom2dom, 3:4-dom2mix, 5:6-mix2dom
    % trans_to_eval = 3:4; % 1:2 - dom2dom, 3:4-dom2mix, 5:6-mix2dom
    % trans_to_eval = 5:6; % 1:2 - dom2dom, 3:4-dom2mix, 5:6-mix2dom
end



if isempty(selectrodes)
    numelect = size(epochs{1}, 1);
    selectrodes = cell(1, numf);
    for nf = 1:numf
        selectrodes{nf} = 1:numelect;
    end
    return_allectrodes = true;
end

for rb = 1:nepochs
    nn = size(epochs{rb}, 2);
    tt = linspace(0, nn/fs, nn);
    
    if exist('artif_mask', 'var') && ~isempty(artif_mask)
        am = single(artif_mask{rb});
        am(~am) = NaN;
    else
        am = [];
    end
    
    %     for nf = 1:numf
    if eeg_or_env
        sm_env = get_envelope(epochs{rb}, fs, ...
            selectrodes, freqs_to_plot, highpass_filt, raw_eeg_trans, ...
            avg_all_chans, smoothdur, am, smooth_env);
    else
        sm_env = {epochs{rb}(1,:), epochs{rb}(2,:)};
    end

    if do_plots
        for nf = 1:numf
            dp = tt>=durplot(1) & tt<=durplot(2);
            subplot(nepochs, 1, rb)
            hold on
            plot(tt(dp), mean(sm_env{nf}(:, dp), 1), plot_cols(nf))
        end
    end
    
    %%% behavioral analysis
    if do_behav
        %     sn = scans(ns);
        psycho = results(rb).psycho;
        params = results(rb).params;
        scan_dur = params.scanDuration;
        if numel(scan_dur)>1, scan_dur = scan_dur(2); end
        
        [ttb, respb] = interp_to_next_response(psycho.tKeyPress, psycho.responseKey, ...
            0, scan_dur, nsamps, [], [mix_btn mixdur_to_exclude]);
        
        if do_plots
            hold on, plot(ttb, 2*(respb==1), 'b', ttb, 2*(respb==2), 'r')
        end
        
        
        
        if nargout==1
            btn_trans2 = get_transition_envelopes(sm_env, fs, respb, ttb, ...
                dur_before_after, trans_to_eval, trans_to_plot, raw_eeg_trans);
            for type_trans = 1:numtrans
                for nf = 1:numf
                    if ~isempty(btn_trans2{type_trans})

                        btn_trans{type_trans}{nf} = [btn_trans{type_trans}{nf}; ...
                            btn_trans2{type_trans}{nf}];
                    end
                end
            end
        elseif nargout==2
            
            %%% needs work
            [btn_trans2, trans_durations2] = get_transition_envelopes(sm_env, fs, respb, ttb, ...
                dur_before_after, trans_to_eval, trans_to_plot, raw_eeg_trans);
            for type_trans = 1:numtrans
                for nf = 1:numf
                    %                     if isempty(btn_trans{type_trans}{nf})
                    %                         btn_trans{type_trans}{nf} = btn_trans2{type_trans}{nf};
                    %                         trans_durations{type_trans}{nf} = trans_durations2{type_trans}{nf};
                    %                     else
                    %                         btn_trans{type_trans}{nf} = [btn_trans{type_trans}{nf}; ...
                    %                             btn_trans2{type_trans}{nf}];
                    %                         trans_durations{type_trans}{nf} = ...
                    %                             [trans_durations{type_trans}{nf}; trans_durations2{type_trans}{nf}];
                    %                     end
                    if ~isempty(btn_trans2{type_trans})
                        btn_trans{type_trans}{nf} = [btn_trans{type_trans}{nf}; ...
                            btn_trans2{type_trans}{nf}];
                        trans_durations{type_trans} = ...
                            [trans_durations{type_trans}; trans_durations2{type_trans}];
                    end
                    
                end
            end
        end
        

        
    end
end

if ~nargout
    sizind = max(size(btn_trans{1}));
    tt3 = linspace(-trans_to_plot(1), trans_to_plot(2), sizind);
    figure
    af1 = squeeze(nanmean(nanmean(btn_trans{1}, 1), 3));
    af2 = squeeze(nanmean(nanmean(btn_trans{2}, 1), 3));
    subplot(131), plot(tt3, af1')
    subplot(132), plot(tt3, af2')
    clear af3
    af3 = [btn_trans{1}(:, [1 2], :, :); btn_trans{2}(:, [2 1], :, :)];
    af3 = squeeze(nanmean(nanmean(af3, 1), 3));
    subplot(133), plot(tt3, af3')
end

varargout{1} = btn_trans;
if nargout>1
    varargout{2} = trans_durations;
end
end