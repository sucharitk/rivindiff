function varargout = get_transition_envelopes(sm_env, fs, respb, ttb, ...
    dur_before_after, trans_to_eval, trans_to_plot, raw_eeg_trans, ...
    align_ssvep_freqs)

if ~exist('align_ssvep_freqs', 'var')||isempty(align_ssvep_freqs), align_ssvep_freqs=0; end
numtrans = numel(trans_to_eval);
btn_trans = cell(1, numtrans); % numbeforeafter x numtranstypes
trans_durations = cell(1, numtrans);
trans_num = zeros(1, numtrans);

if align_ssvep_freqs==1
    durfreqs = 1./freqs_to_plot;
    sampdur = 1/fs;
end

if nargout==1
    trans_times = response_transition_points(respb, ttb, dur_before_after, [1 2 3]);
else
    [trans_times, trans_durs] = ...
        response_transition_points(respb, ttb, dur_before_after, [1 2 3]);
    
end

numf = numel(sm_env);
nn = numel(sm_env{1});
tt = linspace(0, nn/fs, nn);

for type_trans = 1:numtrans
    typeoftrans = trans_to_eval(type_trans);
    switch typeoftrans
        case 1
            % dom 1 to dom 2
            trans_times2 = trans_times{1, 2};
            if nargout>1
                trans_durs2 = trans_durs{1, 2};
            end
        case 2
            % dom 2 to dom 1
            trans_times2 = trans_times{2, 1};
            if nargout>1
                trans_durs2 = trans_durs{2, 1};
            end
        case 3
            % dom 1 to mix
            trans_times2 = trans_times{1, 3};
            if nargout>1
                trans_durs2 = trans_durs{1, 3};
            end
        case 4
            % dom 2 to mix
            trans_times2 = trans_times{2, 3};
            if nargout>1
                trans_durs2 = trans_durs{2, 2};
            end
        case 5
            % mix to dom 1
            trans_times2 = trans_times{3, 1};
            if nargout>1
                trans_durs2 = trans_durs{3, 1};
            end
        case 6
            % mix to dom 2
            trans_times2 = trans_times{3, 2};
            if nargout>1
                trans_durs2 = trans_durs{3, 2};
            end
    end
    nts = size(trans_times2, 2);
    tslen = size(sm_env{1}, 2);
    if nts>0
        for ttn = 1:nts
            [~, ttind] = min(abs(trans_times2(ttn)-tt));
            if raw_eeg_trans==1&&align_ssvep_freqs
                % shift the data by remainder of flicker frequency
                % cycles
                mnind = zeros(1, numf); mxind = mnind;
                for nf = 1:numf
                    rr = tt(ttind)/durfreqs(nf); % how many cycles of that frequency till time of transition
                    rr = rr-floor(rr); % get the remaining fracition of cycles
                    if rr<.5
                        % shift back
                        trr = rr*durfreqs(nf); % duration of the fraction
                        ttind2 = ttind-round(trr/sampdur); % shift however much of that fraction
                    else
                        % shift forward
                        trr = (1-rr)*durfreqs(nf); % duration of the fraction
                        ttind2 = ttind+round(trr/sampdur); % shift however much of that fraction
                    end
                    
                    mnind(nf) = ttind2-round(trans_to_plot(1)*fs);
                    mxind(nf) = ttind2+round(trans_to_plot(2)*fs);
                end
                if all(mnind>0) && all(mxind<tslen)
                    trans_num(type_trans) = trans_num(type_trans)+1;
                    %                         if return_allectrodes
                    for nf = 1:numf
                        btn_trans{type_trans}{nf}(trans_num(type_trans), :, :) ...
                            = sm_env{nf}(:, mnind(nf):mxind(nf));
                        if nargout>1
                            trans_durations{type_trans}(trans_num(type_trans), :) ...
                                = trans_durs2(ttn);
                        end
                    end
                end
            else
                mnind = ttind-round(trans_to_plot(1)*fs);
                mxind = ttind+round(trans_to_plot(2)*fs);
                if mnind>0 && mxind < tslen
                    trans_num(type_trans) = trans_num(type_trans)+1;
                    %                         if return_allectrodes
                    for nf = 1:numf
                        btn_trans{type_trans}{nf}(trans_num(type_trans), :, :) ...
                            = sm_env{nf}(:, mnind:mxind);
                    end
                    if nargout>1
                        trans_durations{type_trans}(trans_num(type_trans), :, :) ...
                            = trans_durs2(:, ttn);
                    end
                end
            end
            
            
        end
    end
end

varargout{1} = btn_trans;
if nargout>1
    varargout{2} = trans_durations;
end

end