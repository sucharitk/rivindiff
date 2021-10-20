function sm_env = get_envelope(epochs, fs, selectrodes, freqs_to_plot, pass_filt, raw_eeg_trans, ...
    avg_all_chans, smoothdur, artif_mask, smooth_env)

if ~exist('pass_filt', 'var') || isempty(pass_filt), pass_filt = 0.1; end
if ~exist('avg_all_chans', 'var')||isempty(avg_all_chans), avg_all_chans = false; end
if ~exist('smooth_env', 'var')||isempty(smooth_env), smooth_env = false; end

if pass_filt
    if numel(pass_filt)==1
        % High pass filter
        if fs>250
            filter_order = 4;
        else
            filter_order = 5;
        end
        hpfilt_cutoff = pass_filt/(fs/2);
        [hpfilt_b, hpfilt_a] = butter(filter_order, hpfilt_cutoff, 'high');
        
    elseif numel(pass_filt)==2
        % Bands pass filter
        
    end
end

if smooth_env
    gsize = fs/4;
    gauss_filt = normpdf(0:gsize-1,gsize/2,gsize/2);
    gauss_filt = gauss_filt/sum(gauss_filt); % In terms of the eeg data sampling rate
end

numf = numel(freqs_to_plot);
sm_env = cell(1, numf);


for nf = 1:numf
    numelect = 1;
    
    switch raw_eeg_trans
        % 1- for aligning trials according to flicker frequency, 2-for
        % taking the hilber transform
        case 1
            if avg_all_chans
                %%% average channels after computing envelope
                h1sm = double(epochs(selectrodes{nf}, :));
                if smooth_env
                    h1sm = filtfilt(gauss_filt,1, double(h1sm));
                end
                if pass_filt
                    h1sm = filtfilt(hpfilt_b, hpfilt_a, h1sm);
                end
                % for erps with frequency sliding, take the raw signal
                sm_env{nf} = h1sm;
            else
                %%% first average channels then compute envelop
                h1sm = double(mean(epochs(selectrodes{nf}, :), 1));
                if pass_filt
                    h1sm = filtfilt(hpfilt_b, hpfilt_a, h1sm);
                end
                % for erps with frequency sliding, take the raw signal
                sm_env{nf} = h1sm;
            end
            
        case 2
            %%% take the hilber transform
            if avg_all_chans
                
                h1sm = double(mean(epochs(selectrodes{nf}, :), 1));
                if pass_filt
                    h1sm = filtfilt(hpfilt_b, hpfilt_a, h1sm);
                end
                senv = abs(hilbert(h1sm));
                if exist('smoothdur', 'var') && ~isempty(smoothdur) && smoothdur
                    senv = smooth(senv, smoothdur*fs)';
                end
                sm_env{nf} = senv;
            else
                numelect = numel(selectrodes{nf});                 
                %%% average channels after computing envelope
                h1sm = double(epochs(selectrodes{nf}, :));
                % for lookinat erps around different frequency bands
                sm_env{nf} = abs(hilbert(h1sm'))';
            end
            
        case 0
            %%% compute the RLS envelop
            epochsize = smoothdur*fs;
            rls_delay = round(.5*smoothdur*fs); %see dumbtest_rls_simulation
            
            % 0- for computing RLS envelop
            if avg_all_chans
                % average all the possible channels before computing envelop
                %             numelect = 1;
                [hcn1, hsn1] = AdaptiveRLSfast(round(epochsize), ...
                    mean(epochs(selectrodes{nf}, :), 1), ...
                    freqs_to_plot(nf), [], 1, fs);
                
                h1cplx = complex(hcn1,hsn1);
                
                %%%%% shift the envelopes BACK IN TIME because RLS filter induces a
                %%%%% delay
                h1cplx = [h1cplx(rls_delay:end); zeros(rls_delay-1,1)];
                % take the magnitude
                h1env = double(abs(h1cplx));
                
                h1sm = h1env;
                if smooth_env
                    h1sm = filtfilt(gauss_filt,1, double(h1sm));
                end
                if pass_filt
                    h1sm = filtfilt(hpfilt_b, hpfilt_a, h1sm);
                end
                h1sm(end-rls_delay+2:end) = NaN;
                sm_env{nf}(1, :) = h1sm;
                %             sm_env{nf}(numelect, :) = h1sm;
            else
                % return the channels separately
                numelect = numel(selectrodes{nf});
                for nel = 1:numelect
                    [hcn1, hsn1] = AdaptiveRLSfast(round(epochsize), ...
                        epochs(selectrodes{nf}(nel), :), ...
                        freqs_to_plot(nf), [], 1, fs);
                    
                    h1cplx = complex(hcn1,hsn1);
                    
                    % take the magnitude
                    h1env = double(abs(h1cplx));
                    
                    h1sm = h1env;
                    if smooth_env
                        h1sm = filtfilt(gauss_filt,1, double(h1sm));
                    end
                    if pass_filt
                        h1sm = filtfilt(hpfilt_b, hpfilt_a, h1sm);
                    end
                    
                    %%%%% shift the envelopes BACK IN TIME because RLS filter induces a
                    %%%%% delay
                    h1sm = [h1sm(rls_delay:end); zeros(rls_delay-1,1)];

                    sm_env{nf}(nel, :) = h1sm;
                end
            end
    end
    if exist('artif_mask', 'var') && ~isempty(artif_mask)
        am2 = repmat(artif_mask, [numelect, 1]);
        sm_env{nf} = am2.*sm_env{nf};
    end
end

end