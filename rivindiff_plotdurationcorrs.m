function rivindiff_plotdurationcorrs(expt, pre_post_subj)
%
%
% plot correlations between sfm and rivalry as well as between sessions of
% sfm and riv

cd(expt.session_dir)

nsubj = expt.nsubj;

subj_ii = 0;
sess_ii = 0;
subj_ii_s = 0;
subj_num = []; subj_num_s = [];

for ns = 1:nsubj
    % subject loop
    
    subj_data = expt.data(ns);
    
    if subj_data.subj_valid
        if  (pre_post_subj==1 && subj_data.nsess>1) || ...
                (~pre_post_subj && subj_data.nsess==1) || ...
                pre_post_subj>=3
            % for subjects that have two sessions
            subj_ii = subj_ii + 1;
            subj_num = [subj_num ns];
            subj_codes{subj_ii} = [subj_data.subj_code];
            
            notempcell_riv = ~isemptycell(subj_data.invalid_riv);
            notempcell_sfm = ~isemptycell(subj_data.invalid_sfm);
            invalid_riv_runs = []; invalid_sfm_runs = [];
            for ii = 1:2
                if notempcell_riv(ii)
                    invalid_riv_runs = [invalid_riv_runs subj_data.invalid_riv{ii}];
                end
                if notempcell_sfm(ii)
                    invalid_sfm_runs = [invalid_sfm_runs subj_data.invalid_sfm{ii}];
                end
            end
            invalid_riv_runs = unique(invalid_riv_runs);
            invalid_sfm_runs = unique(invalid_sfm_runs);
            
            for snum = 1:subj_data.nsess
                sess_ii = sess_ii + 1;
                subj_valid_riv(subj_ii, snum) = subj_data.subj_valid_riv;
                if subj_data.subj_valid_riv
                    subj_median_riv(subj_ii, snum) = subj_data.riv_durations.mean_all(snum);
                    mx = subj_data.riv_durations.means(snum, 3); if isnan(mx), mx=0; end
                    subj_riv_mixed(subj_ii, snum) = mx;
                end
                
                subj_valid_sfm(subj_ii, snum) = subj_data.subj_valid_sfm;
                if subj_data.subj_valid_sfm
                    subj_median_sfm(subj_ii, snum) = subj_data.sfm_durations.mean_all(snum);
                    mx = subj_data.sfm_durations.means(snum, 3); if isnan(mx), mx=0; end
                    subj_sfm_mixed(subj_ii, snum) = mx;
                end
                
                subj_age(subj_ii) = subj_data.subj_age;
                subj_sex(subj_ii) = subj_data.subj_sex;
            end
        end
    end
end

if pre_post_subj==1
    figure('Name', 'compare session 1 and 2 rivalry and sfm dominance durations')
    subplot(121)
    riv_valid_inds = subj_valid_riv(:, 1)&subj_valid_riv(:, 2);
    scatter(subj_median_riv(riv_valid_inds, 1), subj_median_riv(riv_valid_inds, 2))
    text(subj_median_riv(riv_valid_inds, 1)+.025, subj_median_riv(riv_valid_inds, 2), ...
        subj_codes(riv_valid_inds))
    lc = [min(min(subj_median_riv(riv_valid_inds, :))) ...
        max(max(subj_median_riv(riv_valid_inds, :)))];
    line(lc, lc, 'LineWidth', 3)
    [~, pp] = ttest(diff(subj_median_riv(riv_valid_inds, :)'));
    title(sprintf('rivalry: sess2 > sess1, p-val: %.5g', pp))
    xlabel('session 1')
    ylabel('session 2')
    
    subplot(122)
    sfm_valid_inds = subj_valid_sfm(:, 1)&subj_valid_sfm(:, 2);
    scatter(subj_median_sfm(:, 1), subj_median_sfm(:, 2))
    text(subj_median_sfm(:, 1)+.2, subj_median_sfm(:, 2), subj_codes)
    lc = [min(subj_median_sfm(:)) max(subj_median_sfm(:))];
    [~, pp] = ttest(diff(subj_median_sfm(sfm_valid_inds, :)'));
    title(sprintf('sfm: sess2 > sess1, p-val: %.5g', pp))
    line(lc, lc, 'LineWidth', 3)
    xlabel('session 1')
    ylabel('session 2')
    
    figure('Name', 'compare session 1 and 2 rivalry and sfm mixed durations')
    subplot(121)
    scatter(subj_riv_mixed(riv_valid_inds, 1), subj_riv_mixed(riv_valid_inds, 2))
    text(subj_riv_mixed(riv_valid_inds, 1)+.025, subj_riv_mixed(riv_valid_inds, 2), ...
        subj_codes(riv_valid_inds))
    lc = [min(min(subj_riv_mixed(riv_valid_inds, :))) ...
        max(max(subj_riv_mixed(riv_valid_inds, :)))];
    line(lc, lc, 'LineWidth', 3)
    [~, pp] = ttest(diff(subj_riv_mixed(riv_valid_inds, :)'));
    title(sprintf('rivalry: sess2 > sess1, p-val: %.5g', pp))
    xlabel('session 1')
    ylabel('session 2')
    
    subplot(122)
    scatter(subj_sfm_mixed(:, 1), subj_sfm_mixed(:, 2))
    text(subj_sfm_mixed(:, 1)+.2, subj_sfm_mixed(:, 2), subj_codes)
    lc = [min(subj_sfm_mixed(:)) max(subj_sfm_mixed(:))];
    [~, pp] = ttest(diff(subj_sfm_mixed(sfm_valid_inds, :)'));
    title(sprintf('sfm: sess2 > sess1, p-val: %.5g', pp))
    line(lc, lc, 'LineWidth', 3)
    xlabel('session 1')
    ylabel('session 2')
    
    figure
    subplot(121)
    riv_valid_inds = subj_valid_riv(:, 1)&subj_valid_sfm(:, 1);
    riv1 = subj_median_riv(riv_valid_inds, 1);
    sfm1 = subj_median_sfm(riv_valid_inds, 1);
    valsubj = logical(riv1); valsubj = valsubj&logical(sfm1);
    riv1 = riv1(valsubj); sfm1 = sfm1(valsubj);
    scatter(riv1, sfm1)
    [rr, pp] = corrcoef(riv1, sfm1);
    title(sprintf('correlation coeff: %0.2g, p-val: %0.4g', rr(2), pp(2)))
    p = polyfit(riv1, sfm1, 1);
    sfit = polyval(p, riv1);
    hold on, plot(riv1, sfit, 'LineWidth', 3)
    xlabel('rivalry alternation duration');
    ylabel('sfm alternation duration');
    
    subplot(122)
    riv1 = subj_median_riv(riv_valid_inds, 2);
    sfm1 = subj_median_sfm(riv_valid_inds, 2);
    valsubj = logical(riv1); valsubj = valsubj&logical(sfm1);
    riv1 = riv1(valsubj); sfm1 = sfm1(valsubj);
    scatter(riv1, sfm1)
    [rr, pp] = corrcoef(riv1, sfm1);
    title(sprintf('correlation coeff: %0.2g, p-val: %0.4g', rr(2), pp(2)))
    p = polyfit(riv1, sfm1, 1);
    sfit = polyval(p, riv1);
    hold on, plot(riv1, sfit, 'LineWidth', 3)
    xlabel('rivalry alternation duration');
    ylabel('sfm alternation duration');
    
    %%% 
    subj_valid_riv = logical(subj_valid_riv);
    figure('Name', 'correlate mixed durations with mean dominance')
    subplot(121)
    scatter(subj_riv_mixed(subj_valid_riv(:, 1), 1), subj_median_riv(subj_valid_riv(:, 1), 1))
    text(subj_riv_mixed(subj_valid_riv(:, 1), 1)+.025, subj_median_riv(subj_valid_riv(:, 1), 1), ...
        subj_codes(subj_valid_riv(:, 1)))
    lc = [min([subj_riv_mixed(subj_valid_riv(:, 1), 1); subj_median_riv(subj_valid_riv(:, 1), 1)]) ...
        max([subj_riv_mixed(subj_valid_riv(:, 1), 1); subj_median_riv(subj_valid_riv(:, 1), 1)])];
    line(lc, lc, 'LineWidth', 3)
    [rr, pp] = corrcoef(subj_riv_mixed(subj_valid_riv(:, 1), 1), ...
        subj_median_riv(subj_valid_riv(:, 1), 1));
    title(sprintf('correlation coeff: %0.2g, p-val: %0.4g', rr(2), pp(2)))
    xlabel('mixed')
    ylabel('dominance')
    
    subplot(122)
    scatter(subj_riv_mixed(subj_valid_riv(:, 1), 2), subj_median_riv(subj_valid_riv(:, 1), 2))
    text(subj_riv_mixed(subj_valid_riv(:, 1), 2)+.025, subj_median_riv(subj_valid_riv(:, 1), 2), ...
        subj_codes(subj_valid_riv(:, 1)))
    lc = [min([subj_riv_mixed(subj_valid_riv(:, 1), 2); subj_median_riv(subj_valid_riv(:, 1), 2)]) ...
        max([subj_riv_mixed(subj_valid_riv(:, 1), 2); subj_median_riv(subj_valid_riv(:, 1), 2)])];
    line(lc, lc, 'LineWidth', 3)
    [rr, pp] = corrcoef(subj_riv_mixed(subj_valid_riv(:, 1), 2), ...
        subj_median_riv(subj_valid_riv(:, 1), 2));
    title(sprintf('correlation coeff: %0.2g, p-val: %0.4g', rr(2), pp(2)))
    xlabel('mixed')
    ylabel('dominance')
    
else
    riv_valid = [expt.data.subj_valid_riv];
    sfm_valid = [expt.data.subj_valid_sfm];
    sub_valid = [expt.data.subj_valid];
    v1 = riv_valid(subj_num) & sfm_valid(subj_num) & sub_valid(subj_num);
    millness = [expt.data.subj_millness];
    millness = millness(subj_num);
    
    % v2 = riv_valid(subj_num_s) & sfm_valid(subj_num_s) & sub_valid(subj_num_s);
    % corr between riv and sfm for all session
    figure('Name', 'All first session correlations')
    riv1 = subj_median_riv(v1, 1);
    sfm1 = subj_median_sfm(v1, 1);
    scatter(riv1, sfm1)
    [rr, pp] = corrcoef(riv1, sfm1);
    title(sprintf('correlation coeff: %0.2g, p-val: %0.4g', rr(2), pp(2)))
    p = polyfit(riv1, sfm1, 1);
    sfit = polyval(p, riv1);
    hold on, plot(riv1, sfit, 'LineWidth', 3)
    xlabel('rivalry alternation duration');
    ylabel('sfm alternation duration');
    
    hold on
    scatter(riv1(millness(v1)==3), sfm1(millness(v1)==3), 'r+')
    scatter(riv1(millness(v1)==2), sfm1(millness(v1)==2), 'gs')
    
    subj_valid_riv = logical(subj_valid_riv);
    figure('Name', 'correlate mixed durations with mean dominance')
    subplot(121)
    scatter(subj_riv_mixed(subj_valid_riv(:, 1), 1), subj_median_riv(subj_valid_riv(:, 1), 1))
    text(subj_riv_mixed(subj_valid_riv(:, 1), 1)+.025, subj_median_riv(subj_valid_riv(:, 1), 1), ...
        subj_codes(subj_valid_riv(:, 1)))
    lc = [min([subj_riv_mixed(subj_valid_riv(:, 1), 1); subj_median_riv(subj_valid_riv(:, 1), 1)]) ...
        max([subj_riv_mixed(subj_valid_riv(:, 1), 1); subj_median_riv(subj_valid_riv(:, 1), 1)])];
    line(lc, lc, 'LineWidth', 3)
    [rr, pp] = corrcoef(subj_riv_mixed(subj_valid_riv(:, 1), 1), ...
        subj_median_riv(subj_valid_riv(:, 1), 1));
    title(sprintf('correlation coeff: %0.2g, p-val: %0.4g', rr(2), pp(2)))
    xlabel('mixed')
    ylabel('dominance')
    
   %%% figure 1c for the paper
   %    subj_median_riv;
end
end