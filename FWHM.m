%% FWHM calculation
% improile to get line profile manually if needed 
% 2026 2/9

profile_DAS  = a;          %profile_DAS;
profile_DMAS = b;            %  profile_DMAS;
% profile_DAS  = 20*log10(profile_DAS / max(profile_DAS(:)));
% profile_DMAS = 20*log10(profile_DMAS / max(profile_DMAS(:)));

% fwhm_das  = calculate_fwhm(x_img, profile_DAS);
% fwhm_dmas = calculate_fwhm(x_img, profile_DMAS);
[fwhm_das, x_L_das, x_R_das, y_half_das]   = calculate_fwhm_stats(x_img, profile_DAS);
[fwhm_dmas, x_L_dmas, x_R_dmas, y_half_dmas] = calculate_fwhm_stats(x_img, profile_DMAS);
fprintf('A FWHM: %.3f mm | B FWHM: %.3f mm\n', fwhm_das, fwhm_dmas);

fig2 = figure('Name', 'A vs B Profile with FWHM', 'Color', 'w');

yyaxis left
plot(x_img, profile_DAS, 'b--', 'LineWidth', 1.5);
hold on;
plot([x_L_das, x_R_das], [y_half_das, y_half_das], 'b-', 'LineWidth', 2); 
text(x_R_das + 0.5, y_half_das, sprintf('%.2f mm', fwhm_das), 'Color', 'b', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('A Amplitude (Linear)', 'Color', 'b');
ylim([0, max(profile_DAS)*1.2]); 
ax = gca; ax.YColor = 'b';

yyaxis right
plot(x_img, profile_DMAS, 'r-', 'LineWidth', 2);
hold on;
plot([x_L_dmas, x_R_dmas], [y_half_dmas, y_half_dmas], 'r-', 'LineWidth', 2);
text(x_R_dmas + 0.5, y_half_dmas, sprintf('%.2f mm', fwhm_dmas), 'Color', 'r', 'FontSize', 10, 'FontWeight', 'bold');
ylabel('B Amplitude (Linear)', 'Color', 'r');
ylim([0, max(profile_DMAS)*1.2]); 
ax = gca; ax.YColor = 'r';

grid on;
xlabel('Lateral Position x [mm]');
title({['Lateral Profile at z=40',' mm'], ...
       ['Resolution Improvement: ' sprintf('%.1f%%', (fwhm_das - fwhm_dmas)/fwhm_das * 100)]});

legend({'A Signal', 'A FWHM', 'B Signal', 'B FWHM'}, 'Location', 'best');




function [width, x_left, x_right, half_val] = calculate_fwhm_stats(x_axis, y_signal)

    [max_val, max_idx] = max(y_signal);
    half_val = max_val / 2;
    
    idx_left = find(y_signal(1:max_idx) < half_val, 1, 'last');
    if isempty(idx_left)
        x_left = x_axis(1); 
    else
        % interpolation
        x1 = x_axis(idx_left);      y1 = y_signal(idx_left);
        x2 = x_axis(idx_left+1);    y2 = y_signal(idx_left+1);
        x_left = x1 + (half_val - y1) * (x2 - x1) / (y2 - y1);
    end

    idx_right_temp = find(y_signal(max_idx:end) < half_val, 1, 'first');
    if isempty(idx_right_temp)
        x_right = x_axis(end);
    else
        idx_right = max_idx + idx_right_temp - 1;
        x1 = x_axis(idx_right-1);   y1 = y_signal(idx_right-1);
        x2 = x_axis(idx_right);     y2 = y_signal(idx_right);
        x_right = x1 + (half_val - y1) * (x2 - x1) / (y2 - y1);
    end
    
    width = x_right - x_left;
end