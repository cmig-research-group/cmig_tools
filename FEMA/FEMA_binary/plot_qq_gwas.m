function plot_qq_gwas(neg_log10_p_values, varargin)


p = inputParser;
addRequired(p, 'neg_log10_p_values', @isnumeric);
addParameter(p, 'title', 'QQ Plot for GWAS', @ischar);
addParameter(p, 'color', [0.2 0.4 0.8], @(x) isnumeric(x) && length(x)==3);
addParameter(p, 'markersize', 20, @isnumeric);
addParameter(p, 'show_lambda', false, @islogical);
parse(p, neg_log10_p_values, varargin{:});

valid_idx = isfinite(neg_log10_p_values);
neg_log10_p_values = neg_log10_p_values(valid_idx);

% sort the observed value
observed_sorted = sort(neg_log10_p_values);
n = length(observed_sorted);

% calculate theoretical quantiles
theoretical_p = ((1:n) - 0.5)' / n;
expected_sorted = -log10(theoretical_p);
expected_sorted = flipud(expected_sorted); 

% plot
figure;

p_lower_ci = betainv(0.025, (1:n)', n - (1:n)' + 1);
p_upper_ci = betainv(0.975, (1:n)', n - (1:n)' + 1);
    
neg_log10_ci_upper = sort(-log10(p_lower_ci));
neg_log10_ci_lower = sort(-log10(p_upper_ci));

plot(expected_sorted, neg_log10_ci_upper, '--', 'Color', [0.5,0.5,0.5], 'LineWidth', 0.8);
hold on;
plot(expected_sorted, neg_log10_ci_lower, '--', 'Color', [0.5,0.5,0.5], 'LineWidth', 0.8);
hold on;

scatter(expected_sorted, observed_sorted, p.Results.markersize, ...
        'filled', 'MarkerFaceColor', p.Results.color, ...
        'MarkerFaceAlpha', 0.6);
hold on;

% add off-diagonal line
max_val = max([expected_sorted; observed_sorted]);
plot([0 max_val], [0 max_val], 'r--', 'LineWidth', 1.5);

% label and title
xlabel('Expected -log_{10}(P)');
ylabel('Observed -log_{10}(P)');
title(p.Results.title);
grid on;

% axis
axis equal;
xlim([0 max_val]);
ylim([0 max_val]);

% GC value - lambda
if p.Results.show_lambda
    lambda = calculate_lambda(neg_log10_p_values);
    text(0.05*max_val, 0.95*max_val, ...
         sprintf('\\lambda = %.3f', lambda), ...
         'FontSize', 12, 'BackgroundColor', 'white');
end

end

function lambda = calculate_lambda(neg_log10_p_values)
p_values = 10.^(-neg_log10_p_values);
chi_squared = chi2inv(1 - p_values, 1);
lambda = median(chi_squared) / 0.4549364;
end
