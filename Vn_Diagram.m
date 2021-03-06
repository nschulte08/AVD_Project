function [ output_args ] = Vn_Diagram( MTOW, S, altitude, M_cruise, M_max, CL_max )
% Creates the V-n diagram of the vehicle
% ALL INPUT VALUES MUST BE IN ENGLISH UNITS (lbf, ft^2, ft)
% First altitude input should be 0 ft (sea level)
figure;
hold on;
line_styles = {'-','--',':','-.','--',':','-.'};

% Get max positive and negative load factors
n_max = 2.1 + (24000 / (MTOW + 10000)); % FAR 25.337b
if n_max < 2.5
    n_max = 2.5; % must be between 2.5 and 3.8
end
n_min = -1; % FAR 25.337c
curve = [0 : .01 : n_max]';
max_negative_n_index = find(curve(:,1) <= abs(n_min), 1, 'last');
iterations = numel(curve);

for i = 1 : numel(altitude) % create a V-n diagram for each altitude
    [~, ~, ~, rho, a] = ATMO(altitude(i), 'E');
    for j = 1 : iterations
        V_stall = sqrt((2 * curve(j,1) * MTOW) / (rho * S * CL_max));
        curve(j, i+1) = convvel(V_stall, 'ft/s', 'm/s');
    end
    plot(curve(:,i+1), curve(:,1), 'LineStyle', line_styles{i}, 'LineWidth', 2, 'Color', 'k'); % positive curve
    plot(curve(1:max_negative_n_index,i+1), -curve(1:max_negative_n_index,1), 'LineStyle', line_styles{i}, 'LineWidth', 2, 'Color', 'k'); % negative curve
    V_max = convvel(M_max * a, 'ft/s', 'm/s');
    V_cruise = convvel(M_cruise * a, 'ft/s', 'm/s');
    points = [ V_max, n_max; V_max, 0; V_cruise, n_min ]; 
    plot(points(:,1), points(:,2), 'LineStyle', line_styles{i}, 'LineWidth', 2, 'Color', 'k'); % straight border lines
    if i == 1 % plot the top and bottom boundary lines
        V_corner = curve(end,2);
        V_corner_neg = curve(max_negative_n_index, 2);
        plot([V_corner, V_max], [n_max, n_max] ,'LineStyle', line_styles{i}, 'LineWidth', 2, 'Color', 'k'); % max positive load factor (n) line
        plot([V_corner_neg, V_cruise], [n_min, n_min] ,'LineStyle', line_styles{i}, 'LineWidth', 2, 'Color', 'k'); % max negative load factor (-n) line
    end
    % set up legend for each line
    LH(i) = plot(nan, nan, 'LineStyle', line_styles{i}, 'LineWidth', 2, 'Color', 'k');
    L{i} = sprintf('%0.0f %s', convlength(altitude(i), 'ft', 'm'), 'm');
end

grid on;
set(gca, 'fontsize', 12, 'fontweight', 'bold');
title('V-n Diagram', 'FontSize', 16, 'FontWeight', 'bold');
xlabel('True Airspeed, V [m/s]', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Load Factor, n', 'FontSize', 16, 'FontWeight', 'bold');
legend(LH, L, 'Location', 'NW', 'FontSize', 14, 'FontWeight', 'bold');
ylim([n_min - 0.5, n_max + 0.5]);
hold off;
end
