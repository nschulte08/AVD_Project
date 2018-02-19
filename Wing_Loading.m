function [ max_load, min_load ] = AVD_Wing_Loading( b, MTOW, lambda )
% Schrenk's method to determine wing loading, p536 Nicolai
% UNITS IN METRIC (m and N)

span = [0 : .25 : b/2]'; % divide wingspan into segments
loads = zeros(1); % 1st column: trapezoidal, 2nd column: elliptical, 3rd column: average
n = 1; % gust load factor, Nicolai p 535
L = MTOW * n; % total lift applied to wing

for i = 1 : numel(span)
    y = span(i, 1); % wing station
    loads(i, 1) = (2 * L) / (b * (1 + lambda)) * (1 - (2 * y / b) * (1 - lambda)); % trapezoidal
    loads(i, 2) = (4 * L) / (pi * b) * sqrt(1 - (2 * y / b)^2); % elliptical
    loads(i, 3) = (loads(i, 1) + loads(i, 2)) / 2; % average
end

hold on;
plot(span(:,1), loads(:,1), 'k :', 'LineWidth', 2.5); % trapezoidal
plot(span(:,1), loads(:,2), 'k --', 'LineWidth', 2.5); % elliptical
plot(span(:,1), loads(:,3), 'k -', 'LineWidth', 2.5); % average
grid on; grid minor;
xlabel('Wing Half-Span [m]', 'FontSize', 16, 'FontWeight', 'bold');
ylabel('Lift/Span [N/m]', 'FontSize', 16, 'FontWeight', 'bold');
title('Wing Loading Distribution', 'FontSize', 16, 'FontWeight', 'bold');
legend({'Trapezoidal','Elliptical','Average'}, 'Location', 'SW', 'FontSize', 20);
hold off;

max_load = loads(1,3);
min_load = loads(end,3);
end
