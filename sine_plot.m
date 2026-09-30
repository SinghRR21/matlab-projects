% sine_plot.m
% Plots sin(x) and cos(x) from 0 to 2*pi on the same graph.

clear; clc; close all;                  % Clear the workspace, clear the command window, and close open figures

x = linspace(0, 2*pi, 200);             % Create 200 evenly spaced points from 0 to 2*pi
y_sin = sin(x);                         % Compute the sine of each x value
y_cos = cos(x);                         % Compute the cosine of each x value

figure;                                 % Open a new figure window
plot(x, y_sin, 'b-', 'LineWidth', 1.5); % Plot sin(x) as a solid blue line
hold on;                                % Keep the current plot so the next plot is added to it
plot(x, y_cos, 'r--', 'LineWidth', 1.5);% Plot cos(x) as a dashed red line
hold off;                               % Release the plot so future plots replace it

xlim([0 2*pi]);                         % Limit the x-axis to the range 0 to 2*pi
xticks(0:pi/2:2*pi);                    % Place x-axis ticks at multiples of pi/2
xticklabels({'0', '\pi/2', '\pi', '3\pi/2', '2\pi'}); % Label the ticks in terms of pi
grid on;                                % Show grid lines to make values easier to read

xlabel('x (radians)');                  % Label the x-axis
ylabel('y');                            % Label the y-axis
title('sin(x) and cos(x) from 0 to 2\pi'); % Add a title to the graph
legend('sin(x)', 'cos(x)', 'Location', 'best'); % Add a legend identifying each curve
