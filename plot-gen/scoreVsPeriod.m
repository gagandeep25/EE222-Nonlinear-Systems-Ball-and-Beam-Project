% Your data list
data_pid_sin = [1.04622523485385 1.01172646124932 0.985261787725355 0.964247504458132 0.947669450772139 0.932920923826508 0.921065968286952 0.911071001839102 0.902098014256173];

data_pid_square = [5.13730839233362 4.93337884489977 4.76819171987716 4.57298598761860 4.35285865532941 4.24069616487708 4.09006559041723 3.89773306158793 3.73342406841827];


data_lqr_sin = [0.887800873616232 0.870998868294365 0.857508398114077 0.846235029643763 0.836666052926347 0.828234349653315 0.820718254290085 0.813999820406251 0.808079672624737];

data_lqr_square = [5.61465640268846 5.20848819139829 4.89731832059958 4.58368669488369 4.27438989262143 4.10235271231202 3.92435258277520 3.69738953708598 3.50615211286383];


% X-axis labels
per = 6:0.5:10;
% Create the plot
figure;
plot(per, data_pid_sin, '-o', 'LineWidth', 2);
hold on;
plot(per, data_lqr_sin, '-o', 'LineWidth', 2);
hold off;
%xticks(1:length(labels));
%xticklabels(labels);

legend('PID', 'LQR');
% Add labels and title
xlabel('Period (in seconds)');
ylabel('Score');
title('Score vs Period (Sine Wave)');
grid on;

% Create the plot
figure;
plot(per, data_pid_square, '-o', 'LineWidth', 2);
hold on;
plot(per, data_lqr_square, '-o', 'LineWidth', 2);
hold off;
%xticks(1:length(labels));
%xticklabels(labels);

legend('PID', 'LQR');
% Add labels and title
xlabel('Period (in seconds)');
ylabel('Score');
title('Score vs Period (Square Wave)');
grid on;
