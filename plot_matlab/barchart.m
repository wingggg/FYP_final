load('runtime_results');

%% single variable barchart
% figure
% b = bar(num_vector, mem_all_before_opt);
% axis([0.5,33.5,0,30]);    % x-axis min, max; y-axis min, max; (no limit:=inf)
% title('Percentage Memory Usage During Runtime')
% xlabel('Data Subsets');
% ylabel('MEM %');
% b.FaceColor = 'b';
% b.EdgeColor = 'w';

%% two variable barchart
vector = zeros(33,2);
vector(:,1) = mem_all_before_opt(:);
vector(:,2) = mem_all_parallel(:);

figure
b = bar(vector);
axis([0.5,33.5,0,35]);    % x-axis min, max; y-axis min, max; (no limit:=inf)
title('Percentage Memory Usage During Runtime On Pi Before and After Optimisation')
xlabel('Data Subsets');
ylabel('MEM %');
b(1).FaceColor = 'y';
b(2).FaceColor = 'b';
b(1).EdgeColor = 'w';
b(2).EdgeColor = 'w';
l = legend('Before optimisation','Parallelisation');

%% multi variable stacked
% vector = zeros(4,33);
% vector(1,:) = runtime_all_pc(:)';
% vector(2,:) = runtime_all_before_opt(:)';
% vector(3,:) = runtime_all_threading(:)';
% vector(4,:) = runtime_all_parallelising(:,2)';

% figure
% b = bar(vector, 'stacked');
% title('Total Execution Time For Complete Dataset')
% xlabel('Data Subsets');
% ylabel('Execution Time (seconds)');
% for i=1:33
%     b(i).EdgeColor = 'w';
% end
% legend('1: PC','2: Before optimisation','3: Multi-threading','4: Parallelisation');
