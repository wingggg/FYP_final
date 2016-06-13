load('runtime_results');

%% single variable barchart
% figure
% bar(num_vector, runtime_all_pc);
% title('Execution Time on PC and Raspberry Pi')
% xlabel('Data Subsets');
% ylabel('Execution Time (seconds)');

%% two variable barchart
% vector = zeros(33,2);
% vector(:,1) = runtime_all_before_opt(:);
% vector(:,2) = runtime_all_threading(:);
% 
% figure
% b = bar(vector);
% title('Execution Time on Raspberry Pi Before and After Optimisation')
% xlabel('Data Subsets');
% ylabel('Execution Time (seconds)');
% b(1).FaceColor = 'c';
% b(2).FaceColor = 'g';
% b(1).EdgeColor = 'w';
% b(2).EdgeColor = 'w';
% l = legend('Before optimisation','Multi-threading');

%% multi variable stacked
% vector = zeros(4,33);
% vector(1,:) = runtime_all_pc(:)';
% vector(2,:) = runtime_all_before_opt(:)';
% vector(3,:) = runtime_all_threading(:)';
% vector(4,:) = runtime_all_parallelising(:,2)';

figure
b = bar(vector, 'stacked');
title('Total Execution Time For Complete Dataset')
xlabel('Data Subsets');
ylabel('Execution Time (seconds)');
for i=1:33
    b(i).EdgeColor = 'w';
end
legend('1: PC','2: Before optimisation','3: Multi-threading','4: Parallelisation');
