load('runtime_results');

%% single variable barchart
% figure
% bar(num_vector, runtime_all_pc);
% title('Execution Time on PC and Raspberry Pi')
% xlabel('Data Subsets');
% ylabel('Execution Time (seconds)');

%% two variable barchart
vector = zeros(33,2);
vector(:,1) = runtime_all_pc(:);
vector(:,2) = runtime_all_before_opt(:);

figure
b = bar(vector);
title('Execution Time on PC and Raspberry Pi')
xlabel('Data Subsets');
ylabel('Execution Time (seconds)');
b(1).FaceColor = 'r';
b(2).FaceColor = 'c';
b(1).EdgeColor = 'w';
b(2).EdgeColor = 'w';
l = legend('PC','Raspberry Pi');
