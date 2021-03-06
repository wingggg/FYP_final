load('runtime_results');
o = load('acc_results_all_original');
p = load('acc_results_pi_para');

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
% vector = zeros(33,2);
% vector(:,1) = mem_all_before_opt(:);
% vector(:,2) = mem_all_parallel(:);
% 
% figure
% b = bar(vector);
% axis([0.5,33.5,0,35]);    % x-axis min, max; y-axis min, max; (no limit:=inf)
% title('Percentage Memory Usage During Runtime On Pi Before and After Optimisation')
% xlabel('Data Subsets');
% ylabel('MEM %');
% b(1).FaceColor = 'y';
% b(2).FaceColor = 'b';
% b(1).EdgeColor = 'w';
% b(2).EdgeColor = 'w';
% l = legend('Before optimisation','Parallelisation');

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

%% subplots accuracy side by side
% figure
% subplot(1,2,1);
% b = bar(num_vector, validAccuracy);
% axis([0.5,33.5,0,1]);    % x-axis min, max; y-axis min, max; (no limit:=inf)
% title('"Valid" accuracy performace for complete dataset')
% xlabel('Data Subsets');
% ylabel('Accuracy');
% b.FaceColor = 'b';
% b.EdgeColor = 'w';
% 
% subplot(1,2,2);
% b = bar(num_vector, regionsAccuracy);
% axis([0.5,33.5,0,1]);    % x-axis min, max; y-axis min, max; (no limit:=inf)
% title({'"Region" accuracy performace for those image sets','that are classified as "valid" correctly'});
% xlabel('Data Subsets');
% ylabel('Accuracy');
% b.FaceColor = 'g';
% b.EdgeColor = 'w';

%% two variable barchart subplots side by side

vector = zeros(33,2);
vector(:,1) = o.validAccuracy(:);
vector(:,2) = p.validAccuracy(:);

figure
subplot(2,1,1);
b = bar(vector);
axis([0.5,33.5,0,1]);    % x-axis min, max; y-axis min, max; (no limit:=inf)
title('"Valid" accuracy performace for complete dataset')
xlabel('Data Subsets');
ylabel('Accuracy');
b(1).FaceColor = 'c';
b(2).FaceColor = 'b';
b(1).EdgeColor = 'w';
b(2).EdgeColor = 'w';
l = legend('Before optimisation','Parallelisation');


vector = zeros(33,2);
vector(:,1) = o.regionsAccuracy(:);
vector(:,2) = p.regionsAccuracy(:);

hold on
subplot(2,1,2);
b = bar(vector);
axis([0.5,33.5,0,1]);    % x-axis min, max; y-axis min, max; (no limit:=inf)
title({'"Region" accuracy performace for those image sets that are classified as "valid" correctly'});
xlabel('Data Subsets');
ylabel('Accuracy');
b(1).FaceColor = 'y';
b(2).FaceColor = 'g';
b(1).EdgeColor = 'w';
b(2).EdgeColor = 'w';
l = legend('Before optimisation','Parallelisation');

hold off