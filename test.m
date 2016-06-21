% data_size = 33;
% 
% testResults = struct([]);
% for i=1:data_size
%     testResults = [testResults; struct([])];
% end
% validAccuracy = zeros(data_size,1);
% regionsAccuracy = zeros(data_size,1);
% 
% directory = './data/bcode/bin/datasets/';
p = load('./plot_matlab/acc_results_all_original_1');

% % change x_size and y_size
% [testResults{1,1}, validAccuracy(1), regionsAccuracy(1)] = testDataset([directory 'output01_pi_imgResized16']);
% [testResults{2,1}, validAccuracy(2), regionsAccuracy(2)] = testDataset([directory 'output01_pi_imgResized4']);
% [testResults{3,1}, validAccuracy(3), regionsAccuracy(3)] = testDataset([directory 'output01_pi_imgResized916']);
% [testResults{4,1}, validAccuracy(4), regionsAccuracy(4)] = testDataset([directory 'output01_pi_orig']);
% [testResults{5,1}, validAccuracy(5), regionsAccuracy(5)] = testDataset([directory 'output01_pi_imgResized2417']);
% [testResults{6,1}, validAccuracy(6), regionsAccuracy(6)] = testDataset([directory 'output01_pi_imgResized4817']);

% [testResults{1,1}, validAccuracy(1), regionsAccuracy(1)] = testDataset([directory 'K01_out']);
% [testResults{2,1}, validAccuracy(2), regionsAccuracy(2)] = testDataset([directory 'K07_out']);
% [testResults{3,1}, validAccuracy(3), regionsAccuracy(3)] = testDataset([directory 'KDC2-KD04_out']);
% [testResults{4,1}, validAccuracy(4), regionsAccuracy(4)] = testDataset([directory 'Mo02_out']);
% [testResults{5,1}, validAccuracy(5), regionsAccuracy(5)] = testDataset([directory 'Mo03_out']);
% [testResults{6,1}, validAccuracy(6), regionsAccuracy(6)] = testDataset([directory 'Mo06_out']);
% [testResults{7,1}, validAccuracy(7), regionsAccuracy(7)] = testDataset([directory 'Mo12_out']);
% [testResults{8,1}, validAccuracy(8), regionsAccuracy(8)] = testDataset([directory 'Mo13_out']);
% [testResults{9,1}, validAccuracy(9), regionsAccuracy(9)] = testDataset([directory 'MT09_out']);
% [testResults{10,1}, validAccuracy(10), regionsAccuracy(10)] = testDataset([directory 'MT10_out']);
% [testResults{11,1}, validAccuracy(11), regionsAccuracy(11)] = testDataset([directory 'MT30_out']);
% [testResults{12,1}, validAccuracy(12), regionsAccuracy(12)] = testDataset([directory 'MTEXT06_out']);
% [testResults{13,1}, validAccuracy(13), regionsAccuracy(13)] = testDataset([directory 'PI05_out']);
% [testResults{14,1}, validAccuracy(14), regionsAccuracy(14)] = testDataset([directory 'PI10_out']);
% [testResults{15,1}, validAccuracy(15), regionsAccuracy(15)] = testDataset([directory 'PI30_out']);
% [testResults{16,1}, validAccuracy(16), regionsAccuracy(16)] = testDataset([directory 'PI35_out']);
% [testResults{17,1}, validAccuracy(17), regionsAccuracy(17)] = testDataset([directory 'SY09_out']);
% [testResults{18,1}, validAccuracy(18), regionsAccuracy(18)] = testDataset([directory 'SY13_out']);
% [testResults{19,1}, validAccuracy(19), regionsAccuracy(19)] = testDataset([directory 'SY16_out']);
% [testResults{20,1}, validAccuracy(20), regionsAccuracy(20)] = testDataset([directory 'SY31_out']);
% [testResults{21,1}, validAccuracy(21), regionsAccuracy(21)] = testDataset([directory 'T02_out']);
% [testResults{22,1}, validAccuracy(22), regionsAccuracy(22)] = testDataset([directory 'T05_out']);
% [testResults{23,1}, validAccuracy(23), regionsAccuracy(23)] = testDataset([directory 'T08_out']);
% [testResults{24,1}, validAccuracy(24), regionsAccuracy(24)] = testDataset([directory 'T11_out']);
% [testResults{25,1}, validAccuracy(25), regionsAccuracy(25)] = testDataset([directory 'T12_out']);
% [testResults{26,1}, validAccuracy(26), regionsAccuracy(26)] = testDataset([directory 'T15_out']);
% [testResults{27,1}, validAccuracy(27), regionsAccuracy(27)] = testDataset([directory 'W09_out']);
% [testResults{28,1}, validAccuracy(28), regionsAccuracy(28)] = testDataset([directory 'W11_out']);
% [testResults{29,1}, validAccuracy(29), regionsAccuracy(29)] = testDataset([directory 'W11b_out']);
% [testResults{30,1}, validAccuracy(30), regionsAccuracy(30)] = testDataset([directory 'W12_out']);
% [testResults{31,1}, validAccuracy(31), regionsAccuracy(31)] = testDataset([directory 'W13_out']);
% [testResults{32,1}, validAccuracy(32), regionsAccuracy(32)] = testDataset([directory 'W13all_out']);
% [testResults{33,1}, validAccuracy(33), regionsAccuracy(33)] = testDataset([directory 'W17_out']);
