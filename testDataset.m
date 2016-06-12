function [testResults, validAccuracy, regionsAccuracy] = testDataset(outputFileName)

%% Read files
file = fopen(['./data/bcode/bin/' outputFileName '.txt']);
outString = textscan(file, '%s');
outString = outString{1,1};
fclose(file);

S = strsplit(outString{1,1}, '/');
datasetName = S(length(S)-1);
file = fopen(strjoin(['./benchmark/benchmark_' datasetName '.txt'], ''));
benchString = textscan(file, '%s');
benchString = benchString{1,1};
fclose(file);

%% Constructing structures
outStruct = struct([]);
for i=1:8:length(outString)
    outStruct = [outStruct struct(  'filename', outString{i,1}, ...
                                    'valid', str2num(outString{i+1,1}), ...
                                    'start_x', str2num(outString{i+2,1}), ...
                                    'start_y', str2num(outString{i+3,1}), ...
                                    'range_x', str2num(outString{i+4,1}), ...
                                    'range_y', str2num(outString{i+5,1}), ...
                                    'confidence', str2num(outString{i+6,1}), ...
                                    'coverage', str2num(outString{i+7,1}) ...
                                    )];
end

benchStruct = struct([]);
for i=1:length(benchString)
    str = [0 strsplit(benchString{i,1}, '/')];
    if strcmp(str{1, length(str)-1}, datasetName) == 1
        filename = benchString{i,1};
        state = str2num(benchString{i+1,1});
        
        regions = [];
        if state == 1
            count = 2;
            while length(benchString{i+count,1}) == 1
                regions = [regions str2num(benchString{i+count,1})];
                if (i+count) < length(benchString)
                    count = count + 1;
                else
                    break;
                end
            end
        end
        
        benchStruct = [benchStruct struct('filename', filename, 'state', state, 'regions', regions)];
    end
end

%% Image dimensions
    x_size = 340;
    y_size = 180;
    x_blockSize = x_size/3;
    y_blockSize = y_size/3;

%% Comparing results

testResults = struct([]);
for k=1:length(outStruct)
    % make sure the correct files are being compared
    if strcmp(benchStruct(k).filename, outStruct(k).filename)
        output = outStruct(k);
        validMatch = 0;
        regionsMatch = 0;
        regions = [];
        
        % validity match
        if benchStruct(k).state == outStruct(k).valid
            validMatch = 1;
        end
        
        % regions match
        if output.valid == 1
            % accumulate regions
            range_y = output.range_y;
            start_y = output.start_y;
            for i=1:3
                range_x = output.range_x;
                start_x = output.start_x;
                for j=1:3
                    regions = [regions findRegion(start_x, start_y)];
                    if range_x > (x_blockSize - mod(start_x, x_blockSize))
                        range_x = range_x - x_blockSize;
                        start_x = start_x + x_blockSize;
                    end
                end
                if range_y > (y_blockSize - mod(start_y, y_blockSize))
                    range_y = range_y - y_blockSize;
                    start_y = start_y + y_blockSize;
                end
            end
            regions = unique(regions);  % in ascending ordering

            % compare
            if benchStruct(k).state == 1
                inter = intersect(regions, benchStruct(k).regions);
                union = unique(cat(2, regions, benchStruct(k).regions));
                % Accuracy test comparing to the 'union' of the sets instead 
                % of to just one of the sets because error consists of both 
                % cases: +- a region
                regionsMatch = length(inter)/length(union);       % matching result out of 1
            end
        end
        
        % output results
        testResults = [testResults struct('outputIndex', k, ...
                                        'filename', outStruct(k).filename, ...
                                        'validDetected', outStruct(k).valid, ...
                                        'validBenchmark', benchStruct(k).state, ...
                                        'regionsDetected', regions, ...
                                        'regionsBenchmark', benchStruct(k).regions, ...
                                        'validMatch', validMatch, ...
                                        'regionsMatch', regionsMatch ...        % as a percentage
                                        )];
    end
end

%% Accuracy
% Accuracy test for animal detected or not
validAccuracy = sum([testResults(:).validMatch])/length([testResults(:).validMatch]);

% Accuracy test for detection bounding box for those images accurately
% filtered only
regionsAccuracy = 0;
index = intersect(find([testResults(:).validDetected]), find([testResults(:).validBenchmark]));
for i=1:length(index)
    regionsAccuracy = regionsAccuracy + testResults(index(i)).regionsMatch;
end
regionsAccuracy = regionsAccuracy/length(index);


end
