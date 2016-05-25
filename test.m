%% Read files
file = fopen('./data/bcode/bin/output01.txt');
outString = textscan(file, '%s');
outString = outString{1,1};
fclose(file);

S = strsplit(outString{1,1}, '/');
datasetName = S(length(S)-1);
file = fopen(strjoin(['benchmark_' datasetName '.txt'], ''));
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

%% Generate look-up table

% Image dimensions
x_size = 340;
y_size = 180;

% Construct regions
grid = [0 1 2; 3 4 5; 6 7 8];
lut = struct([]);
for j=1:3
    for i=1:3
        lut = [lut struct(  'x_start', (i-1)*floor(x_size/3), ...
                            'y_start', (j-1)*floor(y_size/3), ...
                            'x_index', i, ...
                            'y_index', j ...
                            )];
    end
end

%% Comparing results

cmpResults = struct([]);
for i=1:length(outStruct)
    if strcmp(benchStruct(i).filename, outStruct(i).filename)
        % write code
        
        
        cmpResults = [cmpResults struct('filename', outStruct(i).filename, ...
                                        'binaryMatch', binaryMatch, ...
                                        'regionsMatch', regionsMatch ...        % as a percentage
                                        )];
    end
end


%% Accuracy

