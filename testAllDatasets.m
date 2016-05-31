% Program tests parameters in outputAll.txt against those in
% benchmarkAll.txt and writes results to testResultsAll.txt

%% Image dimensions
    x_size = 340;
    y_size = 180;
    x_blockSize = x_size/3;
    y_blockSize = y_size/3;

%% Read files
file1 = fopen('./data/bcode/bin/outputAll.txt');
outString = textscan(file1, '%s');
outString = outString{1,1};
fclose(file1);

file2 = fopen('./benchmark/benchmarkAll.txt');
benchString = textscan(file2, '%s');
benchString = benchString{1,1};
fclose(file2);

outfile = fopen('./benchmark/testResultsAll.txt', 'wt');

%% Count number of files being processed
numFiles = 0;
for i=1:length(benchString)
    str = strsplit(benchString{i,1}, '.');
    if strcmp(str(length(str)), 'JPG') || strcmp(str(length(str)), 'jpg')
        numFiles = numFiles + 1;
    end
end

%% For each file:
k = 1;
for i=1:numFiles
    % constructing outStruct
    j = (i-1)*8+1;
    outStruct = struct( 'filename', outString{j,1}, ...
                        'valid', str2num(outString{j+1,1}), ...
                        'start_x', str2num(outString{j+2,1}), ...
                        'start_y', str2num(outString{j+3,1}), ...
                        'range_x', str2num(outString{j+4,1}), ...
                        'range_y', str2num(outString{j+5,1}), ...
                        'confidence', str2num(outString{j+6,1}), ...
                        'coverage', str2num(outString{j+7,1}) ...
                        );

    S = strsplit(outStruct.filename, '/');
    datasetName = S(length(S)-1);

    % constructing benchStruct
    str = [0 strsplit(benchString{k,1}, '/')];
    
    if strcmp(str(length(str)-1), datasetName)
        filename = benchString{k,1};
        state = str2num(benchString{k+1,1});
        
        regions = [];
        count = 2;
        if state == 1
            while length(benchString{k+count,1}) == 1
                regions = [regions str2num(benchString{k+count,1})];
                if (k+count) < length(benchString)      % so that index does not exceed matrix dimension
                    count = count + 1;
                else
                    break;
                end
            end
        end
        k = k + count;
        benchStruct = struct('filename', filename, 'state', state, 'regions', regions);
    end

    %% Comparing results
    % make sure the correct files are being compared
    if strcmp(benchStruct.filename, outStruct.filename)
        output = outStruct;
        validMatch = 0;
        regionsMatch = 0;
        regions = [];

        % validity match
        if benchStruct.state == outStruct.valid
            validMatch = 1;
        end

        % regions match
        if output.valid == 1
            % accumulate regions
            range_y = output.range_y;
            start_y = output.start_y;
            for m=1:3
                range_x = output.range_x;
                start_x = output.start_x;
                for n=1:3
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
            if benchStruct.state == 1
                inter = intersect(regions, benchStruct.regions);
                union = unique(cat(2, regions, benchStruct.regions));
                % Accuracy test comparing to the 'union' of the sets instead 
                % of to just one of the sets because error consists of both 
                % cases: +- a region
                regionsMatch = length(inter)/length(union);       % matching result out of 1
            end
        end

        % output results
        testResults = struct(   'outputIndex', i, ...
                                'filename', outStruct.filename, ...
                                'validDetected', outStruct.valid, ...
                                'validBenchmark', benchStruct.state, ...
                                'regionsDetected', regions, ...
                                'regionsBenchmark', benchStruct.regions, ...
                                'validMatch', validMatch, ...
                                'regionsMatch', regionsMatch ...        % as a percentage
                                );
                            
        formatSpec = '%d %s  \n';
        fprintf(outfile, formatSpec, [i testResults.filename]);
    else
        fprintf(outfile, formatSpec, [i 'ERROR']);
    end

end


% %% Accuracy
% % Accuracy test for animal detected or not
% validAccuracy = sum([testResults(:).validMatch])/length([testResults(:).validMatch]);
% 
% % Accuracy test for detection bounding box for those images accurately
% % filtered only
% regionsAccuracy = 0;
% index = intersect(find([testResults(:).validDetected]), find([testResults(:).validBenchmark]));
% for i=1:length(index)
%     regionsAccuracy = regionsAccuracy + testResults(index(i)).regionsMatch;
% end
% regionsAccuracy = regionsAccuracy/length(index);
% 


fclose(outfile);
