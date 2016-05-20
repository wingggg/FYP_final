% M = dlmread('output01.txt',' ',0,1);
file = fopen('output01.txt');
C = textscan(file, '%s');
fclose(file);

% -------------
params = [];
for i=2:length(C{1,1}) 
    params = [params str2num(C{1,1}{i,1})];
end

