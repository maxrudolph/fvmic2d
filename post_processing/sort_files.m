function [sortedfiles ltezero] = sort_files(files)

%This is modified specifically to sort files of the form Markers.0.1000 
%http://www.mathworks.com/support/solutions/en/data/1-OGU22/index.html?product=SL&solution=1-OGU22

%second returned object is a list of files whose number is less than zero

numfiles = size(files,1); % Find number of files
numdelim = 3; % Number of delimiters
delims = ['.']; % Delimiters used

% Create a matrix of the file list index and the delimiters
filenums=[ [1:numfiles]' zeros(numfiles,3) ];
for i=1:numfiles % Cycles through list of files
    rem=files(i).name;
    for j=1:numdelim % Cycles through the filename delimiters
        [token,rem] = strtok(rem,delims)
        if(j>1)
        filenums(i,j) = str2num(token); 
        end
    end
end

% Sort the matrix by rows
filenums = sortrows(filenums,[2:numdelim+1]);

% Show the filenames in the new order
files(filenums(:,1)).name

% Create a cell array of the filenames in new order
for i=1:numfiles
    sortedfiles{i,1} = files(filenums(i,1)).name;
end

ltezero = any((filenums<0)');