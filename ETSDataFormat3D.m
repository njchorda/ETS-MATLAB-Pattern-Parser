function varargout = ETSDataFormat3D(filename)
% If you only care about total power or only one polarization is measured:
%   [D, phi, theta, freqs] = ETSDataFormat3D(filename)
% If you want all phi, theta, and total polarizations:
% - [D_total, D_phi, D_theta, phi, theta, freqs] = ETSDataFormat3D(filename)
% - If only one vertical/horizontal/total is present, that is the only data
% array that will be output in the three-argout case
% - D is a numPhi x numTheta x numFreqs array. The function will extract all of the data as an array, as well as
% the angular step in phi, angular step in theta, and measured frequencies as separate vectors.
% The function will save a .mat file with the same name as the input file
% and process that if it exists.
% clear;clc;close all
% filename = '251211_3DPrintedParasiticPatch_Hcut.csv';

[filePath, fname_noExt, ~] = fileparts(filename);
matFilename = [filePath fname_noExt '.mat'];
if exist(matFilename, 'file')
    matFileExists = true;
    disp('mat file already exists. Loading...')
    varsInFile = who('-file', matFilename);
    load(matFilename);
    switch numel(varsInFile)
        case 4
            if exist('D_total', 'var')
                varargout{1} = D_total;
            elseif exist('D_vert', 'var')
                varargout{1} = D_vert;
            elseif exist('D_horiz', 'var')
                varargout{1} = D_horiz;
            end
            varargout{2} = phi;
            varargout{3} = theta;
            varargout{4} = freqs;
        case 6
        if ~(exist('D_total', 'var') && exist('D_phi', 'var') && exist('D_theta', 'var'))
            error('Too many output arguments. One of the vertical/horizontal/total polarizations does not exist');
        end
        varargout{1} = D_total;
        varargout{2} = D_phi;
        varargout{3} = D_theta;
        varargout{4} = phi;
        varargout{5} = theta;
        varargout{6} = freqs;
        otherwise
            error('Number of variables in saved mat file is incorrect.')
    end
    return % Exit the function if the mat file already exists
else
    matFileExists = false;
    disp('mat file does not exist. Loading data from csv...')
end

D = readcell(filename);
D = D(3:end, :);
disp('Processed input file')


%Find indices of polarization "type" (vertical/horizontal/total)
idxTypeCell = cellfun(@ismissing, D(:, 1), 'UniformOutput', false);
% idxTypeCell = ismissing(D(:,1));
% idxsType = find(~idxTypeCell);
idxsType = [];
for i = 1:numel(idxTypeCell)
    if ~idxTypeCell{i}
        idxsType(end + 1) = i;
    end
end


% disp('Num tables')
numTables = numel(idxsType); %Determine the number of tables (phi/theta/total)
%Potentially an issue if the table includes the antenna parameters or not
% disp('Num rows per table')
numRowsPerTable = idxsType(2) - idxsType(1) - 1; % determine the length of each table (num theta plus header)

%For each table, determine if it is for the horizontal/vertical/total power
for i = 1:numTables
    if contains(lower(D{idxsType(i), 1}), 'theta')
        D_thetaTemp = D(idxsType(i):idxsType(i)+numRowsPerTable, 2:end);
        if ~(exist('theta', 'var') && exist('phi', 'var') && exist('freqs', 'var')) %Only do this processing if hasn't been done yet
            theta = D_thetaTemp(2, 3:end);
            thetaMask = ~cell2mat(cellfun(@ismissing, D_thetaTemp(2, 3:end), 'UniformOutput', false));
            % disp('theta--------------------------------------------------')
            theta = cell2mat(theta(thetaMask)); %Extract theta values and put in vector
            phi = D_thetaTemp(:, 2);
            phimask = cellfun('isclass', phi, 'double');
            % disp('phi--------------------------------------------------')
            phi = sort(unique(cell2mat(phi(phimask)))); %Extract phi values and put in vector
            freqs = D_thetaTemp(2:end, 1)';
            freqMask = ~cell2mat(cellfun(@ismissing, freqs, 'UniformOutput', false));
            freqs = cell2mat(freqs(freqMask));
            freqs = freqs(freqs>0); %Extract all frequencies, determine which are above zero, and put in vector
        end
        D_thetaTemp = D_thetaTemp(2:end, 3:3+numel(theta)-1); %trim data
        D_theta = zeros(numel(phi), numel(theta), numel(freqs));
        for j = 0:numel(freqs)-1
            % D_thetaTemp((3+j*(3+numel(phi))):(3+j*(3+numel(phi)) + numel(phi)-1),:)
            % D_thetaTemp((3+j*(3+numel(phi))):(3+j*(3+numel(phi)) + stopCond),:)
            % cell2mat(D_thetaTemp((3+j*(3+numel(phi))):(3+j*(3+numel(phi)) + numel(phi)-1),:))
            if numel(phi) == 1
                D_theta(:, :, j + 1) = cell2mat(D_thetaTemp(3*(j + 1), :));
            else
                blockStart = 3 + j * (2 + numel(phi));
                blockStop = blockStart + numel(phi) - 1;
                D_theta(:, :, j+1) = cell2mat(D_thetaTemp(blockStart:blockStop,:));
            end
        end
        disp('Processed theta polarization')
    elseif contains(lower(D{idxsType(i), 1}), 'phi')
        D_phiTemp = D(idxsType(i):idxsType(i)+numRowsPerTable, 2:end);
        if ~(exist('theta', 'var') && exist('phi', 'var') && exist('freqs', 'var'))
            theta = D_phiTemp(2, 3:end);
            thetaMask = ~cell2mat(cellfun(@ismissing, D_phiTemp(2, 3:end), 'UniformOutput', false));
            theta = cell2mat(theta(thetaMask));
            phi = D_phiTemp(:, 2);
            phimask = cellfun('isclass', phi, 'double');
            phi = sort(unique(cell2mat(phi(phimask))));
            freqs = D_phiTemp(2:end, 1)';
            freqMask = ~cell2mat(cellfun(@ismissing, freqs, 'UniformOutput', false));
            freqs = cell2mat(freqs(freqMask));
            freqs = freqs(freqs>0);
        end
        D_phiTemp = D_phiTemp(2:end, 3:3+numel(theta)-1); %trim data
        D_phi = zeros(numel(phi), numel(theta), numel(freqs));
        for j = 0:numel(freqs)-1
            % D_phi(:, :, j+1) = cell2mat(D_phiTemp((3+j*(3+numel(phi))):(3+j*(3+numel(phi)) + numel(phi)-1),:));
            if numel(phi) == 1
                D_phi(:, :, j + 1) = cell2mat(D_phiTemp(3*(j + 1), :));
            else
                blockStart = 3 + j * (2 + numel(phi));
                blockStop = blockStart + numel(phi) - 1;
                D_phi(:, :, j+1) = cell2mat(D_phiTemp(blockStart:blockStop,:));
            end
        end
        disp('Processed phi polarization')
    elseif contains(lower(D{idxsType(i), 1}), 'total power')
        D_totalTemp = D(idxsType(i):idxsType(i)+numRowsPerTable, 2:end);
        if ~(exist('theta', 'var') && exist('phi', 'var') && exist('freqs', 'var'))
            theta = D_totalTemp(2, 3:end);
            thetaMask = ~cell2mat(cellfun(@ismissing, D_totalTemp(2, 3:end), 'UniformOutput', false));
            theta = cell2mat(theta(thetaMask));
            phi = D_totalTemp(:, 2);
            phimask = cellfun('isclass', phi, 'double');
            phi = sort(unique(cell2mat(phi(phimask))));
            freqs = D_totalTemp(2:end, 1)';
            freqMask = ~cell2mat(cellfun(@ismissing, freqs, 'UniformOutput', false));
            freqs = cell2mat(freqs(freqMask));
            freqs = freqs(freqs>0);
        end
        D_totalTemp = D_totalTemp(2:end, 3:3+numel(theta)-1); %trim data
        D_total = zeros(numel(phi), numel(theta), numel(freqs));
        for j = 0:numel(freqs)-1
            % D_total(:, :, j+1) = cell2mat(D_totalTemp((3+j*(3+numel(phi))):(3+j*(3+numel(phi)) + numel(phi)-1),:));
            if numel(phi) == 1
                D_total(:, :, j + 1) = cell2mat(D_totalTemp(3*(j + 1), :));
            else
                blockStart = 3 + j * (2 + numel(phi));
                blockStop = blockStart + numel(phi) - 1;
                D_total(:, :, j+1) = cell2mat(D_totalTemp(blockStart:blockStop,:));
            end
        end
        disp('Processed total power')
    end
end

freqs = freqs.*1e6; % Format frequencies in Hz

switch nargout
    case 4
        if exist('D_total', 'var')
            varargout{1} = D_total;
            D = D_total;
        elseif exist('D_vert', 'var')
            varargout{1} = D_vert;
            D = D_vert;
        elseif exist('D_horiz', 'var')
            varargout{1} = D_horiz;
            D = D_horiz;
        end
        varargout{2} = phi;
        varargout{3} = theta;
        varargout{4} = freqs;
        if ~matFileExists
            save(matFilename, 'D', 'phi', 'theta', 'freqs');
        end
    case 6
        if ~(exist('D_total', 'var') && exist('D_phi', 'var') && exist('D_theta', 'var'))
            error('Too many output arguments. One of the vertical/horizontal/total polarizations does not exist');
        end
        varargout{1} = D_total;
        varargout{2} = D_phi;
        varargout{3} = D_theta;
        varargout{4} = phi;
        varargout{5} = theta;
        varargout{6} = freqs;
        if ~matFileExists
            save(matFilename, 'D_total', 'D_phi', 'D_theta', 'phi', 'theta', 'freqs');
        end
    otherwise
        error('Output should be in the format [Data, phi, theta, freqs] or [Data_total, Data_phi, Data_theta, phi, theta, freqs]');
end

end