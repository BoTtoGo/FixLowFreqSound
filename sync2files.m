% Define the input folders for the two sets of files
inputFolder1 = 'C:\Users\oror8\Downloads\ForSync\ForSync\הקלטות מצלמה\לאחר הוספת קצב'; % Folder for dpos1_upsampled
inputFolder2 = 'C:\Users\oror8\Downloads\ForSync\ForSync\הקלטות מיקרופון\לאחר הפחתת קצב';       % Folder for Recording files

% Define the output folder for the aligned files
outputFolder = 'C:\Users\oror8\Downloads\ForSync\ForSync\אחרי סנכרון';   % Set the folder for output files

% Ensure the output folder exists
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

% Loop through indices from 1 to 900
startIndex = 558;  % First index
endIndex = 566; % Last index

% Go through all indices and process each pair
for i = startIndex:endIndex
    % Construct the filenames for the input files
    file1Name = sprintf('%d_pos2_upsampled.wav', i); % First WAV file name
    file2Name = sprintf('Recording (%d)_downsampled.wav', i); % Second WAV file name

    % Full paths for the input files
    file1Path = fullfile(inputFolder1, file1Name); % Path for dpos1_upsampled files
    file2Path = fullfile(inputFolder2, file2Name); % Path for Recording files

    % Check if both files exist
    if isfile(file1Path) && isfile(file2Path)
        % Load the first WAV file and ensure it's a single channel (mono)
        [file1, fs1] = audioread(file1Path);
        if size(file1, 2) > 1
            file1 = mean(file1, 2); % Convert to single channel if stereo
        end

        % Load the second WAV file and ensure it's a single channel (mono)
        [file2, fs2] = audioread(file2Path);
        if size(file2, 2) > 1
            file2 = mean(file2, 2); % Convert to single channel if stereo
        end

        % Ensure the sample rates are the same
        if fs1 ~= fs2
            warning('Sample rates of %s and %s do not match. Skipping alignment.', file1Name, file2Name);
            continue; % Move on to the next pair
        end

        % Calculate the cross-correlation and find the maximum lag
        [crossCorr, lags] = xcorr(file1, file2); % Find cross-correlation and lags
        [~, maxIdx] = max(crossCorr); % Index of the maximum correlation
        lag = lags(maxIdx); % Lag with the highest correlation
        timeShift = lag / fs1; % Convert the lag to time (in seconds)

        % Align the signals based on the calculated lag
        if lag > 0
            % Shift the second signal to align with the first
            alignedFile2 = [zeros(lag, 1); file2]; % Add zeros at the beginning
            alignedFile1 = file1; % First signal stays unchanged
        else
            % Shift the first signal to align with the second
            alignedFile1 = [zeros(abs(lag), 1); file1]; % Add zeros at the beginning
            alignedFile2 = file2; % Second signal stays unchanged
        end

        % Ensure both aligned signals have the same length
        minLength = min(length(alignedFile1), length(alignedFile2));
        alignedFile1 = alignedFile1(1:minLength); % Truncate to the same length
        alignedFile2 = alignedFile2(1:minLength); % Truncate to the same length

        % Construct output filenames
        alignedFile1Name = sprintf('%d_pos2_upsampled_aligned.wav', i); % Output name for the first aligned file
        alignedFile2Name = sprintf('Recording_%d_downsampled_aligned.wav', i); % Output name for the second aligned file

        % Full paths for the output files
        alignedFile1Path = fullfile(outputFolder, alignedFile1Name);
        alignedFile2Path = fullfile(outputFolder, alignedFile2Name);

        % Save the aligned files
        audiowrite(alignedFile1Path, alignedFile1, fs1); % Save aligned first WAV file
        audiowrite(alignedFile2Path, alignedFile2, fs2); % Save aligned second WAV file

        % Output a message indicating success
        fprintf('Files %s and %s have been aligned. Time shift: %.4f seconds.\n', ...
                file1Name, file2Name, timeShift);

    else
        % Warn if one or both files are missing and continue
        warning('Files %s and/or %s not found. Skipping alignment.', file1Path, file2Path);
    end
end
