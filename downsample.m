%מוריד 6

% Load the original WAV file
[inputAudio, originalFs] = audioread('Recording (13).wav'); % Ensure valid input

% Calculate the new sample rate after downsampling
downsamplingFactor = 6;
newFs = originalFs / downsamplingFactor; % New sample rate

% Original and downsampled time vectors
tOriginal = (0:length(inputAudio) - 1) / originalFs;
tDownsampled = linspace(0, tOriginal(end), length(inputAudio) / downsamplingFactor);

% Use interpolation to downsample
outputAudio = interp1(tOriginal, inputAudio, tDownsampled, 'linear'); % Linear interpolation

% Save the downsampled audio to a new WAV file
audiowrite('Recording (13)_downsampled.wav', outputAudio, newFs);

disp('Audio file has been downsampled using interpolation and saved as "output_downsampled_interp.wav".');


% --------------------------------------------------------------
%מפה הקוד

% Set the folder where the recordings are located
folderPath = 'C:\Users\oror8\Downloads\ForSync\ForSync'; % Change this to your actual folder path

% Downsampling factor
downsamplingFactor = 6;

% Loop through the specified range of recordings
for i = 558:567
    % Construct the filename
    fileName = sprintf('Recording (%d).wav', i);
    fullFilePath = fullfile(folderPath, fileName);
    
    % Check if the file exists
    if exist(fullFilePath, 'file')
        % Load the original WAV file
        [inputAudio, originalFs] = audioread(fullFilePath);

        % Calculate the new sample rate after downsampling
        newFs = originalFs / downsamplingFactor; % New sample rate

        % Create the original and downsampled time vectors
        tOriginal = (0:length(inputAudio) - 1) / originalFs;
        tDownsampled = linspace(0, tOriginal(end), ceil(length(inputAudio) / downsamplingFactor)); % Ensure integer length

        % Use linear interpolation to downsample
        outputAudio = interp1(tOriginal, inputAudio, tDownsampled, 'linear');

        % Construct the output filename
        outputFileName = sprintf('Recording (%d)_downsampled.wav', i);
        outputFullFilePath = fullfile(folderPath, outputFileName);

        % Save the downsampled audio to a new WAV file
        audiowrite(outputFullFilePath, outputAudio, newFs);

        % Inform that the file has been processed
        disp(['File "', fileName, '" has been downsampled and saved as "', outputFileName, '".']);
    else
        % Inform if the file does not exist
        disp(['File "', fileName, '" does not exist in the specified folder.']);
    end
end
