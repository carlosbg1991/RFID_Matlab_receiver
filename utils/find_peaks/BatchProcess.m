% Script for applying a signal processing operation (such as curve fitting)
% to an entire folder (directory) of data files.
%  "C:\Users\Tom\Downloads\DataFiles" is the directory where the data are
%  stored. Change this for your computer.
% After running this script, you can import the dairy file so produced into
% Excel by opening an Excel worksheet, click on a cell, click Data > From
% Text, select the diary file and click Import.  This will put all the
% collected terminal output into that spreadsheet.
format compact
format short g
DataDirectory='C:\Users\Tom\Downloads\DataFiles'; % DataDirectory is where data are stored
FunctionsDirectory='C:\Users\Tom\Dropbox\MATLAB6p5'; % FunctionsDirectory is functions are located
fn=ls(DataDirectory); % generates a character matrix of all files in data directory
Sizefn=size(fn);
NumFiles=Sizefn(1); % Number of files in directory
NameLength=Sizefn(2); % Maximum length of file names in that directory
disp([num2str(NumFiles) ' files found in this directory.'])
disp(' ')
% Open a "diary" file that collects all the terminial window output in one
% text file. Change the directory and file name as desired.
diary([FunctionsDirectory '\BatchProcess_' date '.txt']) % Creates "BatchProcess.txt" and adds date to file name
% DataDirectory is the directory where the Matlab
% functions are stored.  Change this for your computer.
cd(FunctionsDirectory);
% Now go through all the files in that folder and apply signal processing.
% Skip any files that generate and error and keep going.
for FileNum=1:NumFiles
    try
        datamatrix=load([DataDirectory '\' fn(FileNum,1:NameLength)]);
        DataArray=[datamatrix(:,1) datamatrix(:,2)]; % Change to suit format of data
        disp([num2str(FileNum) ': ' fn(FileNum,1:NameLength)])% Prints out the file number and name
        disp('         Peak#   Position       Height       Width         Area')
        % Curve fitting settings:
        windowcenter=0; % Change this to fit other regions of the spectrum
        windowwidth=0; % Change t50 o zoom in or out to smaller or larger region.
        NumPeaks=2; % Change this to the expected number of peaks.
        Shape=1;  % 1-Gaussian, 2-Lorentzian, etc. Type "help peakfit" for list.
        NumTrials=10; % Usually 1 to 10. Lower numbers faster but less stable.
        startvector=0; % Difficult cases may require start vector, see help file.
        BaselineMode=0; % Can be 0, 1, 2,or 3 (1=linear baseline subtraction)
        [FitResults,GoodnessOfFit]=peakfit(DataArray,windowcenter,windowwidth,NumPeaks,Shape,1,NumTrials,startvector,BaselineMode);
        disp(FitResults)
        disp('   % fitting error    R2')
        disp(GoodnessOfFit)
        drawnow
        disp(' ')
        % Add code here if you want to perform another set of 
        % operations on the same data file with different settings, i.e.,
        % different data region, number of peaks, peak shape, etc.
    catch me
        disp([ num2str(FileNum) ': Error encounted with file: ' me.identifier])
        disp(' ')
    end
end
diary off
disp('--------------------------------------------------------------- ')