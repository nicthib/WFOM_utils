function spectrum = importOceanOptics(pattern)
%IMPORTOCEANOPTICS  Import an Ocean Optics spectrum file.
%   S = IMPORTOCEANOPTICS(FILENAME) loads data from FILENAME into S.
%
%   IMPORTOCEANOPTICS can import the spectral data contained in Ocean
%   Optics processed spectrum ('ProcSpec') files in the OOI binary format.
%   These files can be saved from the SpectraSuite spectrometer operating
%   software which can be purchased from Ocean Optics.
%
%   FILENAME is the name of the file to import, e.g. 'Spectrum.ProcSpec'.
%   Multiple files can be imported at once by using wildcards, e.g.
%   '*.ProcSpec'.
%
%   The function returns a structure containing metadata and data. If
%   multiple files are imported at once, a structure array is returned.
%   The fields of the structure are:
%
%      path is the directory the file was imported from.
%      filename is the name of the file that was imported.
%      datetime is the date and time the acquisition was made.
%      serial is the serial number of the spectrometer.
%      wavelength_calibration_coefficients is an array of four values:
%        1. The intercept (also called the zeroth order coefficient)
%        2. The first coefficient
%        3. The second coefficient
%        4. The third coefficient
%      integration_time is the integration time.
%      integration_time_unit is 'ms' (milliseconds).
%      saturated is 'yes' or 'no' depending if the specturm is saturated.
%      wavelength is an array of wavelength values.
%      wavelength_unit is 'nm' (nanometres).
%      spectrum is the unprocessed spectral data (as a column vector).
%      dark is the corresponding dark spectrum (if recorded).
%      referecnce is the corresponding reference spectrum (if recorded).
%
%   Note: this function has been tested with files created by Ocean Optics
%   SpectraSuite version 2.0.154.
%
%   Note: this function only reads in only a selection of the metadata in
%   the file.
%
%   Note: this funciton cannot import files saved in the 'tab delimited' or
%   'tab delimited, no header' format. To import these files use MATLAB's
%   IMPORTDATA function or use the Import Wizard by selecting File->Import
%   Data...
%
%   Author: Iain Robinson
%   Contact: iain@physics.org
%   Requirements: MATLAB R2009b or later
%   Date: 2015-12-01
%   Version: 1.0

% The argument filePattern may be a single file name, or a pattern with
% wildcards (e.g. '*.sig') and may or may not include a path.
[directory, wildcard, extension] = fileparts(pattern); % Separate out the parts.
list = dir(fullfile(directory, [wildcard, extension])); % Get a list of the files and directories that match the specified pattern.
listOfFiles = list([list.isdir]==0); % Exclude directories (such as . and ..) from the list, to give a list of files only.
listOfFileNames = {listOfFiles.name};  % Copy the names of the files into a cell array.

% Check that the list of file names is not empty.
if isempty(listOfFileNames)
    error('No files were found which match the name or pattern:\n\t%s', pattern);
end

% Import each file in the list, one at a time.
for j = 1:length(listOfFileNames)
    fprintf('Importing file %s...', listOfFileNames{j});
    fileName = fullfile(directory, listOfFileNames{j});
    
    % Check that the file exists.
    if ~exist(fileName, 'file')
        error('The file:\n\t%s\nwas not found.', fileName);
    end
    
    % Check that the file is a zip file.
    fid = fopen(fileName);
    if ~( fread(fid, 1, 'char') == 'P' && fread(fid, 1, 'char') == 'K' )
        error('The file:\n\t%s\ndoes not appear to be an Ocean Optics binary processed spectrum file. Note that this function cannot read Ocean Optics text files such as the "tab-delimited" format exported by SpectraSuite. If you have saved files in this format, try using MATLAB''s built-in IMPORTDATA function.', fileName);
    end
    fclose(fid);
	
	% Unzip the ProcSpec file to the temporary directory.
	% This may fail if no temporary directory is available, or if the user
	% does not have write permission for the temporary directory. If that
	% happens, change 'tempdir' to the name of some other directory.
	unzippedFileNames = unzip(fileName, tempdir);
	xmlFileName = '';
	for n = 1:numel(unzippedFileNames)
		if strfind(unzippedFileNames{n}, 'ps_')
			xmlFileName = unzippedFileNames{n};
			break
		end
	end
	if isempty(xmlFileName)
		error('The XML data file could not be found in the zip archive:\n\t%s', fileName);
	end

	% Read the XML file.
	rawXML = fileread(xmlFileName);
	
	% Ocean Optics SpectraSuite produces files which are not valid XML as
	% they contain invalid characters and do not specify a character
	% encoding. Fix these problems before parsing the file.

	% Open a file to save the data. This will overwrite the existing
	% (temporary) file.
	fid = fopen(xmlFileName, 'w');

	% Save an XML declaration to specify an encoding.
	fwrite(fid, '<?xml version="1.0" encoding="ISO-8859-1" ?>', 'uchar');

	% Remove invalid characters from the XML.
	cleanXML = char(regexprep(rawXML, '\o020', '')); % Removes any invalid "data link escape" characters.

	% Write the cleaned-up XML.
	fwrite(fid, cleanXML, 'uchar');

	% Close the file.
	fclose(fid);
	
	% Parse the cleaned-up file.
	document = xmlread(xmlFileName);

	% Delete the temporary files which were unzipped from the archive.
	for fileToDelete = unzippedFileNames
		delete(char(fileToDelete))
	end
	
    % Process the parsed document.
    rootNode = document.getDocumentElement;
    
    % Read the "source spectrum".
    sourceSpectraNode = rootNode.getElementsByTagName('sourceSpectra').item(0);
    
    % Read whether the spectrum is saturated.
    saturatedChar = char(sourceSpectraNode.getElementsByTagName('saturated').item(0).getTextContent);
    if strcmp(saturatedChar, 'false')
        saturated = 'no';
    elseif strcmp(saturatedChar, 'true')
        saturated = 'yes';
    else
        saturated = 'unknown';
    end
    
    % Read the integration time (which seems to be in microseconds) and
    % convert to milliseconds.
    integrationTimeInMicroseconds = str2double(sourceSpectraNode.getElementsByTagName('integrationTime').item(0).getTextContent);
    integrationTimeInMilliseconds = integrationTimeInMicroseconds / 1000;

    % Read the data.
    pixelValuesNode = sourceSpectraNode.getElementsByTagName('pixelValues').item(0);
    doubleNodeList = pixelValuesNode.getElementsByTagName('double');
    spectrumNumberOfPixels = doubleNodeList.getLength;
    metadataNumberOfPixels = str2double(sourceSpectraNode.getElementsByTagName('numberOfPixels').item(0).getTextContent);
    if spectrumNumberOfPixels ~= metadataNumberOfPixels
        error('Wrong number of pixels in file. The metadata states that the spectrum contains %d pixels, but %d pixels were found in the spectrum.', metadataNumberOfPixels, spectrumNumberOfPixels);
    end
    numberOfPixels = spectrumNumberOfPixels;
    pixels = (0:1:numberOfPixels-1)';
    
    % Read the "source" spectrum.
    data = zeros(numberOfPixels, 1); % Preallocate array for speed.
    for i=1:numberOfPixels
        data(i) = str2double(doubleNodeList.item(i-1).getTextContent);
    end
    
    % Read the acquisition date and time. The datetime is in Java
    % milliseconds. Convert it to a date string.
    acquisitionTimeNode = sourceSpectraNode.getElementsByTagName('acquisitionTime').item(0);
    milliTime = str2double(acquisitionTimeNode.getElementsByTagName('milliTime').item(0).getTextContent);
    acquisitionDatetime = datestr(milliTime/86400/1000 + datenum(1970,1,1));
    
    % Read the wavelengths
    channelWavelengthsNode = sourceSpectraNode.getElementsByTagName('channelWavelengths').item(0);
    doubleNodeList = channelWavelengthsNode.getElementsByTagName('double');
    if doubleNodeList.getLength ~= numberOfPixels
        error('Number of wavelength values not equal to number of data values in file:\n\t%s', fileName);
    end
    wavelength = zeros(numberOfPixels,1);
    for i=1:numberOfPixels
        wavelength(i) = str2double(doubleNodeList.item(i-1).getTextContent);
    end
    
    % Read the wavelength calibration coefficients.
    channelCoefficientsNode = sourceSpectraNode.getElementsByTagName('channelCoefficients').item(0);
    intercept = str2double(channelCoefficientsNode.getElementsByTagName('__WlIntercept').item(0).getTextContent);
    firstCoefficient = str2double(channelCoefficientsNode.getElementsByTagName('__WlFirst').item(0).getTextContent);
    secondCoefficient = str2double(channelCoefficientsNode.getElementsByTagName('__WlSecond').item(0).getTextContent);
    thirdCoefficient = str2double(channelCoefficientsNode.getElementsByTagName('__WlThird').item(0).getTextContent);
    % Combine wavelength calibration cofefficients into a single array.
    wavelengthCalibrationCoefficients = [ intercept firstCoefficient secondCoefficient thirdCoefficient ];
    
    % Read the spectrometer's serial number.
    serialNumber = char(sourceSpectraNode.getElementsByTagName('spectrometerSerialNumber').item(0).getTextContent);
    
    % Look for a dark spectrum. If there is a dark spectrum in the file
    % then one of <darkSpectrum> nodes will contain the actual data, the
    % others will just be references to it. If there is no dark spectrum
    % then set the variable dark to be an empty matrix.
    dark = [];
    darkSpectrumNodeList = rootNode.getElementsByTagName('darkSpectrum');
    for d = 1:darkSpectrumNodeList.getLength
        darkSpectrumNode = darkSpectrumNodeList.item(d-1);
        if ~darkSpectrumNode.hasAttribute('reference') % If it's not a reference...
            % Read the spectrum data.
            pixelValuesNode = darkSpectrumNode.getElementsByTagName('pixelValues').item(0);
            doubleNodeList = pixelValuesNode.getElementsByTagName('double');
            if doubleNodeList.getLength ~= numberOfPixels
                error('Dark spectrum has different number of data values to raw spectrum.');
            end
            dark = zeros(numberOfPixels,1);
            for i=1:numberOfPixels
                dark(i) = str2double(doubleNodeList.item(i-1).getTextContent);
            end
        end
    end
    if isempty(dark)
        warning('No dark spectrum in file:\n\t%s', fileName);
    end
    
    % Read the reference spectrum (if there is one).
    reference = [];
    referenceSpectrumNodeList = rootNode.getElementsByTagName('referenceSpectrum');
    for r = 1:referenceSpectrumNodeList.getLength
        referenceSpectrumNode = referenceSpectrumNodeList.item(r-1);
        if ~ referenceSpectrumNode.hasAttribute('reference') % If it's not a reference to a reference spectrum.
            % Read the data.
            pixelValuesNode = referenceSpectrumNode.getElementsByTagName('pixelValues').item(0);
            doubleNodeList = pixelValuesNode.getElementsByTagName('double');
            if doubleNodeList.getLength ~= numberOfPixels
                error('Reference spectrum has different number of data values to raw spectrum.');
            end
            reference = zeros(numberOfPixels,1);
            for i=1:numberOfPixels
                reference(i) = str2double(doubleNodeList.item(i-1).getTextContent);
            end
        end
    end
    if isempty(reference)
        warning('No reference spectrum in file:\n\t%s', fileName);
    end
    
    %
    % COPY DATA INTO A STRUCTURE
    %
    % Give the spectrum a name. The name could probably be taken from
    % inside the file, but here the fileName (without the path and
    % extension) is used.
    %
    % Target spectra have an extra field called 'pair' which matches it
    % up with the corresponding reference spectrum. These values can be
    % edited in the workspace if need be.
    [pathstr, name, ext] = fileparts(fileName);
    
    spectrum(j).path = pathstr;
    spectrum(j).filename = [name, ext];
    spectrum(j).datetime = acquisitionDatetime;
    spectrum(j).serial = serialNumber;
    spectrum(j).wavelength_calibration_coefficients = wavelengthCalibrationCoefficients;
    spectrum(j).integration_time = integrationTimeInMilliseconds;
    spectrum(j).integration_time_unit = 'ms';
    spectrum(j).saturated = saturated;
    spectrum(j).pixels = pixels;
    spectrum(j).wavelength = wavelength;
    spectrum(j).wavelength_unit = 'nm';
    spectrum(j).spectrum = data;
    spectrum(j).dark = dark;
    spectrum(j).reference = reference;
end
end