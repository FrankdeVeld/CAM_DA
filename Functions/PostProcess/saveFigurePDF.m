function saveFigurePDF(fig, fileName, width, height, units)
% saveFigurePDF  Resize a figure window and save it as a PDF
%
%   saveFigurePDF(fig, fileName, width, height)
%   saveFigurePDF(fig, fileName, width, height, units)
%
% Inputs:
%   fig      - figure handle
%   fileName - output file name, with or without ".pdf"
%   width    - figure width
%   height   - figure height
%   units    - units for width/height: 'pixels', 'inches', 'centimeters'
%              default: 'inches'
%
% Example:
%   f = figure;
%   plot(rand(10,1))
%   saveFigurePDF(f, 'myplot', 6, 4, 'inches')

    if nargin < 5 || isempty(units)
        units = 'inches';
    end

    if ~isa(fig, 'matlab.ui.Figure')
        error('First input must be a valid figure handle.');
    end

    if ~ischar(fileName) && ~isstring(fileName)
        error('fileName must be a string or character vector.');
    end

    fileName = char(fileName);
    if ~endsWith(lower(fileName), '.pdf')
        fileName = [fileName '.pdf'];
    end

    oldUnits = fig.Units;
    fig.Units = units;

    pos = fig.Position;
    fig.Position = [pos(1), pos(2), width, height];

    exportgraphics(fig, fileName, 'ContentType', 'vector');

    fig.Units = oldUnits;
end