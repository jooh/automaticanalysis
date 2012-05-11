% Create a movie

function aas_image_avi(imageFN, axisDim, frameSize, rotations, outlineFN)
if ischar(imageFN)
    imageFN = {imageFN};
end
if nargin < 2
    axisDim = 1;
end
if nargin < 3
    frameSize = [400 500];
end
if nargin < 4
    rotations = 0;
end
if nargin < 5
    outlineFN = [];
end

% Save graphical output to common diagnostics directory
if ~exist(fullfile(getHome, 'diagnostics'), 'dir')
    mkdir(fullfile(getHome, 'diagnostics'))
end

[~, movieFN] = fileparts(imageFN{1});

% Make a movie file
movieFilename = fullfile(getHome, 'diagnostics', ...
    [movieFN '.avi']);

% Create movie file by defining aviObject
if exist(movieFilename,'file')
    delete(movieFilename);
end

aviObject = avifile(movieFilename,'compression','none');

% Get the figure!
figure(2)
try
    % Try to set figure 1 to be on top!
    fJframe = get(2, 'JavaFrame');
    fJframe.fFigureClient.getWindow.setAlwaysOnTop(true)
catch
end

% This does not work if it's larger than the window, conservative...
windowSize = [1 1 frameSize(1) frameSize(2)];
set(2,'Position', windowSize)

Y = cell(size(imageFN));
% Load the image
for f = 1:length(imageFN)
    Y{f} = spm_read_vols(spm_vol(imageFN{f}));
    limsY{f} = [min(Y{f}(:)) max(Y{f}(:))];
end
% Load the outline
if ~isempty(outlineFN)
    oY = spm_read_vols(spm_vol(outlineFN));
end

colormap gray

for d = 1:size(Y{1},axisDim)
    for f = 1:length(imageFN)
        h = subplot(1, length(imageFN), f);
        if axisDim == 1
            imageSlice = squeeze(Y{f}(d,:,:));
        elseif axisDim == 2
            imageSlice = squeeze(Y{f}(:,d,:));
        elseif axisDim == 3
            imageSlice = squeeze(Y{f}(:,:,d));
        end
        
        if ~isempty(outlineFN)
            if axisDim == 1
                outlineSlice = edge(squeeze(oY(d,:,:)),'canny');
            elseif axisDim == 2
                outlineSlice = edge(squeeze(oY(:,d,:)),'canny');
            elseif axisDim == 3
                outlineSlice = edge(squeeze(oY(:,:,d)),'canny');
            end
            
            %imageSlice(logical(outlineSlice)) = NaN;
        end
        
        for r = 1:rotations
            imageSlice = rot90(imageSlice);
        end
        
        imagescnan(imageSlice, 'NanColor', [1 0 0])
        if ~isempty(outlineFN)
           hold on
           [x y] = find(flipdim(outlineSlice,2));
           scatter(x,y,3,'r', 'd') 
           hold off
        end
        
        caxis([limsY{f}(1), limsY{f}(2)])
        axis equal off
        zoomSubplot(h, 1.2)
    end
    pause(0.01)
    drawnow
    
    % Capture frame and store in aviObject
    F = getframe(2,windowSize);
        
    aviObject = addframe(aviObject,F);
end

aviObject = close(aviObject);

try
    % Return figure 1 to not be on top!
    fJframe.fFigureClient.getWindow.setAlwaysOnTop(false)
catch
end