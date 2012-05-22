% Create a movie
% aas_image_avi(imageFN, outlineFN, movieFN, axisDim, frameSize, rotations)
function aas_image_avi(imageFN, outlineFN, movieFN, axisDim, frameSize, rotations)

if ischar(imageFN)
    imageFN = {imageFN};
end
if nargin < 2
    outlineFN = [];
end
if nargin < 3
    [~, movieFN] = fileparts(imageFN{1});
    
    % Make a movie file
    movieFN = fullfile(getHome, ...
        [movieFN '.avi']);
end
if nargin < 4 || isempty(axisDim)
    axisDim = 1;
end
if nargin < 5 || isempty(frameSize)
    frameSize = [400 500];
end
if nargin < 6 || isempty(rotations)
    rotations = 0;
end

% Create movie file by defining aviObject
if exist(movieFN,'file')
    delete(movieFN);
end

aviObject = avifile(movieFN,'compression','none');

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
    % Get a good threshold
    thresh = zeros(size(Y{1},axisDim),2);
    for d = 1:size(Y{1},axisDim)
        if axisDim == 1
                outlineSlice = squeeze(oY(d,:,:));
            elseif axisDim == 2
                outlineSlice = squeeze(oY(:,d,:));
            elseif axisDim == 3
                outlineSlice = squeeze(oY(:,:,d));
            end
            [outlineSlice thresh(d,:)] = edge(outlineSlice, 'canny');
    end
    thresh = mean(thresh);
end

colormap gray

for d = 1:size(Y{1},axisDim)
    for f = 1:length(imageFN)
        h = subplot(1, length(imageFN), f);
        
        % Get image slice to draw
        if axisDim == 1
            imageSlice = squeeze(Y{f}(d,:,:));
        elseif axisDim == 2
            imageSlice = squeeze(Y{f}(:,d,:));
        elseif axisDim == 3
            imageSlice = squeeze(Y{f}(:,:,d));
        end
        
        % If present, get outline slice to draw
        if ~isempty(outlineFN)
            if axisDim == 1
                outlineSlice = squeeze(oY(d,:,:));
            elseif axisDim == 2
                outlineSlice = squeeze(oY(:,d,:));
            elseif axisDim == 3
                outlineSlice = squeeze(oY(:,:,d));
            end
            outlineSlice = edge(outlineSlice, 'canny', thresh);
        end
        
        % Rotate slices
        
        imageSlice = rot90(imageSlice, rotations);
        if ~isempty(outlineFN)
            outlineSlice = rot90(outlineSlice,rotations - 1);
        end
        
        % Draw slices
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
% Save video
aviObject = close(aviObject);

try
    % Return figure 1 to not be on top!
    fJframe.fFigureClient.getWindow.setAlwaysOnTop(false)
catch
end