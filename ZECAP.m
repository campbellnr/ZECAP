% 
% ************************************************************************
% (c) 2017-2020, Nathaniel R. Campbell
% Xavier Lab & White Lab, MSKCC
% nathaniel.r.campbell@gmail.com
% 
% Zebrafish Embryo Cancer Analysis Pipeline (ZECAP) v1.0
% ************************************************************************
%
% Requirements:
% - Bio-Formats for MATLAB
%       https://docs.openmicroscopy.org/bio-formats/5.7.3/users/matlab/index.html
% - user input to find orientation of fish where segmentation of eyes and
%       yolk does not work correctly
%
% Notes:
% - For speed, best to run with showIm and showKeyIm both set to 0
%
%
% MATLAB Version Info:
% ----------------------------------------------------------------------------------------------------
% MATLAB Version: 9.0.0.341360 (R2016a)
% MATLAB License Number: 922760
% Operating System: Mac OS X  Version: 10.13.6 Build: 17G11023 
% Java Version: Java 1.7.0_75-b13 with Oracle Corporation Java HotSpot(TM) 64-Bit Server VM mixed mode
% ----------------------------------------------------------------------------------------------------
% MATLAB                                                Version 9.0         (R2016a)
% Bioinformatics Toolbox                                Version 4.6         (R2016a)
% Curve Fitting Toolbox                                 Version 3.5.3       (R2016a)
% Fuzzy Logic Toolbox                                   Version 2.2.23      (R2016a)
% Image Processing Toolbox                              Version 9.4         (R2016a)
% Neural Network Toolbox                                Version 9.0         (R2016a)
% Optimization Toolbox                                  Version 7.4         (R2016a)
% Parallel Computing Toolbox                            Version 6.8         (R2016a)
% Partial Differential Equation Toolbox                 Version 2.2         (R2016a)
% Signal Processing Toolbox                             Version 7.2         (R2016a)
% Statistics and Machine Learning Toolbox               Version 10.2        (R2016a)

%% setup

clear all
close all

testing = 0; % run on a subset of images (1) or on all images (1)
showIm = 1; % generate figures (1) or not (0)
showKeyIm = 1; % generate key figures (1) or not (0)
smallMem = 0; % reduce memory usage (1) or not (0)
saveRes = 1; % save figures at end (1) or not (0)
doManual = 1; % do manual scoring (1) or not (0)

datadir = 'Example_Images/'; % folder with images
% store images in nested folders:
% datadir \ day \ group \ image.czi

fileID = '*.czi';

date = ['ZECAP_' char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))];
mkdir(['~/Desktop/' date]); % place to save outputs

cd(datadir)

datadir = [pwd '/'];

%% experiment parameters

% days fish were imaged
days = {'1dpt'};
realDays = [1]; % dpt of imaging

% group names
plates = {'A1', 'E1'};

% microscope channel information
BF = 1;
mCh1 = 2;
mCh2 = 3;
tdT1 = 4;
tdT2 = 5;
GFP1 = 6;
GFP2 = 7;
GFP3 = 8;

chName = {'BF', 'mCh1', 'mCh2', 'tdTomato1', 'tdTomato2', 'GFP1', 'GFP2', 'GFP3'};
chUse = GFP1; % data channel to use
chUseBkg = tdT1; % background channel to use

% segmentation thresholds for after background correction
thresh = NaN(1,length(chName));
thresh(GFP1) = 10;

% threshold for displaying images in figures
adjThresh = [thresh(GFP1)/(2^16-1) 0.01]; % will be reset later if rawGFP == 1

% image to use if testing == 1
testDay = 1;
testPlate = 1;
testFish = 1;

% get code version information
version = mfilename('fullpath');

%% 


% ************************************************************************
%
%
%
%           Load Files & Background Correction
%
%
%
% ************************************************************************



%% list of files

data = struct;

for j = 1:length(days)
    for k = 1:length(plates)
        cd([datadir days{j} '/' plates{k}])
        data(j).plate(k).fish = dir(fileID);
    end
end

%% load files into memory
tic

% use all files if not testing
if testing == 0
    testDay = 1:length(days);
    testPlate = 1:length(plates);
end

% read images into memory
for j = testDay
    for k = testPlate
        cd([datadir days{j} '/' plates{k}])
        for m = 1:length(data(j).plate(k).fish)
            % read image into memory
            data(j).plate(k).fish(m).Im = bfopen(data(j).plate(k).fish(m).name);
        end
    end
end

display('Load Files')
toc

%% Detect edge of each fish
tic

% segment BF images
for j = testDay
    for k = testPlate
        for m = 1:length(data(j).plate(k).fish)
            BW = edge(data(j).plate(k).fish(m).Im{1,1}{BF,1}, 'LoG');
            SE = strel('disk', 1);
            BW1 = imdilate(BW, SE);
            SE2 = strel('disk', 2);
            BW2 = bwareaopen(BW1, 140);
            BW3 = imdilate(BW2, SE2);
            BW6 = bwareaopen(~BW3, 8000);
            BW6 = bwareaopen(~BW6, 70000);
            SE3 = strel('octagon', 60);
            data(j).plate(k).fish(m).seg = imclose(BW6, SE3);
        end
    end
end

% Show Segmentation Results
if showIm == 1
    for j = testDay
        for k = testPlate
            figure
            for m = 1:length(data(j).plate(k).fish)
                n = ceil(sqrt(length(k)*length(data(j).plate(k).fish)));
                subplot(n,n,m)
                imshowpair(data(j).plate(k).fish(m).Im{1,1}{BF,1}, data(j).plate(k).fish(m).seg)
                title(['Day ' num2str(j) ', Plate ' plates{k} ', Image ' num2str(m)])
            end
        end
    end
end

display('Segment BF')
toc

%% crop images using BF segmentation
tic

% apply mask from BF image to all channels and store in "crop"
for j = testDay
    for k = testPlate
        for m = 1:length(data(j).plate(k).fish)
            data(j).plate(k).fish(m).crop = data(j).plate(k).fish(m).Im;
            for i = 1:length(chName)
                data(j).plate(k).fish(m).crop{1,1}{i,1}(data(j).plate(k).fish(m).seg == 0) = 0;
            end
        end
    end
end

% show results of mask for GFP1 channel
if showIm == 1
    for j = testDay
        for k = testPlate
            figure
            for m = 1:length(data(j).plate(k).fish)
                n = ceil(sqrt(length(k)*length(data(j).plate(k).fish)));
                subplot(n,n,m)
                imshow(imadjust(data(j).plate(k).fish(m).crop{1,1}{GFP1,1}))
                title(['Day ' num2str(j) ', Plate ' plates{k} ', Image ' num2str(m)])
            end
        end
    end
end

display('Crop Images')
toc

%% plot ratio vs. background
tic

All = struct;

for j = testDay
    % initialize variables
    for i = find(~strcmp(chName, 'BF'))
        All(j).ch(i).all = [];
    end
    for k = testPlate
        for m = 1:length(data(j).plate(k).fish)
            for i = find(~strcmp(chName, 'BF'))
                % compile all pixel values within area of fish
                All(j).ch(i).all = [All(j).ch(i).all; data(j).plate(k).fish(m).Im{1,1}{i,1}(data(j).plate(k).fish(m).seg ~= 0)];
            end
        end
    end
end

% make figure
if showIm == 1    
    figure
    scatter(All(j).ch(tdT1).all, (double(All(j).ch(GFP1).all) ./ double(All(j).ch(tdT1).all)), 'b.')
    title('GFP1/tdT1 vs. tdT1')
    ylabel('GFP1/tdT1')
    xlabel('tdT1')

    figure
    scatter(All(j).ch(mCh1).all, (double(All(j).ch(GFP1).all) ./ double(All(j).ch(mCh1).all)), 'b.')
    title('GFP1/mCh1 vs. mCh1')
    ylabel('GFP1/mCh1')
    xlabel('mCh1')

    figure
    scatter(All(j).ch(tdT1).all, (double(All(j).ch(GFP2).all) ./ double(All(j).ch(tdT1).all)), 'b.')
    title('GFP2/tdT1 vs. tdT1')
    ylabel('GFP2/tdT1')
    xlabel('tdT1')

    figure
    scatter(All(j).ch(tdT2).all, (double(All(j).ch(GFP2).all) ./ double(All(j).ch(tdT2).all)), 'b.')
    title('GFP2/tdT2 vs. tdT2')
    ylabel('GFP2/tdT2')
    xlabel('tdT2')
end


display('Compile Data & Plot')
toc

% will use channels from parameters for rest of code
% users can switch to different combinations of channels based on needs

%% Regression on ratio vs. tdT1
tic

mdl = struct;

% run a model for each experimental day
for j = testDay
    mdl(j).mdl = fitlm(double(All(j).ch(chUseBkg).all), double(All(j).ch(chUse).all) ./ double(All(j).ch(chUseBkg).all),'constant','RobustOpts','on');
end

% put model coefficients into matrix and then clear mdl
% arrangement of coeffs in table:
% mdl_coeff(data_ch, background_ch, day)

mdl_coeff = NaN(length(chName), length(chName), length(days));
for j = testDay
    mdl_coeff(chUse, chUseBkg, j) = mdl(j).mdl.Coefficients.Estimate;
end

clear mdl
display('Clear mdl')


% scale results of regression to be more conservative in background corr.
regScale = 1.7;

% make figure
if showIm == 1
    for j = testDay
        figure
        plot(All(j).ch(chUseBkg).all, (double(All(j).ch(chUse).all) ./ double(All(j).ch(chUseBkg).all)), 'b.')
        title([chName{chUse} '/' chName{chUseBkg} ' vs. ' chName{chUseBkg} ' || ' days{j}])
        ylabel([chName{chUse} '/' chName{chUseBkg}])
        xlabel(chName{chUseBkg})
        line([min(All(j).ch(chUseBkg).all(:)), max(All(j).ch(chUseBkg).all(:))], [mdl_coeff(chUse, chUseBkg, j), mdl_coeff(chUse, chUseBkg, j)], 'color', 'r', 'linewidth', 1.5)
        line([min(All(j).ch(chUseBkg).all(:)), max(All(j).ch(chUseBkg).all(:))], regScale.*[mdl_coeff(chUse, chUseBkg, j), mdl_coeff(chUse, chUseBkg, j)], 'color', 'm', 'linewidth', 1.5)
        legend({'Pixel Ratio', 'Fit', [num2str(regScale) ' * Fit']})
        set(gca, 'Fontsize', 16)
    end
end

% Free up some memory
if smallMem == 1
    clear('All')
end

display('Regression')
toc

%% Background Correction using regScale*mdl
tic

for j = testDay
    for k = testPlate
        for m = 1:length(data(j).plate(k).fish)
            data(j).plate(k).fish(m).cor = uint16(double(data(j).plate(k).fish(m).Im{1,1}{chUse,1}) - (regScale .* mdl_coeff(chUse, chUseBkg, j) .* double(data(j).plate(k).fish(m).Im{1,1}{chUseBkg,1})));
        end
    end
end

% show results
if showIm == 1
    for j = testDay
        for k = testPlate
            figure
            for m = 1:length(data(j).plate(k).fish)
                n = ceil(sqrt(length(k)*length(data(j).plate(k).fish)));
                subplot(n,n,m)
                imagesc(data(j).plate(k).fish(m).cor, [0,500])
                axis image
                colorbar
                title(['Day ' num2str(j) ', Plate ' plates{k} ', Image ' num2str(m)])
            end
        end
    end
end

display('Background Correction')
toc

%% Show example linearized image plots from last dish
if testing == 1
    j = testDay;
    k = testPlate;
end

mdl = mdl_coeff(chUse, chUseBkg, j);
range = [0, 5000];

if showIm == 1
    for m = 1:length(data(j).plate(k).fish)
        figure, 
        subplot(3,1,1), plot(data(j).plate(k).fish(m).Im{1,1}{chUse,1}(:), 'color', 'g'),
            title(['Day ' num2str(j) ', Plate ' plates{k} ', Image ' num2str(m) '  ||  ' chName{chUse}])
            set(gca, 'ylim', range)
        subplot(3,1,2), plot(regScale .* mdl .* double(data(j).plate(k).fish(m).Im{1,1}{chUseBkg,1}(:)), 'color', 'r'),
            title([chName{chUseBkg} ' * ' num2str(regScale * mdl)])
            set(gca, 'ylim', range)
        subplot(3,1,3), plot(data(j).plate(k).fish(m).cor(:), 'color', 'c'),
            title('Corrected')
            set(gca, 'ylim', range)
    end
end


%%


% ***************************************************************
%
%
%               Background Correction Complete
%
%               Extract data from images & compile
%               Make basic plots
%
%
% ***************************************************************



%% Total area above threshold

% count pixels above threshold within fish outline
for j = testDay
    for k = testPlate
        for m = 1:length(data(j).plate(k).fish)
            data(j).plate(k).fish(m).pxCountSeg = sum(sum(data(j).plate(k).fish(m).cor > thresh(chUse) & data(j).plate(k).fish(m).seg == 1));
        end
    end
end

%% Integrate intensity over area

% only within fish outline
for j = testDay
    for k = testPlate
        for m = 1:length(data(j).plate(k).fish)
            data(j).plate(k).fish(m).IntegSeg = sum(data(j).plate(k).fish(m).cor(data(j).plate(k).fish(m).cor > thresh(chUse) & data(j).plate(k).fish(m).seg == 1));
        end
    end
end

%% Compile into Table

% results = table;
r_day = [];
r_lin = [];
r_plates = [];
r_areaseg = [];
r_iaSeg = [];
r_fish = [];
r_ch = [];

for j = testDay
    for k = testPlate
        for n = 1:length(data(j).plate(k).fish)
            r_day = [r_day; j];
            r_lin = [r_lin; strncmp(plates(k), 'E',1)];
            r_plates = [r_plates; plates(k)];
            r_areaseg = [r_areaseg; data(j).plate(k).fish(n).pxCountSeg];
            r_iaSeg = [r_iaSeg; data(j).plate(k).fish(n).pxCountSeg];
            r_fish = [r_fish; n];
            r_ch = [r_ch; chName(chUse)];
        end
    end
end

results = table(r_day, r_plates, r_lin, r_ch, r_fish, r_areaseg, r_iaSeg);
results.Properties.VariableNames = {'Day', 'Plate', 'LIN', 'Ch', 'Fish', 'AreaSeg', 'I_Aseg'};

clear('r_*')

% add RealDay variable that has actual dpt, not indexed day
results.RealDay = results.Day;
for i = 1:length(realDays)
    results.RealDay(results.Day == i) = realDays(i);
end

results(1:5,:)


%% Histograms of tumor burden

if showKeyIm == 1
    % 1 -- GFP1 Area
    figure('position', [100,100,900,500])
    for j = testDay
        subplot(1, length(testDay), j)
        h1 = histogram(results.AreaSeg(results.Day == j & results.LIN == 0));
        hold on
        h2 = histogram(results.AreaSeg(results.Day == j & results.LIN == 1));

        h1.BinWidth = max(results.AreaSeg(results.Day == j))/50;
        h2.BinWidth = h1.BinWidth;

        legend('A', 'E', 'location', 'northeast')
    end
    title('AreaSeg Histogram')


    % 2 -- GFP1 I*A
    figure('position', [100,100,900,500])
    for j = testDay
        subplot(1, length(testDay), j)
        h1 = histogram(results.I_Aseg(results.Day == j & results.LIN == 0));
        hold on
        h2 = histogram(results.I_Aseg(results.Day == j & results.LIN == 1));

        h1.BinWidth = max(results.I_Aseg(results.Day == j))/50;
        h2.BinWidth = h1.BinWidth;

        legend('A', 'E', 'location', 'northeast')
    end
    title('I*Aseg Histogram')

end

%% Boxplots
if showKeyIm == 1
    chPlot = chName{chUse};
    
    % 1 -- AreaSeg by day, LIN
    figure
    boxplot(results.AreaSeg(strcmp(results.Ch, chPlot)), {results.Day(strcmp(results.Ch, chPlot)), results.LIN(strcmp(results.Ch, chPlot))})
    title([chPlot ' AreaSeg'], 'interpreter', 'none')
    ylabel('AreaSeg')
    set(gca, 'yscale', 'log')
    

    % 2 -- I_Aseg by day, LIN
    figure
    chPlot = 'GFP1';
    boxplot(results.I_Aseg(strcmp(results.Ch, chPlot)), {results.Day(strcmp(results.Ch, chPlot)), results.LIN(strcmp(results.Ch, chPlot))})
    title([chPlot ' I_Aseg'], 'interpreter', 'none')
    ylabel('AreaSeg * Intensity')
    set(gca, 'yscale', 'log')

    
end

%%


% ***************************************************************
%
%
%
%       spatial analysis -- change of basis
%
%
%
% ***************************************************************


%% find eyes
tic

% set parameters
eyeThresh = 0.0024; % threshold for segmenting eyes
ovrSeg = 4; % oversegmentation correction for watershed


for j = testDay
    for k = testPlate
        for m = 1:length(data(j).plate(k).fish)
            bw = im2bw(data(j).plate(k).fish(m).Im{1,1}{BF,1}, eyeThresh);
            SE_eye2 = strel('disk', 4);
            SE_eye = strel('disk', 5);
            SE_eyeE = strel('line', 10, 5);
            SE_eyeE2 = strel('line', 10, -5);

            bwc = imclose(imcomplement(bw), SE_eye2);
            bwE = imerode(bwc, SE_eyeE);
            bwE = imerode(bwE, SE_eyeE2);
            bwo = imclose(bwE, SE_eye);

            bwoo = ~bwareaopen(~bwo, 500);
            % keep it within fish
            bwoo_fish = logical(bwoo .* data(j).plate(k).fish(m).seg);
            bwoBig_eye = bwareaopen(bwoo_fish, 250);

            % watershed the two eyes (slower but better performance)
            BWdist = -bwdist(~double(bwoBig_eye));
            % watershed transform
            Ld = watershed(BWdist);
            % change image based on watershed
            bwoBig_eye2 = bwoBig_eye;
            bwoBig_eye2(Ld == 0) = 0;
            % correct oversegmentation
            mask = imextendedmin(BWdist,ovrSeg,4); % handle to adjust (Im, scalar, 4 or 8)
            % overlay masks onto size threshold BW image
            bwoBig_eye2 = imimposemin(double(bwoBig_eye), mask);
            Ld2 = watershed(bwoBig_eye2);
            bwoBig_eye3 = bwoBig_eye;
            bwoBig_eye3(Ld2 == 0) = 0;

            bwoBig_eye4 = bwareaopen(bwoBig_eye3, 250);

            s = regionprops(bwoBig_eye4,'centroid');
            data(j).plate(k).fish(m).centroids_eye = cat(1, s.Centroid);
            data(j).plate(k).fish(m).eyeMid = [mean(data(j).plate(k).fish(m).centroids_eye(:,1)), mean(data(j).plate(k).fish(m).centroids_eye(:,2))];
        end
    end
end

% show results
if showIm == 1
    for j = testDay
        for k = testPlate
            figure
            for m = 1:length(data(j).plate(k).fish)
                if ~isnan(data(j).plate(k).fish(m).Im{1,1}{BF,1})
                    n = ceil(sqrt(length(k)*length(data(j).plate(k).fish)));
                    subplot(n,n,m), 
                    imshow(imadjust(data(j).plate(k).fish(m).Im{1,1}{BF,1}))
                    hold on
                    % add eyes
                    if ~isempty(data(j).plate(k).fish(m).centroids_eye)
                        scatter(data(j).plate(k).fish(m).centroids_eye(:,1), data(j).plate(k).fish(m).centroids_eye(:,2), '*r')
                        scatter(data(j).plate(k).fish(m).eyeMid(:,1), data(j).plate(k).fish(m).eyeMid(:,2), '*c')
                    end
                    hold off
                    title(['D' num2str(j) ', P' num2str(k) ', F' num2str(m)])
                end
            end
        end
    end
end

display('Find Eyes')
toc

%% find yolk sac
tic

% set parameters
% threshold for segmenting yolk
yolkThresh = 26;

% convert yolkThresh to fraction
yolkThresh = yolkThresh / (2^16-1);

% background correction & segment
for j = testDay
    for k = testPlate
        for m = 1:length(data(j).plate(k).fish)
            
            % background correction
            gauss = imgaussfilt(data(j).plate(k).fish(m).Im{1,1}{tdT1,1}, 200);
            gcor = imsubtract(data(j).plate(k).fish(m).Im{1,1}{tdT1,1}, gauss);
            gcor_seg = uint16(double(gcor) .* double(data(j).plate(k).fish(m).seg));

            % find yolk
            bw = im2bw(gcor_seg, yolkThresh);
            SE_yolk = strel('disk', 14);
            bwc = imclose(bw, SE_yolk);
            bwcBig = ~bwareaopen(~bwc, 500);
            bwo = imopen(bwcBig, SE_yolk);
            bwoBig_yolk = bwareaopen(bwo, 4000);

            bwoBig_yolk2 = logical(data(j).plate(k).fish(m).seg .* bwoBig_yolk);

            s = regionprops(bwoBig_yolk2, 'centroid');
            data(j).plate(k).fish(m).centroids_yolk = cat(1, s.Centroid); % yolk centroid(s)

            data(j).plate(k).fish(m).yolkHull = regionprops(bwoBig_yolk2, 'ConvexHull'); % yolk hull
            data(j).plate(k).fish(m).yolkEdge = regionprops(edge(bwoBig_yolk2), 'PixelList'); % yolk edge
        end
    end
end

% show results
if showIm == 1
    for j = testDay
        for k = testPlate
            figure
            for m = 1:length(data(j).plate(k).fish)
                n = ceil(sqrt(length(k)*length(data(j).plate(k).fish)));
                subplot(n,n,m),
                imshow(imadjust(data(j).plate(k).fish(m).Im{1,1}{tdT1,1}))
                hold on
                % add yolk
                if ~isempty(data(j).plate(k).fish(m).centroids_yolk)
                    % yolk centroid
                    scatter(data(j).plate(k).fish(m).centroids_yolk(:,1), data(j).plate(k).fish(m).centroids_yolk(:,2), '*r')
                    % yolk hull
                    if length(data(j).plate(k).fish(m).yolkHull) == 1
                        plot(data(j).plate(k).fish(m).yolkHull.ConvexHull(:,1), data(j).plate(k).fish(m).yolkHull.ConvexHull(:,2), 'm')
                    end
                end
                % fish hull
%                     plot(data(j).plate(k).fish(m).hull.ConvexHull(:,1), data(j).plate(k).fish(m).hull.ConvexHull(:,2), 'b')
                % eyes
                if ~isempty(data(j).plate(k).fish(m).centroids_eye)
                    plot(data(j).plate(k).fish(m).eyeMid(1), data(j).plate(k).fish(m).eyeMid(2), '*c')
                    plot(data(j).plate(k).fish(m).centroids_eye(:,1), data(j).plate(k).fish(m).centroids_eye(:,2), '*y')
                end
                hold off
                title(['D' num2str(j) ', P' num2str(k) ', F' num2str(m)])
            end
        end
    end
end

display('Find Yolk')
toc

%% find orientation of fish (A1 and A2 only)
tic

for j = testDay
    for k = testPlate
        for m = 1:length(data(j).plate(k).fish)
            % find hull for display purposes
            data(j).plate(k).fish(m).hull = regionprops(data(j).plate(k).fish(m).seg, 'ConvexHull');


            % make rotation matrix "A1"
            angle = regionprops(data(j).plate(k).fish(m).seg, 'orientation'); % angle of major axis from horizontal
            data(j).plate(k).fish(m).angle = angle; % assign angle to data structure
            data(j).plate(k).fish(m).A1 = [cos(deg2rad(angle.Orientation)), sin(deg2rad(angle.Orientation));...
                                                -sin(deg2rad(angle.Orientation)), cos(deg2rad(angle.Orientation))]'; % http://mathworld.wolfram.com/RotationMatrix.html


            % Find which end eyes are on and make matrix "A2"
            box = regionprops(data(j).plate(k).fish(m).seg, 'BoundingBox');
            dist1 = [box.BoundingBox(1), box.BoundingBox(2)]; % left corner of box
            dist2 = [box.BoundingBox(1) + box.BoundingBox(3), box.BoundingBox(2)]; % right corner of box

            % make matrix
            if pdist2(data(j).plate(k).fish(m).eyeMid, dist1, 'euclidean') < pdist2(data(j).plate(k).fish(m).eyeMid, dist2, 'euclidean')
                % do nothing
                display('Eyes on Left')
                data(j).plate(k).fish(m).A2 = [1 0; 0 1];
            else
                % flip fish horizontally
                display('Eyes on Right')
                data(j).plate(k).fish(m).A2 = [-1 0; 0 1];
            end
        end
    end
end


% show anchor points
if showIm == 1
    for j = testDay
        for k = testPlate
            figure
            for m = 1:length(data(j).plate(k).fish)
                n = ceil(sqrt(length(k)*length(data(j).plate(k).fish)));
                subplot(n,n,m),
                imshow(imadjust(data(j).plate(k).fish(m).Im{1,1}{BF,1}))
                hold on
                % add eyes
                if ~isempty(data(j).plate(k).fish(m).eyeMid)
                    plot(data(j).plate(k).fish(m).eyeMid(1), data(j).plate(k).fish(m).eyeMid(2), '*c')
                    plot(data(j).plate(k).fish(m).centroids_eye(:,1), data(j).plate(k).fish(m).centroids_eye(:,2), '*y')
                end
                % add yolk
                if ~isempty(data(j).plate(k).fish(m).centroids_yolk)
                    plot(data(j).plate(k).fish(m).centroids_yolk(:,1), data(j).plate(k).fish(m).centroids_yolk(:,2), '*m')
                end
                hold off
                title(['D' num2str(j) ', P' num2str(k) ', F' num2str(m)])
            end
        end
    end
end


% show rotation results
if showIm == 1
    for j = testDay
        for k = testPlate
            figure
            for m = 1:length(data(j).plate(k).fish)
                n = ceil(sqrt(length(k)*length(data(j).plate(k).fish)));
                subplot(n,n,m),
                % show hull
                plot(data(j).plate(k).fish(m).hull.ConvexHull(:,1), data(j).plate(k).fish(m).hull.ConvexHull(:,2), 'r')
                hold on
                % show rotated hull
                A1_hull = (data(j).plate(k).fish(m).A1 * data(j).plate(k).fish(m).hull.ConvexHull')';
                plot(A1_hull(:,1), A1_hull(:,2), 'b')
                hold off
                legend({'hull', 'A1_hull'}, 'location', 'northeast', 'interpreter', 'none')
                title(['D' num2str(j) ', P' num2str(k) ', F' num2str(m)])
            end
        end
    end
end

display('Find Fish Orientation')
toc

%% Change of Basis

for j = testDay
    for k = testPlate
        for m = 1:length(data(j).plate(k).fish)

            % apply rotation & eye reflection (A2*A1*)
            A1 = data(j).plate(k).fish(m).A1;
            A2 = data(j).plate(k).fish(m).A2;
            eyeMid = data(j).plate(k).fish(m).eyeMid;

            % COB of landmarks
            % hull
            data(j).plate(k).fish(m).d_hull = (A2*A1 * ...
                [data(j).plate(k).fish(m).hull.ConvexHull(:,1) - eyeMid(1,1), ...
                data(j).plate(k).fish(m).hull.ConvexHull(:,2) - eyeMid(1,2)]')';
            % hull edge (make then COB)
            % make hull edge
            data(j).plate(k).fish(m).hullEdge = regionprops(edge(data(j).plate(k).fish(m).seg), 'PixelList');
            % COB hull edge
            data(j).plate(k).fish(m).d_hullEdge = (A2*A1 * ...
                [data(j).plate(k).fish(m).hullEdge.PixelList(:,1) - eyeMid(1,1), ...
                data(j).plate(k).fish(m).hullEdge.PixelList(:,2) - eyeMid(1,2)]')';
            % eye centroids
            data(j).plate(k).fish(m).d_eye = (A2*A1 * ...
                [data(j).plate(k).fish(m).centroids_eye(:,1) - eyeMid(1,1), ...
                data(j).plate(k).fish(m).centroids_eye(:,2) - eyeMid(1,2)]')';
            % eyeMid
            data(j).plate(k).fish(m).d_eyeMid = (A2*A1 * ...
                [data(j).plate(k).fish(m).eyeMid(:,1) - eyeMid(1,1), ...
                data(j).plate(k).fish(m).eyeMid(:,2) - eyeMid(1,2)]')';
            % yolk
            if ~isempty(data(j).plate(k).fish(m).centroids_yolk)
                % yolk centroid
                data(j).plate(k).fish(m).d_yolk = (A2*A1 * ...
                    [data(j).plate(k).fish(m).centroids_yolk(:,1) - eyeMid(1,1), ...
                    data(j).plate(k).fish(m).centroids_yolk(:,2) - eyeMid(1,2)]')';
                % yolk hull
                data(j).plate(k).fish(m).d_yolkHull = (A2*A1 * ...
                    [data(j).plate(k).fish(m).yolkHull(1).ConvexHull(:,1) - eyeMid(1,1), ...
                    data(j).plate(k).fish(m).yolkHull(1).ConvexHull(:,2) - eyeMid(1,2)]')';
                % yolk edge
                data(j).plate(k).fish(m).d_yolkEdge = (A2*A1 * ...
                    [data(j).plate(k).fish(m).yolkEdge(1).PixelList(:,1) - eyeMid(1,1), ...
                    data(j).plate(k).fish(m).yolkEdge(1).PixelList(:,2) - eyeMid(1,2)]')';
            end

            % find which fish edge the yolk edge is closer to
            % want to have yolk on positive side of fish

            % if not 2 eyes and 1 yolk, request user input
            if size(data(j).plate(k).fish(m).centroids_eye, 1) ~= 2 || size(data(j).plate(k).fish(m).centroids_yolk, 1) ~= 1
                % request user input for ventral/dorsal
                figure(999), imshow(imadjust(data(j).plate(k).fish(m).Im{1,1}{BF,1}))
                title('Click Top/Bottom For Ventral (Yolk) Side and L/R For Eye Side')
                [x,y] = ginput(1);

                % Manual Yolk Direction
                if y < 1088/2
                    display('Yolk Up = Yolk (-)')
                    data(j).plate(k).fish(m).d_A3 = [1 0; 0 -1];
                elseif y > 1088/2
                    display('Yolk Down = Yolk (+)')
                    data(j).plate(k).fish(m).d_A3 = [1 0; 0 1];
                end

                % Manual Eye Direction
                if x < 1388/2
                    display('Eyes Left')
                    data(j).plate(k).fish(m).A2 = [1 0; 0 1];
                elseif x> 1388/2
                    display('Eyes Right')
                    data(j).plate(k).fish(m).A2 = [-1 0; 0 1];
                end

                close figure 999
            else

                hullLim = 20;

                if min(min(pdist2(data(j).plate(k).fish(m).d_hullEdge(data(j).plate(k).fish(m).d_hullEdge(:,2) < -hullLim,:),...
                    data(j).plate(k).fish(m).d_yolkEdge)))...
                    < min(min(pdist2(data(j).plate(k).fish(m).d_hullEdge(data(j).plate(k).fish(m).d_hullEdge(:,2) > hullLim,:),...
                    data(j).plate(k).fish(m).d_yolkEdge)))

                    % yolk is closer to negative side of fish
                    % need to flip fish vertically
                    display('Yolk (-) = Yolk Up')
                    data(j).plate(k).fish(m).d_A3 = [1 0; 0 -1];

                else
                    % yolk is closer to positive side of fish
                    % don't need to flip fish (identity matrix)
                    display('Yolk (+) = Yolk Down')
                    data(j).plate(k).fish(m).d_A3 = [1 0; 0 1];
                end
            end


            % apply yolk reflection (A3*)

            A3 = data(j).plate(k).fish(m).d_A3;

            % COB of landmarks
            % hull
            data(j).plate(k).fish(m).dd_hull = (A3 * data(j).plate(k).fish(m).d_hull')';
            % fish edge
            data(j).plate(k).fish(m).dd_hullEdge = (A3 * data(j).plate(k).fish(m).d_hullEdge')';
            % eye centroids
            data(j).plate(k).fish(m).dd_eye = (A3 * data(j).plate(k).fish(m).d_eye')';
            % eyeMid
            data(j).plate(k).fish(m).dd_eyeMid = (A3 * data(j).plate(k).fish(m).d_eyeMid')';
            % yolk
            if ~isempty(data(j).plate(k).fish(m).centroids_yolk)
                % yolk centroid
                data(j).plate(k).fish(m).dd_yolk = (A3 * data(j).plate(k).fish(m).d_yolk')';
                % yolk hull
                data(j).plate(k).fish(m).dd_yolkHull = (A3 * data(j).plate(k).fish(m).d_yolkHull')';
                % yolk edge
                data(j).plate(k).fish(m).dd_yolkEdge = (A3 * data(j).plate(k).fish(m).d_yolkEdge')';
            end

        end
    end
end


% show A2*A1
if showIm == 1
    for j = testDay
        for k = testPlate
            figure
            for m = 1:length(data(j).plate(k).fish)
                n = ceil(sqrt(length(k)*length(data(j).plate(k).fish)));
                subplot(n,n,m),
                hold on
                % add features
                if ~isempty(data(j).plate(k).fish(m).d_yolk)
                    % yolk centroid
                    scatter(data(j).plate(k).fish(m).d_yolk(:,1), data(j).plate(k).fish(m).d_yolk(:,2), '*r')
                    % yolk hull
                    plot(data(j).plate(k).fish(m).d_yolkHull(:,1), data(j).plate(k).fish(m).d_yolkHull(:,2), 'm', 'linewidth', 1.5)
                    % fish hull
                    plot(data(j).plate(k).fish(m).d_hull(:,1), data(j).plate(k).fish(m).d_hull(:,2), 'b', 'linewidth', 1.5)
                    % eyes
                    plot(data(j).plate(k).fish(m).d_eyeMid(1), data(j).plate(k).fish(m).d_eyeMid(2), '*c')
                    plot(data(j).plate(k).fish(m).d_eye(:,1), data(j).plate(k).fish(m).d_eye(:,2), '*y')
                end
                hold off
                title(['A2*A1 - D' num2str(j) ', P' num2str(k) ', F' num2str(m)])
            end
        end
    end
end


% show A3*A2*A1
if showIm == 1
    for j = testDay
        for k = testPlate
            figure
            for m = 1:length(data(j).plate(k).fish)
                n = ceil(sqrt(length(k)*length(data(j).plate(k).fish)));
                subplot(n,n,m),
                hold on
                % add features
                if ~isempty(data(j).plate(k).fish(m).dd_yolk)
                    % yolk centroid
                    scatter(data(j).plate(k).fish(m).dd_yolk(:,1), data(j).plate(k).fish(m).dd_yolk(:,2), '*r')
                    % yolk hull
                    plot(data(j).plate(k).fish(m).dd_yolkHull(:,1), data(j).plate(k).fish(m).dd_yolkHull(:,2), 'm', 'linewidth', 1.5)
                    % yolk edge
                    plot(data(j).plate(k).fish(m).dd_yolkEdge(:,1), data(j).plate(k).fish(m).dd_yolkEdge(:,2), '*k')
                    % fish hull
                    plot(data(j).plate(k).fish(m).dd_hull(:,1), data(j).plate(k).fish(m).dd_hull(:,2), 'b', 'linewidth', 1.5)
                    % fish edge
                    plot(data(j).plate(k).fish(m).dd_hullEdge(:,1), data(j).plate(k).fish(m).dd_hullEdge(:,2), '*k')
                    % eyes
                    plot(data(j).plate(k).fish(m).dd_eyeMid(1), data(j).plate(k).fish(m).dd_eyeMid(2), '*c')
                    plot(data(j).plate(k).fish(m).dd_eye(:,1), data(j).plate(k).fish(m).dd_eye(:,2), '*y')
                end
                hold off
                title(['A3*A2*A1 - D' num2str(j) ', P' num2str(k) ', F' num2str(m)])
            end
        end
    end
end


%% Eye & Yolk QC
display(' ')
display(' ')
display('The following fish were flagged:')

% eyes
display(' ')
for j = testDay
    for k = testPlate
        for m = 1:length(data(j).plate(k).fish)
            if size(data(j).plate(k).fish(m).centroids_eye,1) ~= 2
                display(['D' num2str(j) ', P' num2str(k) ', F' num2str(m) ' Eyes ' num2str(size(data(j).plate(k).fish(m).centroids_eye,1))])
            end
        end
    end
end


% yolk
display(' ')
for j = testDay
    for k = testPlate
        for m = 1:length(data(j).plate(k).fish)
            if size(data(j).plate(k).fish(m).centroids_yolk,1) ~= 1
                display(['D' num2str(j) ', P' num2str(k) ', F' num2str(m) ' Yolks ' num2str(size(data(j).plate(k).fish(m).centroids_yolk,1))])
            end
        end
    end
end


%% find centroids of tumors in autofluorescence-corrected images (within fish outline)
tic

% Centroid for Area
% WeightedCentroid for I_A

for j = testDay
    for k = testPlate
        for m = 1:length(data(j).plate(k).fish)

            % make BW using threshold from AreaSeg measurement
            bw = (im2bw(data(j).plate(k).fish(m).cor, (thresh(chUse) / ...
                (2^16-1)))) .* data(j).plate(k).fish(m).seg;
            L = bwlabel(bw);

            % concatenate properties into table
            % apply mask from BW to corrected image
            data(j).plate(k).fish(m).tumor_stats = regionprops('table', L, data(j).plate(k).fish(m).cor,...
                'Centroid', 'WeightedCentroid', 'Area', 'MeanIntensity', 'Solidity', 'Perimeter');
        end
    end
end


% show AreaSeg results
if showIm == 1
    for j = testDay
        for k = testPlate
            figure
            for m = 1:length(data(j).plate(k).fish)
                n = ceil(sqrt(length(k)*length(data(j).plate(k).fish)));
                subplot(n,n,m),
                imshow(imadjust(data(j).plate(k).fish(m).cor, adjThresh, []))
                hold on
                % show centroids
                if ~isempty(data(j).plate(k).fish(m).tumor_stats)
                    scatter(data(j).plate(k).fish(m).tumor_stats.Centroid(:,1), ...
                        data(j).plate(k).fish(m).tumor_stats.Centroid(:,2), '.r')
                end
                hold off
                title(['GFP D' num2str(j) ', P' num2str(k) ', F' num2str(m)])
            end
        end
    end
end


% show I_Aseg results
if showIm == 1
    for j = testDay
        for k = testPlate
            figure
            for m = 1:length(data(j).plate(k).fish)
                n = ceil(sqrt(length(k)*length(data(j).plate(k).fish)));
                subplot(n,n,m),
                imshow(imadjust(data(j).plate(k).fish(m).cor, adjThresh, []))
                hold on
                % show weighted centroids
                if ~isempty(data(j).plate(k).fish(m).tumor_stats)
                    scatter(data(j).plate(k).fish(m).tumor_stats.WeightedCentroid(:,1), ...
                        data(j).plate(k).fish(m).tumor_stats.WeightedCentroid(:,2), '.r')
                end
                hold off
                title(['GFP D' num2str(j) ', P' num2str(k) ', F' num2str(m)])
            end
        end
    end
end

display('Find Tumor Centroids')
toc

%% Change of basis to fish space using A1, A2, d_A3, & eyeMid
tic
% b = A2*d_A3*A1 * (data - eye)

for j = testDay
    for k = testPlate
        for m = 1:length(data(j).plate(k).fish)

            % put COB info into temp variables
            A1 = data(j).plate(k).fish(m).A1;
            A2 = data(j).plate(k).fish(m).A2;
            A3 = data(j).plate(k).fish(m).d_A3;
            eyeMid = data(j).plate(k).fish(m).eyeMid;

            % COB of tumor centroids
            if ~isempty(data(j).plate(k).fish(m).tumor_stats)
                % Centroids (AreaSeg)
                % re-center using eyeMid
                data(j).plate(k).fish(m).tumor_stats.Centroid_A = ...
                    [data(j).plate(k).fish(m).tumor_stats.Centroid(:,1) - eyeMid(1,1), ...
                    data(j).plate(k).fish(m).tumor_stats.Centroid(:,2) - eyeMid(1,2)];
                % rotate/reflect using A1,A2,A3
                data(j).plate(k).fish(m).tumor_stats.Centroid_A = ...
                    (A2*A3*A1 * data(j).plate(k).fish(m).tumor_stats.Centroid_A')';

                % WeightedCentroids (I_Aseg)
                % re-center using eyeMid
                data(j).plate(k).fish(m).tumor_stats.WeightedCentroid_A = ...
                    [data(j).plate(k).fish(m).tumor_stats.WeightedCentroid(:,1) - eyeMid(1,1), ...
                    data(j).plate(k).fish(m).tumor_stats.WeightedCentroid(:,2) - eyeMid(1,2)];
                % rotate/reflect using A1,A2,A3
                data(j).plate(k).fish(m).tumor_stats.WeightedCentroid_A = ...
                    (A2*A3*A1 * data(j).plate(k).fish(m).tumor_stats.WeightedCentroid_A')';

                % make I*A
                data(j).plate(k).fish(m).tumor_stats.I_A = ...
                    data(j).plate(k).fish(m).tumor_stats.Area .* data(j).plate(k).fish(m).tumor_stats.MeanIntensity;
            end
        end
    end
end


% show GFP1 (centroid, AreaSeg) results
if showIm == 1
    plotScale1 = 0.05;
    for j = testDay
        for k = testPlate
            figure
            for m = 1:length(data(j).plate(k).fish)
                n = ceil(sqrt(length(k)*length(data(j).plate(k).fish)));
                subplot(n,n,m),
                % add hull
                plot(data(j).plate(k).fish(m).dd_hull(:,1), data(j).plate(k).fish(m).dd_hull(:,2), 'b', 'linewidth', 2)
                hold on
                % add tumor centroids (scaled b area)
                if ~isempty(data(j).plate(k).fish(m).tumor_stats)
                    scatter(data(j).plate(k).fish(m).tumor_stats.Centroid_A(:,1), ...
                        data(j).plate(k).fish(m).tumor_stats.Centroid_A(:,2),...
                        plotScale1.*data(j).plate(k).fish(m).tumor_stats.Area, 'filled', 'og')
                end
                % add yolk
                if ~isempty(data(j).plate(k).fish(m).centroids_yolk)
                    % centroid
                    plot(data(j).plate(k).fish(m).dd_yolk(:,1), data(j).plate(k).fish(m).dd_yolk(:,2), '*r')
                    % hull
                    plot(data(j).plate(k).fish(m).dd_yolkHull(:,1), data(j).plate(k).fish(m).dd_yolkHull(:,2), 'r', 'linewidth', 2)
                end
                % add eye and eye mid
                plot(data(j).plate(k).fish(m).dd_eye(:,1), data(j).plate(k).fish(m).dd_eye(:,2), '*k')
                plot(data(j).plate(k).fish(m).dd_eyeMid(:,1), data(j).plate(k).fish(m).dd_eyeMid(:,2), '*c')
                hold off
                xlabel('FS1 (px)')
                ylabel('FS2 (px)')
                set(gca, 'xlim', [-200, 1200])
                set(gca, 'ylim', [-250, 250])
                title(['D' num2str(j) ', ' plates{k} ', F' num2str(m) ' AreaSeg'])
            end
        end
    end
end

% show GFP1 (WeightedCentroid, I_Aseg) results
if showIm == 1
    plotScale2 = 0.0005;
    for j = testDay
        for k = testPlate
            figure
            for m = 1:length(data(j).plate(k).fish)
                n = ceil(sqrt(length(k)*length(data(j).plate(k).fish)));
                subplot(n,n,m),
                % add hull
                plot(data(j).plate(k).fish(m).dd_hull(:,1), data(j).plate(k).fish(m).dd_hull(:,2), 'b', 'linewidth', 2)
                hold on
                % add tumor weighted centroids (scaled by I_A)
                if ~isempty(data(j).plate(k).fish(m).tumor_stats)
                    scatter(data(j).plate(k).fish(m).tumor_stats.WeightedCentroid_A(:,1), ...
                            data(j).plate(k).fish(m).tumor_stats.WeightedCentroid_A(:,2),...
                        plotScale2.*data(j).plate(k).fish(m).tumor_stats.I_A, 'filled', 'og')
                end
                % add yolk
                if ~isempty(data(j).plate(k).fish(m).centroids_yolk)
                    % centroid
                    plot(data(j).plate(k).fish(m).dd_yolk(:,1), data(j).plate(k).fish(m).dd_yolk(:,2), '*r')
                    plot(data(j).plate(k).fish(m).dd_yolkHull(:,1), data(j).plate(k).fish(m).dd_yolkHull(:,2), 'r', 'linewidth', 2)
                end
                % add eye and eye mid
                plot(data(j).plate(k).fish(m).dd_eye(:,1), data(j).plate(k).fish(m).dd_eye(:,2), '*k')
                plot(data(j).plate(k).fish(m).dd_eyeMid(:,1), data(j).plate(k).fish(m).dd_eyeMid(:,2), '*c')
                hold off
                xlabel('FS1 (px)')
                ylabel('FS2 (px)')
                set(gca, 'xlim', [-200, 1200])
                set(gca, 'ylim', [-250, 250])
                title(['D' num2str(j) ', ' plates{k} ', F' num2str(m) ' I_Aseg'])
            end
        end
    end
end

display('Change of Basis')
toc


%% Compile tumor_stats into Table
tic

r_day = [];
r_LIN = [];
r_plates = [];
r_fish = [];
r_ch = [];
r_area = [];
r_I_A = [];
r_centroid_x = [];
r_centroid_y = [];
r_wcentroid_x = [];
r_wcentroid_y = [];
r_solidity = [];
r_perim = [];
r_centroid_Ax = [];
r_centroid_Ay = [];
r_wcentroid_Ax = [];
r_wcentroid_Ay = [];


for j = testDay
    for k = testPlate
        for n = 1:length(data(j).plate(k).fish)
            if ~isempty(data(j).plate(k).fish(n).tumor_stats)
                for q = 1:length(data(j).plate(k).fish(n).tumor_stats.Area)
                    r_day = [r_day; j];
                    r_LIN = [r_LIN; isempty(strfind(plates{k}, 'E'))]; % E(1) or A(0)
                    r_plates = [r_plates; plates(k)];
                    r_fish = [r_fish; n];
                    r_ch = [r_ch; chName(chUse)];
                    r_area = [r_area; data(j).plate(k).fish(n).tumor_stats.Area(q)];
                    r_I_A = [r_I_A; data(j).plate(k).fish(n).tumor_stats.I_A(q)];
                    r_centroid_x = [r_centroid_x; data(j).plate(k).fish(n).tumor_stats.Centroid(q,1)];
                    r_centroid_y = [r_centroid_y; data(j).plate(k).fish(n).tumor_stats.Centroid(q,2)];
                    r_wcentroid_x = [r_wcentroid_x; data(j).plate(k).fish(n).tumor_stats.WeightedCentroid(q,1)];
                    r_wcentroid_y = [r_wcentroid_y; data(j).plate(k).fish(n).tumor_stats.WeightedCentroid(q,2)];
                    r_solidity = [r_solidity; data(j).plate(k).fish(n).tumor_stats.Solidity(q)];
                    r_perim = [r_perim; data(j).plate(k).fish(n).tumor_stats.Perimeter(q)];
                    r_centroid_Ax = [r_centroid_Ax; data(j).plate(k).fish(n).tumor_stats.Centroid_A(q,1)];
                    r_centroid_Ay = [r_centroid_Ay; data(j).plate(k).fish(n).tumor_stats.Centroid_A(q,2)];
                    r_wcentroid_Ax = [r_wcentroid_Ax; data(j).plate(k).fish(n).tumor_stats.WeightedCentroid_A(q,1)];
                    r_wcentroid_Ay = [r_wcentroid_Ay; data(j).plate(k).fish(n).tumor_stats.WeightedCentroid_A(q,2)];
                end
            else
                r_day = [r_day; j];
                r_LIN = [r_LIN; isempty(strfind(plates{k}, 'E'))]; % E(1) or A(0)
                r_plates = [r_plates; plates(k)];
                r_fish = [r_fish; n];
                r_ch = [r_ch; chName(chUse)];
                r_area = [r_area; NaN];
                r_I_A = [r_I_A; NaN];
                r_centroid_x = [r_centroid_x; NaN];
                r_centroid_y = [r_centroid_y; NaN];
                r_wcentroid_x = [r_wcentroid_x; NaN];
                r_wcentroid_y = [r_wcentroid_y; NaN];
                r_solidity = [r_solidity; NaN];
                r_perim = [r_perim; NaN];
                r_centroid_Ax = [r_centroid_Ax; NaN];
                r_centroid_Ay = [r_centroid_Ay; NaN];
                r_wcentroid_Ax = [r_wcentroid_Ax; NaN];
                r_wcentroid_Ay = [r_wcentroid_Ay; NaN];
            end
        end
    end
end

tumor_stats = table(r_day, r_LIN, r_plates, r_fish, r_ch, r_area, r_I_A, r_centroid_x, r_centroid_y, r_wcentroid_x, r_wcentroid_y, ...
    r_solidity, r_perim, r_centroid_Ax, r_centroid_Ay, r_wcentroid_Ax, r_wcentroid_Ay);
tumor_stats.Properties.VariableNames = {'Day', 'LIN', 'Plate', 'Fish', 'Ch', 'Area', 'I_A', 'Centroid_x', 'Centroid_y', 'WeightedCentroid_x', 'WeightedCentroid_y',...
    'Solidity', 'Perimeter', 'Centroid_Ax', 'Centroid_Ay', 'WeightedCentroid_Ax', 'WeightedCentroid_Ay'};

clear('r_*')

tumor_stats(1:10,:)


%% Calculate magnitude of centroid vectors (eyeDist)

tumor_stats.EyeDist = sqrt(((tumor_stats.Centroid_Ax).^2) + ((tumor_stats.Centroid_Ay).^2));
tumor_stats.WeightedEyeDist = sqrt(((tumor_stats.WeightedCentroid_Ax).^2) + ((tumor_stats.WeightedCentroid_Ay).^2));


%% Calculate Center of mass

summaryStats = grpstats(results, {'Ch', 'Day', 'LIN'}, {'mean', 'meanci', 'std', 'median'}, 'DataVars', {'AreaSeg', 'I_Aseg'});

summaryStats.xCOM = NaN(length(summaryStats.Day),1);
summaryStats.yCOM = NaN(length(summaryStats.Day),1);
summaryStats.xCOMweighted = NaN(length(summaryStats.Day),1);
summaryStats.yCOMweighted = NaN(length(summaryStats.Day),1);

for j = testDay
    for k = unique(tumor_stats.LIN)'
        % x center of mass
        summaryStats.xCOM(summaryStats.Day == j & summaryStats.LIN == k & strcmp(summaryStats.Ch, chName{chUse})) = ...
            nansum(tumor_stats.Centroid_Ax(tumor_stats.Day == j & tumor_stats.LIN == k & strcmp(tumor_stats.Ch, chName{chUse})) .* ...
            tumor_stats.Area(tumor_stats.Day == j & tumor_stats.LIN == k & strcmp(tumor_stats.Ch, chName{chUse})) / ...
            nansum(tumor_stats.Area(tumor_stats.Day == j & tumor_stats.LIN == k & strcmp(tumor_stats.Ch, chName{chUse}))));
        % y center of mass
        summaryStats.yCOM(summaryStats.Day == j & summaryStats.LIN == k & strcmp(summaryStats.Ch, chName{chUse})) = ...
            nansum(tumor_stats.Centroid_Ay(tumor_stats.Day == j & tumor_stats.LIN == k & strcmp(tumor_stats.Ch, chName{chUse})) .* ...
            tumor_stats.Area(tumor_stats.Day == j & tumor_stats.LIN == k & strcmp(tumor_stats.Ch, chName{chUse})) / ...
            nansum(tumor_stats.Area(tumor_stats.Day == j & tumor_stats.LIN == k & strcmp(tumor_stats.Ch, chName{chUse}))));
        % weighted x center of mass
        summaryStats.xCOMweighted(summaryStats.Day == j & summaryStats.LIN == k & strcmp(summaryStats.Ch, chName{chUse})) = ...
            nansum(tumor_stats.WeightedCentroid_Ax(tumor_stats.Day == j & tumor_stats.LIN == k & strcmp(tumor_stats.Ch, chName{chUse})) .* ...
            tumor_stats.I_A(tumor_stats.Day == j & tumor_stats.LIN == k & strcmp(tumor_stats.Ch, chName{chUse})) / ...
            nansum(tumor_stats.I_A(tumor_stats.Day == j & tumor_stats.LIN == k & strcmp(tumor_stats.Ch, chName{chUse}))));
        % weighted y center of mass
        summaryStats.yCOMweighted(summaryStats.Day == j & summaryStats.LIN == k & strcmp(summaryStats.Ch, chName{chUse})) = ...
            nansum(tumor_stats.WeightedCentroid_Ay(tumor_stats.Day == j & tumor_stats.LIN == k & strcmp(tumor_stats.Ch, chName{chUse})) .* ...
            tumor_stats.I_A(tumor_stats.Day == j & tumor_stats.LIN == k & strcmp(tumor_stats.Ch, chName{chUse})) / ...
            nansum(tumor_stats.I_A(tumor_stats.Day == j & tumor_stats.LIN == k & strcmp(tumor_stats.Ch, chName{chUse}))));
    end
end


%%


% ***********************************************************************
%
%
%
%                   Plots of Change of Basis
%
%
% ***********************************************************************


%% plot overlaid tumors by group

if showKeyIm == 1
    % Area
    plotScale3 = 0.4;
    transparency = 0.7;
    color = [0, 0.5, 0];
    figure('position', [100, 100, 1200, 900])
    for j = testDay
        for k = unique(tumor_stats.LIN)'
            % add outline
            subplot(length(testDay),length(unique(tumor_stats.LIN)),(j-1)*length(unique(tumor_stats.LIN))+k+1),
            plot(data(j).plate(1).fish(1).dd_hull(:,1), data(j).plate(1).fish(1).dd_hull(:,2), 'b')
            hold on
            % add tumors
            t = scatter(tumor_stats.Centroid_Ax(tumor_stats.Day == j & tumor_stats.LIN == k & strcmp(tumor_stats.Ch, chName{chUse})),...
                tumor_stats.Centroid_Ay(tumor_stats.Day == j & tumor_stats.LIN == k & strcmp(tumor_stats.Ch, chName{chUse})), ...
                plotScale3 .* tumor_stats.Area(tumor_stats.Day == j & tumor_stats.LIN == k & strcmp(tumor_stats.Ch, chName{chUse})), color, 'filled');
            alpha(t, transparency)
            title([chName{chUse} ' Area -- D' num2str(j) ', LIN ' num2str(k)])
            set(gca, 'xlim', [-200, 1200])
            set(gca, 'ylim', [-250, 250])
            % add COM
            scatter(summaryStats.xCOM(summaryStats.Day == j & summaryStats.LIN == k & strcmp(summaryStats.Ch, chName{chUse})) ,...
                summaryStats.yCOM(summaryStats.Day == j & summaryStats.LIN == k & strcmp(summaryStats.Ch, chName{chUse})), '*r')
        end
    end
    
    % I_A GFP2
    plotScale4 = 1e-3;
    color = [0, 0.5, 0];
    figure('position', [100, 100, 1200, 900])
    for j = testDay
        for k = unique(tumor_stats.LIN)'
            % add outline
            subplot(length(testDay),length(unique(tumor_stats.LIN)),(j-1)*length(unique(tumor_stats.LIN))+k+1),
            plot(data(j).plate(1).fish(1).dd_hull(:,1), data(j).plate(1).fish(1).dd_hull(:,2), 'b')
            hold on
            % add tumors
            t = scatter(tumor_stats.WeightedCentroid_Ax(tumor_stats.Day == j & tumor_stats.LIN == k & strcmp(tumor_stats.Ch, chName{chUse})),...
                tumor_stats.WeightedCentroid_Ay(tumor_stats.Day == j & tumor_stats.LIN == k & strcmp(tumor_stats.Ch, chName{chUse})), ...
                plotScale4 .* tumor_stats.I_A(tumor_stats.Day == j & tumor_stats.LIN == k & strcmp(tumor_stats.Ch, chName{chUse})), color, 'filled');
            alpha(t, 2/length(unique(tumor_stats.Fish(tumor_stats.Day == j & tumor_stats.LIN == k & strcmp(tumor_stats.Ch, chName{chUse})))))
            title([chName{chUse} ' I_A -- D' num2str(j) ', LIN ' num2str(k)])
            set(gca, 'xlim', [-200, 1200])
            set(gca, 'ylim', [-250, 250])
            % add COM
            scatter(summaryStats.xCOMweighted(summaryStats.Day == j & summaryStats.LIN == k & strcmp(summaryStats.Ch, chName{chUse})) ,...
                summaryStats.yCOMweighted(summaryStats.Day == j & summaryStats.LIN == k & strcmp(summaryStats.Ch, chName{chUse})), '*r')
        end
    end
end


%%


% ***************************************************************
%
%
%
%       Tumor Burden By Area
%
%
%
% ***************************************************************



%% Transform Images to Standardized Orientation (with cropping)
tic

for j = testDay
    for k = testPlate
        for m = 1:length(data(j).plate(k).fish)

            % find box to crop image
            seg_rot = imrotate(data(j).plate(k).fish(m).seg, -data(j).plate(k).fish(m).angle.Orientation);
            seg_box = regionprops(seg_rot, 'BoundingBox');

            % rotate (A1) and crop image in each channel
            Im_rot = [];
            Im_rot_crop = [];
            seg_rot = [];

            % rotate all channels
            Im_rot = data(j).plate(k).fish(m).Im;
            Im_rot_crop = Im_rot;
            for n = 1:length(chName)
                Im_rot{1,1}{n,1} = imrotate(data(j).plate(k).fish(m).Im{1,1}{n,1}, -data(j).plate(k).fish(m).angle.Orientation, 'nearest');
                Im_rot_crop{1,1}{n,1} = imcrop(Im_rot{1,1}{n,1}, seg_box.BoundingBox);
            end
            
            
            % cor
            cor_rot.Im = imrotate(data(j).plate(k).fish(m).cor, -data(j).plate(k).fish(m).angle.Orientation, 'nearest');
            cor_rot_crop.Im = imcrop(cor_rot.Im, seg_box.BoundingBox);
            
            
            % BF seg
            seg_rot.Im = imrotate(data(j).plate(k).fish(m).seg, -data(j).plate(k).fish(m).angle.Orientation, 'nearest');
            seg_rot_crop.Im = imcrop(seg_rot.Im, seg_box.BoundingBox);


            % flip images horizontally (A2)
            Im_rot_crop_h = Im_rot_crop;
            if data(j).plate(k).fish(m).A2(1,1) == -1
                for n = 1:length(chName)
                    Im_rot_crop_h{1,1}{n,1} = flip(Im_rot_crop{1,1}{n,1}, 2);
                end
                % cor
                cor_rot_crop_h.Im = flip(cor_rot_crop.Im, 2);
                % BF seg
                seg_rot_crop_h.Im = flip(seg_rot_crop.Im, 2);
            else
                % cor
                cor_rot_crop_h.Im = cor_rot_crop.Im;
                % BF seg
                seg_rot_crop_h.Im = seg_rot_crop.Im;
            end


            % flip images vertically (A3)
            Im_rot_crop_hv = Im_rot_crop_h;
            if data(j).plate(k).fish(m).d_A3(2,2) == -1
                for n = 1:length(chName)
                    Im_rot_crop_hv{1,1}{n,1} = flip(Im_rot_crop_h{1,1}{n,1}, 1);
                end
                % cor
                cor_rot_crop_hv.Im = flip(cor_rot_crop_h.Im, 1);
                % BF seg
                seg_rot_crop_hv.Im = flip(seg_rot_crop_h.Im, 1);
            else
                % cor
                cor_rot_crop_hv.Im = cor_rot_crop_h.Im;
                % BF seg
                seg_rot_crop_hv.Im = seg_rot_crop_h.Im;
            end


            % save to data structure
            data(j).plate(k).fish(m).rchv = Im_rot_crop_hv;
            % cor
            data(j).plate(k).fish(m).rchv_cor = cor_rot_crop_hv.Im;
            % BF seg
            data(j).plate(k).fish(m).rchv_seg = seg_rot_crop_hv.Im;

        end
    end
end

if showKeyIm == 1
    for j = testDay
        for k = testPlate
            figure
            for m = 1:length(data(j).plate(k).fish)
                n = ceil(sqrt(length(k)*length(data(j).plate(k).fish)));
                subplot(n,n,m)
                imshow(imadjust(data(j).plate(k).fish(m).rchv{1,1}{BF,1}))
                title(['Day ' num2str(j) ', Plate ' plates{k} ', Image ' num2str(m)])
            end
        end
    end
end


display('Transform images')
toc


%% show transformed images with tumor overlay

if showKeyIm == 1
    for j = testDay
        for k = testPlate
            figure('Position', [100, 100, 1500, 600])
            for m = 1:length(data(j).plate(k).fish)
                n = ceil(sqrt(length(k)*length(data(j).plate(k).fish)));
                subplot(n,n,m)
                i1 = imshow(imadjust(data(j).plate(k).fish(m).rchv{1,1}{BF,1}));
                hold on
                % show GFP image with transparency
                
                % convert image to RGB
                tmpIm = data(j).plate(k).fish(m).rchv_cor;
                
                zzz = zeros(size(tmpIm));
                hhh = [];
                hhh(:,:,1) = uint8(zzz);
                hhh(:,:,2) = uint8(zzz);
                hhh(:,:,3) = uint8(zzz);

                hhh(:,:,2) = uint8((imadjust(tmpIm, adjThresh, [])./(2^16-1)).*255);
                
                i2 = imshow(hhh);
                alpha(i2, 0.6);
                hold off
                
                title(['Day ' num2str(j) ', Plate ' plates{k} ', Fish ' num2str(m)])
            end
        end
    end
end


%% calculate tumor burden binned by region
tic
% number of regions
nReg = 2;

% count pixels above threshold
for j = testDay
    for k = testPlate
        for m = 1:length(data(j).plate(k).fish)
            stepSize = ceil(size(data(j).plate(k).fish(m).rchv_cor,2)/nReg);
            p = 1;
            tmp = [];
            while p < nReg
                tmp(p).Im = data(j).plate(k).fish(m).rchv_cor(:,((p-1)*stepSize+1):(p*stepSize));
                % pixel counts
                data(j).plate(k).fish(m).rchv_pxCount(p).pxCount = length(tmp(p).Im(tmp(p).Im > thresh(chUse)));
                % pixel integration
                data(j).plate(k).fish(m).rchv_pxCount(p).Integ = sum(tmp(p).Im(tmp(p).Im > thresh(chUse)));

                p = p+1;
            end

                % do final region
                tmp(p).Im = data(j).plate(k).fish(m).rchv_cor(:,((p-1)*stepSize+1):end);
                % pixel counts
                data(j).plate(k).fish(m).rchv_pxCount(p).pxCount = length(tmp(p).Im(tmp(p).Im > thresh(chUse)));
                % pixel integration
                data(j).plate(k).fish(m).rchv_pxCount(p).Integ = sum(tmp(p).Im(tmp(p).Im > thresh(chUse)));
        end
    end
end

display('Tumor Burden by Region')
toc

%% Compile into Table
% new column for each region of fish (wide)

results3 = results;

for p = 1:nReg
    % results = table;
    r_rchv_count = [];
    r_rchv_integ = [];

    for j = testDay
        for k = testPlate
            for m = 1:length(data(j).plate(k).fish)
                r_rchv_count = [r_rchv_count; data(j).plate(k).fish(m).rchv_pxCount(p).pxCount];
                r_rchv_integ = [r_rchv_integ; data(j).plate(k).fish(m).rchv_pxCount(p).Integ];
            end
        end
    end

    results3 = [results3, array2table(r_rchv_count), array2table(r_rchv_integ)];
    results3.Properties.VariableNames{end-1} = ['rchv_count' num2str(p)];
    results3.Properties.VariableNames{end} = ['rchv_integ' num2str(p)];

    clear('r_*')

end

%% Compile into Table
% column describing which region of fish (long)

r_day = [];
r_type = [];
r_plates = [];
r_ch = [];
r_fish = [];
r_fishID = [];
r_reg = [];
r_area = [];
r_ia = [];

for j = testDay
    for k = testPlate
        for n = 1:length(data(j).plate(k).fish)
            for p = 1:nReg
                r_day = [r_day; j];
                r_type = [r_type; strncmp(plates(k), 'E', 1)]; % E (1) or A (0)
                r_plates = [r_plates; plates(k)];
                r_ch = [r_ch; chName(chUse)];
                r_fish = [r_fish; n];
                r_fishID = [r_fishID; (k-1)*12+n];
                r_reg = [r_reg; p];
                r_area = [r_area; data(j).plate(k).fish(n).rchv_pxCount(p).pxCount];
                r_ia = [r_ia; data(j).plate(k).fish(n).rchv_pxCount(p).Integ];
            end
        end
    end
end

results2 = table(r_day, r_type, r_plates, r_ch, r_fish, r_fishID, r_reg, r_area, r_ia);
results2.Properties.VariableNames = {'Day', 'LIN', 'Plate', 'Ch', 'Fish', 'FishID', 'Reg', 'Area', 'I_A'};

clear('r_*')

results2.RealDay = results2.Day;
for i = 1:length(realDays)
    results2.RealDay(results2.Day == i) = realDays(i);
end

results2(1:10,:)


%% Boxplots

if showKeyIm == 1
    chPlot = chName(chUse);
    
    % 1 -- GFP1 Area by day, region
    figure('position', [100, 100, 1600, 900])
    boxplot(results2.Area(strcmp(results2.Ch, chPlot)), {results2.Day(strcmp(results2.Ch, chPlot)), results2.Reg(strcmp(results2.Ch, chPlot))})
    title([chPlot ' Area'], 'interpreter', 'none')
    ylabel([chPlot ' Area'])
    xlabel('Day | Region')
    set(gca, 'yscale', 'log')

    % 2 -- GFP1 Integ by day, region
    figure('position', [100, 100, 1600, 900])
    boxplot(results2.I_A(strcmp(results2.Ch, chPlot)), {results2.Day(strcmp(results2.Ch, chPlot)), results2.Reg(strcmp(results2.Ch, chPlot))})
    title([chPlot ' I_A'], 'interpreter', 'none')
    ylabel([chPlot ' I_A'])
    xlabel('Day | Region')
    set(gca, 'yscale', 'log')

    % 3 -- GFP1 Area by day, LIN, region
    figure('position', [100, 100, 1600, 900])
    boxplot(results2.Area(strcmp(results2.Ch, chPlot)), ...
        {results2.Day(strcmp(results2.Ch, chPlot)), results2.Reg(strcmp(results2.Ch, chPlot)), results2.LIN(strcmp(results2.Ch, chPlot))})
    title([chPlot ' Area'], 'interpreter', 'none')
    ylabel([chPlot ' Area'])
    xlabel('Day | Region | LIN')
    set(gca, 'yscale', 'log')

    % 4 -- GFP1 Integ by day, LIN, region
    figure('position', [100, 100, 1600, 900])
    boxplot(results2.I_A(strcmp(results2.Ch, chPlot)), ...
        {results2.Day(strcmp(results2.Ch, chPlot)), results2.Reg(strcmp(results2.Ch, chPlot)), results2.LIN(strcmp(results2.Ch, chPlot))})
    title([chPlot ' I_A'], 'interpreter', 'none')
    ylabel([chPlot ' I_A'])
    xlabel('Day | Region | LIN')
    set(gca, 'yscale', 'log')
    
end


%%


% ***************************************************************
%
%
%
%       Save Outputs
%
%
%
% ***************************************************************


%% save version info

if saveRes == 1
    verFile = fopen(['~/Desktop/' date '/version.log'], 'w');
    fprintf(verFile, version);
    fclose(verFile);
end


%% Save all figures to files

if saveRes == 1
    h = get(0,'children');
    
    for i = 1:length(h)
        figure(h(i).Number)
        print(['~/desktop/' date '/figure_' num2str(h(i).Number, '%03i') '.png'], '-dpng')
    end

end
