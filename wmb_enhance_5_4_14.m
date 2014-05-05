%Finds WM regions in entire image (no strio/matrix distinction)
%
%INPUT: folder of cropped MOR1 images
%OUTPUT: folder of: grayscale small WM, BW all WM, BW large WM, and BW
%small WM images

clear all
close all

N = 20;
low_in_increment = 0.001;
high_in_increment = 0.001;
nonblackfrac_target = 0.2;
nonblack_areathresh_percent = 0.1;
medPixVal_target = 0.5;
foldername = 'WMB Processed test 2 originals';
outputfolder = 'WMB Processed test 2';
wmbimages = {};
names = {};

files = dir(fullfile(foldername,'*.tif'));

mkdir(outputfolder); %creates a new folder for processed images (used later)

for k = 1:length(files)
    img = imread(fullfile([foldername,'/',files(k).name]));
    wmbimages{end + 1} = im2double(img);
    filename = files(k).name;
    [path, name, ext] = fileparts(filename);
    names{end+1} = name;
end

for k = 1:length(wmbimages)
    current = wmbimages{k};
    [a,b] = size(current);
    
    black = current == 0; 
    
    pixlists = regionprops(black, 'PixelList');
    [numRegions, ~] = size(pixlists);
    
    blackedge = {};
    
    %finds black borders
    for m = 1:numRegions
        edge = 0;
        currentPixList = pixlists(m).PixelList;
        [currentUnwantedArea, ~] = size(currentPixList);
        for n = 1:currentUnwantedArea
            curY = currentPixList(n,1);
            curX = currentPixList(n,2);
            if curY == 1 || curY == b || curX == 1 || curX == a
                edge = 1;                
            end            
        end
        if edge == 1
            blackedge{end + 1} = currentPixList;
        end
    end
    
    currentinv = imcomplement(current);
    
    numignoredpx_edges = 0;
    
    for m = 1:length(blackedge)
        currentPixList = blackedge{m};
        [currentUnwantedArea, ~] = size(currentPixList);
        numignoredpx_edges = numignoredpx_edges + currentUnwantedArea;
        for n = 1:currentUnwantedArea
            curY = currentPixList(n,1);
            curX = currentPixList(n,2);
            currentinv(curX,curY) = 0;
        end
    end

    low_in = 1;
    high_in = 1;
    nonblackfrac = 0;
    nonblack_areathresh = nonblack_areathresh_percent*(a*b - numignoredpx);

    %decreases low_in until nonblackfrac >= nonblackfrac_target
    while nonblackfrac < nonblackfrac_target
        low_in = low_in - low_in_increment;
        current2 = imadjust(currentinv,[low_in;high_in], [0;1]);
        nonblackimg = current2 > 0;
        pixlists = regionprops(nonblackimg, 'PixelList');
        [numRegions, ~] = size(pixlists);
        numignoredpx2 = numignoredpx_edges;
        for m = 1:numRegions
            edge = 0;
            currentPixList = pixlists(m).PixelList;
            [currentUnwantedArea, ~] = size(currentPixList);
            if currentUnwantedArea > nonblack_areathresh
                for n = 1:currentUnwantedArea
                    curY = currentPixList(n,1);
                    curX = currentPixList(n,2);
                    if curY == 1 || curY == b || curX == 1 || curX == a
                        edge = 1;                
                    end            
                end
                if edge == 1
                    current2(curX,curY) = 0;                    
                    numignoredpx2 = numignoredpx2 + currentUnwantedArea;
                end
            end
        end
        nonblackimg = current2 > 0;
        low_in_previous = low_in;
        nonblackfrac_previous = nonblackfrac;
        nonblackfrac = sum(nonblackimg(:))/(a*b - numignoredpx2);
    end
    if abs(nonblackfrac_target - nonblackfrac_previous) < abs(nonblackfrac_target - nonblackfrac)
        nonblackfrac = nonblackfrac_previous;
        low_in = low_in_previous;
    end
    
    %saves all wm image before large regions are removed
    newname3 = [names{k},'_wm_bw.png'];
    imwrite(nonblackimg,[outputfolder,'/', newname3], 'png');
    
    %makes large WM/ventricles black in both images to be saved
    largewm = zeros(size(current2));
    nonblackimg = current2 > 0;
    pixlists = regionprops(nonblackimg, 'PixelList');
    [numRegions, ~] = size(pixlists);
    for m = 1:numRegions
        edge = 0;
        currentPixList = pixlists(m).PixelList;
        [currentUnwantedArea, ~] = size(currentPixList);
        if currentUnwantedArea > nonblack_areathresh
            for n = 1:currentUnwantedArea
                curY = currentPixList(n,1);
                curX = currentPixList(n,2);
                if curY == 1 || curY == b || curX == 1 || curX == a
                    edge = 1;                
                end            
            end
            if edge == 1
                current2(curX,curY) = 0;
                nonblackimg(curX,curY) = 0;
                largewm(curX,curY) = 1;
            end
        end
    end
    
    %decreases high_in until medPixVal >= medPixVal_target
    medPixVal = median(current2(current2 > 0));
    if medPixVal < medPixVal_target
        while medPixVal < medPixVal_target
            high_in = high_in - high_in_increment;
            if low_in >= high_in
                break
            end
            current2 = imadjust(currentinv,[low_in;high_in], [0;1]);
            high_in_previous = high_in;
            medPixVal_previous = medPixVal;
            medPixVal = median(current2(current2 > 0));
        end
        if abs(medPixVal_target - medPixVal_previous) < abs(medPixVal_target - medPixVal)
            medPixVal = medPixVal_previous;
            high_in = high_in_previous;
        end
    end   
    
    %saves output images
    newname4 = [names{k},'_large_wm_bw.png'];
    imwrite(largewm,[outputfolder,'/', newname4], 'png');

    newname = [names{k},'_wmb.png'];
    newname2 = [names{k},'_small_wm_bw.png'];   
    imwrite(current2,[outputfolder,'/', newname], 'png');  
    imwrite(nonblackimg,[outputfolder,'/', newname2], 'png');  
    
    display(nonblackfrac);
    display(medPixVal);
end