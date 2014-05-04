%finds WM regions, differentiating between matrix and striosome regions in
%calculations
%
%INPUT: folder of cropped MOR1 images
%OUTPUT: folder of: grayscale small WM, BW all WM, BW large WM, and BW
%small WM images

close all
clear all
N = 20;
low_in_increment = 0.001;
high_in_increment = 0.0001;
nonblackfrac_target = 0.2;
nonblack_areathresh_percent = 0.1;
medPixVal_target = 0.5;
foldername = 'WMB Processed test 2 originals';
outputfolder = 'WMB Processed test 2';
wmbimages = {};
names = {};

files = dir(fullfile(foldername,'*.tif'));

mkdir(outputfolder); 

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
    
    level = graythresh(current);
    currentbw = im2bw(current,level);
    currentbw = imfill(currentbw, 'holes');
    currentbw = bwareaopen(currentbw, round(0.985*pi*5^2));
    Nstriopx = sum(sum(currentbw));
    Nmatpx = sum(sum(~currentbw));
    
    black = current == 0;
    
    pixlists = regionprops(black, 'PixelList');
    [numRegions, ~] = size(pixlists);
    
    blackedge = {};
    
    %finding black borders
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
    
    %separating strio & matrix
    striowm = zeros(a,b);
    matwm = zeros(a,b);
    for i = 1:a
        for j = 1:b
            if currentbw(i,j) == 1
                striowm(i,j) = currentinv(i,j);
            else
                matwm(i,j) = currentinv(i,j);
            end
        end
    end
    invs = {striowm, matwm};

    results = cell(1,8);
    
    for t = 1:length(invs) %loop through strio and matrix separately 
        
        currentinv = invs{t};
        numignoredpx_edges = 0;
        
        %remove black border regions
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
        nonblack_areathresh = nonblack_areathresh_percent*(a*b - numignoredpx_edges);
        
        %decreases low_in until nonblackfrac >= nonblackfrac_target
        while nonblackfrac < nonblackfrac_target
            low_in = low_in - low_in_increment;
            if low_in < 0
                break
            end
            current2 = imadjust(currentinv,[low_in;high_in], [0;1]);
            nonblackimg = current2 > 0;
            pixlists = regionprops(nonblackimg, 'PixelList');
            [numRegions, ~] = size(pixlists);
            numignoredpx2 = numignoredpx_edges;
            
            for m = 1:numRegions
                currentPixList = pixlists(m).PixelList;
                [currentUnwantedArea, ~] = size(currentPixList);                
                if currentUnwantedArea > nonblack_areathresh               
                    numignoredpx2 = numignoredpx2 + currentUnwantedArea;
                end              
            end
            
            nonblackimg = current2 > 0;
            low_in_previous = low_in;
            nonblackfrac_previous = nonblackfrac;
            if t == 1
                Ntotalpx = Nstriopx;
            elseif t == 2
                Ntotalpx = Nmatpx;
            end
            nonblackfrac = sum(nonblackimg(:))/(Ntotalpx - numignoredpx2);
        end        
        if abs(nonblackfrac_target - nonblackfrac_previous) < abs(nonblackfrac_target - nonblackfrac)
            nonblackfrac = nonblackfrac_previous;
            low_in = low_in_previous;
        end
        
        results{1,4*t - 1} = nonblackimg; %all wm bw
        
        %make large WM/ventricles black in both images to be saved
        largewm = zeros(size(current2));
        nonblackimg = current2 > 0;
        pixlists = regionprops(nonblackimg, 'PixelList');
        [numRegions, ~] = size(pixlists);
        for m = 1:numRegions
            currentPixList = pixlists(m).PixelList;
            [currentUnwantedArea, ~] = size(currentPixList);
            
            if currentUnwantedArea > nonblack_areathresh
                for n = 1:currentUnwantedArea
                    curY = currentPixList(n,1);
                    curX = currentPixList(n,2);
                    current2(curX, curY) = 0;
                    nonblackimg(curX, curY) = 0;   
                    largewm(curX,curY) = 1;
                end
            end            
        end
        
        results{1,4*t} = largewm; %large wm bw
        
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
        
        results{1,4*t - 3} = current2; %wmb gray
        results{1,4*t - 2} = nonblackimg; %small wm bw
        
        display(nonblackfrac);
        display(medPixVal);        
    end
    
    finalwmb_gray = imadd(results{1,1},results{1,5});
    finalwmb_bw = results{1,2} | results{1,6};
    finalwm_bw = results{1,3} | results{1,7};
    finallarge_wm = results{1,4} | results{1,8};
    
    newname = [names{k},'_wmb.png'];
    newname2 = [names{k},'_small_wm_bw.png'];
    newname3 = [names{k},'_wm_bw.png'];
    newname4 = [names{k},'_large_wm_bw.png'];
    imwrite(finalwmb_gray,[outputfolder,'/', newname], 'png');    
    imwrite(finalwmb_bw,[outputfolder,'/', newname2], 'png');  
    imwrite(finalwm_bw,[outputfolder,'/', newname3], 'png'); 
    imwrite(finallarge_wm,[outputfolder,'/', newname4], 'png'); 
end