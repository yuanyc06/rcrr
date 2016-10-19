function salOutput = rcrr_saliency(imDir, imNameIn)
% Function salOutput = saliency_RCRR(imDir, imNameIn) executes the 
% proposed regularized random walk ranking algorithm and output the saliency map
% result. See demo.m for an example of usage.
% 
% Inputs:       imDir - path of image input image
%                   imNameIn - name of the input image (without path and
%                   suffix); only *.jpg files are supported
% Outputs:    salOutput - the pixel-wise saliency map with the same size as
%                   the input image
% 
% 27/04/15 - by Yuchen Yuan
% Based on the paper:
% Yuchen Yuan, Changyang Li, Jinman Kim, Weidong Cai, and David Dagan Feng. 
% "Saliency Detection via Reversion Correction and Regularized Random Walks
% Ranking". TPAMI

% Parameter initialization & superpixel segmentation    
    spNum = 200; 
    theta = 10;
    alpha = 0.99;    
    beta = 90;
    reverseTh = 1.5;
    imName = fullfile(imDir, [imNameIn(1:end-4) '.jpg']);
    imgBmpName = [imName(1:(end-4)), '.bmp'];
    [img, wid] = removeFrame(imName);
    img = uint8(img*255);
    [m, n, ~] = size(img);
    comm = [fullfile('support','SLIC.exe') ' ' imgBmpName ' ' int2str(20) ' ' int2str(spNum) ' '];
    if ispc
        evalc('system(comm);');
    elseif isunix
        comm = ['wine ',comm];
        evalc('unix(comm);');
    else
        error('Only Windows and Linux systems are currently supported.');
    end
    spName = [imgBmpName(1:end-4)  '.dat'];
    superpixels = readDat([m,n], spName);
    spNum = max(superpixels(:));
    delete(imgBmpName);
    delete(spName);
    delete([spName(1:end-4) '_SLIC.bmp']);
    
% Step 1: background saliency approximation
    W = calcWeights4MR(img, superpixels, spNum, theta);
    d = sum(W); 
    D = sparse(1:spNum,1:spNum,d); 
    clear d;
    A = (D-alpha*W)\eye(spNum); 
    A = A.*(~diag(ones(spNum,1)));
    bgdRow = {1, m, 1:m, 1:m};
    bgdCol = {1:n, 1:n, 1, n};
    boundPixels = zeros(size(superpixels));  
    
    bgdSal = zeros(spNum, 4);
    for i = 1:4
        y = zeros(spNum,1);
        b = unique(superpixels(bgdRow{i}, bgdCol{i}));
        y(b) = 1;
        salTmp = A*y;
        salTmp = (salTmp-min(salTmp(:)))/(max(salTmp(:))-min(salTmp(:)));  
        bgdSal(:,i) = 1 - salTmp;        
        boundPixels(ismember(superpixels, b)) = 1;
    end
    salStep1 = bgdSal(:,1).*bgdSal(:,2).*bgdSal(:,3).*bgdSal(:,4);
	salStep1 = (salStep1-min(salStep1(:)))/(max(salStep1(:))-min(salStep1(:)));

% Step 2: foreground saliency approximation
	th = mean(salStep1);
	salStep1(salStep1<th) = 0;
	salStep1(salStep1>=th) = 1;
	salStep2 = A*salStep1;  
    salStep2 = (salStep2-min(salStep2(:)))/(max(salStep2(:))-min(salStep2(:)));
    %Detect if the saliency estimation result is reversed
    salMapStep2 = superpixels;
    for i = 1:spNum
        salMapStep2(salMapStep2 == i) = salStep2(i);
    end
    [labelKmeans, ~] = kMeans(uint8(salMapStep2*255), 2);
    reverse = mean(labelKmeans(boundPixels==1));
   if (reverse > reverseTh)
        for i = 1:4
            y = zeros(spNum,1);
            b = superpixels(bgdRow{i}, bgdCol{i});
            boundKmeansLabels = labelKmeans(bgdRow{i}, bgdCol{i});
            spToDiscard = unique(b(boundKmeansLabels == 1));
            b = unique(b);
            b(ismember(b, spToDiscard)) = [];
            if isempty(b)
                bgdSal(:,i) = ones(spNum, 1);
            else
                y(b) = 1;
                salTmp = A*y;
                salTmp = (salTmp-min(salTmp(:)))/(max(salTmp(:))-min(salTmp(:)));  
                bgdSal(:,i) = 1 - salTmp;
            end
        end
        salStep1 = bgdSal(:,1).*bgdSal(:,2).*bgdSal(:,3).*bgdSal(:,4);
        salStep1 = (salStep1-min(salStep1(:)))/(max(salStep1(:))-min(salStep1(:)));
        th = mean(salStep1);
        salStep1(salStep1<th) = 0;
        salStep1(salStep1>=th) = 1;
        salStep2 = A*salStep1;  
        salStep2 = (salStep2-min(salStep2(:)))/(max(salStep2(:))-min(salStep2(:)));
   end

% Step 3: regularized random walk ranking
	th1 = (mean(salStep2) + max(salStep2)) / 2;
	th2 = mean(salStep2);
	mu = (1-alpha) / alpha;
	[seedAll, label] = seed4RW(salStep2, th1, th2);

    img = im2double(img);
    for i = 1:length(seedAll)
        [seedY, seedX] = find(superpixels == seedAll(i));
        seedXM = round(mean(seedX));
        seedYM = round(mean(seedY));
        seedAll(i) = (seedXM - 1) * m + seedYM;
    end
    seedAllLabel = sortrows([seedAll,label]);
    diffSeedAllLabel = diff(seedAllLabel(:,1));
    seedAllLabel(diffSeedAllLabel == 0, :) = [];
    seedAll = seedAllLabel(:,1);
    label = seedAllLabel(:,2);    

    %Build graph
    N=m*n;
    edges=[(1:N)',((1:N)+1)'];
    edges=[[edges(:,1);(1:N)'],[edges(:,2);(1:N)'+m]];
    excluded=find((edges(:,1)>N)|(edges(:,1)<1)|(edges(:,2)>N)| ...
        (edges(:,2)<1));
    edges([excluded;(m:m:((n-1)*m))'],:)=[]; 

    %Generate weights and Laplacian matrix
    tmp = img(:,:,1);
    imgVals = tmp(:);
    tmp = img(:,:,2);
    imgVals(:,2) = tmp(:);
    tmp = img(:,:,3);
    imgVals(:,3) = tmp(:);
    
    dis = sqrt(sum((imgVals(edges(:,1),:) - imgVals(edges(:,2),:)).^2, 2));
    dis = normalize(dis);
    weights = exp(-beta*dis);
    N=max(max(edges));
    W=sparse([edges(:,1);edges(:,2)],[edges(:,2);edges(:,1)],[weights;weights],N,N);
    D = diag(sum(W));
    L=D-W;

    %Determine which label values have been used
    labelAdjust = min(label); 
    label = label - labelAdjust + 1; %Adjust labels to be > 0
    labelRecord(label) = 1;
    labelPresent = find(labelRecord);
    labelNum = length(labelPresent);

    %Set up Dirichlet problem
    bound = zeros(length(seedAll), labelNum);
    for k=1:labelNum
        bound(:,k) = (label(:) == labelPresent(k));
    end

    %Solve the combinatorial Dirichlet problem
    saliencyFull = zeros(size(superpixels));
    for i = 1:length(salStep2)
        saliencyFull(superpixels == i) = salStep2(i);
    end
    y = saliencyFull(:);
    y = (y - min(y)) / (max(y) - min(y));
    y = [y, (1-y)];
    index = seedAll(:);
    N = length(L);
    antiIndex = 1:N;
    antiIndex(index) = [];
    antiY = y;
    antiY(index,:) = [];
    b = -L(antiIndex,index)*(bound);
    muI = sparse((1:length(antiIndex)), (1:length(antiIndex)), mu*ones(length(antiIndex),1));
    x = (L(antiIndex,antiIndex) + muI)\(mu * antiY + b);
    probabilities = zeros(size(bound));
    probabilities(index,:) = bound;
    probabilities(antiIndex,:) = x;
    probabilities = reshape(probabilities,[m n labelNum]);

    % Wrap up result for output
	sal = probabilities(:,:,1);
    salOutput = zeros(wid(1),wid(2));
	salOutput(wid(3):wid(4),wid(5):wid(6)) = sal;
end