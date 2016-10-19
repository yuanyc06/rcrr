
function bSPSal = contextAware(img, superpixels)

[m,n,z] = size(img);
img_pixels = reshape(img, m*n, z);
% c = 3;

bRow = {1, m, 1:m, 1:m};
bCol = {1:n, 1:n, 1, n};
bSP = [];
    
for i = 1:4
    bSP = union(bSP, superpixels(bRow{i}, bCol{i}));
end

nBSP = length(bSP);

fColor = zeros(nBSP, 3);
% fPosition = zeros(nBSP, 2);

for i = 1:nBSP
    fColor(i, :) = mean(img_pixels((superpixels == bSP(i)), :), 1);
%     [row, col] = find(superpixels == bSP(i));
%     fPosition(i, :) = [mean(col), mean(row)]; % x, y
end
% fColor = colorspace('Lab<-', fColor);
% fPosition = fPosition/max(m,n);
% 
% dColor = zeros(nBSP, nBSP);
% dPosition = zeros(nBSP, nBSP);
% 
% for i = 1:nBSP
%     for j = 1:nBSP
%         dColor(i, j) = sqrt(sum((fColor(i,:)-fColor(j,:)).^2));
%         dPosition(i, j) = sqrt(sum((fPosition(i,:)-fPosition(j,:)).^2));
% %         d(i,j) = dColor(i,j)/(1+c*dPosition(i,j));
%     end
% end
% dColor = normalize(dColor);
% d = dColor./(1+c*dPosition);
% D = sum(d);
% S = 1 - exp(-D/nBSP);
% S = round((S - min(S))/(max(S) - min(S)) * 255);

% bSPSalShow = zeros(size(superpixels));
% for i = 1:length(bSP)
%     bSPSalShow(superpixels == bSP(i)) = S(i);
% end
% figure, imshow(bSPSalShow, [0 255]);

% if max(S) - min(S) > 0.5*max(S)
    fColor = reshape(fColor, nBSP, 1, 3);
    errKmeans = Inf;
    for clusterNum = 2:3
        [label, vec] = kMeans(fColor, clusterNum);
        err = 0;
        for j = 1:nBSP
            err = err + sqrt(sum(squeeze(fColor(j,1,:)) - vec(:,label(j))).^2);
        end
        if err < errKmeans;
            errKmeans = err;
            clusterNumFinal = clusterNum;
            labelFinal = label;
        end
    end
    
    for j = 1:clusterNumFinal
        numOf(j) = sum(labelFinal == j);
    end
    
    [val, ind] = min(numOf);
    if val < (nBSP / clusterNumFinal / 2)
        labelInd = 1:clusterNumFinal;
        labelInd(ind) = [];
        bSPSal = bSP(ismember(labelFinal, labelInd));
    else
        bSPSal = bSP;
    end
% 
% else
%     bSPSal = [];
% end
    