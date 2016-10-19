function W = calcWeights4MR(inputImg, superpixels, spNum, theta)
% Function W = calcWeights4MR(input_img, superpixels, spno, theta)
% calculates the weight matrix for manifold-ranking-based background and
% foreground saliency approximations.
% 
% Inputs:       input_img - the input RGB image
%                   superpixels - the superpixel segmentation result
%                   spno - total superpixel number
%                   theta - controlling parameter
% Outputs:    W - the weight matrix
% 
% 27/04/15 - by Yuchen Yuan
% Based on the paper:
% Yuchen Yuan, Changyang Li, Jinman Kim, Weidong Cai, and David Dagan Feng. 
% "Saliency Detection via Reversion Correction and Regularized Random Walks
% Ranking". TPAMI

    if nargin < 4
        theta = 10;
    end
    
    [m,n,k] = size(inputImg);
    
% Calculate average Lab value of each superpixel
    inputPixels=reshape(inputImg, m*n, k);
    spRGB=zeros(spNum,1,3);
    for i=1:spNum
        spRGB(i,1,:)=mean(inputPixels((superpixels==i),:),1);
    end  
    spLab = colorSpace('Lab<-', spRGB); 
    spLab=reshape(spLab,spNum,3);
 
% Calculate edges    
    neighb = zeros(spNum);
    for i = 1:m-1
        for j = 1:n-1
            if(superpixels(i,j)~=superpixels(i,j+1))
                neighb(superpixels(i,j),superpixels(i,j+1)) = 1;
                neighb(superpixels(i,j+1),superpixels(i,j)) = 1;
            end;
            if(superpixels(i,j)~=superpixels(i+1,j))
                neighb(superpixels(i,j),superpixels(i+1,j)) = 1;
                neighb(superpixels(i+1,j),superpixels(i,j)) = 1;
            end;
            if(superpixels(i,j)~=superpixels(i+1,j+1))
                neighb(superpixels(i,j),superpixels(i+1,j+1)) = 1;
                neighb(superpixels(i+1,j+1),superpixels(i,j)) = 1;
            end;
            if(superpixels(i+1,j)~=superpixels(i,j+1))
                neighb(superpixels(i+1,j),superpixels(i,j+1)) = 1;
                neighb(superpixels(i,j+1),superpixels(i+1,j)) = 1;
            end;
        end;
    end;    
    bd=unique([superpixels(1,:),superpixels(m,:),superpixels(:,1)',superpixels(:,n)']);
    for i=1:length(bd)
        for j=i+1:length(bd)
            neighb(bd(i),bd(j))=1;
            neighb(bd(j),bd(i))=1;
        end
    end

    edges=[];
    for i=1:spNum
        indexT=[];
        ind=find(neighb(i,:)==1);
        for j=1:length(ind)
            indj=find(neighb(ind(j),:)==1);
            indexT=[indexT,indj];
        end
        indexT=[indexT,ind];
        indexT=indexT((indexT>i));
        indexT=unique(indexT);
        if(~isempty(indexT))
            ed=ones(length(indexT),2);
            ed(:,2)=i*ed(:,2);
            ed(:,1)=indexT;
            edges=[edges;ed];
        end
    end

% Calculate weight matrix
    weights=exp(-theta*normalize(sqrt(sum((spLab(edges(:,1),:)-spLab(edges(:,2),:)).^2,2))));
    W=sparse([edges(:,1);edges(:,2)],[edges(:,2);edges(:,1)],[weights;weights],spNum,spNum);
end
