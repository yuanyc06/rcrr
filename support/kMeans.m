%% Fast K means Algorithm for clustering a Gray Image or Color Image
% It uses Preallocation and parallel operations to optimize algorithm time.
function [labelIm,vecMean] = kMeans(img,clusterNum,varargin)
% IM input Image. NO_OF_CLUSTER is number of cluster.
% VARARGIN will define Colorspace if it is empty RGB colorspace will be
% choose for clustering, or if its not empty than L*a*b* Colorspace will be
% used for clustering the image.
% LABEL_IM out put clustered Image. VEC_MEAN are Centers of
% Corresponding Clusters.

% Written By Ankit Dixit.
% ankitdx@gmail.com.

[m,n,p] = size(img);

%% ================= For Gray Images ====================================
if p<3
    gray = img;
    
    minimum = min(gray(:));
    vector=double((gray(:)-minimum)+1);% 1
    vector = repmat(vector,[1,clusterNum]);
    vecMean=(1:clusterNum).*max((vector))/(clusterNum+1);
    num = length(vector);
    itr = 0;
    %================ for gray image ==========================
    while(true)
        itr = itr+1;
        old_mean=vecMean;
        vecMean = repmat(vecMean,[num,1]);
        distance=(((vector-vecMean)).^2);
        vecMean(2:end,:)=[];
        [~,label_vector] = min(distance,[],2);
        for i=1:clusterNum
            index=(label_vector==i);
            vecMean(:,i)=sum(vector(index))/nnz(index);
        end
        if (vecMean==old_mean | itr>25)% You can change it accordingly
            break;
        end
    end
    labelIm = reshape(label_vector,size(gray));
else
    if isempty(varargin)
%% ============================== For RGB colorspace ==================
        img = double(img);
        
        red = img(:,:,1);
        green = img(:,:,2);
        blue = img(:,:,3);
        
        vecRed= (red(:)-min(red(:))+1);
        vecGreen= (green(:)-min(green(:))+1);
        vecBlue= (blue(:)-min(blue(:))+1);
        
        vecRed = repmat(vecRed,[1,clusterNum]);
        vecGreen = repmat(vecGreen,[1,clusterNum]);
        vecBlue = repmat(vecBlue,[1,clusterNum]);
        
        meanRed1=(1:clusterNum).*max((vecRed(:)))/(clusterNum+1);
        meanGreen1=(1:clusterNum).*max((vecGreen(:)))/(clusterNum+1);
        meanBlue1=(1:clusterNum).*max((vecBlue(:)))/(clusterNum+1);
        
        
        num = length(vecRed);
        itr = 0;
        while(true)
            
            itr = itr+1;
            oldRed=meanRed1;
            oldGreen=meanGreen1;
            oldBlue=meanBlue1;
            
            meanRed = repmat(meanRed1,[num,1]);
            meanGreen = repmat(meanGreen1,[num,1]);
            meanBlue = repmat(meanBlue1,[num,1]);
            
            distRed=(((vecRed-meanRed)).^2);
            distGreen=(((vecGreen-meanGreen)).^2);
            distBlue=(((vecBlue-meanBlue)).^2);
            
            clear meanRed meanGreen meanBlue
            
            distance = sqrt(distRed+distGreen+distBlue);
            
            [~,label_vector] = min(distance,[],2);
            
            for i=1:clusterNum
                index=(label_vector==i);
                meanRed1(:,i)=ceil(mean(vecRed(index)));
                meanGreen1(:,i)=ceil(mean(vecGreen(index)));
                meanBlue1(:,i)=ceil(mean(vecBlue(index)));
            end
            
            if ((meanRed1==oldRed & meanGreen1==oldGreen & meanBlue1==oldBlue | itr>25))% You can change it accordingly
                break;
            end
        end
        
        labelIm = reshape(label_vector,[m,n]);
        vecMean = [meanRed1;meanGreen1;meanBlue1];
        
    else
%% ================= For L*a*b* Color Space==============================
        cForm = makecform('srgb2lab');
        lab = double(applycform(img,cForm));
        a = lab(:,:,2);
        b = lab(:,:,3);
        
        veca= (a(:)-min(a(:))+1);
        vecb= (b(:)-min(b(:))+1);
        
        
        veca = repmat(veca,[1,clusterNum]);
        vecb = repmat(vecb,[1,clusterNum]);
        
        meana1=(1:clusterNum).*max((veca(:)))/(clusterNum+1);
        meanb1=(1:clusterNum).*max((vecb(:)))/(clusterNum+1);
         
        num = length(veca);
        itr = 0;
        while(true)
            
            itr = itr+1;
            olda=meana1;
            oldb=meanb1;
            
            meana = repmat(meana1,[num,1]);
            meanb= repmat(meanb1,[num,1]);
            
            dista=(((veca-meana)).^2);
            distb=(((vecb-meanb)).^2);
            
            clear meana meanb
            
            distance = sqrt(dista+distb);
            
            [~,label_vector] = min(distance,[],2);
            
            for i=1:clusterNum
                index=(label_vector==i);
             
                meana1(:,i)=ceil(mean(veca(index)));
                meanb1(:,i)=ceil(mean(vecb(index)));
                
            end
            
            if ((meana1==olda & meanb1==oldb))% You can change it accordingly
                break;
            end
        end
        
        labelIm = reshape(label_vector,[m,n]);
        vecMean = [meana1;meanb1];
    end
    
end
end

