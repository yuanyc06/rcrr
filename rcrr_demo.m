% Demo for publication "Reversion Correction and Regularized 
% Random Walks Ranking for Saliency Detection" 
% by Yuchen Yuan
% The BMIT Group, The University of Sydney 2015

function rcrr_demo()

%% Initialization
addpath('support');
imDir = 'image';
salDir = 'result';% Output path of the saliency map
if ~exist(salDir, 'dir')
    mkdir(salDir);
end
imFiles = dir(fullfile(imDir, '*.jpg'));
imFiles = {imFiles.name};
imNum = length(imFiles);

%% Algorithm start
for i = 1:imNum
    fprintf('%s: processing image %d of %d\n', mfilename, i, imNum);
% Calculate saliency
    sal = rcrr_saliency(imDir, imFiles{i});
% Output saliency map to file
    salName = fullfile(salDir, [imFiles{i}(1:end-4), '_rcrr.png']);
    imwrite(sal, salName);
end

