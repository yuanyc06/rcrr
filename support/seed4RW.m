function [seeds, label] = seed4RW(sal, th1, th2)
% Function [seeds, label] = seed4RW(sal, th1, th2) extracts the seeds and
% corresponding labels from the former saliency approximation result.
% 
% Inputs:       sal - the k*1 saliency approximation result, where k is the
%                           total superpixel number
%                   th1 - higher threshold used to extract foreground seeds
%                   th2 - lower threshold used to extract background seeds
% Outputs:    seeds - a vector of the foreground/background seeds, each
%                   element correspond to the location of the seed
%                   label - the corresponding labels of the seeds, with 1
%                   as foreground, and 2 as background
% 
% 27/04/15 - by Yuchen Yuan
% Based on the paper:
% Yuchen Yuan, Changyang Li, Jinman Kim, Weidong Cai, and David Dagan Feng. 
% "Saliency Detection via Reversion Correction and Regularized Random Walks
% Ranking". TPAMI

bsalc_copy = sal;
sal(sal<th1)=0;
sal(sal>=th1)=1;

fgdseeds = find(sal);
fgdlabel = ones(size(fgdseeds));
bgdseeds = find(bsalc_copy < th2);
bgdlabel = 2*ones(size(bgdseeds));

seeds = [fgdseeds; bgdseeds];
label = [fgdlabel; bgdlabel];

seedslabel = sortrows([seeds,label]);
diffseedslabel = diff(seedslabel(:,1));
seedslabel(diffseedslabel == 0, :) = [];

seeds = seedslabel(:,1);
label = seedslabel(:,2);

end