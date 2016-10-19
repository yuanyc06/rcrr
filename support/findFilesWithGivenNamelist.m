% Find and copy files in the given namelist from source folder to destination folder

source = 'E:\1.Work\Program\shc_rrwr_matlab\evaluation\rrwr_reverse\MSRA10K\';
destination = 'E:\1.Work\Program\shc_rrwr_matlab\evaluation\rrwr\MSRA10K\';
namelist = [27609
    46662
    52445
    62922
    77023
    83877
    92826
    94250
    100991
    105010
    108590
    129898
    142406
    153280
    153969
    154872
    164434
    165356
    173675
    202152
    202819];
filetype = '.mat';

% imnames=dir(['./precisionrecall/DUT-OMRON/Ours/' '*' 'mat']);

for i = 1:length(namelist)
    filename = strcat(source, num2str(namelist(i)), filetype);
%     filename = strcat(source, imnames(namelist(i)).name);
    copyfile(filename, destination);
end
