clear all
close all

tor = [];
for i = 1:16
    rng(i,'twister');
    r = randperm(4,4);
    tor = [tor;r];
end
taskO = repmat({'nC-nR'}, size(tor));
taskO(tor==2) = {'nC-R'};
taskO(tor==3) = {'C-nR'};
taskO(tor==4) = {'C-R'};
taskO

writecell(taskO,'subjects/taskOrder_02232022.txt','Delimiter','tab')
type 'subjects/taskOrder_02232022.txt'
