load('C:/Users/frank/Documents/GitHub/CAM_DA/data/conjunctions_leo.mat','data');
subset_Table = data(:, [3:8, 15:20]); 
writetable(subset_Table, 'xpxsdata.csv');