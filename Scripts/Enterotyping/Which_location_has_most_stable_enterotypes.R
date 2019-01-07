# which location is the most stable in enterotyping by different amplicons?
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Which_location_has_most_stable_enterotypes.R')
# get the result for separate enterotyping
enterotype.stability.diff.across.locat_ent.sep <- whichLocMostStabl(N=100000)
# Now, get result for all-together enterotyping
work_dir_1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_06_enterotyping/Sample-enterotype/Split_by_locat'
enterotype.stability.diff.across.locat_ent.togeth <- whichLocMostStabl(N=100000,work_dir=work_dir_1)
