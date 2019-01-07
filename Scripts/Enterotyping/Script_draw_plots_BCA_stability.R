
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Graphs_enterotype_stability_v5.R')
# First, without removing outliers
# V1V2 - 3 locations
plotEntStabil(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V1V2_3e_IL_sample-enterotype.csv',
               file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V3V4_3e_IL_sample-enterotype.csv', 
               file3 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V5V6_3e_IL_sample-enterotype.csv',
               qiime_out = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
               outdir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs/BCA',
               location='.IL')
plotEntStabil(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V1V2_3e_TR_sample-enterotype.csv',
               file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V3V4_3e_TR_sample-enterotype.csv', 
               file3 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V5V6_3e_TR_sample-enterotype.csv',
               qiime_out = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
               outdir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs/BCA',
               locat='.TR')
plotEntStabil(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V1V2_3e_RE_sample-enterotype.csv',
               file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V3V4_3e_RE_sample-enterotype.csv', 
               file3 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V5V6_3e_RE_sample-enterotype.csv',
               qiime_out = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
               outdir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs/BCA',
               locat='.RE')
# V3V4 - 3 locations
plotEntStabil(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V1V2_3e_IL_sample-enterotype.csv',
              file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V3V4_3e_IL_sample-enterotype.csv', 
              file3 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V5V6_3e_IL_sample-enterotype.csv',
              qiime_out = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
              outdir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs/BCA',
              location='.IL')
plotEntStabil(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V1V2_3e_TR_sample-enterotype.csv',
              file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V3V4_3e_TR_sample-enterotype.csv', 
              file3 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V5V6_3e_TR_sample-enterotype.csv',
              qiime_out = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
              outdir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs/BCA',
              locat='.TR')
plotEntStabil(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V1V2_3e_RE_sample-enterotype.csv',
              file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V3V4_3e_RE_sample-enterotype.csv', 
              file3 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V5V6_3e_RE_sample-enterotype.csv',
              qiime_out = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
              outdir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs/BCA',
              locat='.RE')
# V5V6 - 3 locations
plotEntStabil(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V1V2_3e_IL_sample-enterotype.csv',
              file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V3V4_3e_IL_sample-enterotype.csv', 
              file3 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V5V6_3e_IL_sample-enterotype.csv',
              qiime_out = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
              outdir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs/BCA',
              location='.IL')
plotEntStabil(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V1V2_3e_TR_sample-enterotype.csv',
              file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V3V4_3e_TR_sample-enterotype.csv', 
              file3 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V5V6_3e_TR_sample-enterotype.csv',
              qiime_out = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
              outdir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs/BCA',
              locat='.TR')
plotEntStabil(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V1V2_3e_RE_sample-enterotype.csv',
              file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V3V4_3e_RE_sample-enterotype.csv', 
              file3 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Sample-enterotype/V5V6_3e_RE_sample-enterotype.csv',
              qiime_out = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
              outdir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_09_26_sep_ent_by_location/Graphs/BCA',
              locat='.RE')


# Than, with removing outliers
source('/home/lev-genetik/Desktop/Projects/liege/src/Lev/Graphs_enterotype_stability_v5.R')
# V1V2 - 3 locations
plotEntStabil(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V1V2_3e_IL_sample-enterotype.csv',
              file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V3V4_3e_IL_sample-enterotype.csv', 
              file3 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V5V6_3e_IL_sample-enterotype.csv',
              qiime_out = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
              outdir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Graphs/BCA',
              location='.IL',
              outliers.removed =TRUE)
plotEntStabil(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V1V2_3e_TR_sample-enterotype.csv',
              file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V3V4_3e_TR_sample-enterotype.csv', 
              file3 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V5V6_3e_TR_sample-enterotype.csv',
              qiime_out = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
              outdir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Graphs/BCA',
              location='.TR',
              outliers.removed =TRUE)
plotEntStabil(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V1V2_3e_RE_sample-enterotype.csv',
              file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V3V4_3e_RE_sample-enterotype.csv', 
              file3 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V5V6_3e_RE_sample-enterotype.csv',
              qiime_out = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V1V2/otu_table_L6.txt',
              outdir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Graphs/BCA',
              location='.RE',
              outliers.removed =TRUE)
# V3V4 - 3 locations
plotEntStabil(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V1V2_3e_IL_sample-enterotype.csv',
              file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V3V4_3e_IL_sample-enterotype.csv', 
              file3 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V5V6_3e_IL_sample-enterotype.csv',
              qiime_out = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
              outdir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Graphs/BCA',
              location='.IL',
              outliers.removed =TRUE)
plotEntStabil(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V1V2_3e_TR_sample-enterotype.csv',
              file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V3V4_3e_TR_sample-enterotype.csv', 
              file3 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V5V6_3e_TR_sample-enterotype.csv',
              qiime_out = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
              outdir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Graphs/BCA',
              location='.TR',
              outliers.removed =TRUE)
plotEntStabil(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V1V2_3e_RE_sample-enterotype.csv',
              file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V3V4_3e_RE_sample-enterotype.csv', 
              file3 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V5V6_3e_RE_sample-enterotype.csv',
              qiime_out = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V3V4/otu_table_L6.txt',
              outdir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Graphs/BCA',
              location='.RE',
              outliers.removed =TRUE)
# V5V6 - 3 locations
plotEntStabil(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V1V2_3e_IL_sample-enterotype.csv',
              file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V3V4_3e_IL_sample-enterotype.csv', 
              file3 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V5V6_3e_IL_sample-enterotype.csv',
              qiime_out = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
              outdir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Graphs/BCA',
              location='.IL',
              outliers.removed =TRUE)
plotEntStabil(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V1V2_3e_TR_sample-enterotype.csv',
              file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V3V4_3e_TR_sample-enterotype.csv', 
              file3 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V5V6_3e_TR_sample-enterotype.csv',
              qiime_out = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
              outdir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Graphs/BCA',
              location='.TR',
              outliers.removed =TRUE)
plotEntStabil(file1 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V1V2_3e_RE_sample-enterotype.csv',
              file2 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V3V4_3e_RE_sample-enterotype.csv', 
              file3 = '/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Sample-enterotype/V5V6_3e_RE_sample-enterotype.csv',
              qiime_out = '/media/lev-genetik/980E73270E72FD96/Liege/Work_dir/DATA/summarize_taxa_counts/V5V6/otu_table_L6.txt',
              outdir='/media/lev-genetik/980E73270E72FD96/Liege/Results/Enterotyping/2018_11_09_enterotyping/Graphs/BCA',
              location='.RE',
              outliers.removed =TRUE)