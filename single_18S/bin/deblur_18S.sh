#import data to qiime2
qiime tools import \
    --type 'SampleData[PairedEndSequencesWithQuality]' \
    --input-path /SAN/Susanas_den/EimeriaMicrobiome/data/18S_2018_22_Mb_1 \
    --input-format CasavaOneEightSingleLanePerSampleDirFmt \
    --output-path single_18S/tmp/demux-18S_201822Mb1.qza

# join pair-ends

qiime vsearch join-pairs \
--i-demultiplexed-seqs tmp/demux-18S_201822Mb1.qza \
--o-joined-sequences tmp/demux-18S_201822Mb1-JOINED.qza

# quality filter
qiime quality-filter q-score \
    --i-demux  tmp/demux-18S_201822Mb1-JOINED.qza \
    --o-filtered-sequences  tmp/demux-18S_201822Mb1-filtered.qza \
    --o-filter-stats  tmp/demux-18S_201822Mb1-filtered-STATS.qza

# denoise with deblur
qiime deblur denoise-other \
    --i-reference-seqs /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/SILVA_138.1_SSURef_NR99_tax_silva.qza \
    --i-demultiplexed-seqs tmp/demux-18S_201822Mb1-filtered.qza \
    --p-trim-length 200 \
    --p-sample-stats \
    --p-jobs-to-start 10 \
    --o-representative-sequences tmp/req-seqs-18S.qza \
    --o-table tmp/table-18S.qza \
    --o-stats tmp/deblur-stats-18S.qza

# De novo clustering to create 99% OTUs
qiime vsearch cluster-features-de-novo \
    --i-table tmp/table-18S.qza \
    --i-sequences tmp/req-seqs-18S.qza \
    --p-perc-identity 0.99 \
    --o-clustered-table tmp/table-dn-99-18S.qza \
    --o-clustered-sequences tmp/rep-seqs-dn-99-18S.qza

# De novo chimera checking with denovo clustering
qiime vsearch uchime-denovo \
    --i-table tmp/table-dn-99-18S.qza \
    --i-sequences tmp/rep-seqs-dn-99-18S.qza \
    --output-dir tmp/uchime-dn-out-18S

# Exclude chimeras and borderline chimeras for denovo cluster
qiime feature-table filter-features \
--i-table tmp/table-dn-99-18S.qza \
--m-metadata-file tmp/uchime-dn-out-18S/nonchimeras.qza \
--o-filtered-table tmp/uchime-dn-out-18S/table-nonchimeric-wo-borderline.qza

qiime feature-table filter-seqs \
--i-data tmp/rep-seqs-dn-99-18S.qza \
--m-metadata-file tmp/uchime-dn-out-18S/nonchimeras.qza \
--o-filtered-data tmp/uchime-dn-out-18S/rep-seqs-nonchimeric-wo-borderline.qza

qiime feature-table summarize \
--i-table tmp/uchime-dn-out-18S/table-nonchimeric-wo-borderline.qza \
--o-visualization tmp/uchime-dn-out-18S/table-nonchimeric-wo-borderline.qzv


#assign taxonomy
qiime feature-classifier classify-sklearn \
--i-classifier /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/classifier-SILVA-138.1.qza \
--i-reads tmp//uchime-dn-out-18S/rep-seqs-nonchimeric-wo-borderline.qza \
--p-n-jobs -2 \
--o-classification tmp/taxonomy-18S.qza

qiime metadata tabulate \
--m-input-file tmp/taxonomy-18S.qza \
--o-visualization tmp/taxonomy-18S.qzv


#export qiime to phyloseq
qiime tools export \
--input-path tmp/uchime-dn-out-18S/table-nonchimeric-wo-borderline.qza \
--output-path phyloseq

biom convert \
-i phyloseq/feature-table.biom \
-o phyloseq/otu_table.txt \
--to-tsv

qiime tools export \
--input-path tmp/taxonomy-18S.qza \
--output-path phyloseq/

# manually change #OTUID to OTUID in otu_table.txt
#manually vhange "feature ID" to "OTUID"
