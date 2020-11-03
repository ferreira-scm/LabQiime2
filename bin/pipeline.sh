# mv *I_* ../path remove the 'I' files

#create artifact for qiime
qiime tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path /SAN/Susanas_den/EimeriaMicrobiome/data/2018_22_Eie_FullRun_1 \
      --input-format CasavaOneEightSingleLanePerSampleDirFmt \
      --output-path /SAN/Susanas_den/gitProj/LabQiime2/tmp/demux-paired-end_FullRun_1.qza

#join pair-ends
qiime vsearch join-pairs \
      --i-demultiplexed-seqs /SAN/Susanas_den/gitProj/LabQiime2/tmp/demux-paired-end_FullRun_1.qza \
      --o-joined-sequences /SAN/Susanas_den/gitProj/LabQiime2/tmp/demux-Fullrun-JOINED.qza

#make summary of reads
qiime demux summarize   \
      --i-data /SAN/Susanas_den/gitProj/LabQiime2/tmp/demux-Fullrun-JOINED.qza \
      --o-visualization /SAN/Susanas_den/gitProj/LabQiime2/tmp/demux-Fullrun-JOINED.qzv

#qualily filter
qiime quality-filter q-score \
      --i-demux  /SAN/Susanas_den/gitProj/LabQiime2/tmp/demux-Fullrun-JOINED.qza \
      --o-filtered-sequences  /SAN/Susanas_den/gitProj/LabQiime2/tmp/demux-Fullrun-Joined-FILTERED.qza \
      --o-filter-stats  /SAN/Susanas_den/gitProj/LabQiime2/tmp/demux-Fullrun-Joined-filter-STATS.qza

# visualize first quality filter

qiime metadata tabulate \
      --m-input-file tmp/demux-Fullrun-Joined-filter-STATS.qza \
      --o-visualization tmp/demux-Fullrun-Joined-filter-STATS.qzv

# denoise with deblur
qiime deblur denoise-other \
      --i-reference-seqs /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/SILVA_138.1_SSURef_NR99_tax_silva.qza \
      --i-demultiplexed-seqs tmp/demux-Fullrun-Joined-FILTERED.qza \
      --p-trim-length 420 \
      --p-sample-stats \
      --o-representative-sequences tmp/req-seqs-Fullrun.qza \
      --o-table tmp/table-Fullrun.qza \
      --o-stats tmp/deblur-stats-Fullrun.qza
