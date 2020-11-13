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
      --p-jobs-to-start 10 \
      --o-representative-sequences tmp/req-seqs-Fullrun.qza \
      --o-table tmp/table-Fullrun.qza \
      --o-stats tmp/deblur-stats-Fullrun.qza

# view summary of deblur feature table
qiime feature-table summarize \
      --i-table tmp/table-Fullrun.qza \
      --o-visualization tmp/table-Fullrun.qzv

# view stats from deblur
qiime deblur visualize-stats \
      --i-deblur-stats tmp/deblur-stats-Fullrun.qza \
      --o-visualization tmp/deblur-stats.qzv

# De novo clustering to create 99% OTUs
qiime vsearch cluster-features-de-novo \
      --i-table tmp/table-Fullrun.qza \
      --i-sequences tmp/req-seqs-Fullrun.qza \
      --p-perc-identity 0.99 \
      --o-clustered-table tmp/table-dn-99-Fullrun.qza \
      --o-clustered-sequences tmp/rep-seqs-dn-99-Fullrun.qza

# alternative Closed-reference clustering with 99% identity against the Silva SSU RefNR99 database

qiime vsearch cluster-features-closed-reference \
      --i-table tmp/table-Fullrun.qza \
      --i-sequences tmp/req-seqs-Fullrun.qza \
      --i-reference-sequences /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/SILVA_138.1_SSURef_NR99_tax_silva.qza \
      --p-perc-identity 0.99 \
      --o-clustered-table tmp/table-cr-99-Fullrun.qza \
      --o-clustered-sequences tmp/rep-seqs-cr-99-Fullrun.qza \
      --o-unmatched-sequences tmp/unmatched-cr-Fullrun.qza

# alternative Open-reference clustering with 99% identity agains the same reference database above
qiime vsearch cluster-features-open-reference \
      --i-table tmp/table-Fullrun.qza \
      --i-sequences tmp/req-seqs-Fullrun.qza \
      --i-reference-sequences /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/SILVA_138.1_SSURed_NR99_tax_dna.qza \
      --p-perc-identity 0.99 \
      --o-clustered-table tmp/table-or-99-Fullrun.qza \
      --o-clustered-sequences tmp/rep-seqs-or-99-Fullrun.qza \
      --o-new-reference-sequences tmp/new-ref-seqs-or-Fullrun.qza

# De novo chimera checking with denovo clustering
qiime vsearch uchime-denovo \
      --i-table tmp/table-dn-99-Fullrun.qza \
      --i-sequences tmp/rep-seqs-dn-99-Fullrun.qza \
      --output-dir tmp/uchime-dn-out-Fullrun

# Visualize summary stats
qiime metadata tabulate \
      --m-input-file tmp/uchime-dn-out-Fullrun/stats.qza \
      --o-visualization tmp/uchime-dn-out-Fullrun/stats.qzv

# Exclude chimeras and borderline chimeras for denovo cluster
qiime feature-table filter-features \
      --i-table tmp/table-dn-99-Fullrun.qza \
      --m-metadata-file tmp/uchime-dn-out-Fullrun/nonchimeras.qza \
      --o-filtered-table tmp/uchime-dn-out-Fullrun/table-nonchimeric-wo-borderline.qza

qiime feature-table filter-seqs \
      --i-data tmp/rep-seqs-dn-99-Fullrun.qza \
      --m-metadata-file tmp/uchime-dn-out-Fullrun/nonchimeras.qza \
      --o-filtered-data tmp/uchime-dn-out-Fullrun/rep-seqs-nonchimeric-wo-borderline.qza

qiime feature-table summarize \
      --i-table tmp/uchime-dn-out-Fullrun/table-nonchimeric-wo-borderline.qza \
      --o-visualization tmp/uchime-dn-out-Fullrun/table-nonchimeric-wo-borderline.qzv

# alternative chimera checking
qiime vsearch uchime-ref \
      --i-sequences tmp/rep-seqs-or-99-Fullrun.qza \
      --i-table tmp/table-or-99-Fullrun.qza \
      --i-reference-sequences /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/SILVA_138.1_SSURef_NR99_tax_silva.qza \
      --output-dir tmp/uchime-ref-Fullrun
      
qiime metadata tabulate \
      --m-input-file tmp/uchime-ref-Fullrun/stats.qza \
      --o-visualization tmp/uchime-ref-Fullrun/stats.qzv

qiime feature-table filter-features \
      --i-table tmp/table-or-99-Fullrun.qza \
      --m-metadata-file tmp/uchime-ref-Fullrun/nonchimeras.qza \
      --o-filtered-table tmp/uchime-ref-Fullrun/table-nonchimeric-wo-borderline.qza

qiime feature-table filter-seqs \
      --i-data tmp/rep-seqs-or-99-Fullrun.qza \
      --m-metadata-file tmp/uchime-ref-Fullrun/nonchimeras.qza \
      --o-filtered-data tmp/uchime-ref-Fullrun/rep-seqs-nonchimeric-wo-borderline.qza

qiime feature-table summarize \
      --i-table tmp/uchime-ref-Fullrun/table-nonchimeric-wo-borderline.qza \
      --o-visualization tmp/uchime-ref-Fullrun/table-nonchimeric-wo-borderline.qzv

### Assign taxonomy de novo
qiime feature-classifier classify-sklearn \
      --i-classifier /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/classifier-SILVA-138.1.qza \
      --i-reads tmp//uchime-dn-out-Fullrun/rep-seqs-nonchimeric-wo-borderline.qza \
      --p-n-jobs -2 \
      --o-classification tmp/taxonomy-Fullrun.qza

qiime metadata tabulate \
      --m-input-file tmp/taxonomy-Fullrun.qza \
      --o-visualization tmp/taxonomy-Fullrun.qzv

#assign taxonomy open reference chimera
qiime feature-classifier classify-sklearn \
      --i-classifier /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/classifier-SILVA-138.1.qza \
      --i-reads tmp//uchime-ref-Fullrun/rep-seqs-nonchimeric-wo-borderline.qza \
      --p-n-jobs -2 \
      --o-classification tmp/taxonomy-OR-Fullrun.qza

qiime metadata tabulate \
      --m-input-file tmp/taxonomy-OR-Fullrun.qza \
      --o-visualization tmp/taxonomy-OR-Fullrun.qzv

# desperate attempt - taxonomy annotation of non clustered, non checked for chimera table

qiime feature-classifier classify-sklearn \
      --i-classifier /SAN/db/RDP_Silva/Silva_138.1/SSU_RefNR99/classifier-SILVA-138.1.qza \
      --i-reads tmp/req-seqs-Fullrun.qza \
      --p-n-jobs -2 \
      --o-classification tmp/taxonomy-notClutesterd-Fullrun.qza

qiime metadata tabulate \
      --m-input-file tmp/taxonomy-notClutesterd-Fullrun.qza \
      --o-visualization tmp/taxonomy-notClusterd-Fullrun.qzv

# did the same but filtered for 200bp, I will use that from now

qiime taxa filter-table \
      --i-table uchime-200-dn-out-Fullrun/table-nonchimeric-wo-borderline.qza \
      --i-taxonomy taxonomy-200-Fullrun.qza \
      --p-mode exact \
      --p-exclude "f__Eimeriorina" \
      --o-filtered-table table-200-Fullrun-Eimeriorina.qza


# export qiime to phyloseq
qiime tools export \
      --input-path table-no-mitochondria-no-chloroplast.qza \
      --output-dir phyloseq

biom convert \
     -i phyloseq/feature-table.biom \
     -o phyloseq/otu_table.txt \
     --to-tsv

# Manually change #OTUID to OTUID

# 2 Export taxonomy table
qiime tools export \
      --input-path taxonomy.qza \
      --output-dir phyloseq

# Manually change "feature ID" to "OTUID"
