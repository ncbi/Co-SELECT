#!/bin/env bash

dest="/home/pals2/Work/Co-SELECT-paper"
src=`readlink -f ../results`

declare -a normal_files=(
"fig_heatmap_shapemers_other1_cycle4_d1enriched.pdf"
"fig_heatmap_shapemers_publish_cycle4_d1enriched.pdf"
"fig_pca_shapemers_other1_cycle4_d1enriched.pdf"
"fig_pca_shapemers_publish_cycle4_d1enriched.pdf"
"fig_promiscuous_shapemers_other1_cycle4_d1enriched.pdf"
"fig_promiscuous_shapemers_publish_cycle4_d1enriched.pdf"
"fig_qvalue_combined_family_cycle4_d1enriched.pdf"
"fig_qvalue_selected_cycle4_d1enriched.pdf"
"fig_qvalue_separate_family_cycle4_d1enriched.pdf"
"table_significant_tfs_at_fdr_cycle4_d1enriched.pdf"
"detailed_results_cycle4_d1enriched.xlsx"
"enriched_shapemers_other1_cycle4_d1enriched.xlsx"
"enriched_shapemers_publish_cycle4_d1enriched.xlsx"
)


declare -a cycle3_files=(
"fig_heatmap_shapemers_other1_cycle3_d0.pdf"
"fig_heatmap_shapemers_publish_cycle3_d0.pdf"
"fig_pca_shapemers_other1_cycle3_d0.pdf"
"fig_pca_shapemers_publish_cycle3_d0.pdf"
"fig_promiscuous_shapemers_other1_cycle3_d0.pdf"
"fig_promiscuous_shapemers_publish_cycle3_d0.pdf"
"fig_qvalue_combined_family_cycle3_d0.pdf"
"fig_qvalue_selected_cycle3_d0.pdf"
"fig_qvalue_separate_family_cycle3_d0.pdf"
"table_significant_tfs_at_fdr_cycle3_d0.pdf"
"detailed_results_cycle3_d0.xlsx"
"enriched_shapemers_other1_cycle3_d0.xlsx"
"enriched_shapemers_publish_cycle3_d0.xlsx"
)


declare -a cycle4_files=(
"fig_heatmap_shapemers_other1_cycle4_d0.pdf"
"fig_heatmap_shapemers_publish_cycle4_d0.pdf"
"fig_pca_shapemers_other1_cycle4_d0.pdf"
"fig_pca_shapemers_publish_cycle4_d0.pdf"
"fig_promiscuous_shapemers_other1_cycle4_d0.pdf"
"fig_promiscuous_shapemers_publish_cycle4_d0.pdf"
"fig_qvalue_combined_family_cycle4_d0.pdf"
"fig_qvalue_selected_cycle4_d0.pdf"
"fig_qvalue_separate_family_cycle4_d0.pdf"
"table_significant_tfs_at_fdr_cycle4_d0.pdf"
"detailed_results_cycle4_d0.xlsx"
"enriched_shapemers_other1_cycle4_d0.xlsx"
"enriched_shapemers_publish_cycle4_d0.xlsx"
"fig_seqlogo_promiscuous_shapemers_publish_d0.pdf"
"fig_seqlogo_enriched_shapemers_publish_d0_th1.20.pdf"
#"fig_seqlogo_enriched_shapemers_publish_d0_th1.10.pdf"
)

cd $dest

for i in "${cycle4_files[@]}"; do
   echo "$i"
   cp "$src/d0/$i" "$dest/NewResults/NormalForeground"
   git add "NewResults/NormalForeground/$i"
done

for i in "${cycle3_files[@]}"; do
   echo "$i"
   cp "$src/d0/$i" "$dest/NewResults/NormalForegroundCycle3"
   git add "NewResults/NormalForegroundCycle3/$i"
done

for i in "${normal_files[@]}"; do
   echo "$i"
   cp "$src/d1enriched/$i" "$dest/NewResults/ForegroundWithEnrichedPartialMotifs"
   git add "NewResults/ForegroundWithEnrichedPartialMotifs/$i"
done



git commit -am "updated figures"
git pull
git push
