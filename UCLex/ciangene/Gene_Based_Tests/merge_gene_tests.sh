ldak="/cluster/project8/vyp/cian/support/ldak/ldak"
for dir in */
do
    echo $dir
	$ldak --join-genes-reml $dir	
done

R CMD BATCH plot_results.R
