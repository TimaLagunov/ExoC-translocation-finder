# ExoC-translocation-founder
The code for translocation search on mcool maps after Exo-C experiment

## Data preparation
Make the mcool files for your samples with [cooler](https://github.com/open2c/cooler).
The mcool files with your samples and the mcool with control sample should have the same resolutions set.
For ExoC control sample usage, the resolution set is standard cooler set (powers of 2).

## Example usage
1. Create your working directory:

'''
mkdir example_exoc_dir
cd example_exoc_dir
'''

3. Link all scripts to the working directory:
> scripts_path = _path to ExoC-translocation-finder folder_
'''
for t in $(find ${scripts_path} -mindepth 1 -maxdepth 1 -name "*.py"); do 
	n=$(echo $t | rev | cut -d '/' -f 1 | rev) 
	ln -s -T $t $n
done
'''
3. Download the [sample](https://genedev.bionet.nsc.ru/ftp/by_Project/ExoC/mcools/hg19/s132_P10.mcool) and the [control](https://genedev.bionet.nsc.ru/ftp/by_Project/ExoC/mcools/merged_hg19/sup_pat.XX.mcool) files to your working directory from URLs, or via _wget_:
'''
wget https://genedev.bionet.nsc.ru/ftp/by_Project/ExoC/mcools/hg19/s132_P10.mcool
wget https://genedev.bionet.nsc.ru/ftp/by_Project/ExoC/mcools/merged_hg19/sup_pat.XX.mcool
'''
4. Inside your cluster queue run _trans_founder_complete_new.py_ script to make csv files with all artifacts. Then run _bedpe_maker_complete.py_ on the csv files to make bedpe files with filtered translocations.
'''
n_proc = _your number of valid threads (max threads minus 2 is highly recomended)_
python trans_founder_complete_new.py -d . -p s132_P10.mcool -c sup_pat.XX.mcool -n ${n_proc} --logs
python bedpe_maker_complete.py -f s132_P10-all_artifacts.csv
python bedpe_maker_complete.py -f s132_P10-all_artifacts_line.csv
'''
5. If you want to perform deeper analysis of artifacts, you can use jupyter notebook _Complete method analysis .ipynb_
