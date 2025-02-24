# ExoC-translocation-finder
The code for translocation search on mcool maps after Exo-C experiment

# Requirements
1. Python >= 3.9
2. Install [cooltools](https://github.com/open2c/cooltools) in empty conda environment (this tool is really picky and may cause dependency probmlems).

## Example usage
1. Create your working directory:

```
mkdir example_exoc_dir
cd example_exoc_dir
```

2. Link all scripts to the working directory:
> [!NOTE]
> scripts_path = _path to ExoC-translocation-finder folder_

```
for t in $(find ${scripts_path} -mindepth 1 -maxdepth 1 -name "*.py"); do 
	n=$(echo $t | rev | cut -d '/' -f 1 | rev) 
	ln -s -T $t $n
done
```

3. Download the [sample](https://genedev.bionet.nsc.ru/ftp/by_Project/ExoC/hg38_files/s132_P10_exoc.mcool) and the [control](https://genedev.bionet.nsc.ru/ftp/by_Project/ExoC/hg38_files/sup_pat.XX.exoc.mcool) files to your working directory from URLs, or via _wget_:

```
wget https://genedev.bionet.nsc.ru/ftp/by_Project/ExoC/hg38_files/s132_P10_exoc.mcool
wget https://genedev.bionet.nsc.ru/ftp/by_Project/ExoC/hg38_files/sup_pat.XX.exoc.mcool
```

4. Inside your ***cluster queue*** run _trans_founder_complete_new.py_ script to make csv files with all artifacts. Then run _bedpe_maker_complete.py_ on the csv files to make bedpe files with filtered translocations.

>[!WARNING]
>**!Using cluster is highly recomended because of high RAM usage!** (~250gb for example)

>[!NOTE]
>n_proc = _your number of valid threads_ (max threads minus 2 is recomended)

```
python trans_founder_complete_new.py -d . -p s132_P10_exoc.mcool -c sup_pat.XX.mcool -n ${n_proc} --logs
python bedpe_maker_complete.py -f s132_P10_exoc-all_artifacts.csv
python bedpe_maker_complete.py -f s132_P10_exoc-all_artifacts_line.csv
```

5. Check that your result is the same with _s132_P10_exoc-all_artifacts.example.bedpe_ and _s132_P10_exoc-all_artifacts_line.example.bedpe_
6. If you want to perform deeper analysis of artifacts, you can use jupyter notebook _Complete method analysis .ipynb_

## Data and working directory preparation
1. Align your **fastq** files with [juicer](https://github.com/aidenlab/juicer)
2. From *merged_nodups.txt* file, make the mcool file for your ***sample*** with [cooler](https://github.com/open2c/cooler).
	The mcool file with your ***sample*** and the mcool with ***control*** sample should have the same resolutions set.
	For ExoC ***control*** sample usage, the resolution set is "1000,10000,16000,100000,250000,256000,1000000,4096000".
>[!NOTE]
>mnd_path = _path to merged_nodups.txt file_
>
>chr_size_file = _path to chrsize file_

```
awk 'BEGIN {OFS="\t"; FS=" "} ($9>=30)&&($12>=30) {print $2,$3,$6,$7}' ${mnd_path} | cooler cload pairs -c1 1 -p1 2 -c2 3 -p2 4 ${chr_size_file}:1000 - sample.cool
cooler zoomify -r 1000,10000,16000,100000,250000,256000,1000000,4096000 sample.cool
rm sample.cool
```

3. Follow steps **1** and **2** from the example.
4. Link your _sample.mcool_ file to working directory:

```
ln -s -T ../sample.mcool sample.mcool
```

5. Download the [XX control](https://genedev.bionet.nsc.ru/ftp/by_Project/ExoC/hg38_files/sup_pat.XX.exoc.mcool) or [XY control](https://genedev.bionet.nsc.ru/ftp/by_Project/ExoC/hg38_files/sup_pat.XY.exoc.mcool) (depends on your ***sample*** sex):

```
wget https://genedev.bionet.nsc.ru/ftp/by_Project/ExoC/hg38_files/sup_pat.XX.exoc.mcool
wget https://genedev.bionet.nsc.ru/ftp/by_Project/ExoC/hg38_files/sup_pat.XY.exoc.mcool
```

**OR** make your own control by merging the *merged_nodups.txt* files from your experiment and making *mcool* file like in step **2** from data preparation.
## Run translocation finder for your experiment
1. Inside your ***cluster queue*** run _trans_founder_complete_new.py_ script to make csv files with all artifacts. Then run _bedpe_maker_complete.py_ on the csv files to make bedpe files with filtered translocations.
>[!WARNING]
>**!Using cluster is highly recomended because of high RAM usage!** (~250gb for example)

>[!NOTE]
>n_proc = _your number of valid threads_ (max threads minus 2 is recomended)

```
python trans_founder_complete_new.py -d . -p sample.mcool -c control.mcool -n ${n_proc} --logs
python bedpe_maker_complete.py -f sample-all_artifacts.csv
python bedpe_maker_complete.py -f sample-all_artifacts_line.csv
```

If you need to test only specific chromosome pairs you can use *--chrpairs* flag to optimase calculations.
