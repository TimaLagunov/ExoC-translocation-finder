import pandas as pd
import numpy as np
import argparse
import warnings
import os
warnings.filterwarnings("ignore")

def parser(command_line):
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-f", "--file", help="work file", type=str)
    parser.add_argument("-a", "--all", help="output all results", action='store_true')
    args = parser.parse_args(command_line)
    return args.file, args.all

def bedpe_line(file: str, all_results: bool):
    data_out = pd.DataFrame()
    data = pd.read_csv(file, sep=',')

    if not all_results:
        k = 1.73
        data['answer'] = ((np.exp(k*np.log10(data['counts_pat'])) < data['from_stat_binom']) | (data['from_stat_binom'] > 330))
        data = data[data['answer']]

    data['from_start'] = np.array(data['from_start']*data['resolution'], dtype=int)
    data['from_end'] = np.array((data['from_end']+1)*data['resolution'], dtype=int)
    data['to_start'] = np.array(data['to_start']*data['resolution'], dtype=int)
    data['to_end'] = np.array((data['to_end']+1)*data['resolution'], dtype=int)
    data['size_new'] = 1
    data['color'] = "0,0,0"

    data_out = data[['from_chr', 'from_start', 'from_end', 'to_chr', 'to_start', 'to_end', 'from_stat_binom', 'color', 'size']]
    #data_out.set_axis(["#chr1\tx1\tx2\tchr2\ty1\ty2\tscore\tcolor\tsize"], axis=1)
    data_out.to_csv(file[:-3]+'bedpe', sep="\t", header="#chr1\tx1\tx2\tchr2\ty1\ty2\tscore\tcolor\tsize".split("\t"), index=False)

def bedpe_base(file, all_results: bool):
    data_out = pd.DataFrame()
    data = pd.read_csv(file, sep=',')

    if not all_results:
        k_base = (2.4 - 0.74)/(3.62 - 1.81)
        b_base = 0.74 - k_base*1.81
        data['answer'] = (data['resolution'] >= 16e3) & (np.log10(data['prob']) >= k_base*np.log10(data['pat_sum'])+b_base) & (data['pat_sum'] > 10)
        data = data[data['answer']]

    data['x_start'] = 1
    data['x_end'] = 1
    data['y_start'] = 1
    data['y_end'] = 1
    data['size_new'] = 1
    data['color'] = "0,0,0"
    
    mask = (data['sx']==1) & (data['sy']==1)
    pp = data[mask][['x_start','x_end','y_start','y_end']].copy()
    pp[['x_start','x_end','y_start','y_end']] = np.vstack([(data['x'][mask]*data['resolution'][mask]).values, 
                                                         ((data['x'][mask]+data['x_lim'][mask])*data['resolution'][mask]).values,
                                                         (data['y'][mask]*data['resolution'][mask]).values, 
                                                         ((data['y'][mask]+data['y_lim'][mask])*data['resolution'][mask]).values,
                                                        ]).T
    
    mask = (data['sx']==-1) & (data['sy']==-1)
    mm = data[mask][['x_start','x_end','y_start','y_end']].copy()
    mm[['x_start','x_end','y_start','y_end']] = np.vstack([((data['x'][mask]+1-data['x_lim'][mask])*data['resolution'][mask]).values, 
                                                         ((data['x'][mask]+1)*data['resolution'][mask]).values,
                                                         ((data['y'][mask]+1-data['y_lim'][mask])*data['resolution'][mask]).values, 
                                                         ((data['y'][mask]+1)*data['resolution'][mask]).values,
                                                        ]).T

    mask = (data['sx']==1) & (data['sy']==-1)
    pm = data[mask][['x_start','x_end','y_start','y_end']].copy()
    pm[['x_start','x_end','y_start','y_end']] = np.vstack([(data['x'][mask]*data['resolution'][mask]).values, 
                                                         ((data['x'][mask]+data['x_lim'][mask])*data['resolution'][mask]).values,
                                                         ((data['y'][mask]+1-data['y_lim'][mask])*data['resolution'][mask]).values, 
                                                         ((data['y'][mask]+1)*data['resolution'][mask]).values,
                                                        ]).T

    mask = (data['sx']==-1) & (data['sy']==1)
    mp = data[mask][['x_start','x_end','y_start','y_end']].copy()
    mp[['x_start','x_end','y_start','y_end']] = np.vstack([((data['x'][mask]+1-data['x_lim'][mask])*data['resolution'][mask]).values, 
                                                         ((data['x'][mask]+1)*data['resolution'][mask]).values,
                                                         (data['y'][mask]*data['resolution'][mask]).values, 
                                                         ((data['y'][mask]+data['y_lim'][mask])*data['resolution'][mask]).values,
                                                        ]).T

    data[['x_start','x_end','y_start','y_end']] =  pp.append(mm).append(pm).append(mp)
    data['size_new'] = data['size']*data['resolution']

    data_out = data[['chrX', 'x_start', 'x_end', 'chrY', 'y_start', 'y_end', 'type', 'prob', 'sx', 'sy', 'color', 'size_new', 's_init']]
    #data_out.set_axis(["#chr1\tx1\tx2\tchr2\ty1\ty2\tname\tscore\tstrand1\tstrand2\tcolor\tsize\ts_init"], axis=1)
    data_out.to_csv(file[:-3]+'bedpe', sep="\t", header="#chr1\tx1\tx2\tchr2\ty1\ty2\tname\tscore\tstrand1\tstrand2\tcolor\tsize\ts_init".split("\t"), index=False)

def main(command_line=None):
    file, all_results = parser(command_line)
    if 'line' in file:
        bedpe_line(os.path.realpath(file), all_results)
    else:
        bedpe_base(os.path.realpath(file), all_results)

if __name__ == '__main__':
    main()