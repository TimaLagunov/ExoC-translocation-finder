import numpy as np
import pandas as pd
import scipy.sparse as sparse
from scipy.stats import binom
import cooltools
import os
from exoc_funcs_new import data_loader, interval, interval_maker, verbose_timedelta
from time import time
from multiprocessing import Pool
import random
import argparse
import warnings
warnings.filterwarnings("ignore")

print("All libraries are loaded")

def parse_line(command_line):
    global logs_dir, leave_logs, calc_cis
    # Command line parser block
    parser = argparse.ArgumentParser(description='')
    parser.add_argument("-d", "--dir", help="work direktory", type=str, required=True)
    parser.add_argument("-p", "--patient", help="patient file", type=str, required=True)
    parser.add_argument("-c", "--control", help="control file", type=str, required=True)
    parser.add_argument("-n", "--nproc", help="number of threads", type=int, required=True, default=1)
    parser.add_argument("-m", "--mosaik", help="list of mosaik percent", nargs='*', type=float, default=[1])
    parser.add_argument("-w", "--wsize", help="maximum window size", type=int, default=32)
    parser.add_argument("-r", "--resolutions", help="list of resolutions for translocations search", nargs='*', type=int, default=[int(4096e3), int(256e3), int(16e3), int(1e3)])
    parser.add_argument("--ttrans", help="list of thresholds for TRANS contacts to get them into analisys (one for each resolution)", nargs='*', type=float, default=[-4.0, -2.0, -1.0, 0.0])
    parser.add_argument("--tcis", help="list of thresholds for CIS contacts to get them into analisys (one for each resolution)", nargs='*', type=float, default=[-20.0, -10.0, -5.0, -2.5])
    parser.add_argument("--logs", help="leave all logs after complite", action='store_true')
    parser.add_argument("--cis", help="cmeasure artifact level for cis contacts (not recommended)", action='store_true')
    parser.add_argument("--chrpairs", help="cchromosome pairs for statistics calculation (format chr1-chr2)", nargs='*', type=str, default=[])
    args = parser.parse_args(command_line)

    # Log directory maker
    logs_dir = os.path.join(os.path.realpath(args.dir), 'tmp_logs')
    if not os.path.exists(logs_dir):
        os.mkdir(logs_dir)

    leave_logs = args.logs
    calc_cis = args.cis

    return (os.path.realpath(args.dir),
            os.path.realpath(args.patient), 
            os.path.realpath(args.control), 
            args.mosaik, 
            args.nproc, 
            args.wsize,
            args.resolutions,
            args.ttrans,
            args.tcis,
            args.chrpairs)

def calculations_trans(chr_name1, chr_name2):
    global t0, chrshape, pat_slice, control_slice, pat_cov, expected_prob_slice, same_coords_inds, control_same_dots_cov, ir, tmp_artifacts,\
        res, resolutions, tmp_old_dois, doi_threshold_trans, mosaik_array, prob_window, control_same_dots_slice, p_s, window_size, \
        pat_id, control_cov

    chr_time = time() - t0 # Time of the beginning of calculations for current cromosome 
    chr1, chr2 = interval(*chrshape[chr_name1],'left'), interval(*chrshape[chr_name2],'left')

    # Cutting current chromosome coords and data from sparce matrixes
    tmp_pat_slice = pat_slice[chr1.to_slice(), chr2.to_slice()].copy()
    tmp_control_slice = control_slice[chr1.to_slice(), chr2.to_slice()].copy()
    tmp_control_same_dots_slice = control_same_dots_slice[chr1.to_slice(), chr2.to_slice()].copy()
    tmp_expected_prob_slice = expected_prob_slice[chr1.to_slice(), chr2.to_slice()].copy()

    # Calculating "probabilities" for pixels to be artifacts
    prob = binom.pmf(tmp_pat_slice.data, pat_cov, tmp_expected_prob_slice.data.astype(float))
    prob2 = binom.pmf(tmp_control_same_dots_slice.data, control_same_dots_cov, tmp_expected_prob_slice.data.astype(float))
    prob[prob<=0] = prob[prob>0].min() if len(prob[prob>0])>0 else -1e-300
    prob = np.log10(prob)
    prob2[prob2<=0] = prob2[prob2>0].min() if len(prob2[prob2>0])>0 else -1e-300
    prob2 = np.log10(prob2)
    prob = np.nan_to_num(prob, nan=0.0, posinf=0.0, neginf=0.0)
    prob2 = np.nan_to_num(prob2, nan=0.0, posinf=0.0, neginf=0.0)
    del(tmp_control_same_dots_slice, tmp_expected_prob_slice)

    r,c = tmp_pat_slice.tocoo().row, tmp_pat_slice.tocoo().col

    # Preparing masks for pixels that already parts of artifacts
    if ir == 0:
        artifacts_doi = False
        old_dois_doi = True
    else:
        artifacts_chr_mask = (tmp_artifacts['chrX']==chr_name1) & (tmp_artifacts['chrY']==chr_name2)
        artifacts_doi = (tmp_artifacts['x_interval'][artifacts_chr_mask].apply(lambda x: x(r)).sum(axis=0)) & \
                        (tmp_artifacts['y_interval'][artifacts_chr_mask].apply(lambda x: x(c)).sum(axis=0))
        old_dois_doi = True

    if res == resolutions[-1]:
        # At last resolution we try only dots, that were parts of DOIs at previous step
        old_dois_chr_mask = (tmp_old_dois['chrX']==chr_name1) & (tmp_old_dois['chrY']==chr_name2)
        old_dois_doi =  (tmp_old_dois['x_interval'][old_dois_chr_mask].apply(lambda x: x(r)).sum(axis=0)) & \
                        (tmp_old_dois['y_interval'][old_dois_chr_mask].apply(lambda x: x(c)).sum(axis=0))
        doi_mask = np.array(((~artifacts_doi) & (old_dois_doi)), dtype=bool)
    else:
        # Making a mask of Dots-Of-Interest (DOIs)
        doi_mask = np.array(((prob < prob2 + doi_threshold_trans) & (~artifacts_doi) & (old_dois_doi)), dtype=bool)
    del(prob, prob2, artifacts_doi, old_dois_doi)

    best_prob = np.ones((doi_mask.sum(), 1), dtype = float)*1e4
    best_params = np.zeros((doi_mask.sum(), 10), dtype = int)
    best_mosaik = np.zeros((doi_mask.sum(), 1), dtype = float)

    for ibp, (x,y) in enumerate(zip(r[doi_mask], c[doi_mask])):
        for sx in [1,-1]:
            for sy in [1,-1]:
                for mos in mosaik_array:
                    m_k = mos/2
                    roi_window = interval_maker(prob_window.shape[2],sx,sy) + (x,y)
                    roi_pat = tmp_pat_slice[roi_window[0].to_slice(), roi_window[1].to_slice()].toarray()
                    roi_control = tmp_control_slice[roi_window[0].to_slice(), roi_window[1].to_slice()].toarray()
                    
                    trans_prob = 1/2 - np.sqrt(1/4 + 
                                                   1/4*np.square(roi_control/control_cov) - 
                                                   1/4*roi_control/control_cov
                                                  )

                    tmp_answer_prob = binom.pmf(roi_pat, pat_cov, 
                                                (roi_control/control_cov).astype(float)
                                               )
                    tmp_answer_prob[tmp_answer_prob<=0] = tmp_answer_prob[tmp_answer_prob>0].min() if len(tmp_answer_prob[tmp_answer_prob>0])>0 else -1e-300
                    tmp_answer_prob = np.log10(tmp_answer_prob)
                    tmp_answer_prob = np.nan_to_num(tmp_answer_prob, nan=0.0, posinf=0.0, neginf=0.0)

                    answer_arm = np.zeros((roi_pat.shape[0], roi_pat.shape[1], prob_window.shape[2]))
                    for i,ep in enumerate(p_s[:prob_window.shape[2]]):
                        answer_arm[:,:,i] = binom.pmf(roi_pat, pat_cov, (1/2 - np.sqrt(1/4 + 
                        4*(m_k/2*np.square(ep) + (1-m_k/2)*np.square(trans_prob) 
                             - m_k/2*ep - (1-m_k/2)*(trans_prob)))).astype(float))
                    answer_arm[answer_arm<=0] = answer_arm[answer_arm>0].min() if len(answer_arm[answer_arm>0])>0 else -1e-300
                    answer_arm = np.log10(answer_arm)
                    answer_arm = np.nan_to_num(answer_arm, nan=0.0, posinf=0.0, neginf=0.0)

                    if ir > 0:
                        for art in tmp_artifacts[artifacts_chr_mask].loc[:,['x_interval','y_interval']].values:
                            tmp_x_slice = art[0](roi_window[0])
                            tmp_y_slice = art[1](roi_window[1])
                            
                            if tmp_x_slice is not None and tmp_y_slice is not None: 
                                x_slice = (tmp_x_slice-(x+(window_size-1)*(sx-1)/2)).to_slice()
                                y_slice = (tmp_y_slice-(y+(window_size-1)*(sy-1)/2)).to_slice()
                                tmp_answer_prob[x_slice,y_slice] = 1
                                answer_arm[x_slice,y_slice] = 1
                    
                    tmp_answer_prob = tmp_answer_prob[::sx,::sy]
                    answer_arm = answer_arm[::sx,::sy, :]
                    roi_pat = roi_pat[::sx,::sy]
                    roi_control = roi_control[::sx,::sy]

                    for s_init in range(5):
                        tmp_answer_arm = answer_arm[:,:,s_init:]
                        tmp_prob_window = prob_window[:tmp_answer_arm.shape[0], :tmp_answer_arm.shape[1], :tmp_answer_arm.shape[2]]

                        tmp_r_slice = min(tmp_answer_arm.shape[0],tmp_answer_prob.shape[0],tmp_prob_window.shape[0])
                        tmp_c_slice = min(tmp_answer_arm.shape[1],tmp_answer_prob.shape[1],tmp_prob_window.shape[1])
                        tmp_s_slice = min(tmp_answer_arm.shape[2],tmp_prob_window.shape[2])

                        tmp_prob_window = tmp_prob_window[:tmp_r_slice, :tmp_c_slice, :tmp_s_slice]
                        tmp_answer_arm = tmp_answer_arm[:tmp_r_slice, :tmp_c_slice, :tmp_s_slice] * tmp_prob_window
                        tmp_answer_prob = tmp_answer_prob[:tmp_r_slice, :tmp_c_slice]
                        roi_pat = roi_pat[:tmp_r_slice, :tmp_c_slice]
                        roi_control = roi_control[:tmp_r_slice, :tmp_c_slice]

                        for x_lim in range(1,tmp_answer_arm.shape[0]+1):
                            for y_lim in range(1,tmp_answer_arm.shape[1]+1):
                                s_lim = max(x_lim,y_lim)
                                tmp_prob = (tmp_answer_prob[:x_lim, :y_lim] * tmp_prob_window[:x_lim, :y_lim, :s_lim].sum(axis=2)).sum() - \
                                            (tmp_answer_arm[:x_lim, :y_lim, :s_lim]).sum()
                                        
                                pat_sum = (roi_pat[:x_lim, :y_lim] * tmp_prob_window[:x_lim, :y_lim, :s_lim].sum(axis=2)).sum()
                                control_sum = (roi_control[:x_lim, :y_lim] * tmp_prob_window[:x_lim, :y_lim, :s_lim].sum(axis=2)).sum()

                                if best_prob[ibp,0] >= tmp_prob:
                                    best_params[ibp,:] = [x_lim,y_lim,s_init,sx,sy,0,pat_sum,control_sum,x,y]
                                    best_prob[ibp] = tmp_prob
                                    best_mosaik[ibp] = mos

    print(f'{chr_name1}-{chr_name2}: {best_params.shape[0]} done!')
        
    end_time = time() - t0
    with open(os.path.join(logs_dir, f'{pat_id}.{chr_name1}-{chr_name2}.tmp_logs'), 'w') as file:
        file.write(f'\n{chr_name1}-{chr_name2} ({best_params.shape[0]}): Finish at {verbose_timedelta(end_time)}, delta: {verbose_timedelta(end_time-chr_time)}\n')

    tmp_data_frame = pd.DataFrame(columns=['x','y','chrX','chrY','prob','x_lim','y_lim','s_init','sx','sy'])
    tmp_data_frame['x'] = best_params[:,-2].reshape(-1)
    tmp_data_frame['y'] = best_params[:,-1].reshape(-1)
    tmp_data_frame['chrX'] = chr_name1
    tmp_data_frame['chrY'] = chr_name2
    tmp_data_frame['prob'] = -best_prob[:,0].reshape(-1) # Artifact level was negative (smaller value means greater artifact level). Change it to normal
    tmp_data_frame['x_lim'] = best_params[:,0].reshape(-1)
    tmp_data_frame['y_lim'] = best_params[:,1].reshape(-1)
    tmp_data_frame['s_init'] = best_params[:,2].reshape(-1)
    tmp_data_frame['sx'] = best_params[:,3].reshape(-1)
    tmp_data_frame['sy'] = best_params[:,4].reshape(-1)
    tmp_data_frame['type'] = best_params[:,5].reshape(-1)
    tmp_data_frame['type'] = tmp_data_frame['type'].apply(lambda x: {0: 'cis', 1: 'trans'}[x])
    tmp_data_frame['pat_sum'] = best_params[:,6].reshape(-1)
    tmp_data_frame['control_sum'] = best_params[:,7].reshape(-1)
    tmp_data_frame['resolution'] = res
    tmp_data_frame['mosaik'] = best_mosaik

    return tmp_data_frame
                    
def calculations_cis(chr_name1, chr_name2):
    global t0, chrshape, pat_slice, control_slice, pat_cov, expected_prob_slice, same_coords_inds, control_same_dots_cov, ir, tmp_artifacts,\
        res, resolutions, tmp_old_dois, doi_threshold_trans, mosaik_array, prob_window, control_same_dots_slice, p_s, window_size, \
        pat_id, expected_trans, control_cov

    chr_time = time() - t0 # Time of the beginning of calculations for current cromosome 
    chr1, chr2 = interval(*chrshape[chr_name1],'left'), interval(*chrshape[chr_name2],'left')

    # Cutting current chromosome coords and data from sparce matrixes
    tmp_pat_slice = pat_slice[chr1.to_slice(), chr2.to_slice()].copy()
    tmp_control_slice = control_slice[chr1.to_slice(), chr2.to_slice()].copy()
    tmp_control_same_dots_slice = control_same_dots_slice[chr1.to_slice(), chr2.to_slice()].copy()
    tmp_expected_prob_slice = expected_prob_slice[chr1.to_slice(), chr2.to_slice()].copy()

    # Calculating "probabilities" for pixels to be artifacts
    prob = binom.pmf(tmp_pat_slice.data, pat_cov, tmp_expected_prob_slice.data.astype(float))
    prob2 = binom.pmf(tmp_control_same_dots_slice.data, control_same_dots_cov, tmp_expected_prob_slice.data.astype(float))
    prob[prob<=0] = prob[prob>0].min() if len(prob[prob>0])>0 else -1e-300
    prob = np.log10(prob)
    prob2[prob2<=0] = prob2[prob2>0].min() if len(prob2[prob2>0])>0 else -1e-300
    prob2 = np.log10(prob2)
    prob = np.nan_to_num(prob, nan=0.0, posinf=0.0, neginf=0.0)
    prob2 = np.nan_to_num(prob2, nan=0.0, posinf=0.0, neginf=0.0)
    del(tmp_control_same_dots_slice, tmp_expected_prob_slice)

    r,c = tmp_pat_slice.tocoo().row, tmp_pat_slice.tocoo().col
    
    if res == resolutions[0]:
        diag_mask = (c-r) > 4
    else:
        diag_mask = ((c-r) > 4) & ((c-r) < 4+window_size)

    # Preparing masks for pixels that already parts of artifacts
    if ir == 0:
        artifacts_doi = False
        old_dois_doi = True
    else:
        artifacts_chr_mask = (tmp_artifacts['chrX']==chr_name1) & (tmp_artifacts['chrY']==chr_name2)
        artifacts_doi = (tmp_artifacts['x_interval'][artifacts_chr_mask].apply(lambda x: x(r)).sum(axis=0)) & \
                        (tmp_artifacts['y_interval'][artifacts_chr_mask].apply(lambda x: x(c)).sum(axis=0))
        old_dois_doi = True

    if res == resolutions[-1]:
        # At last resolution we try only dots, that were parts of DOIs at previous step
        old_dois_chr_mask = (tmp_old_dois['chrX']==chr_name1) & (tmp_old_dois['chrY']==chr_name2)
        old_dois_doi =  (tmp_old_dois['x_interval'][old_dois_chr_mask].apply(lambda x: x(r)).sum(axis=0)) & \
                        (tmp_old_dois['y_interval'][old_dois_chr_mask].apply(lambda x: x(c)).sum(axis=0))
        doi_mask = np.array(((~artifacts_doi) & (old_dois_doi)), dtype=bool)
    else:
        # Making a mask of Dots-Of-Interests (DOIs)
        doi_mask = np.array(((prob < prob2 + doi_threshold_cis) & (~artifacts_doi) & (old_dois_doi)), dtype=bool)
    del(prob, prob2, artifacts_doi, old_dois_doi)
    doi_mask = doi_mask & diag_mask

    best_prob = np.ones((doi_mask.sum(), 1), dtype = float)*1e4
    best_params = np.zeros((doi_mask.sum(), 10), dtype = int)
    best_mosaik = np.zeros((doi_mask.sum(), 1), dtype = float)

    for ibp, (x,y) in enumerate(zip(r[doi_mask], c[doi_mask])):
        for sx in [1,-1]:
            for sy in [1,-1]:
                for mos in mosaik_array:
                    m_k = mos/2
                    roi_window = interval_maker(prob_window.shape[2],sx,sy) + (x,y)
                    roi_pat = tmp_pat_slice[roi_window[0].to_slice(), roi_window[1].to_slice()].toarray()
                    roi_control = tmp_control_slice[roi_window[0].to_slice(), roi_window[1].to_slice()].toarray()

                    ref_prob = 1/2 - np.sqrt(1/4 + 
                                                   1/2*np.square(roi_control/control_cov) - 
                                                   1/2*roi_control/control_cov + 
                                                   expected_trans - 
                                                   np.square(expected_trans)
                                                  )
                    tmp_answer_prob = binom.pmf(roi_pat, pat_cov, 
                                                (1/2 - np.sqrt(1/4 + 2*(np.square(ref_prob) - (ref_prob) + 
                                                                                  np.square(expected_trans) - (expected_trans)))).astype(float)
                                               )
                    tmp_answer_prob[tmp_answer_prob<=0] = tmp_answer_prob[tmp_answer_prob>0].min() if len(tmp_answer_prob[tmp_answer_prob>0])>0 else -1e-300
                    tmp_answer_prob = np.log10(tmp_answer_prob)
                    tmp_answer_prob = np.nan_to_num(tmp_answer_prob, nan=0.0, posinf=0.0, neginf=0.0)

                    answer_arm = np.zeros((roi_pat.shape[0], roi_pat.shape[1], prob_window.shape[2]))
                    for i,ep in enumerate(p_s[:prob_window.shape[2]]):
                        answer_arm[:,:,i] = binom.pmf(roi_pat, pat_cov, (1/2 - np.sqrt(1/4 + 
                                                                                    (1*m_k*np.square(ep) - 1*m_k*ep + 
                                                                                    1*m_k*np.square(ref_prob) - 1*m_k*ref_prob + 
                                                                                    (4-2*m_k)*np.square(expected_trans) - (4-2*m_k)*(expected_trans))
                                                                                    )).astype(float))
                    answer_arm[answer_arm<=0] = answer_arm[answer_arm>0].min() if len(answer_arm[answer_arm>0])>0 else -1e-300
                    answer_arm = np.log10(answer_arm)
                    answer_arm = np.nan_to_num(answer_arm, nan=0.0, posinf=0.0, neginf=0.0)

                    answer_trans = binom.pmf(roi_pat, pat_cov, (1/2 - np.sqrt(1/4 + 
                                                                                    ((4-3*m_k)*np.square(expected_trans) - (4-3*m_k)*(expected_trans) + 
                                                                                     1*m_k*np.square(ref_prob) - 1*m_k*ref_prob)
                                                                                    )).astype(float))
                    answer_trans[answer_trans<=0] = answer_trans[answer_trans>0].min() if len(answer_trans[answer_trans>0])>0 else -1e-300
                    answer_trans = np.log10(answer_trans)
                    answer_trans = np.nan_to_num(answer_trans, nan=0.0, posinf=0.0, neginf=0.0)

                    if ir > 0:
                        for art in tmp_artifacts[artifacts_chr_mask].loc[:,['x_interval','y_interval']].values:
                            tmp_x_slice = art[0](roi_window[0])
                            tmp_y_slice = art[1](roi_window[1])
                            
                            if tmp_x_slice is not None and tmp_y_slice is not None: 
                                x_slice = (tmp_x_slice-(x+(window_size-1)*(sx-1)/2)).to_slice()
                                y_slice = (tmp_y_slice-(y+(window_size-1)*(sy-1)/2)).to_slice()
                                tmp_answer_prob[x_slice,y_slice] = 1
                                answer_arm[x_slice,y_slice] = 1
                                answer_trans[x_slice,y_slice] = 1
                    
                    tmp_answer_prob = tmp_answer_prob[::sx,::sy]
                    answer_arm = answer_arm[::sx,::sy, :]
                    answer_trans = answer_trans[::sx,::sy]
                    roi_pat = roi_pat[::sx,::sy]
                    roi_control = roi_control[::sx,::sy]

                    for s_init in range(5):
                        tmp_answer_arm = answer_arm[:,:,s_init:]
                        tmp_prob_window = prob_window[:tmp_answer_arm.shape[0], :tmp_answer_arm.shape[1], :tmp_answer_arm.shape[2]]

                        tmp_r_slice = min(tmp_answer_arm.shape[0],tmp_answer_prob.shape[0],tmp_prob_window.shape[0])
                        tmp_c_slice = min(tmp_answer_arm.shape[1],tmp_answer_prob.shape[1],tmp_prob_window.shape[1])
                        tmp_s_slice = min(tmp_answer_arm.shape[2],tmp_prob_window.shape[2])

                        tmp_prob_window = tmp_prob_window[:tmp_r_slice, :tmp_c_slice, :tmp_s_slice]
                        tmp_answer_arm = tmp_answer_arm[:tmp_r_slice, :tmp_c_slice, :tmp_s_slice] * tmp_prob_window
                        tmp_answer_prob = tmp_answer_prob[:tmp_r_slice, :tmp_c_slice]
                        answer_trans = answer_trans[:tmp_r_slice, :tmp_c_slice]
                        roi_pat = roi_pat[:tmp_r_slice, :tmp_c_slice]
                        roi_control = roi_control[:tmp_r_slice, :tmp_c_slice]

                        for x_lim in range(1,tmp_answer_arm.shape[0]+1):
                            for y_lim in range(1,tmp_answer_arm.shape[1]+1):
                                s_lim = max(x_lim,y_lim)
                                tmp_prob = (tmp_answer_prob[:x_lim, :y_lim] * tmp_prob_window[:x_lim, :y_lim, :s_lim].sum(axis=2)).sum() - \
                                            (tmp_answer_arm[:x_lim, :y_lim, :s_lim]).sum()
                                        
                                tmp_prob_trans =(tmp_answer_prob[:x_lim, :y_lim] * tmp_prob_window[:x_lim, :y_lim, :s_lim].sum(axis=2)).sum() - \
                                                (answer_trans[:x_lim, :y_lim] * tmp_prob_window[:x_lim, :y_lim, :s_lim].sum(axis=2)).sum()
                                        
                                pat_sum = (roi_pat[:x_lim, :y_lim] * tmp_prob_window[:x_lim, :y_lim, :s_lim].sum(axis=2)).sum()
                                control_sum = (roi_control[:x_lim, :y_lim] * tmp_prob_window[:x_lim, :y_lim, :s_lim].sum(axis=2)).sum()

                                if best_prob[ibp,0] >= tmp_prob:
                                    best_params[ibp,:] = [x_lim,y_lim,s_init,sx,sy,0,pat_sum,control_sum,x,y]
                                    best_prob[ibp] = tmp_prob
                                    best_mosaik[ibp] = mos

                                if best_prob[ibp,0] >= tmp_prob_trans:
                                    best_params[ibp,:] = [x_lim,y_lim,s_init,sx,sy,1,pat_sum,control_sum,x,y]
                                    best_prob[ibp] = tmp_prob_trans
                                    best_mosaik[ibp] = mos

    print(f'{chr_name1}-{chr_name2}: {best_params.shape[0]} done!'+' '*100, end='\r')
        
    end_time = time() - t0
    with open(os.path.join(logs_dir, f'{pat_id}.{chr_name1}-{chr_name2}.tmp_logs'), 'w') as file:
        file.write(f'\n{chr_name1}-{chr_name2} ({best_params.shape[0]}): Finish at {verbose_timedelta(end_time)}, delta: {verbose_timedelta(end_time-chr_time)}\n')

    tmp_data_frame = pd.DataFrame(columns=['x','y','chrX','chrY','prob','x_lim','y_lim','s_init','sx','sy'])
    tmp_data_frame['x'] = best_params[:,-2].reshape(-1)
    tmp_data_frame['y'] = best_params[:,-1].reshape(-1)
    tmp_data_frame['chrX'] = chr_name1
    tmp_data_frame['chrY'] = chr_name2
    tmp_data_frame['prob'] = -best_prob[:,0].reshape(-1)
    tmp_data_frame['x_lim'] = best_params[:,0].reshape(-1)
    tmp_data_frame['y_lim'] = best_params[:,1].reshape(-1)
    tmp_data_frame['s_init'] = best_params[:,2].reshape(-1)
    tmp_data_frame['sx'] = best_params[:,3].reshape(-1)
    tmp_data_frame['sy'] = best_params[:,4].reshape(-1)
    tmp_data_frame['type'] = best_params[:,5].reshape(-1)
    tmp_data_frame['type'] = tmp_data_frame['type'].apply(lambda x: {0: 'cis', 1: 'trans'}[x])
    tmp_data_frame['pat_sum'] = best_params[:,6].reshape(-1)
    tmp_data_frame['control_sum'] = best_params[:,7].reshape(-1)
    tmp_data_frame['resolution'] = res
    tmp_data_frame['mosaik'] = best_mosaik

    return tmp_data_frame

def calculations_trans_line(chr_name1, chr_name2):
    global t0, chrshape

    chr_time = time() - t0 # Time of the beginning of calculations for current cromosome 
    chr1, chr2 = interval(*chrshape[chr_name1],'left'), interval(*chrshape[chr_name2],'left')

    tmp_pat_slice = pat_slice[chr1.to_slice(), chr2.to_slice()].copy()
    tmp_control_slice = control_slice[chr1.to_slice(), chr2.to_slice()].copy()

    tmp_df_ans = pd.DataFrame()
    for ax in [0,1]:
        tmp_df = pd.DataFrame(columns=['chr_from', 'chr_to', 'x', 'stat'])
        
        tmp_df['counts_pat'] = list(np.array(tmp_pat_slice.sum(axis=ax)).reshape(-1))
        tmp_df['counts_control'] = list(np.array(tmp_control_slice.sum(axis=ax)).reshape(-1))
        tmp_df['chr_from'] = chr_name1 if ax==1 else chr_name2
        tmp_df['chr_to'] = chr_name2 if ax==1 else chr_name1
        tmp_df['x'] = list(range(chr1.maxint()+1-chr1.minint())) if ax==1 else list(range(chr2.maxint()+1-chr2.minint()))

        tmp_answers = list(
            -np.nan_to_num(np.log10(binom.pmf(
                np.array(tmp_pat_slice.sum(axis=ax)).reshape(-1), 
                tmp_pat_slice.sum(), 
                np.array(tmp_control_slice.sum(axis=ax)).reshape(-1)/tmp_control_slice.sum()
            )), nan=0.0, posinf=350.0, neginf=-350.0) \
            + np.nan_to_num(np.log10(binom.pmf(
                np.array(tmp_control_slice.sum(axis=ax)).reshape(-1), 
                tmp_control_slice.sum(), 
                np.array(tmp_control_slice.sum(axis=ax)).reshape(-1)/tmp_control_slice.sum()
            )), nan=0.0, posinf=350.0, neginf=-350.0)
        )
        tmp_df['stat_binom'] = tmp_answers

        tmp_df = tmp_df[(tmp_df['stat_binom']>6) & (tmp_df['counts_control']>0)]

        if len(tmp_df) >= 1:
            
            if len(tmp_df) > 1:
                x_minmax_pairs = []
                x_min_tmp = np.sort(tmp_df['x'].values)[0]
                x_max_tmp = x_min_tmp
                x_prev = x_min_tmp
                for x_tmp in np.sort(tmp_df['x'].values)[1:]:
                    if x_tmp - x_prev <= 4:
                        x_max_tmp = x_tmp
                        x_prev = x_tmp
                    else:
                        x_minmax_pairs += [(x_min_tmp, x_max_tmp)]
                        x_min_tmp = x_tmp
                        x_max_tmp = x_min_tmp
                        x_prev = x_min_tmp
                x_minmax_pairs += [(x_min_tmp, x_max_tmp)]
            else:
                x_minmax_pairs = [(tmp_df['x'].values[0],tmp_df['x'].values[0])]

            for x_pair in x_minmax_pairs:
                tmp_ans = pd.DataFrame()
                roi_slice = interval(x_pair[0],x_pair[1]+1,'left')
                
                tmp_ans.loc[0,'counts_pat'] = np.array(tmp_pat_slice.sum(axis=ax), dtype=int).reshape(-1)[roi_slice.to_slice()].sum()
                tmp_ans['counts_control'] = np.array(tmp_control_slice.sum(axis=ax), dtype=int).reshape(-1)[roi_slice.to_slice()].sum()
                tmp_ans['from_chr'] = chr_name1 if ax==1 else chr_name2
                tmp_ans['to_chr'] = chr_name2 if ax==1 else chr_name1
                tmp_ans['from_start'] = roi_slice.minint()
                tmp_ans['from_end'] = roi_slice.maxint()
                tmp_ans['size'] = (roi_slice.maxint()+1-roi_slice.minint())
                tmp_ans['x_interval'] = roi_slice*res
                tmp_ans['from_stat_binom'] = -np.nan_to_num(np.log10(binom.pmf(
                        np.array(tmp_pat_slice.sum(axis=ax)).reshape(-1)[roi_slice.to_slice()].sum(), 
                        tmp_pat_slice.sum(), 
                        np.array(tmp_control_slice.sum(axis=ax)).reshape(-1)[roi_slice.to_slice()].sum()/tmp_control_slice.sum()
                    )), nan=0.0, posinf=350.0, neginf=-350.0) + np.nan_to_num(np.log10(binom.pmf(
                        np.array(tmp_control_slice.sum(axis=ax)).reshape(-1)[roi_slice.to_slice()].sum(), 
                        tmp_control_slice.sum(), 
                        np.array(tmp_control_slice.sum(axis=ax)).reshape(-1)[roi_slice.to_slice()].sum()/tmp_control_slice.sum()
                    )), nan=0.0, posinf=350.0, neginf=-350.0)
                tmp_slice_y_binom = -np.nan_to_num(np.log10(binom.pmf(
                                        np.array(tmp_pat_slice[roi_slice.to_slice(),:].sum(axis=0)).reshape(-1), 
                                        tmp_pat_slice.sum(), 
                                        np.array(tmp_control_slice[roi_slice.to_slice(),:].sum(axis=0)).reshape(-1)/tmp_control_slice.sum()
                                    )), nan=0.0, posinf=350.0, neginf=-350.0) \
                                    + np.nan_to_num(np.log10(binom.pmf(
                                        np.array(tmp_control_slice[roi_slice.to_slice(),:].sum(axis=0)).reshape(-1), 
                                        tmp_control_slice.sum(), 
                                        np.array(tmp_control_slice[roi_slice.to_slice(),:].sum(axis=0)).reshape(-1)/tmp_control_slice.sum()
                                    )), nan=0.0, posinf=350.0, neginf=-350.0) \
                                    if ax==1 else \
                                    -np.nan_to_num(np.log10(binom.pmf(
                                        np.array(tmp_pat_slice[:,roi_slice.to_slice()].sum(axis=1)).reshape(-1), 
                                        tmp_pat_slice.sum(), 
                                        np.array(tmp_control_slice[:,roi_slice.to_slice()].sum(axis=1)).reshape(-1)/tmp_control_slice.sum()
                                    )), nan=0.0, posinf=350.0, neginf=-350.0) \
                                    + np.nan_to_num(np.log10(binom.pmf(
                                        np.array(tmp_control_slice[:,roi_slice.to_slice()].sum(axis=1)).reshape(-1), 
                                        tmp_control_slice.sum(), 
                                        np.array(tmp_control_slice[:,roi_slice.to_slice()].sum(axis=1)).reshape(-1)/tmp_control_slice.sum()
                                    )), nan=0.0, posinf=350.0, neginf=-350.0)

                tmp_y_max = np.arange(len(tmp_slice_y_binom))[tmp_slice_y_binom==tmp_slice_y_binom.max()]
                if len(tmp_y_max) > 1: tmp_y_max = tmp_y_max.mean()
                tmp_ans['to_start'] = int(max(tmp_y_max - 10*(roi_slice.maxint()+1-roi_slice.minint()), 0))
                tmp_ans['to_end'] = int(min(tmp_y_max + 10*(roi_slice.maxint()+1-roi_slice.minint()), chr1.maxint())) \
                    if ax==0 else \
                    int(min(tmp_y_max + 10*(roi_slice.maxint()+1-roi_slice.minint()), chr2.maxint()))
                tmp_ans['y_interval'] = interval(tmp_ans['to_start'][0]*res, (tmp_ans['to_end'][0]+1)*res, 'left')
                tmp_ans['to_stat_binom'] = tmp_slice_y_binom.max()
                tmp_df_ans = tmp_df_ans.append(tmp_ans, ignore_index=True)

    tmp_df_ans['resolution'] = res
    
    end_time = time() - t0
    with open(os.path.join(logs_dir, f'{pat_id}.{chr_name1}-{chr_name2}.tmp_logs'), 'w') as file:
        file.write(f'\n{chr_name1}-{chr_name2} ({tmp_df_ans.shape[0]}): Finish at {verbose_timedelta(end_time)}, delta: {verbose_timedelta(end_time-chr_time)}\n')

    return tmp_df_ans

def main(command_line=None):
    global log_file, t0, chrshape, pat_slice, control_slice, pat_cov, expected_prob_slice, same_coords_inds, control_same_dots_cov, ir, tmp_artifacts,\
        res, resolutions, tmp_old_dois, doi_threshold_trans, doi_threshold_cis, mosaik_array, prob_window, control_same_dots_slice, p_s, window_size,\
        pat_id, expected_trans, control_cov

    # Save start time for logs
    t0 = time()
    end_time = t0

    # Parse line in dif variables
    (workdir,
     pat_file, 
     control_file, 
     mosaik_array, 
     nproc,
     max_w_size,
     resolutions,
     doi_thrs_trans,
     doi_thrs_cis,
     chrpairs) = parse_line(command_line)

    # DataFrames initialization
    cols_list = ['x', 'y', 'chrX', 'chrY', 'prob', 'x_lim', 'y_lim', 's_init', 'sx', 'sy', 'type', 'pat_sum', 'control_sum', 'resolution', 'mosaik', 'size', 'x_interval', 'y_interval']
    artifacts = pd.DataFrame(columns=cols_list)
    old_dois = pd.DataFrame(columns=cols_list)
    artifacts_line = pd.DataFrame()

    # Few constants that will help us
    chunksize = 10**7
    max_coord_const = int(10**7)
    ir = 0 # Valuable iteration counter

    pat_id = os.path.splitext(os.path.split(pat_file)[-1])[0]
    log_file = os.path.join(workdir, f'./{pat_id}.logs')
    all_start_time = time() - t0
    with open(log_file, 'w') as file:
        file.write(f'Start whole process at {verbose_timedelta(all_start_time)}\n')

    for (res, doi_threshold_trans, doi_threshold_cis) in zip(resolutions, doi_thrs_trans, doi_thrs_cis):
        print("Loading patient and control maps")
        pat, control, chrshape = data_loader(pat_file, control_file, res)

        print("All maps are loaded")

        start_time = time() - t0
        with open(log_file, 'a') as file:
            file.write(f'\nStart work for resolution {res} at {verbose_timedelta(start_time)}\n')

        # Loading sparse matrixes without main diogonal and sum all counts
        pat_slice = sparse.triu(pat.matrix(balance=False, sparse=True)[:,:], k=1)
        control_slice = sparse.triu(control.matrix(balance=False, sparse=True)[:,:], k=1)
        pat_cov = pat_slice.data.sum()

        # Cut only nonzero pixels of patient from control
        _,_, same_coords_inds = np.intersect1d(
            pat_slice.row.astype(int)*max_coord_const + pat_slice.col.astype(int),
            control_slice.row.astype(int)*max_coord_const + control_slice.col.astype(int),
            assume_unique=True, return_indices=True)

        # If not all nonzero pixels of patient are nonzero in control, we add patient map to control map and cut the pixels again
        if len(same_coords_inds) != len(pat_slice.data):
            control_slice = (control_slice+pat_slice).tocoo()
            _,_, same_coords_inds = np.intersect1d(
            pat_slice.row.astype(int)*max_coord_const + pat_slice.col.astype(int),
            control_slice.row.astype(int)*max_coord_const + control_slice.col.astype(int),
            assume_unique=True, return_indices=True)

        control_same_dots_cov = control_slice.data[same_coords_inds].sum()
        control_cov = control_slice.data.sum()
        expected_prob = control_slice.data[same_coords_inds]/control_same_dots_cov # Simple contact probability estimation
        expected_prob_slice = sparse.csr_matrix((expected_prob, (control_slice.row[same_coords_inds], control_slice.col[same_coords_inds])), 
                                                shape = control_slice.shape)
        control_same_dots_slice = sparse.csr_matrix((control_slice.data[same_coords_inds], (control_slice.row[same_coords_inds], control_slice.col[same_coords_inds])), 
                                                shape = control_slice.shape)
        del(expected_prob)

        pat_slice = pat_slice.tocsr()
        control_slice = control_slice.tocsr()

        # Calc contacts vs distance probability for fearther calculations
        cvd_control = cooltools.expected_cis(
            clr=control,
            view_df=None,
            smooth=False,
            aggregate_smoothed=False,
            nproc=nproc,
            chunksize=chunksize,
            clr_weight_name = None,
            ignore_diags=1,
        )
        p_s = ((cvd_control.groupby('dist').sum()['count.sum']/cvd_control.groupby('dist').sum()['n_valid'])/control_same_dots_cov).values[1:]
        del(cvd_control)

        window_size = min(max_w_size,len(p_s))
        r = np.arange(window_size)[:, np.newaxis, np.newaxis]
        c = np.arange(window_size)[np.newaxis, :, np.newaxis]
        s = np.arange(window_size)[np.newaxis, np.newaxis, :]
        prob_window = (r + c) == s # Make mask for fearther calculations
        del(r,c,s)

        # Calc expected trans contact probability for fearther calculations
        cvd_trans_control = cooltools.expected_trans(
            clr=control,
            view_df=None,
            nproc=nproc,
            chunksize=chunksize,
            clr_weight_name = None,
        )
        cvd_trans_control = cvd_trans_control[cvd_trans_control['region2']!='chrM']
        expected_trans = 1/2 - np.sqrt(1/4 + 
                               1/4*np.square((cvd_trans_control['count.sum'].sum()/cvd_trans_control['n_valid'].sum())/control_same_dots_cov) - 
                               1/4*(cvd_trans_control['count.sum'].sum()/cvd_trans_control['n_valid'].sum())/control_same_dots_cov
                              )
        del(cvd_trans_control)

        # Recalculating contacts vs distance probability for diploid chroms
        p_s = 1/2 - np.sqrt(1/4 + 
                           1/2*np.square(p_s) - 
                           1/2*p_s + 
                            expected_trans - 
                            np.square(expected_trans)
                          )
        
        if ir > 0:
            tmp_artifacts = artifacts.copy()
            if len(tmp_artifacts.values) > 0:
                tmp_artifacts['x_interval'] = tmp_artifacts[['x_interval','resolution']].apply(lambda x: interval(x[0].minint()*x[1]/res, (x[0].maxint()+1)*x[1]/res-1, 'both'), axis=1)
                tmp_artifacts['y_interval'] = tmp_artifacts[['y_interval','resolution']].apply(lambda x: interval(x[0].minint()*x[1]/res, (x[0].maxint()+1)*x[1]/res-1, 'both'), axis=1)
            else:
                ir = 0
            
            tmp_old_dois = old_dois.copy()
            if len(tmp_old_dois.values) > 0:
                tmp_old_dois['x_interval'] = tmp_old_dois[['x_interval','resolution']].apply(lambda x: interval(x[0].minint()*x[1]/res, (x[0].maxint()+1)*x[1]/res-1, 'both'), axis=1)
                tmp_old_dois['y_interval'] = tmp_old_dois[['y_interval','resolution']].apply(lambda x: interval(x[0].minint()*x[1]/res, (x[0].maxint()+1)*x[1]/res-1, 'both'), axis=1)

        # Calcs for trans contacts
        if len(chrpairs) == 0:
            all_chrs_combs = []
            for ichr1,chr_name1 in enumerate(chrshape):
                for ichr2,chr_name2 in enumerate(chrshape):
                    if ichr1 < ichr2:
                        all_chrs_combs += [(chr_name1, chr_name2)]
            random.shuffle(all_chrs_combs)
        else:
            all_chrs_combs = [tuple(x.split('-')) for x in chrpairs]
        
        print("Start parallel calculations")
        
        if nproc == 1:
            best_params_data_frame = calculations_trans(*all_chrs_combs[0]).copy()
            for chr_comb in all_chrs_combs[1:]:
                best_params_data_frame = best_params_data_frame.append(calculations_trans(*chr_comb), ignore_index=True).copy()
        else:
            pool = Pool(processes=nproc)
            all_stats = pool.starmap(calculations_trans, all_chrs_combs)
            pool.close()
            best_params_data_frame = all_stats[0].copy()
            for df in all_stats[1:]:
                best_params_data_frame = best_params_data_frame.append(df, ignore_index=True).copy()
            del(all_stats)
        
        if calc_cis:
            #calculuus for cis contacts
            all_chrs_combs = []
            for chr_name in list(chrshape.keys()):
                all_chrs_combs += [(chr_name, chr_name)]
            random.shuffle(all_chrs_combs)

            if nproc == 1:
                for chr_comb in all_chrs_combs[1:]:
                    best_params_data_frame = best_params_data_frame.append(calculations_cis(*chr_comb), ignore_index=True).copy()
            else: 
                pool = Pool(processes=nproc)
                all_stats = pool.starmap(calculations_cis, all_chrs_combs)
                pool.close()
                for df in all_stats[:]:
                    best_params_data_frame = best_params_data_frame.append(df, ignore_index=True).copy()
                del(all_stats)

        all_logs = ''
        for (chr_name1, chr_name2) in all_chrs_combs:
            try:
                with open(os.path.join(logs_dir, f'{pat_id}.{chr_name1}-{chr_name2}.tmp_logs'), 'r') as file:
                    all_logs += file.read()
            except:
                pass
            else:
                os.remove(os.path.join(logs_dir, f'{pat_id}.{chr_name1}-{chr_name2}.tmp_logs'))
        with open(log_file, 'a') as file:
            file.write(all_logs)
        del(all_logs)

        #Calculating with another method
        if res < 500e3:
            if len(chrpairs) == 0:
                all_chrs_combs = []
                for ichr1,chr_name1 in enumerate(chrshape):
                    for ichr2,chr_name2 in enumerate(chrshape):
                        if ichr1 < ichr2:
                            all_chrs_combs += [(chr_name1, chr_name2)]
                random.shuffle(all_chrs_combs)
            else:
                all_chrs_combs = [tuple(x.split('-')) for x in chrpairs]
            
            if nproc == 1:
                answers_line = calculations_trans_line(*all_chrs_combs[0]).copy()
                for chr_comb in all_chrs_combs[1:]:
                    answers_line = answers_line.append(calculations_trans_line(*chr_comb), ignore_index=True).copy()
            else:
                pool = Pool(processes=nproc)
                all_stats = pool.starmap(calculations_trans_line, all_chrs_combs)
                pool.close()
                answers_line = all_stats[0].copy()
                for df in all_stats[1:]:
                    answers_line = answers_line.append(df, ignore_index=True).copy()
                del(all_stats)
            
            all_logs = ''
            for (chr_name1, chr_name2) in all_chrs_combs:
                try:
                    with open(os.path.join(logs_dir, f'{pat_id}.{chr_name1}-{chr_name2}_line.tmp_logs'), 'r') as file:
                        all_logs += file.read()
                except:
                    pass
                else:
                    os.remove(os.path.join(logs_dir, f'{pat_id}.{chr_name1}-{chr_name2}_line.tmp_logs'))
            with open(log_file, 'a') as file:
                file.write(all_logs)
            del(all_logs)

        
        # Some magic with outputs
        if len(best_params_data_frame.values) > 0:
            best_params_data_frame.sort_values('prob', ascending=False, inplace=True)
            
            if ir > 0: artifacts_old = artifacts.copy()
            
            best_params_data_frame['size'] = best_params_data_frame[['x_lim','y_lim']].min(axis=1)
            best_params_data_frame['x_interval'] = best_params_data_frame[['x','x_lim','sx']].apply(lambda x: 
                                                                                                          interval(x[0], x[0] + x[1], closed = 'left') 
                                                                                                          if x[2]>0 else 
                                                                                                          interval(x[0] - x[1], x[0], closed = 'right')
                                                                                                          , axis=1)
            best_params_data_frame['y_interval'] = best_params_data_frame[['y','y_lim','sy']].apply(lambda x: 
                                                                                                          interval(x[0], x[0] + x[1], closed = 'left') 
                                                                                                          if x[2]>0 else 
                                                                                                          interval(x[0] - x[1], x[0], closed = 'right')
                                                                                                          , axis=1)
            artifacts = best_params_data_frame.copy()
            i = 0
            while i < best_params_data_frame.shape[0]:
                if i >= artifacts.shape[0]: break
                    
                tmp_dot = artifacts.iloc[i,:]
                tmp_interval = {'x': (tmp_dot['x'], tmp_dot['x'] + tmp_dot['x_lim']) if tmp_dot['sx']>0 else (tmp_dot['x'] - tmp_dot['x_lim'], tmp_dot['x']),
                                'y': (tmp_dot['y'], tmp_dot['y'] + tmp_dot['y_lim']) if tmp_dot['sy']>0 else (tmp_dot['y'] - tmp_dot['y_lim'], tmp_dot['y']),}
                
                tmp_mask =~((artifacts['chrX'] == tmp_dot['chrX']) *
                            (artifacts['chrY'] == tmp_dot['chrY']) *
                            (artifacts['x'] >= tmp_interval['x'][0]) *
                            (artifacts['x'] <= tmp_interval['x'][1]) *
                            (artifacts['y'] >= tmp_interval['y'][0]) *
                            (artifacts['y'] <= tmp_interval['y'][1]) *
                            (artifacts['prob'] < tmp_dot['prob']))
                self_cheker =((artifacts['chrX'] == tmp_dot['chrX']) *
                              (artifacts['chrY'] == tmp_dot['chrY']) *
                              (artifacts['x'] >= tmp_interval['x'][0]) *
                              (artifacts['x'] <= tmp_interval['x'][1]) *
                              (artifacts['y'] >= tmp_interval['y'][0]) *
                              (artifacts['y'] <= tmp_interval['y'][1]) *
                              (artifacts['prob'] > tmp_dot['prob'])).sum() > 0
                
                if self_cheker > 0:
                    i-=1
                    tmp_mask[tmp_dot.name] = False
                    
                artifacts = artifacts[tmp_mask]
                i+=1
                
            tmp_mask = (((artifacts['size'] >= 3) & (artifacts['prob'] > 0) & (artifacts['type']=='cis')) |
                        ((artifacts['size'] >= 3) & (artifacts['prob'] > 0) & (artifacts['type']=='trans') & (artifacts['chrX']==artifacts['chrY'])))
            if res == resolutions[-1]:
                tmp_mask = (((artifacts[['x_lim','y_lim']].max(axis=1) >= 3) & (artifacts['prob'] > 0) & (artifacts['type']=='cis')) |
                        ((artifacts[['x_lim','y_lim']].max(axis=1) >= 3) & (artifacts['prob'] > 0) & (artifacts['type']=='trans') & (artifacts['chrX']==artifacts['chrY'])))
            artifacts = artifacts[tmp_mask]

            tmp_mask = (((best_params_data_frame['size'] < 3) & (best_params_data_frame['prob'] > 0) & (best_params_data_frame['type']=='cis')) |
                        ((best_params_data_frame['size'] < 3) & (best_params_data_frame['prob'] > 0) & (best_params_data_frame['type']=='trans') & (best_params_data_frame['chrX']==best_params_data_frame['chrY'])))
            old_dois = best_params_data_frame[tmp_mask].copy()
            
            if ir > 0: artifacts = artifacts_old.append(artifacts, ignore_index=True)
        
        if res < 500e3:
            if len(answers_line.values) > 0:
                tmp_mask = (answers_line['size'] > 3) & (answers_line['from_stat_binom'] > 0)
                artifacts_line = artifacts_line.append(answers_line[tmp_mask], ignore_index=True)

        end_time = time() - t0
        with open(log_file, 'a') as file:
            file.write(f'\nFinish at {verbose_timedelta(end_time)}, delta: {verbose_timedelta(end_time-start_time)}\n')
        
        ir += 1
        
    artifacts.to_csv(f'./{pat_id}-all_artifacts.csv', index=False)
    artifacts_line.to_csv(f'./{pat_id}-all_artifacts_line.csv', index=False)
    all_end_time = time() - t0
    with open(log_file, 'a') as file:
        file.write(f'\nFinish whole process at {verbose_timedelta(all_end_time)}, sum delta: {verbose_timedelta(all_end_time-all_start_time)}\n')


if __name__ == '__main__':
    main()
    if not leave_logs:
        os.removedirs(logs_dir)
        os.remove(os.path.realpath(log_file))