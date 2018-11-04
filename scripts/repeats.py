import numpy as np
import scipy.stats
from scipy.integrate import nquad
from joblib import Parallel, delayed
import params

def f_norm(y,z, mu, sigma, n):
    return n*(n-1) * scipy.stats.norm(mu, sigma).pdf(z) * scipy.stats.norm(mu, sigma).pdf(y+z) * ((scipy.stats.norm(mu, sigma).cdf(y+z) - scipy.stats.norm(mu, sigma).cdf(z))**(n-2))

def find_variance(markers, w_size):
    # estimating variance
    samples = []
    for i in range(0, len(markers)-w_size+1):
        window = []
        for j in range(i, i+w_size):
            window.append(markers[j])
        for k in [1, 2]:
            flag, diff_list = check_window(window, k)
            if flag == True:
                samples.append(diff_list)
    N = 0 
    s = 0
    for sample in samples:
        N += len(sample)
        mu = np.mean(sample)
        square_array = np.square(np.asarray(sample) - mu)
        s += np.sum(square_array)
    return s, N
    

def get_var(markers):
    w_size = 10
    return find_variance(markers, w_size)
    

def check_window(window, k):
    diff_list = []
    for i in range(0, len(window)):
        if i-k < 0:
            continue
        diff = window[i] - window[i-k] 
        diff_list.append(diff)
    diff_max = max(diff_list)
    diff_min = min(diff_list)
    if diff_max - diff_min < 1500:
        return True, diff_list
    else:
        return False, None

def check_window_new(window, k, var):
    diff_list = []
    for i in range(0, len(window)):
        if i-k < 0:
            continue
        diff = window[i] - window[i-k] 
        diff_list.append(diff)
    diff_max = max(diff_list)
    diff_min = min(diff_list)
    
    mu = np.mean(diff_list)
    n = len(diff_list)

    test_stat =  diff_max - diff_min
    if test_stat >= 1500:
        return False
    else:
        if params.repeat_naive == True:
            return True
    integ_result = nquad(f_norm, [[test_stat/np.sqrt(var), np.inf], [-np.inf, np.inf]], args=(0, 1, n))
    p_value = integ_result[0]
    if p_value > 0.005:
        return True
    else:
        return False


def parallel_unit(i, markers, w_size, var):
    window = []
    for j in range(i, i+w_size):
        window.append(markers[j])
    flag = False
    for k in [1, 2]:
        flag = check_window_new(window, k, var)
        if flag == True:
            break
    return flag
    

def find_repeats(markers, w_size, var, num_threads):
    if params.repeat_naive == True:
        num_threads = 1 
    # find repeats
    ifrepeats = []
    for i in range(0, len(markers)):
        ifrepeats.append(False)

    flags = Parallel(n_jobs=num_threads)(delayed(parallel_unit)(i, markers, w_size, var) for i in range(0, len(markers)-w_size+1))
    
    for i in range(0, len(markers)-w_size+1):
        if flags[i] == True:
            for j in range(i, i+w_size):
                ifrepeats[j] = True
    repeats = []
    status = False
    for i in range(len(ifrepeats)):
        if ifrepeats[i] == True:
            if status == False:
                status = True
                begin = markers[i]
            end = markers[i]
        else:
            if status == True:
                status = False
                if begin != end:
                    interval = (begin, end)
                    repeats.append(interval)
    if status == True:
        if begin != end:
            interval = (begin, end)
            repeats.append(interval)
    return repeats   
 
def get_repeats(markers, var, num_threads):
    w_size = 10
    return find_repeats(markers, w_size, var, num_threads)

