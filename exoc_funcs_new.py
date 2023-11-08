import numpy as np
import pandas as pd
import scipy.sparse as sparse
from scipy.special import gamma, gammaln
import cooler


def data_loader(pat_path: str, control_path: str, resolution: int) -> tuple:
    pat = cooler.Cooler(f'{pat_path}::resolutions/{resolution}')
    control = cooler.Cooler(f'{control_path}::resolutions/{resolution}')
    pat_chrom_shape = {}
    control_chrom_shape = {}
    
    for i in pat.chromnames:
        pat_chrom_shape[i] = pat.extent(i)

    for i in control.chromnames:
        control_chrom_shape[i] = control.extent(i)
    
    assert pat_chrom_shape == control_chrom_shape, 'Chromosome names, shapes and order have to be the same'

    return pat, control, pat_chrom_shape


class interval:
    def __init__(self,start,end,closed):
        self.start = start if start <= end else end
        self.end = end if start <= end else start
        self.closed = closed
        self.shape = end - start - 1 + ((self.end in self) + (self.start in self))
        
    def __call__(self, other):
        if type(other) is interval:
            return self.intersect(other) 
        elif type(other) is str:
            return self.from_str(other)
        else:
            return self.contains(other)
    
    def __add__(self, other):
        start = self.start + other
        end = self.end + other
        return interval(start, end, self.closed)
    
    def __radd__(self, other):
        return self.__add__(other)
    
    def __sub__(self, other):
        start = self.start - other
        end = self.end - other
        return interval(start, end, self.closed)
    
    def __rsub__(self, other):
        return self.__sub__(other)
    
    def __mul__(self, other):
        start = self.start * other
        end = self.end * other
        return interval(start, end, self.closed)
    
    def __rmul__(self, other):
        return self.__mul__(other)
    
    def __div__(self, other):
        start = self.start / other
        end = self.end / other
        return interval(start, end, self.closed)
    
    def __rdiv__(self, other):
        return self.__div__(other)
    
    
    def __contains__(self, value):
        return self.contains(value)
    
    def __str__(self):
        if self.closed == 'left':
            return f'[{self.start},{self.end})'
        if self.closed == 'right':
            return f'({self.start},{self.end}]'
        if self.closed == 'both':
            return f'[{self.start},{self.end}]'
        if self.closed == 'neither':
            return f'({self.start},{self.end})'
    
    def __repr__(self):
        if self.closed == 'left':
            return f'[{self.start},{self.end})'
        if self.closed == 'right':
            return f'({self.start},{self.end}]'
        if self.closed == 'both':
            return f'[{self.start},{self.end}]'
        if self.closed == 'neither':
            return f'({self.start},{self.end})'
    
    def contains(self, other):
        if self.closed == 'left':
            return (self.start <= other)&(other < self.end)
        if self.closed == 'right':
            return (self.start < other)&(other <= self.end)
        if self.closed == 'both':
            return (self.start <= other)&(other <= self.end)
        if self.closed == 'neither':
            return (self.start < other)&(other < self.end)
    
    def intersect(self, other):
        closed_type = {0: 'neither', 1: 'left', 2: 'right', 3: 'both'}
        if (self.start in other or self.start==other.start) and (self.end in other or self.end==other.end):
            closed = closed_type[(self.end in self and self.end in other)*2 + (self.start in other and self.start in self)]
            return interval(self.start, self.end, closed)
        if (other.start in self or self.start==other.start) and (other.end in self or self.end==other.end):
            closed = closed_type[(other.end in self and other.end in other)*2 + (other.start in other and other.start in self)]
            return interval(other.start, other.end, closed)
        if other.start in self and self.end in other:
            start = other.start
            end = self.end
            closed = closed_type[(self.end in self)*2 + (other.start in other)]
            return interval(start,end,closed)
        if other.end in self and self.start in other:
            start = self.start
            end = other.end
            closed = closed_type[(self.start in self) + (other.end in other)*2]
            return interval(start,end,closed)
        return None
    
    def minint(self):
        x = str(self.start).split('.')
        if len(x) == 0:
            return None
        elif len(x) == 1:
            if self.start in self:
                return int(self.start)
            elif self.start+1 in self:
                return int(self.start+1)
            else:
                return None
        else:
            if int(x[1][0]) == 0:
                if self.start in self:
                    return int(self.start)
                elif self.start+1 in self:
                    return int(self.start+1)
                else:
                    return None
            else:
                if self.start+1 in self:
                    return int(self.start+1)
                else:
                    return None
    
    def maxint(self):
        x = str(self.end).split('.')
        if len(x) == 0:
            return None
        elif len(x) == 1:
            if self.end in self:
                return int(self.end)
            elif self.end-1 in self:
                return int(self.end-1)
            else:
                return None
        else:
            if int(x[1][0]) == 0:
                if self.end in self:
                    return int(self.end)
                elif self.end-1 in self:
                    return int(self.end-1)
                else:
                    return None
            else:
                if self.end-1 in self:
                    return int(self.end-1)
                else:
                    return None
    
    def to_slice(self, step = 1):
        return slice(int(self.minint()), int(self.maxint()+1), int(step)) if step > 0 else slice(int(self.maxint()), int(self.minint()-1), int(step))
    
    def from_str(data: str):
        closed_type = {'()': 'neither', '[)': 'left', '(]': 'right', '[]': 'both'}
        tmp = data.split(',')
        start = float(tmp[0][1:])
        end = float(tmp[-1][:-1])
        tp = tmp[0][0] + tmp[-1][-1]
        return interval(start,end,closed_type[tp])

def interval_maker(size,sx,sy):
    ans = []
    for s in [sx,sy]:
        if s > 0:
            ans += [interval(0,size,'left')]
        else:
            ans += [interval(-size,0,'right')]
    return np.array(ans, dtype=interval)

def verbose_timedelta(delta):
    d, s = divmod(delta, 60*60*24)
    h, s = divmod(s, 60*60)
    m, s = divmod(s, 60)
    s, ms = divmod(s, 1)
    ms = ms*1000
    labels = ['d', 'h', 'm', 's', 'ms']   
    dhms = ['%s%s' % (i, lbl) for i, lbl in zip([int(d), int(h), int(m), int(s), int(ms)], labels)]
    for start in range(len(dhms)):
        if not dhms[start].startswith('0'):
            break
    for end in range(len(dhms)-1, -1, -1):
        if not dhms[end].startswith('0'):
            break  
    return ':'.join(dhms[start:end+1])

def my_logpmf(x, depth, prob):
    return np.nan_to_num(gammaln(depth+1) - gammaln(x+1) - gammaln(depth-x+1) + x*np.log(prob) + (depth-x)*np.log(1-prob), nan=0.0)

def my_log10pmf(x, depth, prob):
    return my_logpmf(x, depth, prob)*np.log10(np.exp(1))
