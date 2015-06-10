# -*- coding: utf-8 -*-
"""
Created on Sun Jun  7 14:27:21 2015

@author: dima
"""
def close(series, ele):
    series = dilate(series, ele)
    series = erode(series, ele)
    return series

def dilate(series, ele):
    new_series = np.copy(series)
    for n in range(len(series)):
        new_series[n] = max(series[max(0, n - ele):min(len(series), n + ele + 1)])
    return new_series

def erode(series, ele):
    new_series = np.copy(series)
    for n in range(len(series)):
        new_series[n] = min(series[max(0, n - ele):min(len(series), n + ele + 1)])
    return new_series
###############################################################################

def find_dense_snps_grid(x, chr_length, snp_number, scale=1, shift=0):
    "characteristic length scale of the Poisson distribution"
    mu = scale * int(chr_length / snp_number)
    x_grid = np.arange(0 - shift * mu, chr_length, mu)
    snps_on_grid = np.empty(x_grid.shape, dtype=int)

    for ii in range(0, len(x_grid)):
        snps_on_grid[ii] = sum(abs(x - x_grid[ii]) < mu)

    return (x_grid, snps_on_grid)

def find_dense_snps(self, scale=1, shift=0):
    "characteristic length scale of the Poisson distribution"
    mu = scale * int(chr_length / self.data.shape[0])
    x_grid = np.copy(self['x'])
    snps_on_grid = np.empty(x_grid.shape, dtype=int)

    for ii in range(0, len(x_grid)):
        snps_on_grid[ii] = sum(abs(self['x'] - x_grid[ii]) < mu)

    return (snps_on_grid)

def threshold_dense_snps(x_, scale=1, shift=0, threshold=3, inv_subscale=10, grid=True):
    x = [None] * inv_subscale;
    y = [None] * inv_subscale
    snp_number = len(x)
    if grid:
        for nn in range(inv_subscale):
            x[nn], y[nn] = find_dense_snps_grid(x_, chr_length, snp_number, scale, nn / inv_subscale)
        x = np.fliplr(np.array(x).T).ravel()
        y = np.fliplr(np.array(y).T).ravel()
    else:
        y = np.array(self.find_dense_snps(scale, 0))
#        x = np.copy(self['x'])
        inv_subscale = 1

    crowded_region = close(y > threshold * scale, 1 * inv_subscale)
    return (crowded_region, y, x)

def plot_snp_density(chr_length, scale=1, shift=0, threshold=3, inv_subscale=10, grid=True):
    crowded_region, y, x = threshold_dense_snps(scale=1, shift=0, threshold=3, inv_subscale=10, grid=True)

    import matplotlib.pyplot as plt

    y_max = np.percentile(y, 99);
    fig = plt.figure()
    fig.suptitle('snp density')
    plt.plot(x * 1e-6, y, 'r.-', label = 'density')
    plt.plot([0, chr_length * 1e-6], np.array([1, 1]) * threshold * scale, 'g-', label = 'threshold')
    plt.plot(x * 1e-6, y_max * self.crowded_region, 'b-', label = 'crowded region label')
    plt.legend()        
    plt.ylim((0, np.ceil(1.1*y_max))) 
    plt.show()

    return crowded_region, x