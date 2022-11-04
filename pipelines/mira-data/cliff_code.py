import numpy as np
import warnings
from scipy import sparse

def _residual_transform(X, pi_j_hat, n_i):
    
    assert(isinstance(X, np.ndarray))
    assert(isinstance(pi_j_hat, np.ndarray))
    assert(isinstance(n_i, np.ndarray))
    pi_j_hat = np.squeeze(pi_j_hat)[np.newaxis, :]
    n_i = np.squeeze(n_i)[:, np.newaxis]

    mu_ij_hat = n_i * pi_j_hat

    count_dif = n_i - X
    expected_count_dif = n_i - mu_ij_hat

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        r_ij = np.multiply(
            np.sign(X - mu_ij_hat), 
            np.sqrt(
            np.where(X > 0, 2 * np.multiply(X, np.log(X / mu_ij_hat)), 0) + \
            2 * np.multiply(count_dif, np.log(count_dif / expected_count_dif))
            )
        )

    return np.clip(np.nan_to_num(r_ij), -10, 10)

def _get_pi(X):
    return np.array(X.sum(0)).reshape(-1)/X.sum()

def _get_n(X):
    return np.array(X.sum(-1)).reshape(-1)


def deviance_transform(X, subset_mask = None):
    '''
    Parameters
    ----------
    X : Scipy sparse matrix (N genes x N cells)
        Scipy sparse matrix of ATAC-seq counts (N cells x N peaks). The values may be binary or positive integers.
    subset_mask : np.ndarray (N cells) or none
        Calculating a dense version of the whole matrix is memory intensive, but the parameters of the transformation can 
        be estimated from the whole matrix, then applied to just a subset. If "subset_mask" is *none*, the transformation
        will be applied to the whole matrix. If subset mask is a boolean mask of shape (N peaks), then the transformed 
        accessibility of only those peaks will be returned.

    Returns
    -------
    residuals : np.ndarray (N genes x N peaks -- subset)

    Example
    -------

    If you are working with an AnnData object called "atac_data",
    to get the transformed values for only the first chromosome:

    >>> import cliff_code
    >>> import numpy as np
    >>> residuals = cliff_code.deviance_transform(atac_data.X, subset_mask = atac_data.var['chrom'].values == 'chr1')

    Then, for very simple KNN smoothing you can use the "connectivities" matrix in the AnnData object, which is calculated
    during the nearest neighbors step of the analysis:

    >>> smoothed = atac_data.obsp['connectivities'].dot(residuals)

    '''

    assert sparse.isspmatrix(X)
    X = X.tocsr()

    if subset_mask is None:
        subset_mask = np.ones(X.shape[-1]).astype(bool)

    n_i = _get_n(X)
    p_i = _get_pi(X)
    X = X[:, subset_mask].toarray()

    residuals = _residual_transform(X, p_i[subset_mask], n_i)

    return residuals

    