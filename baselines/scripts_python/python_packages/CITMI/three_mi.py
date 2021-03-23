import numpy as np
import pandas as pd
import math

from baselines.scripts_python.python_packages.CITMI.tigramite.tigramite.independence_tests import CMIknn
from baselines.scripts_python.python_packages.CITMI.tigramite.tigramite.independence_tests import ParCorr
from baselines.scripts_python.python_packages.CITMI.ctmi_new import aligned_dict_to_df

from scipy import special, stats, spatial

from tigramite import tigramite_cython_code

# try:
#     from baselines.scripts_python.python_packages.CITMI.tigramite.tigramite import tigramite_cython_code
# except:
#     print("Could not import packages for CMIknn and GPDC estimation")


def indep_test(x, y, z=None, sig_samples=10000, p_value=True, measure="cmiknn", k=10, test_indep= True):
    if measure == "cmiknn":
        cd = CMIknn(mask_type=None, significance='shuffle_test', fixed_thres=None, sig_samples=sig_samples,
                    sig_blocklength=3, knn=k, confidence='bootstrap', conf_lev=0.9, conf_samples=10000,
                    conf_blocklength=1, verbosity=0)
    elif measure == "parcorr":
        cd = ParCorr(mask_type=None, significance='shuffle_test', fixed_thres=None, sig_samples=sig_samples,
                     sig_blocklength=3, confidence='bootstrap', conf_lev=0.9, conf_samples=10000, conf_blocklength=1,
                     verbosity=0)
    else:
        cd = None
        print("Independence measure '" + str(measure) + "' do not exist.")
        exit(0)
    dim_x = x.shape[1]
    dim_y = y.shape[1]

    ws = x.shape[1]
    if z is not None:
        # todo validate
        x_past = x[:-ws]
        y_past = y[:-ws]
        x = x[ws:].reset_index(drop=True)
        y = y[ws:].reset_index(drop=True)
        z = z[ws:].reset_index(drop=True)
        # x_past = x_past[x_past.columns[0]]
        # y_past = y_past[y_past.columns[0]]
        z = pd.concat([z, x_past, y_past], axis=1)
        # z = pd.concat([z, y_past], axis=1)

        dim_z = z.shape[1]
        X = np.concatenate((x.values, y.values, z.values), axis=1)
        xyz = np.array([0] * dim_x + [1] * dim_y + [2] * dim_z)
    else:
        # todo validate
        x_past = x[:-ws]
        y_past = y[:-ws]
        x = x[ws:].reset_index(drop=True)
        y = y[ws:].reset_index(drop=True)
        # x_past = x_past[x_past.columns[0]]
        # y_past = y_past[y_past.columns[0]]
        z = pd.concat([x_past, y_past], axis=1)
        # z = y_past.to_frame()

        dim_z = z.shape[1]
        X = np.concatenate((x.values, y.values, z.values), axis=1)
        xyz = np.array([0] * dim_x + [1] * dim_y + [2] * dim_z)

    Xx = np.concatenate((x.values, y.values), axis=1)
    xy = xyz[np.where(xyz != 2)]
    print(np.abs(cd.get_dependence_measure(Xx.T, xy)), np.abs(cd.get_dependence_measure(X.T, xyz)))
    # value = -np.abs(cd.get_dependence_measure(Xx.T, xy)) + np.abs(cd.get_dependence_measure(X.T, xyz))
    c1 = np.abs(cd.get_dependence_measure(Xx.T, xy))
    c2 = np.abs(cd.get_dependence_measure(X.T, xyz))
    value = (c1 - c2) / c1
    if p_value:
        # pvalue = cd.get_shuffle_significance(X.T, xyz, value)
        pvalue = get_shuffle_significance_3point(cd, X.T, xyz, value)
    else:
        pvalue = value
    return pvalue, value

def ctmi_gamma(x,y,z, p_value=True):
    gxy_list = []
    grx_list = []
    gry_list = []
    val_list = []
    pval_list = []
    for gxy in range(-5, 5):
        if gxy >= 0:
            for gry in range(0, 5 + gxy):
                grx = gry - gxy
                pval, val = ctmi_gamma_fix(x, y, z, [gxy, grx, gry], p_value=p_value)
                val_list.append(pval)
                pval_list.append(pval)
                gxy_list.append(gxy)
                grx_list.append(grx)
                gry_list.append(gry)
        else:
            for grx in range(0, 5 + gxy):
                gry = gxy + grx
                pval, val = ctmi_gamma_fix(x, y, z, [gxy, grx, gry], p_value=p_value)
                pval_list.append(pval)
                val_list.append(val)
                gxy_list.append(gxy)
                grx_list.append(grx)
                gry_list.append(gry)

    idx = np.argmax(val_list)
    return gxy_list[idx], grx_list[idx], gry_list[idx], pval_list[idx], val_list[idx]

def ctmi_gamma_fix(x, y, z, gamma_vector, p_value=True):
    gxy, grx, gry = gamma_vector
    names_z = [*z.keys()]
    if len(names_z) > 0:
        z = aligned_dict_to_df(z)
    else:
        z = None

    if gxy >= 0:
        if grx >= 0:
            y = y[gry:]
            if gry>0:
                z = z[:-gry]
            x = x[grx:]

            if gxy > 0:
                x = x[:-gxy]

        else:
            y = y[gxy:]
            x = x[:-gxy]

            z = z[-grx:]
            # x = x[grx:]
            if gry > 0:
                z = z[:-gry]

    elif gxy < 0:
        if gry >= 0:
            x = x[grx:]
            z = z[:-grx]

            if gry > 0:
                y = y[gry:]
            # y = y[gry:]

            y = y[:gxy]
        else:
            x = x[-gxy:]
            y = y[:gxy]

            z = z[-gry:]
            # x = x[grx:]
            if grx > 0:
                z = z[:-grx]
    else:
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True)
        # x = x[x.index % (iter1 + 1) == 0]
        # y = y[y.index % (iter2 + 1) == 0]
        x = x.reset_index(drop=True)
        y = y.reset_index(drop=True)

    x = x.reset_index(drop=True)
    y = y.reset_index(drop=True)
    z = z.reset_index(drop=True)

    m = min(x.shape[0], y.shape[0])
    x = x[:m]
    y = y[:m]
    return indep_test(x, y, z=z, p_value=p_value)



# def get_shuffle_significance_3point(ci, array, xyz, value, return_null_dist=False):
#     """Returns p-value for shuffle significance test.
#
#     For residual-based test statistics only the residuals are shuffled.
#
#     Parameters
#     ----------
#     array : array-like
#         data array with X, Y, Z in rows and observations in columns
#
#     xyz : array of ints
#         XYZ identifier array of shape (dim,).
#
#     value : number
#         Value of test statistic for unshuffled estimate.
#
#     Returns
#     -------
#     pval : float
#         p-value
#     """
#     def get_single_residuals(array, target_var,
#                               standardize=True,
#                               return_means=False):
#         """Returns residuals of linear multiple regression.
#
#         Performs a OLS regression of the variable indexed by target_var on the
#         conditions Z. Here array is assumed to contain X and Y as the first two
#         rows with the remaining rows (if present) containing the conditions Z.
#         Optionally returns the estimated regression line.
#
#         Parameters
#         ----------
#         array : array-like
#             data array with X, Y, Z in rows and observations in columns
#
#         target_var : {0, 1}
#             Variable to regress out conditions from.
#
#         standardize : bool, optional (default: True)
#             Whether to standardize the array beforehand. Must be used for
#             partial correlation.
#
#         return_means : bool, optional (default: False)
#             Whether to return the estimated regression line.
#
#         Returns
#         -------
#         resid [, mean] : array-like
#             The residual of the regression and optionally the estimated line.
#         """
#
#         dim, T = array.shape
#         dim_z = dim - 2
#
#         # Standardize
#         if standardize:
#             array -= array.mean(axis=1).reshape(dim, 1)
#             array /= array.std(axis=1).reshape(dim, 1)
#             if np.isnan(array).sum() != 0:
#                 raise ValueError("nans after standardizing, "
#                                  "possibly constant array!")
#
#         y = array[target_var, :]
#
#         if dim_z > 0:
#             z = np.fastCopyAndTranspose(array[2:, :])
#             beta_hat = np.linalg.lstsq(z, y, rcond=None)[0]
#             mean = np.dot(z, beta_hat)
#             resid = y - mean
#         else:
#             resid = y
#             mean = None
#
#         if return_means:
#             return (resid, mean)
#         return resid
#
#     cd = CMIknn(mask_type=None, significance='shuffle_test', fixed_thres=None, sig_samples=10000,
#                 sig_blocklength=3, knn=10, confidence='bootstrap', conf_lev=0.9, conf_samples=10000,
#                 conf_blocklength=1, verbosity=0)
#
#     x_vals = get_single_residuals(array, target_var=0)
#     y_vals = get_single_residuals(array, target_var=1)
#     array_resid = np.array([x_vals, y_vals])
#     xyz_resid = np.array([0, 1])
#
#     null_dist = _get_shuffle_dist(cd, array_resid, xyz_resid, cd.get_dependence_measure,
#                                        sig_samples=10000,
#                                        sig_blocklength=3,
#                                        verbosity=0)
#
#     # test independence or test dependence
#     pval = (null_dist <= np.abs(value)).mean()
#
#     # Adjust p-value for two-sided measures
#     # if pval < 1.:
#     #     pval *= 2.
#
#     if return_null_dist:
#         return pval, null_dist
#     return pval
#
#
#
# def _get_shuffle_dist(ci, array_cond, array, xyz, dependence_measure,
#                       sig_samples, sig_blocklength=None,
#                       verbosity=0):
#     """Returns shuffle distribution of test statistic.
#
#     The rows in array corresponding to the X-variable are shuffled using
#     a block-shuffle approach.
#
#     Parameters
#     ----------
#     array : array-like
#         data array with X, Y, Z in rows and observations in columns
#
#     xyz : array of ints
#         XYZ identifier array of shape (dim,).
#
#    dependence_measure : object
#        Dependence measure function must be of form
#        dependence_measure(array, xyz) and return a numeric value
#
#     sig_samples : int, optional (default: 100)
#         Number of samples for shuffle significance test.
#
#     sig_blocklength : int, optional (default: None)
#         Block length for block-shuffle significance test. If None, the
#         block length is determined from the decay of the autocovariance as
#         explained in [1]_.
#
#     verbosity : int, optional (default: 0)
#         Level of verbosity.
#
#     Returns
#     -------
#     null_dist : array of shape (sig_samples,)
#         Contains the sorted test statistic values estimated from the
#         shuffled arrays.
#     """
#     xy = xyz[np.where(xyz != 2)]
#
#     dim, T = array_cond.shape
#
#     x_indices = np.where(xyz == 0)[0]
#     dim_x = len(x_indices)
#
#     if sig_blocklength is None:
#         sig_blocklength = ci._get_block_length(array_cond, xyz,
#                                                  mode='significance')
#
#     n_blks = int(math.floor(float(T)/sig_blocklength))
#     # print 'n_blks ', n_blks
#     if verbosity > 2:
#         print("            Significance test with block-length = %d "
#               "..." % (sig_blocklength))
#
#     array_shuffled_cond = np.copy(array_cond)
#     array_shuffled = np.copy(array)
#     block_starts = np.arange(0, T - sig_blocklength + 1, sig_blocklength)
#
#     # Dividing the array up into n_blks of length sig_blocklength may
#     # leave a tail. This tail is later randomly inserted
#     tail_cond = array_cond[x_indices, n_blks*sig_blocklength:]
#     tail = array[x_indices, n_blks*sig_blocklength:]
#
#     null_dist = np.zeros(sig_samples)
#     for sam in range(sig_samples):
#
#         blk_starts = np.random.permutation(block_starts)[:n_blks]
#
#         x_shuffled_cond = np.zeros((dim_x, n_blks*sig_blocklength),
#                               dtype=array_cond.dtype)
#         x_shuffled = np.zeros((dim_x, n_blks*sig_blocklength),
#                               dtype=array.dtype)
#
#         for i, index in enumerate(x_indices):
#             for blk in range(sig_blocklength):
#                 x_shuffled_cond[i, blk::sig_blocklength] = \
#                     array_cond[index, blk_starts + blk]
#                 x_shuffled[i, blk::sig_blocklength] = \
#                     array[index, blk_starts + blk]
#
#
#         # Insert tail randomly somewhere
#         if tail_cond.shape[1] > 0:
#             insert_tail_at = np.random.choice(block_starts)
#             x_shuffled_cond = np.insert(x_shuffled_cond, insert_tail_at,
#                                    tail_cond.T, axis=1)
#             x_shuffled = np.insert(x_shuffled, insert_tail_at,
#                                    tail.T, axis=1)
#
#         for i, index in enumerate(x_indices):
#             array_shuffled_cond[index] = x_shuffled_cond[i]
#             array_shuffled[index] = x_shuffled[i]
#
#         null_dist[sam] = -dependence_measure(array=array_shuffled, xyz=xy) + dependence_measure(array=array_shuffled_cond, xyz=xyz)
#
#     null_dist.sort()
#
#     return null_dist


def get_shuffle_significance_3point(ci, array, xyz, value,
                             return_null_dist=False):
    """Returns p-value for nearest-neighbor shuffle significance test.

    For non-empty Z, overwrites get_shuffle_significance from the parent
    class  which is a block shuffle test, which does not preserve
    dependencies of X and Y with Z. Here the parameter shuffle_neighbors is
    used to permute only those values :math:`x_i` and :math:`x_j` for which
    :math:`z_j` is among the nearest niehgbors of :math:`z_i`. If Z is
    empty, the block-shuffle test is used.

    Parameters
    ----------
    array : array-like
        data array with X, Y, Z in rows and observations in columns

    xyz : array of ints
        XYZ identifier array of shape (dim,).

    value : number
        Value of test statistic for unshuffled estimate.

    Returns
    -------
    pval : float
        p-value
    """
    dim, T = array.shape

    # Skip shuffle test if value is above threshold
    # if value > self.minimum threshold:
    #     if return_null_dist:
    #         return 0., None
    #     else:
    #         return 0.

    # max_neighbors = max(1, int(max_neighbor_ratio*T))
    x_indices = np.where(xyz == 0)[0]
    z_indices = np.where(xyz == 2)[0]

    if len(z_indices) > 0 and ci.shuffle_neighbors < T:
        if ci.verbosity > 2:
            print("            nearest-neighbor shuffle significance "
                  "test with n = %d and %d surrogates" % (
                ci.shuffle_neighbors, ci.sig_samples))

        # Get nearest neighbors around each sample point in Z
        z_array = np.fastCopyAndTranspose(array[z_indices, :])
        tree_xyz = spatial.cKDTree(z_array)
        neighbors = tree_xyz.query(z_array,
                                   k=ci.shuffle_neighbors,
                                   p=np.inf,
                                   eps=0.)[1].astype('int32')

        null_dist = np.zeros(ci.sig_samples)
        for sam in range(ci.sig_samples):

            # Generate random order in which to go through indices loop in
            # next step
            order = np.random.permutation(T).astype('int32')
            # print(order[:5])
            # Select a series of neighbor indices that contains as few as
            # possible duplicates
            restricted_permutation = \
                tigramite_cython_code._get_restricted_permutation_cython(
                    T=T,
                    shuffle_neighbors=ci.shuffle_neighbors,
                    neighbors=neighbors,
                    order=order)

            array_shuffled = np.copy(array)
            for i in x_indices:
                array_shuffled[i] = array[i, restricted_permutation]

            xy_idx = np.where(xyz != 2)
            xy = xyz[xy_idx]
            array_xy_shuffled = array_shuffled[xy_idx[0], :]
            c1 = np.abs(ci.get_dependence_measure(array=array_xy_shuffled, xyz=xy))
            c2 = np.abs(ci.get_dependence_measure(array=array_shuffled, xyz=xyz))
            null_dist[sam] = (c1 - c2)/c1
            # ci.get_dependence_measure(array_shuffled,xyz)
    else:
        null_dist = \
                ci._get_shuffle_dist(array, xyz,
                                       ci.get_dependence_measure,
                                       sig_samples=ci.sig_samples,
                                       sig_blocklength=ci.sig_blocklength,
                                       verbosity=ci.verbosity)

    # Sort
    null_dist.sort()
    pval = (null_dist <= value).mean()

    if return_null_dist:
        return pval, null_dist
    return pval
