import numpy as np


def cgpercentiles(aData_in,
                  aPercentiles_in,
                  missing_value_in=None):
    """
    Calculates the percentiles of a given dataset.

    Args:
      aData_in (array-like): The input data.
      aPercentiles_in (array-like): The percentiles to calculate.
      missing_value_in (float, optional): The value to treat as missing. Defaults to None.

    Returns:
      list: The percentiles of the input data.
    """

    # Assume the worst.
    aData_in0 = np.array(aData_in)
    aData_in = aData_in0.flatten()
    aData_in = np.reshape(aData_in,  (aData_in.size))
    result = -1
    if aData_in.size == 0:
        print('Input data is required.')
        return -1

    aPercentiles_in = np.array(aPercentiles_in)
    aPercentiles_in = aPercentiles_in.flatten()
    if aPercentiles_in.size == 0:
        aPercentiles_in = np.array([0.25, 0.50, 0.75])

    if missing_value_in is not None:
        iFlag_missing_value = 1
        missing_value = missing_value_in
        dummy_index = np.where(aData_in == missing_value)
        aData_in[dummy_index] = np.nan
        good_index = np.where(np.isfinite(aData_in))
        aData_copy = aData_in[good_index]
    else:
        iFlag_missing_value = 0
        good_index = np.where(np.isfinite(aData_in))
        aData_copy = aData_in[good_index]

    num = aData_copy.size
    num_per = aPercentiles_in.size
    min_count = int(num/num_per)

    min_value = np.min(aData_copy)
    nan_index1 = np.where(aData_copy == min_value)

    max_value = np.max(aData_copy)
    nan_index2 = np.where(aData_copy == max_value)

    # np.count_nonzero(aData_copy == min_value)
    nan_count1 = nan_index1[0].size
    # np.count_nonzero(aData_copy == max_value)
    nan_count2 = nan_index2[0].size

    scenario = 4
    # there is no zero problem
    if ((min_value >= 0.0) or (max_value <= 0.0)):
        if ((nan_count1 > min_count) and (nan_count2 > min_count)):
            # both sides have problems
            scenario = 1
        else:
            if ((nan_count1 > min_count) and nan_count2 <= min_count):
                # min side has problem
                scenario = 2
            else:
                if ((nan_count2 > min_count) and (nan_count1 <= min_count)):
                    # max side has problem
                    scenario = 3
                else:
                    # everything is OK
                    scenario = 4

    else:

        nan_index3 = np.where(aData_copy == 0.0)
        nan_count3 = nan_index3[0].size

        if (nan_count3 >= min_count):
            # there is a probelm now

            if ((nan_count1 >= min_count) and (nan_count2 >= min_count)):
                # we have 3 problematic zones
                scenario = 5
            else:
                if ((nan_count1 >= min_count) and (nan_count2 <= min_count)):
                    # min side has problem and zero
                    scenario = 6
                else:
                    if ((nan_count2 >= min_count) and (nan_count1 <= min_count)):
                        # max side has problem and zero
                        scenario = 7
                    else:
                        # only zero problem
                        scenario = 8

        else:
            # not a 0.0 problem
            if ((nan_count1 >= min_count) and (nan_count2 >= min_count)):
                # we have 2 problematic zones
                scenario = 1
            else:
                if ((nan_count1 >= min_count) and (nan_count2 <= min_count)):
                    # min side has problem
                    scenario = 2
                else:
                    if ((nan_count2 >= min_count) and (nan_count1 <= min_count)):
                        # max side has problem
                        scenario = 3
                    else:
                        # only zero problem
                        scenario = 4

    # now

    if scenario == 1:
        good_index = np.where((aData_copy > min_value) &
                              (aData_copy < max_value))
        good_count = good_index[0].size
        if good_count <= num_per:
            return 0
        else:
            sample_count1 = np.ceil(
                float(good_count) / (num_per-1)) + 1  # why plus one?
            bad_copy1 = np.full(sample_count1, fill_value=min_value)
            sample_count2 = np.floor(float(good_count) / (num_per-1)) - 1
            bad_copy2 = np.full(sample_count2, fill_value=max_value)
            data_copy2 = [bad_copy1, bad_copy2, aData_copy[good_index]]

    else:
        if scenario == 2:
            good_index = np.where(aData_copy > min_value)
            good_count = good_index[0].size
            if good_count <= num_per:
                return 0
            else:
                sample_count = int(np.ceil(float(good_count) / (num_per)))
                bad_copy = np.full(sample_count, fill_value=min_value)
                # data_copy2 = [bad_copy,aData_copy[good_index]]
                data_copy2 = np.concatenate((bad_copy, aData_copy[good_index]))

        else:
            if scenario == 3:
                good_index = np.where(aData_copy < max_value)
                good_count = good_index[0].size
                if good_count <= num_per:
                    return 0
                else:
                    sample_count = int(np.ceil(float(good_count) / (num_per)))
                    bad_copy = np.full(sample_count, fill_value=max_value)
                    dummy = aData_copy[good_index]
                    data_copy2 = np.concatenate((bad_copy, dummy))

            else:
                if scenario == 4:
                    data_copy2 = aData_copy
                else:
                    if scenario == 5:
                        good_index = np.where((aData_copy > min_value) & (
                            aData_copy != 0.0) & (aData_copy < max_value))
                        good_count = good_index[0].size
                        if good_count <= num_per:
                            return 0
                        else:
                            sample_count1 = int(
                                np.ceil(float(good_count) / (num_per-2)) + 1)
                            bad_copy1 = np.full(
                                sample_count1, fill_value=min_value)
                            sample_count2 = int(
                                np.floor(float(good_count) / (num_per-2)) - 1)
                            bad_copy2 = np.full(sample_count2, fill_value=0.0)

                            sample_count3 = np.floor(
                                float(good_count) / (num_per-2)) - 1
                            bad_copy2 = np.full(sample_count3, fill_value=0.0)

                            # data_copy2 = [bad_copy1,bad_copy2, aData_copy[good_index]]
                            data_copy2 = np.concatenate(
                                (bad_copy1, bad_copy2, aData_copy[good_index]))

                    else:
                        if scenario == 6:
                            good_index = np.where(
                                (aData_copy > min_value) & (aData_copy != 0.0))
                            good_count = good_index[0].size
                            if good_count <= num_per:
                                return 0
                            else:
                                sample_count1 = int(
                                    np.ceil(float(good_count) / (num_per-1)) + 1)
                                bad_copy1 = np.full(
                                    sample_count1, fill_value=min_value)
                                sample_count2 = int(
                                    np.floor(float(good_count) / (num_per-1)) - 1)
                                bad_copy2 = np.full(
                                    sample_count2, fill_value=0.0)
                                # data_copy2 = [bad_copy1,bad_copy2, aData_copy[good_index]]
                                data_copy2 = np.concatenate(
                                    (bad_copy1, bad_copy2, aData_copy[good_index]))

                        else:
                            if scenario == 7:
                                good_index = np.where(
                                    (aData_copy != 0.0) & (aData_copy < max_value))
                                good_count = good_index[0].size
                                if good_count <= num_per:
                                    return 0
                                else:
                                    sample_count1 = int(
                                        np.ceil(float(good_count) / (num_per-1)) + 1)
                                    bad_copy1 = np.full(
                                        sample_count1, fill_value=0.0)
                                    sample_count2 = int(
                                        np.floor(float(good_count) / (num_per-1)) - 1)
                                    bad_copy2 = np.full(
                                        sample_count2, fill_value=max_value)
                                    # data_copy2 = [bad_copy1,bad_copy2, aData_copy[good_index]]
                                    data_copy2 = np.concatenate(
                                        (bad_copy1, bad_copy2, aData_copy[good_index]))

                            else:  # 8
                                good_index = np.where(aData_copy != 0.0)
                                good_count = good_index[0].size
                                if good_count <= num_per:
                                    return 0
                                else:
                                    sample_count = int(
                                        np.ceil(float(good_count) / (num_per)) - 1)
                                    bad_copy = np.full(
                                        sample_count, fill_value=0.0)
                                    # data_copy2 = [bad_copy,aData_copy[good_index]]
                                    data_copy2 = np.concatenate(
                                        (bad_copy, aData_copy[good_index]))

                                pass

    num = (data_copy2).size
    result = list()
    # Sort the data and find percentiles.
    for p in aPercentiles_in:
        a = np.percentile(data_copy2, p)
        result.append(a)

    return result
