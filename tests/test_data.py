import numpy as np
import pytest

from pyearth.toolbox.data import (
    cgpercentiles,
    check_if_duplicates,
    convert_time_series_daily_to_monthly,
    list_max,
    list_min,
    remove_outliers,
)


class TestListUtils:
    def test_check_if_duplicates(self):
        assert check_if_duplicates([1, 2, 3, 4]) == 1  # Unique -> 1
        assert check_if_duplicates([1, 2, 2, 3]) == 0  # Duplicates -> 0
        assert check_if_duplicates([]) == 1
        assert check_if_duplicates([42]) == 1

    def test_list_min_max(self):
        pList = [[10.0, 20.0, np.nan], [5.0, 30.0], [np.nan, np.nan], []]
        assert list_min(pList) == 5.0
        assert list_max(pList) == 30.0

    def test_list_min_max_all_nan(self):
        pList = [[np.nan, np.nan], []]
        assert np.isnan(list_min(pList))
        assert np.isnan(list_max(pList))


class TestRemoveOutliers:
    def test_remove_outliers(self):
        data = np.array([10.0, 11.0, 10.5, 9.8, 10.2, 100.0, -50.0])
        cleaned = remove_outliers(data, dOutlier_percentage=0.1)
        assert 100.0 not in cleaned
        assert -50.0 not in cleaned
        assert len(cleaned) >= 4


class TestCgPercentiles:
    def test_basic_percentiles(self):
        data = np.linspace(1, 100, 100)
        pcts = cgpercentiles(data, [25, 50, 75])
        assert isinstance(pcts, list)
        assert len(pcts) == 3
        assert np.isclose(pcts[1], 50.5, atol=1.0)


class TestConvertDailyToMonthly:
    def test_daily_to_monthly_mean(self):
        # 31 days of January 2024
        daily_jan = np.full(31, 10.0)
        monthly = convert_time_series_daily_to_monthly(
            daily_jan, 2024, 1, 1, 2024, 1, 31, sType_in="mean"
        )
        assert len(monthly) == 1
        assert np.isclose(monthly[0], 10.0)

    def test_daily_to_monthly_sum(self):
        # 29 days of February 2024 (leap year)
        daily_feb = np.full(29, 2.0)
        monthly = convert_time_series_daily_to_monthly(
            daily_feb, 2024, 2, 1, 2024, 2, 29, sType_in="sum"
        )
        assert len(monthly) == 1
        assert np.isclose(monthly[0], 58.0)
