import datetime
import numpy as np
import pytest

from pyearth.toolbox.date import (
    TimerError,
    day_in_month,
    day_of_year,
    dt2cal,
    from_jd,
    leap_year,
    pytimer,
    to_jd,
)


class TestLeapYear:
    def test_leap_years(self):
        assert leap_year(2024) is True
        assert leap_year(2000) is True
        assert leap_year(2016) is True

    def test_common_years(self):
        assert leap_year(2023) is False
        assert leap_year(1900) is False
        assert leap_year(2100) is False

    def test_invalid_inputs(self):
        with pytest.raises(ValueError):
            leap_year(0)
        with pytest.raises(ValueError):
            leap_year(None)
        with pytest.raises(TypeError):
            leap_year("2024")


class TestDayInMonth:
    def test_standard_months(self):
        assert day_in_month(2024, 1) == 31
        assert day_in_month(2024, 4) == 30
        assert day_in_month(2024, 2) == 29
        assert day_in_month(2023, 2) == 28

    def test_leap_year_override(self):
        # Force 2024 February to 28 days
        assert day_in_month(2024, 2, iFlag_leap_year_in=0) == 28
        # Force 2023 February to 29 days
        assert day_in_month(2023, 2, iFlag_leap_year_in=1) == 29

    def test_invalid_inputs(self):
        with pytest.raises(ValueError):
            day_in_month(2024, 0)
        with pytest.raises(ValueError):
            day_in_month(2024, 13)
        with pytest.raises(ValueError):
            day_in_month(0, 5)
        with pytest.raises(ValueError):
            day_in_month(2024, 2, iFlag_leap_year_in=2)


class TestDayOfYear:
    def test_day_of_year(self):
        assert day_of_year(2024, 1, 1) == 1
        assert day_of_year(2024, 2, 29) == 60
        assert day_of_year(2024, 3, 1) == 61
        assert day_of_year(2023, 3, 1) == 60
        assert day_of_year(2024, 12, 31) == 366
        assert day_of_year(2023, 12, 31) == 365

    def test_invalid_dates(self):
        with pytest.raises(ValueError):
            day_of_year(2023, 2, 29)  # Not a leap year
        with pytest.raises(ValueError):
            day_of_year(2024, 4, 31)  # April has 30 days
        with pytest.raises(TypeError):
            day_of_year("2024", 1, 1)


class TestDt2Cal:
    def test_1d_array(self):
        dates = np.array(["2024-01-15T10:30:45.123456", "2024-02-29T00:00:00"], dtype="datetime64")
        cal = dt2cal(dates)
        assert cal.shape == (2, 7)
        assert np.array_equal(cal[0], [2024, 1, 15, 10, 30, 45, 123456])
        assert np.array_equal(cal[1], [2024, 2, 29, 0, 0, 0, 0])

    def test_scalar_datetime(self):
        dt = np.datetime64("2024-02-29T12:00:00")
        cal = dt2cal(dt)
        assert cal.shape == (7,)
        assert np.array_equal(cal, [2024, 2, 29, 12, 0, 0, 0])

    def test_invalid_type(self):
        with pytest.raises(TypeError):
            dt2cal(np.array([1, 2, 3]))
        with pytest.raises(ValueError):
            dt2cal(np.array([], dtype="datetime64"))


class TestJulian:
    def test_to_and_from_jd(self):
        dt = datetime.datetime(2024, 2, 29, 12, 30, 0)
        jd = to_jd(dt, fmt="jd")
        dt_rec = from_jd(jd, fmt="jd")
        assert dt_rec.year == dt.year
        assert dt_rec.month == dt.month
        assert dt_rec.day == dt.day
        assert dt_rec.hour == dt.hour
        assert dt_rec.minute == dt.minute
        assert abs(dt_rec.second - dt.second) <= 1

    def test_formats(self):
        dt = datetime.datetime(2024, 1, 1, 0, 0, 0)
        jd = to_jd(dt, fmt="jd")
        mjd = to_jd(dt, fmt="mjd")
        rjd = to_jd(dt, fmt="rjd")
        assert np.isclose(mjd, jd - 2400000.5)
        assert np.isclose(rjd, jd - 2400000.0)

        dt_from_mjd = from_jd(mjd, fmt="mjd")
        assert dt_from_mjd.year == 2024
        assert dt_from_mjd.month == 1
        assert dt_from_mjd.day == 1


class TestPyTimer:
    def test_basic_timer(self):
        timer = pytimer(name="test")
        timer.start()
        elapsed = timer.elapsed()
        assert elapsed >= 0.0
        timer.lap("lap1")
        assert len(timer.get_lap_times()) == 1
        stopped_elapsed = timer.stop()
        assert stopped_elapsed >= elapsed

    def test_context_manager(self):
        with pytimer(name="context") as t:
            x = sum(i for i in range(1000))
        assert t.get_total_time() >= 0.0

    def test_error_handling(self):
        timer = pytimer()
        with pytest.raises(TimerError):
            timer.stop()  # Not started
        with pytest.raises(TimerError):
            timer.elapsed()  # Not started
        timer.start()
        with pytest.raises(TimerError):
            timer.start()  # Already started
