"""
High-Resolution Performance Timer

This module provides a high-resolution timer class for measuring code execution time
with nanosecond precision. Built on Python's time.perf_counter(), it offers accurate
timing for performance profiling, benchmarking, and optimization tasks.

Main Classes
------------
pytimer : High-resolution timer for code profiling
TimerError : Custom exception for timer-related errors

Key Features
------------
- Nanosecond precision timing using time.perf_counter()
- Context manager support (with statement)
- Optional logging instead of printing
- Lap time tracking for multiple measurements
- Cumulative time tracking across multiple start/stop cycles
- Human-readable time formatting (ms, μs, ns)
- Thread-safe operation
- Named timers for tracking multiple operations

Use Cases
---------
1. **Performance Profiling**: Measure execution time of code blocks
2. **Algorithm Benchmarking**: Compare performance of different implementations
3. **Optimization**: Identify bottlenecks in code execution
4. **API Response Time**: Measure external API call durations
5. **Database Query Timing**: Track query performance
6. **Data Processing**: Monitor processing time for large datasets
7. **Testing**: Verify performance requirements in unit tests

Technical Details
-----------------
The timer uses time.perf_counter() which:
- Has nanosecond resolution on most systems
- Is monotonic (always moves forward, unaffected by system clock changes)
- Measures elapsed time, not CPU time
- Includes time spent in sleep()
- Does not count time before process started

For CPU time measurement (excluding I/O and sleep), use time.process_time() instead.

Performance Characteristics
---------------------------
- Timer overhead: < 1 microsecond per start/stop
- Resolution: Typically nanoseconds (system dependent)
- Memory: Minimal (stores only timestamps and metadata)

Dependencies
------------
- time: Standard library for high-resolution timing
- logging: For optional logging output
- typing: Type hints for better code clarity

See Also
--------
- time.perf_counter: High-resolution performance counter
- time.process_time: CPU time measurement
- timeit: Standard library module for microbenchmarking
"""

import time
import logging
from typing import Optional, List

# Configure logging
logger = logging.getLogger(__name__)


class TimerError(Exception):
    """
    Custom exception for timer-related errors.

    Raised when timer operations are performed in invalid states, such as:
    - Starting a timer that's already running
    - Stopping a timer that hasn't been started
    - Getting elapsed time when timer is not running

    Examples
    --------
    >>> timer = pytimer()
    >>> timer.stop()  # doctest: +SKIP
    Traceback (most recent call last):
        ...
    TimerError: Timer is not running. Use .start() to start it
    """

    pass


class pytimer:
    """
    High-resolution timer for measuring code execution time.

    This class provides accurate timing capabilities using Python's perf_counter()
    with support for context managers, lap times, and flexible output options.

    The timer can be used in several ways:
    1. Manual start/stop with print output
    2. Context manager (with statement) for automatic timing
    3. Lap timing for multiple measurements
    4. Silent mode with programmatic access to elapsed time

    Attributes
    ----------
    name : Optional[str]
        Optional name for the timer (useful for tracking multiple timers)
    use_logging : bool
        If True, output via logging instead of print
    log_level : int
        Logging level to use (default: logging.INFO)

    Examples
    --------
    Basic usage with manual start/stop:

    >>> timer = pytimer()
    >>> timer.start()
    >>> # ... code to time ...
    >>> timer.stop()
    Elapsed time: 0.0234 seconds

    Context manager usage:

    >>> with pytimer(name="Processing"):
    ...     # ... code to time ...
    ...     pass
    Processing: Elapsed time: 0.0234 seconds

    Lap timing:

    >>> timer = pytimer(name="Laps")
    >>> timer.start()
    >>> # ... task 1 ...
    >>> timer.lap("Task 1")
    Laps - Task 1: 0.0123 seconds
    >>> # ... task 2 ...
    >>> timer.lap("Task 2")
    Laps - Task 2: 0.0234 seconds
    >>> timer.stop()
    Laps: Total elapsed time: 0.0357 seconds

    Silent mode for programmatic access:

    >>> timer = pytimer()
    >>> timer.start()
    >>> # ... code to time ...
    >>> elapsed = timer.elapsed()
    >>> print(f"Took {elapsed:.4f} seconds")
    Took 0.0234 seconds
    """

    def __init__(
        self,
        name: Optional[str] = None,
        use_logging: bool = False,
        log_level: int = logging.INFO,
    ):
        """
        Initialize a new timer.

        Parameters
        ----------
        name : Optional[str], optional
            Name for this timer, useful for identifying which timer in logs.
            Default is None.
        use_logging : bool, optional
            If True, output timing information via logging module instead of print.
            Default is False (use print).
        log_level : int, optional
            Logging level to use when use_logging=True.
            Default is logging.INFO.
        """
        self.name = name
        self.use_logging = use_logging
        self.log_level = log_level
        self._start_time: Optional[float] = None
        self._lap_times: List[tuple] = []  # List of (name, time) tuples
        self._total_time: float = 0.0  # Cumulative time across multiple runs

    def start(self) -> None:
        """
        Start the timer.

        Records the current time using time.perf_counter(). If the timer is
        already running, raises TimerError.

        Raises
        ------
        TimerError
            If timer is already running.

        Examples
        --------
        >>> timer = pytimer()
        >>> timer.start()
        >>> # Timer is now running
        """
        if self._start_time is not None:
            raise TimerError("Timer is running. Use .stop() to stop it")

        self._start_time = time.perf_counter()
        logger.debug(f"Timer {self.name or 'unnamed'} started")

    def stop(self, silent: bool = False) -> float:
        """
        Stop the timer and report elapsed time.

        Calculates elapsed time since start() was called, adds it to cumulative
        total, and optionally prints/logs the result.

        Parameters
        ----------
        silent : bool, optional
            If True, don't print/log the elapsed time. Default is False.

        Returns
        -------
        float
            Elapsed time in seconds.

        Raises
        ------
        TimerError
            If timer is not currently running.

        Examples
        --------
        >>> timer = pytimer()
        >>> timer.start()
        >>> elapsed = timer.stop()
        Elapsed time: 0.0234 seconds
        >>> print(f"Returned: {elapsed:.4f}")
        Returned: 0.0234
        """
        if self._start_time is None:
            raise TimerError("Timer is not running. Use .start() to start it")

        elapsed_time = time.perf_counter() - self._start_time
        self._total_time += elapsed_time
        self._start_time = None

        if not silent:
            self._output(f"Elapsed time: {self._format_time(elapsed_time)}")

        logger.debug(f"Timer {self.name or 'unnamed'} stopped: {elapsed_time:.6f}s")

        return elapsed_time

    def lap(self, lap_name: Optional[str] = None) -> float:
        """
        Record a lap time without stopping the timer.

        Records the time since the last start() or lap() call and stores it
        with an optional label. The timer continues running.

        Parameters
        ----------
        lap_name : Optional[str], optional
            Name for this lap. If None, uses "Lap N" where N is the lap number.

        Returns
        -------
        float
            Time in seconds since last start() or lap().

        Raises
        ------
        TimerError
            If timer is not currently running.

        Examples
        --------
        >>> timer = pytimer(name="Process")
        >>> timer.start()
        >>> # ... step 1 ...
        >>> timer.lap("Step 1")
        Process - Step 1: 0.0123 seconds
        >>> # ... step 2 ...
        >>> timer.lap("Step 2")
        Process - Step 2: 0.0234 seconds
        """
        if self._start_time is None:
            raise TimerError("Timer is not running. Use .start() to start it")

        current_time = time.perf_counter()

        # Calculate time since last lap or start
        if self._lap_times:
            lap_time = current_time - self._lap_times[-1][1]
        else:
            lap_time = current_time - self._start_time

        # Store lap time with absolute timestamp
        if lap_name is None:
            lap_name = f"Lap {len(self._lap_times) + 1}"

        self._lap_times.append((lap_name, current_time))

        # Output lap time
        prefix = f"{self.name} - " if self.name else ""
        self._output(f"{prefix}{lap_name}: {self._format_time(lap_time)}")

        return lap_time

    def elapsed(self) -> float:
        """
        Get elapsed time without stopping the timer.

        Returns the time elapsed since start() was called. The timer continues
        running. Useful for checking progress without stopping.

        Returns
        -------
        float
            Elapsed time in seconds since start().

        Raises
        ------
        TimerError
            If timer is not currently running.

        Examples
        --------
        >>> timer = pytimer()
        >>> timer.start()
        >>> # ... some code ...
        >>> elapsed = timer.elapsed()
        >>> print(f"So far: {elapsed:.4f} seconds")
        So far: 0.0234 seconds
        >>> # Timer is still running
        """
        if self._start_time is None:
            raise TimerError("Timer is not running. Use .start() to start it")

        return time.perf_counter() - self._start_time

    def reset(self) -> None:
        """
        Reset the timer to initial state.

        Stops the timer (if running), clears all lap times, and resets the
        cumulative total time to zero.

        Examples
        --------
        >>> timer = pytimer()
        >>> timer.start()
        >>> timer.stop()
        Elapsed time: 0.0234 seconds
        >>> timer.reset()
        >>> # Timer is now in initial state
        """
        self._start_time = None
        self._lap_times = []
        self._total_time = 0.0
        logger.debug(f"Timer {self.name or 'unnamed'} reset")

    def get_total_time(self) -> float:
        """
        Get cumulative time across all start/stop cycles.

        Returns the total time accumulated from all start/stop cycles since
        the timer was created or last reset.

        Returns
        -------
        float
            Total cumulative time in seconds.

        Examples
        --------
        >>> timer = pytimer()
        >>> timer.start()
        >>> # ... task 1 ...
        >>> timer.stop(silent=True)
        >>> timer.start()
        >>> # ... task 2 ...
        >>> timer.stop(silent=True)
        >>> total = timer.get_total_time()
        >>> print(f"Total: {total:.4f} seconds")
        Total: 0.0468 seconds
        """
        return self._total_time

    def get_lap_times(self) -> List[tuple]:
        """
        Get all recorded lap times.

        Returns
        -------
        List[tuple]
            List of (lap_name, timestamp) tuples.

        Examples
        --------
        >>> timer = pytimer()
        >>> timer.start()
        >>> timer.lap("Task 1")
        Lap Task 1: 0.0123 seconds
        >>> timer.lap("Task 2")
        Lap Task 2: 0.0234 seconds
        >>> laps = timer.get_lap_times()
        >>> print(f"Recorded {len(laps)} laps")
        Recorded 2 laps
        """
        return self._lap_times.copy()

    def is_running(self) -> bool:
        """
        Check if timer is currently running.

        Returns
        -------
        bool
            True if timer is running, False otherwise.

        Examples
        --------
        >>> timer = pytimer()
        >>> print(timer.is_running())
        False
        >>> timer.start()
        >>> print(timer.is_running())
        True
        >>> timer.stop()
        Elapsed time: 0.0234 seconds
        >>> print(timer.is_running())
        False
        """
        return self._start_time is not None

    def _format_time(self, seconds: float) -> str:
        """
        Format time in human-readable form.

        Parameters
        ----------
        seconds : float
            Time in seconds.

        Returns
        -------
        str
            Formatted time string (e.g., "0.0234 seconds", "234 ms", "234 μs").
        """
        if seconds >= 1.0:
            return f"{seconds:.4f} seconds"
        elif seconds >= 0.001:
            return f"{seconds * 1000:.2f} ms"
        elif seconds >= 0.000001:
            return f"{seconds * 1_000_000:.2f} μs"
        else:
            return f"{seconds * 1_000_000_000:.2f} ns"

    def _output(self, message: str) -> None:
        """
        Output a message via print or logging.

        Parameters
        ----------
        message : str
            Message to output.
        """
        if self.name and not message.startswith(self.name):
            message = f"{self.name}: {message}"

        if self.use_logging:
            logger.log(self.log_level, message)
        else:
            print(message)

    def __enter__(self):
        """Start timing when entering context manager."""
        self.start()
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        """Stop timing when exiting context manager."""
        try:
            self.stop()
        except TimerError:
            # If timer wasn't started, don't raise error on exit
            pass
        return False  # Don't suppress exceptions

    def __repr__(self) -> str:
        """String representation of timer."""
        status = "running" if self.is_running() else "stopped"
        name_str = f" '{self.name}'" if self.name else ""
        return f"<pytimer{name_str} ({status})>"


if __name__ == "__main__":
    # Example usage and demonstrations
    print("pytimer Examples:\n")

    # Example 1: Basic timing
    print("1. Basic timing:")
    timer = pytimer()
    timer.start()
    time.sleep(0.1)
    timer.stop()
    print()

    # Example 2: Named timer
    print("2. Named timer:")
    timer = pytimer(name="Processing")
    timer.start()
    time.sleep(0.05)
    timer.stop()
    print()

    # Example 3: Context manager
    print("3. Context manager:")
    with pytimer(name="ContextTimer"):
        time.sleep(0.08)
    print()

    # Example 4: Lap timing
    print("4. Lap timing:")
    timer = pytimer(name="MultiStep")
    timer.start()
    time.sleep(0.03)
    timer.lap("Step 1")
    time.sleep(0.04)
    timer.lap("Step 2")
    time.sleep(0.02)
    timer.stop()
    print()

    # Example 5: Silent mode with elapsed()
    print("5. Silent mode:")
    timer = pytimer()
    timer.start()
    time.sleep(0.05)
    elapsed = timer.elapsed()
    print(f"  Elapsed so far: {elapsed:.4f} seconds")
    time.sleep(0.03)
    total = timer.stop(silent=True)
    print(f"  Total time: {total:.4f} seconds")
    print()

    # Example 6: Cumulative timing
    print("6. Cumulative timing:")
    timer = pytimer(name="Cumulative")
    for i in range(3):
        timer.start()
        time.sleep(0.02)
        timer.stop(silent=True)
    print(f"  Total cumulative time: {timer.get_total_time():.4f} seconds")
