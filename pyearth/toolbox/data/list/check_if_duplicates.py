"""
Duplicate Detection Utilities

This module provides functions for detecting duplicate elements in sequences.
The primary function checks whether a list or sequence contains any duplicate
values, which is useful for data validation, uniqueness verification, and
quality control.

Main Functions
--------------
check_if_duplicates : Check if a sequence contains duplicate elements

Key Features
------------
- Efficient duplicate detection using set comparison
- Support for any hashable data types
- Clear boolean return semantics
- Type hints for better code clarity
- Comprehensive validation and error handling

Use Cases
---------
1. **Data Validation**: Verify uniqueness of IDs or keys
2. **Quality Control**: Check for duplicate entries in datasets
3. **Database Operations**: Validate unique constraints before insertion
4. **File Processing**: Detect duplicate filenames or paths
5. **Configuration Validation**: Ensure unique parameter names
6. **Statistical Analysis**: Verify sample independence
7. **Testing**: Validate test case uniqueness

Technical Details
-----------------
The function uses Python's set data structure for O(n) time complexity
duplicate detection, which is significantly faster than the naive O(n²)
approach of checking each element against all others.

Performance:
- Time Complexity: O(n) where n is the sequence length
- Space Complexity: O(n) for the set conversion
- Optimized for large datasets

Dependencies
------------
- typing: Type hints for better code documentation
- logging: Error reporting and debugging

See Also
--------
- set: Python's built-in set data structure
- collections.Counter: Count element occurrences
"""

import logging
from typing import Sequence, Union

# Configure logging
logger = logging.getLogger(__name__)


def check_if_duplicates(aList_in: Sequence) -> int:
    """
    Check if a sequence contains any duplicate elements.

    This function efficiently determines whether a list, tuple, or other
    sequence contains duplicate values by comparing the sequence length
    with its set representation. Returns 1 if all elements are unique,
    0 if duplicates exist.

    Parameters
    ----------
    aList_in : sequence
        Input sequence to check for duplicates.
        Must contain hashable elements.

        Supported types:
        - list: [1, 2, 3, 4]
        - tuple: (1, 2, 3, 4)
        - str: "abcd" (checks character duplicates)
        - range: range(10)
        - Any sequence with hashable elements

        Elements must be hashable:
        - Supported: int, float, str, tuple, frozenset
        - Supported: Custom objects with __hash__ and __eq__ implemented
        - Not supported: list, dict, set (unhashable types)
        - Not supported: Custom objects without __hash__ method

        **Custom Object Requirements**:
        For custom objects to work with this function, they must:
        1. Implement __hash__() method (returns integer hash)
        2. Implement __eq__() method (defines equality comparison)
        3. Ensure hash consistency: if a == b, then hash(a) == hash(b)

        Example custom object:
        ```python
        class Person:
            def __init__(self, name, id):
                self.name = name
                self.id = id

            def __hash__(self):
                return hash(self.id)  # Hash based on ID

            def __eq__(self, other):
                return isinstance(other, Person) and self.id == other.id
        ```

    Returns
    -------
    int
        Duplicate detection result:

        - 1: No duplicates found (all elements unique)
        - 0: Duplicates found (at least one repeated element)
        - 0: Error occurred (invalid input)

    Raises
    ------
    TypeError
        If input contains unhashable elements (implicit from set conversion).
        Function catches this and returns 0 with error logging.

    Notes
    -----
    1. **Efficiency**: O(n) time complexity using set comparison
    2. **Empty Sequences**: Returns 1 (no duplicates in empty sequence)
    3. **Single Element**: Returns 1 (cannot have duplicates)
    4. **Hashability**: Only works with hashable elements
    5. **Order Preservation**: Does not maintain order (uses set)
    6. **None Values**: None is hashable and handled correctly
    7. **Mixed Types**: Can handle mixed hashable types
    8. **Case Sensitive**: String comparison is case-sensitive
    9. **Numeric Types**: 1 and 1.0 considered same (hash equal)
    10. **Return Type**: Returns integer (not boolean) for legacy compatibility

    Algorithm
    ---------
    1. Convert sequence to set (removes duplicates)
    2. Compare set length with original sequence length
    3. If lengths differ, duplicates exist
    4. If lengths match, all elements unique

    Examples
    --------
    Check list with no duplicates:

    >>> result = check_if_duplicates([1, 2, 3, 4, 5])
    >>> print(result)
    1

    Check list with duplicates:

    >>> result = check_if_duplicates([1, 2, 3, 2, 5])
    >>> print(result)
    0

    Check string for duplicate characters:

    >>> result = check_if_duplicates("hello")
    >>> print(result)  # 'l' appears twice
    0

    Check string with unique characters:

    >>> result = check_if_duplicates("world")
    >>> print(result)
    1

    Empty sequence (no duplicates):

    >>> result = check_if_duplicates([])
    >>> print(result)
    1

    Single element (no duplicates):

    >>> result = check_if_duplicates([42])
    >>> print(result)
    1

    Check tuple with mixed types:

    >>> result = check_if_duplicates((1, 'a', 2.5, 'b', 1))
    >>> print(result)  # 1 appears twice
    0

    Check for duplicate IDs in dataset:

    >>> user_ids = [101, 102, 103, 104, 105]
    >>> if check_if_duplicates(user_ids):
    ...     print("All IDs are unique")
    ... else:
    ...     print("Duplicate IDs found!")
    All IDs are unique

    Validate unique filenames:

    >>> filenames = ['data.csv', 'results.txt', 'data.csv']
    >>> is_unique = check_if_duplicates(filenames)
    >>> print(f"Files unique: {bool(is_unique)}")
    Files unique: False

    Handle None values:

    >>> data = [1, None, 2, None, 3]
    >>> result = check_if_duplicates(data)
    >>> print(result)  # None appears twice
    0

    Check range object:

    >>> result = check_if_duplicates(range(100))
    >>> print(result)
    1

    Numeric type equality (1 == 1.0):

    >>> result = check_if_duplicates([1, 2, 1.0])
    >>> print(result)  # 1 and 1.0 are considered same
    0

    See Also
    --------
    set : Python's set data structure for unique elements
    collections.Counter : Count occurrences of each element
    pandas.Series.duplicated : Pandas duplicate detection

    References
    ----------
    .. [1] Python Data Structures - Sets
           https://docs.python.org/3/tutorial/datastructures.html#sets
    .. [2] Time Complexity of Python Operations
           https://wiki.python.org/moin/TimeComplexity

    Warnings
    --------
    - Input must contain only hashable elements
    - Lists, dicts, sets cannot be elements (unhashable)
    - Large sequences may consume significant memory for set conversion
    - Return is int (1/0) not bool (True/False) for backward compatibility
    """
    # Validate input
    if aList_in is None:
        error_msg = "Input sequence cannot be None"
        logger.error(error_msg)
        print(f"Error: {error_msg}")
        return 0

    try:
        # Convert to list if needed to get length
        # This handles generators and other iterables
        if not hasattr(aList_in, "__len__"):
            aList_in = list(aList_in)
    except Exception as e:
        error_msg = f"Failed to process input sequence: {str(e)}"
        logger.error(error_msg)
        print(f"Error: {error_msg}")
        return 0

    # Check if empty - no duplicates in empty sequence
    if len(aList_in) == 0:
        logger.debug("Empty sequence - no duplicates")
        return 1

    # Check if single element - cannot have duplicates
    if len(aList_in) == 1:
        logger.debug("Single element - no duplicates")
        return 1

    # Efficient duplicate check using set
    try:
        # If set length differs from list length, there are duplicates
        has_duplicates = len(set(aList_in)) != len(aList_in)

        if has_duplicates:
            logger.debug(
                f"Duplicates found: {len(aList_in)} elements, {len(set(aList_in))} unique"
            )
            return 0
        else:
            logger.debug(f"No duplicates: all {len(aList_in)} elements unique")
            return 1

    except TypeError as e:
        # Occurs when sequence contains unhashable elements
        error_msg = f"Input contains unhashable elements: {str(e)}"
        logger.error(error_msg)
        print(f"Error: {error_msg}")
        print(
            f"Hint: Elements must be hashable (int, str, tuple, etc.), not list or dict"
        )
        return 0
    except Exception as e:
        # Catch any other unexpected errors
        error_msg = f"Unexpected error checking duplicates: {str(e)}"
        logger.error(error_msg)
        print(f"Error: {error_msg}")
        return 0
