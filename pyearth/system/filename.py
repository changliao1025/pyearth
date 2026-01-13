import os


def get_extension_from_filename(filename: str) -> str:
    """
    Extract the file extension from a given filename.

    Args:
        filename (str): The name of the file.

    Returns:
        str: The file extension (including the dot), or an empty string if no extension is found.
    """
    _, ext = os.path.splitext(filename)
    return ext


def get_extension_from_path(filepath: str) -> str:
    """
    Extract the file extension from a given filename.

    Args:
        filename (str): The name of the file.

    Returns:
        str: The file extension (including the dot), or an empty string if no extension is found.
    """
    basename = os.path.basename(filepath)
    _, ext = os.path.splitext(basename)
    return ext


def get_filename_from_path_without_extension(filepath: str) -> str:
    """
    Get the filename without its extension.

    Args:
        filename (str): The name of the file.

    Returns:
        str: The filename without the extension.
    """
    basename = os.path.basename(filepath)
    name, _ = os.path.splitext(basename)
    return name


def get_filename_without_extension(filename: str) -> str:
    """
    Get the filename without its extension.

    Args:
        filename (str): The name of the file.

    Returns:
        str: The filename without the extension.
    """

    name, _ = os.path.splitext(filename)
    return name


def get_folder_path(filepath: str) -> str:
    """
    Get the folder path from a full file path.

    Args:
        filepath (str): The full path to the file.

    Returns:
        str: The folder path containing the file.
    """
    return os.path.dirname(filepath)
