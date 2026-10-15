"""Unit tests for pyearth.system package."""

import io
import os
import sys
import tempfile
import unittest
from pathlib import Path

from pyearth.system.create_symlink import create_symlink
from pyearth.system.define_global_variables import (
    cms2cmd,
    earth_radius,
    feet2meter,
    inch2mm,
    missing_value,
    mms2mmd,
    nmonth,
    print_environment_info,
    sMachine,
    sPlatform_os,
)
from pyearth.system.filename import (
    get_extension_from_filename,
    get_extension_from_path,
    get_filename_from_path_without_extension,
    get_filename_without_extension,
    get_folder_path,
)
from pyearth.system.python.get_python_environment import (
    get_conda_environment,
    get_python_environment,
)
from pyearth.system.set_linker_directory import configure_linker_path


class TestFilenameUtils(unittest.TestCase):
    def test_get_extension_from_filename(self):
        self.assertEqual(get_extension_from_filename("data.txt"), ".txt")
        self.assertEqual(get_extension_from_filename("archive.tar.gz"), ".gz")
        self.assertEqual(get_extension_from_filename("no_ext"), "")
        self.assertEqual(get_extension_from_filename(Path("nested/dir/image.png")), ".png")

        with self.assertRaises(ValueError):
            get_extension_from_filename(None)

    def test_get_extension_from_path(self):
        self.assertEqual(get_extension_from_path("/path/to/file.tif"), ".tif")
        self.assertEqual(get_extension_from_path(Path("/path/to/file.shp")), ".shp")
        self.assertEqual(get_extension_from_path("/path/to/noext"), "")

        with self.assertRaises(ValueError):
            get_extension_from_path(None)

    def test_get_filename_from_path_without_extension(self):
        self.assertEqual(get_filename_from_path_without_extension("/path/to/file.tif"), "file")
        self.assertEqual(get_filename_from_path_without_extension("alone.json"), "alone")
        self.assertEqual(get_filename_from_path_without_extension(Path("/dir/sub/test.nc")), "test")

        with self.assertRaises(ValueError):
            get_filename_from_path_without_extension(None)

    def test_get_filename_without_extension(self):
        self.assertEqual(get_filename_without_extension("/path/to/file.tif"), "/path/to/file")
        self.assertEqual(get_filename_without_extension("alone.json"), "alone")
        self.assertEqual(get_filename_without_extension("noext"), "noext")

        with self.assertRaises(ValueError):
            get_filename_without_extension(None)

    def test_get_folder_path(self):
        self.assertEqual(get_folder_path("/path/to/file.tif"), "/path/to")
        self.assertEqual(get_folder_path(Path("/path/to/sub/file.tif")), "/path/to/sub")

        with self.assertRaises(ValueError):
            get_folder_path(None)


class TestCreateSymlink(unittest.TestCase):
    def test_create_symlink_success(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            src_file = Path(tmpdir) / "source.txt"
            src_file.write_text("hello symlink")
            link_file = Path(tmpdir) / "link.txt"

            create_symlink(src_file, link_file)
            self.assertTrue(link_file.is_symlink())
            self.assertEqual(link_file.read_text(), "hello symlink")

            # Replace existing symlink
            src2_file = Path(tmpdir) / "source2.txt"
            src2_file.write_text("hello updated")
            create_symlink(src2_file, link_file)
            self.assertTrue(link_file.is_symlink())
            self.assertEqual(link_file.read_text(), "hello updated")

    def test_create_symlink_validation(self):
        with self.assertRaises(ValueError):
            create_symlink("", "link")
        with self.assertRaises(ValueError):
            create_symlink("source", "")
        with self.assertRaises(ValueError):
            create_symlink(None, "link")

    def test_create_symlink_missing_parent(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            non_existent_dir = Path(tmpdir) / "non_existent" / "link.txt"
            with self.assertRaises(FileNotFoundError):
                create_symlink("source.txt", non_existent_dir)

    def test_create_symlink_target_is_regular_file(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            target = Path(tmpdir) / "regular_file.txt"
            target.write_text("cannot overwrite regular file")
            with self.assertRaises(FileExistsError):
                create_symlink("source.txt", target)


class TestSetLinkerDirectory(unittest.TestCase):
    def test_configure_linker_path(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            lib_dir = Path(tmpdir) / "lib"
            lib_dir.mkdir()

            var_name = "TEST_LINKER_VAR"
            result = configure_linker_path(env_path=tmpdir, prepend=True, var_name=var_name)
            self.assertIn(str(lib_dir.resolve()), result)
            self.assertEqual(os.environ.get(var_name), result)

            # Test append
            extra_dir = Path(tmpdir) / "lib64"
            extra_dir.mkdir()
            result_appended = configure_linker_path(
                env_path=tmpdir, prepend=False, var_name=var_name
            )
            self.assertIn(str(extra_dir.resolve()), result_appended)

    def test_configure_linker_path_invalid(self):
        with self.assertRaises(ValueError):
            configure_linker_path(env_path="/non/existent/path/for/sure")


class TestPythonEnvironment(unittest.TestCase):
    def test_get_python_environment(self):
        env_path, env_name, env_type = get_python_environment()
        self.assertIsInstance(env_path, str)
        self.assertTrue(os.path.isdir(env_path))
        self.assertIsInstance(env_name, str)
        self.assertIn(env_type, ("conda", "venv", "virtualenv", "pyenv", "system"))

    def test_get_conda_environment(self):
        _, _, env_type = get_python_environment()
        if env_type == "conda":
            path, name = get_conda_environment()
            self.assertIsInstance(path, str)
            self.assertIsInstance(name, str)
        else:
            with self.assertRaises(ValueError):
                get_conda_environment()


class TestDefineGlobalVariables(unittest.TestCase):
    def test_global_constants(self):
        self.assertEqual(earth_radius, 6378137.0)
        self.assertEqual(missing_value, -9999.0)
        self.assertEqual(nmonth, 12)
        self.assertEqual(mms2mmd, 86400.0)
        self.assertEqual(feet2meter, 0.3048)
        self.assertEqual(inch2mm, 25.4)
        self.assertEqual(cms2cmd, 86400.0)
        self.assertIn(sPlatform_os, ("Linux", "Darwin", "Windows"))
        self.assertIn(sMachine, ("linux", "mac", "windows", "unix"))

    def test_print_environment_info(self):
        captured_output = io.StringIO()
        old_stdout = sys.stdout
        try:
            sys.stdout = captured_output
            print_environment_info(verbose=False)
        finally:
            sys.stdout = old_stdout
        self.assertIn("PyEarth Environment Configuration", captured_output.getvalue())


if __name__ == "__main__":
    unittest.main()
