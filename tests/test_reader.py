import os
import numpy as np
import pytest

from pyearth.toolbox.reader import (
    line_count,
    parse_xml_file,
    read_configuration_file,
    split_string_into_chunk,
    text_reader_string,
)


class TestLineCount:
    def test_line_count(self, tmp_path):
        f = tmp_path / "test.txt"
        f.write_text("line 1\nline 2\nline 3\n")
        assert line_count(str(f)) == 3

    def test_empty_file(self, tmp_path):
        f = tmp_path / "empty.txt"
        f.write_text("")
        assert line_count(str(f)) == 0

    def test_nonexistent_file(self):
        with pytest.raises(FileNotFoundError):
            line_count("non_existent_file.txt")


class TestSplitStringIntoChunk:
    def test_split_string(self):
        s = "abcdefghij12345"
        chunks = list(split_string_into_chunk(s, 5))
        assert chunks == ["abcde", "fghij", "12345"]

    def test_default_chunk_size(self):
        s = "a" * 25
        chunks = list(split_string_into_chunk(s))
        assert len(chunks) == 3
        assert len(chunks[0]) == 10
        assert len(chunks[1]) == 10
        assert len(chunks[2]) == 5

    def test_invalid_input(self):
        with pytest.raises(TypeError):
            list(split_string_into_chunk(123))
        with pytest.raises(ValueError):
            list(split_string_into_chunk("abc", 0))


class TestReadConfigurationFile:
    def test_read_config(self, tmp_path):
        f = tmp_path / "config.txt"
        f.write_text(
            "# Comment line\n"
            "key1, value1\n"
            "key2 , 12345\n"
            "\n"
            "key3, complex,value\n"
        )
        cfg = read_configuration_file(str(f))
        assert cfg["key1"] == "value1"
        assert cfg["key2"] == "12345"
        assert cfg["key3"] == "complex,value"

    def test_nonexistent(self):
        with pytest.raises(FileNotFoundError):
            read_configuration_file("nonexistent.conf")


class TestParseXmlFile:
    def test_parse_xml(self, tmp_path):
        xml_content = """<?xml version="1.0"?>
<root>
    <group>
        <entry id="CASE" value="my_case">
            <type>char</type>
        </entry>
        <entry id="RUNDIR" value="/scratch/$CASE/run">
            <type>char</type>
        </entry>
        <entry id="NINST" value="4">
            <type>integer</type>
        </entry>
        <entry id="DT" value="1800.0">
            <type>real</type>
        </entry>
        <entry id="FLAG" value="TRUE">
            <type>logical</type>
        </entry>
        <entry id="LEVELS">
            <type>real</type>
            <values>10.0</values>
            <values>20.0</values>
            <values>30.0</values>
        </entry>
    </group>
</root>
"""
        f = tmp_path / "namelist.xml"
        f.write_text(xml_content)
        result = parse_xml_file(str(f))

        assert result["CASE"] == "my_case"
        assert result["RUNDIR"] == "/scratch/my_case/run"
        assert result["NINST"] == 4
        assert result["DT"] == 1800.0
        assert result["FLAG"] is True
        assert isinstance(result["LEVELS"], np.ndarray)
        assert np.allclose(result["LEVELS"], [10.0, 20.0, 30.0])


class TestTextReaderString:
    def test_whitespace_delimited(self, tmp_path):
        f = tmp_path / "data.txt"
        f.write_text("A 10 20.5\nB 30 40.5\n")
        arr = text_reader_string(str(f))
        assert arr.shape == (2, 3)
        assert arr[0, 0] == "A"
        assert arr[0, 1] == "10"
        assert arr[1, 2] == "40.5"

    def test_csv_with_quotes(self, tmp_path):
        f = tmp_path / "data.csv"
        f.write_text('"A",10,"20"\n"B",30,"40"\n')
        arr = text_reader_string(str(f), cDelimiter_in=",", iFlag_remove_quota=True)
        assert arr.shape == (2, 3)
        assert arr[0, 0] == "A"
        assert arr[0, 2] == "20"
        assert arr[1, 1] == "30"

    def test_skip_header_and_limit_rows(self, tmp_path):
        f = tmp_path / "data.txt"
        f.write_text("header line 1\nheader line 2\n1 2\n3 4\n5 6\n")
        arr = text_reader_string(str(f), nrow_in=2, iSkipline_in=2)
        assert arr.shape == (2, 2)
        assert arr[0, 0] == "1"
        assert arr[1, 0] == "3"
