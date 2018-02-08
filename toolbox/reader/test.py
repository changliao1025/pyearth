import sys
print(sys.path)
sys.path.append('/Users/liao313/workspace/python/library/')

print(sys.path)
from eslib_python.toolbox.reader.text import text_reader_string
text_reader_string('/Users/liao313/workspace/python/library/eslib_python/toolbox/reader/data.txt')