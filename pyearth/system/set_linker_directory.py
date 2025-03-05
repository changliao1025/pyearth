import os
from pyearth.system.python.get_python_environment import retrieve_python_environment
sConda_env_path, sConda_env_name = retrieve_python_environment()
os.environ['LD_LIBRARY_PATH'] = f"{sConda_env_path}/lib:{sConda_env_path}/lib64:{os.environ.get('LD_LIBRARY_PATH', '')}"