import sys
import os
# Get the path to the Python interpreter
def get_python_environment():
    slash = os.sep
    sPath_python = sys.executable
    components = sPath_python.split(slash)
    # Find the index of 'envs' in the components
    envs_index = components.index('envs')
    # The environment name is the component right after 'envs'
    sConda_env_name = components[envs_index + 1]

    print(f"Conda environment name: {sConda_env_name}")

    #also retrieve the path to the conda environment
    sConda_env_path = slash.join(components[:envs_index + 2])

    return sConda_env_path , sConda_env_name