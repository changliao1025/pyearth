import os
import errno
def create_symlink(source, target_link):
    """
    """
    try:
        os.symlink(source, target_link)
    except OSError as e:
        if e.errno == errno.EEXIST:
            pass
            os.remove(target_link)
            os.symlink(source, target_link)
        else:
            raise e