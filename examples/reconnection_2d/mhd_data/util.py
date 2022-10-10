"""
Utility functions
"""
import errno
import os

def get_variable_value(variable_name, current_line, content):
    """
    Get the value of one variable from the content of one file

    Args:
        variable_name: the variable name.
        current_line: current line number.
        content: the content of the information file.
    Returns:
        variable_value: the value of the variable.
        line_number: current line number after the operations.
    """
    line_number = current_line
    nline = len(content)
    while not variable_name in content[line_number]:
        line_number += 1
        if line_number >= nline:
            return (0.0, current_line)
    single_line = content[line_number]
    line_splits = single_line.split("=")
    line_splits1 = line_splits[1].split("#")
    variable_value = float(line_splits1[0])
    return (variable_value, line_number)


def mkdir_p(path):
    """Create directory recursively
    """
    try:
        os.makedirs(path)
    except OSError as exc:  # Python >2.5
        if exc.errno == errno.EEXIST and os.path.isdir(path):
            pass
        else:
            raise


if __name__ == "__main__":
    pass
