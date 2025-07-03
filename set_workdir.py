import os
import sys

def set_workdir(py_file=None): # py_path would be getData.py, spectraPlay.py, or PPSD.py
    """Sets the working directory depending on where the Python script or Jupyter notebook is located."""
    # alternative to this function, put 'os.chdir(os.path.dirname(os.path.abspath(getData.__file__)))' before file calls 
    
    try:
        # if running as a script, get py script's directory
        script_dir = os.path.dirname(os.path.abspath(__file__))
    except NameError:
        # if running in a jupyter notebook, get current working directory
        script_dir = os.getcwd()

    # if a specific script name is provided, adjust the path
    if py_file:
        py_path = os.path.join(script_dir, py_file) 
    else:
        py_path = script_dir
    
    # change the working directory 
    os.chdir(os.path.dirname(os.path.abspath(py_path)))

    # if python files and notebook are in different locations, update sys.path
    if script_dir != os.getcwd():
        sys.path.append(script_dir)   
    #print(f"Current working directory set to: {os.getcwd()}")

# call by
#set_workdir('filename.py')

