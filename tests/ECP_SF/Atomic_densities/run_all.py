import os
import subprocess

# Specify the directory
directory = os.path.dirname(os.path.realpath(__file__))

# Loop over all files in the directory
for filename in os.listdir(directory):
    # Check if the file ends with .inp
    if filename.endswith('.inp'):
        # Construct the full file path
        file_path = os.path.join(directory, filename)
        # Execute the system call
        subprocess.run(['C:\\ORCA5\\orca.exe', file_path])