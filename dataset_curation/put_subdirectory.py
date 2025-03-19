import os
import string
import shutil

# Path to the directory containing the files
source_directory = 'a3ms/'  # Replace with your directory path

# Generate the subdirectory names (aa, ab, ..., zz)
subdirs = [f"{first}{second}" for first in string.ascii_lowercase for second in string.ascii_lowercase]

# Create the subdirectories if they don't already exist
for subdir in subdirs:
    os.makedirs(os.path.join(source_directory, subdir), exist_ok=True)

# List all files in the source directory
files = [f for f in os.listdir(source_directory) if os.path.isfile(os.path.join(source_directory, f))]

# Distribute files into subdirectories
for index, file in enumerate(files):
    subdir_name = subdirs[index % len(subdirs)]
    source_file_path = os.path.join(source_directory, file)
    destination_dir_path = os.path.join(source_directory, subdir_name)
    shutil.move(source_file_path, destination_dir_path)

print("Files have been distributed into subdirectories.")

