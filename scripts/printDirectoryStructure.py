import os
import sys

def print_dir_hierarchy_exclude(dir, prefix=""):
    excluded_dirs = {
        '.git',
        'nanoAOD_skim.wiki',
        'JHUGenMELA',
        '__pycache__',
        'condor_logs',
        'external/yaml-cpp',
        'data/__pycache__',
        'scripts/__pycache__',
        'Utils/__pycache__',
        'yaml-cpp',
        'Pdfdata',
        'br.sm1',
        'br.sm2',
        'process.DAT',
        'H4LTools_cc'
    }

    try:
        entries = sorted(os.listdir(dir))  # Sort for consistent output
        entries = [entry for entry in entries if entry not in excluded_dirs]
        entry_count = len(entries)

        for index, entry in enumerate(entries):
            item_path = os.path.join(dir, entry)
            connector = "├──" if index < entry_count - 1 else "└──"
            print(f"{prefix}{connector} {entry}")

            if os.path.isdir(item_path):
                extension = "│   " if index < entry_count - 1 else "    "
                print_dir_hierarchy_exclude(item_path, prefix + extension)
    except PermissionError:
        print(f"{prefix}└── [Permission Denied]")

# If no arguments are provided, print the current directory's hierarchy
if len(sys.argv) == 1:
    print_dir_hierarchy_exclude('.')
else:
    for directory in sys.argv[1:]:
        print_dir_hierarchy_exclude(directory)
