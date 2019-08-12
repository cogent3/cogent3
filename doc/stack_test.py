import glob
import subprocess

base_cmd = 'python -m doctest {}'

relevant_folders = ['cookbook', 'examples']
doctest_files = [rst_file for relevant_folder in relevant_folders for rst_file in glob.glob(relevant_folder + '/*.rst')]


def add_relevant_folder(folder_name):
    relevant_folders.append(folder_name)


for doctest_file in doctest_files:
    cmd = base_cmd.format(doctest_file)
    process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    print('command: ', cmd)
    print('=' * 100)
    print('output:\n')
    print(output.decode('ascii', 'ignore'))
    print('=' * 100)
    print('error:\n')
    print(error)
    print('\n')
