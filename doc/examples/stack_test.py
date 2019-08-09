import glob
import subprocess

base_cmd = 'python -m doctest {}'
for doctest_file in glob.glob("*.rst"):
    cmd = base_cmd.format(doctest_file)
    process = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    print('command: ', cmd)
    print('========================================================')
    print('output:\n')
    print(output.decode('ascii'))
    print('========================================================')
    print('error:\n')
    print(error)
    print('\n')

