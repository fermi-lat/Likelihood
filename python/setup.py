import sys, os

module = sys.argv[1].split('.in')[0]

ifile = open(module + ".in", 'r')
outfile = open(module + ".i", 'w')
for line in ifile:
    if line.find('$(') != -1:
        package = line.split('$(')[1].split(')')[0]
        outfile.write(line.replace('$('+package+')', os.environ[package]))
    else:
        outfile.write(line)
