import sys, os

def replaceEnvVars(target, source, env=None):
    ifile = open(str(source[0]), 'r')
    outfile = open(str(target[0]), 'w')
    for line in ifile:
        if line.find('$(') != -1:
            package = line.split('$(')[1].split(')')[0]
            outfile.write(line.replace('$('+package+')', os.environ[package]))
        else:
            outfile.write(line)

if __name__ == '__main__':
    if len(sys.argv) == 2:
        outfile = sys.argv[1].replace('.in', '.i')
        replaceEnvVars((outfile, ), (sys.argv[1], ))
    elif len(sys.argv) == 3:
        replaceEnvVars((sys.argv[2],), (sys.argv[1],))
