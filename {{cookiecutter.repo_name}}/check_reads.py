import sys
import os

left_reads = sys.argv[1]
right_reads = sys.argv[2]

left_reads = left_reads.split(',')
right_reads = right_reads.split(',')

read_list = zip(left_reads, right_reads)
last_directory = False
checks_passed = True

print 'Checking reads:\n'

for l, r in read_list:
    real_left, real_right = os.path.realpath(l), os.path.realpath(r)
    common_directory = os.path.dirname(os.path.commonprefix([real_left, real_right]))
    rel_left, rel_right = os.path.relpath(real_left, common_directory), os.path.relpath(real_right, common_directory)

    if last_directory != common_directory:
        print 'In directory: %s' % common_directory
    print '%s  <-->  %s' % (rel_left, rel_right)
    if not os.path.exists(real_left):
        print 'Broken symlink detected: %s' % real_left
        checks_passed = False
    if not os.path.exists(real_right):
        print 'Broken symlink detected: %s' % real_right
        checks_passed = False
    if os.path.exists(real_left) and os.path.exists(real_right):
        print 'Symlinks are OK'
    print ""

    last_directory = common_directory

if not checks_passed:
    sys.exit('Checks failed!')

print 'Checks passed'
