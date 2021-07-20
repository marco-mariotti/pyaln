import doctest, os, sys
#unittest, 

thepath = os.path.dirname(os.path.dirname( os.path.abspath( __file__)   )) + '/src'
sys.path.insert(0, thepath)
from pyaln import *

print('='*75)
print(' Testing code examples in docstring of the alignment class '.center(75, '='))
print('='*75)

#doctest.testfile(  thepath+'/alignment.py')
#os.chdir(thepath+'/../docs/')
failure_count1, test_count1 = doctest.testfile('../src/pyaln/alignment.py')
if not failure_count1:
    print('All tests in the alignment class were successful!')

print()
print('='*75)
print(' Testing code examples in the tutorial '.center(75, '='))
print('='*75)


failure_count2, test_count2 = doctest.testfile('../docs/tutorial.rst')
if not failure_count2:
    print('All tests in the tutorial were successful!')

print()

print('=== Summary ==')
print(f'= Number of tests run     = {test_count1+test_count2}')
print(f'= Number of tests failed  = {failure_count1+failure_count2}')
if failure_count1+failure_count2>0:
    sys.exit(1)
    
        
        
