import pickle
import sys

if len(sys.argv) == 2:
    with open(sys.argv[1], 'rb') as f:
        a = pickle.load(f)
        print(a)
else:
    print('invalid number of arguments')
