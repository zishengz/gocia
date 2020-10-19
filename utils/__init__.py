import numpy as np


def conv_str2list(myString):
    tmp = myString
    if myString[0] == '[':
        tmp = tmp[1:-1]
    if ' ,' in tmp:
        delim = ','
    else:
        delim = ' '
    return [eval(i) for i in tmp.split(delim)]
