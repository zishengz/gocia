
import os

def read_frag(fn='fragments'):
    if fn in os.listdir():
        return eval(open(fn, 'r').read())
    else:
        print('Could not find the \"fragments\" file')
        return None

def write_frag(fraglist, fn='fragments'):
    if type(fraglist) is list:
        open(fn, 'w').write(str(fraglist))
    elif type(fraglist) is str:
        open(fn, 'w').write(fraglist)
    else:
        print('fraglist should be either str or list!')

def update_frag_del(list_del, fn='fragments'):
    if len(list_del) > 0:
        fraglist_tmp = read_frag(fn)
        fraglist_new = []
        for f in fraglist_tmp:
            if not set(f) & set(list_del):
                my_f = []
                for i in f:
                    my_f.append(i - len([ii for ii in list_del if ii < i]))
                fraglist_new.append(my_f)
        print(f'New fraglist is: {fraglist_new}')
        write_frag(fraglist_new, fn=fn)