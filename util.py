import operator
from struct import Struct

class FileStruct(Struct):
    def __init__(self, fmt):
        Struct.__init__(self, "<" + fmt)

    def unpack(self, s, buf=None):
        if hasattr(s, 'read'):
            s = s.read(self.size)
        if buf is not None:
            buf.write(s)
        return Struct.unpack(self, s)

def chunk(lst, n):
    ret = []
    items_per = len(lst) // n
    for i in range(n):
        k = i * items_per
        if i < n - 1:
            ret.append(lst[k:(k+items_per)])
        else:
            ret.append(lst[k:])
    return ret

def accumulate(iterable, func=operator.add):
    'Return running totals'
    # accumulate([1,2,3,4,5]) --> 1 3 6 10 15
    # accumulate([1,2,3,4,5], operator.mul) --> 1 2 6 24 120
    it = iter(iterable)
    total = next(it)
    yield total
    for element in it:
        total = func(total, element)
        yield total
