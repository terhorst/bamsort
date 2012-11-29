from collections import namedtuple
import cStringIO as StringIO
import logging
import re

from util import FileStruct
from data_structs import *

ReadHeader = namedtuple("ReadHeader", "block_size refID pos bin_mq_nl "
                                      "flag_nc l_seq next_refID next_pos tlen")
CompressedRead = namedtuple("CompressedRead",  
                            "rid pos worker_num ptr cs bs crc32 rlen")

class RawBAM:
    'An uncompressed BAM file'
    read_hdr = FileStruct('3i2I4i')
    min_read_hdr = FileStruct('3i')

    def __init__(self, f, header):
        self.f = f
        if header:
            self._process_header()
        else:
            self._find_rec_start()

    def _find_rec_start(self):
        rh = self.read_hdr
        s = self.f.read(rh.size)
        hdr = ReadHeader(*rh.unpack(s))
        while not(0 < hdr.block_size < 1000) or not(0 < hdr.l_seq < 200):
            s += self.f.read(1)
            hdr = ReadHeader(*rh.unpack(s[-rh.size:]))
        self.f.seek(len(s)-rh.size)

    def _process_header(self):
        s = FileStruct('4sI')
        f = self.f
        rh = StringIO.StringIO()
        magic,l_text = s.unpack(f)
        assert magic=="BAM\1"
        s = FileStruct('%isI' % l_text)
        hdr, n_ref = s.unpack(f)
        refs = []
        for i in range(n_ref):
            l_name, = int32_t.unpack(f)
            s = FileStruct("%iscI" % (l_name - 1))
            name,_,length = s.unpack(f)
            refs.append((name,length))
        self.header = (hdr,n_ref,refs)
        pos = f.tell()
        f.seek(0)
        self.rawheader = f.read(pos)
        # Edit the header to reflect sorting
        fl = self.rawheader[12:].split('\n')[0]
        fd = dict([t.split(':') for t in fl.split('\t')])
        fd['SO'] = 'coordinate'
        self.rawheader = self.rawheader[:12] + \
                         "\t".join(["%s:%s" % t for t in fd.items()]) + \
                         self.rawheader[(12+len(fl)):]

    def __iter__(self):
        while True:
            ret = StringIO.StringIO()
            ptr = self.f.tell()
            bs = self.f.read(int32_t.size)
            if len(bs) != int32_t.size:
                raise StopIteration
            ret.write(bs)
            bs, = int32_t.unpack(bs)
            dat = self.f.read(bs)
            if len(dat) != bs:
                raise StopIteration
            ret.write(dat)
            yield ptr,ret.getvalue()
