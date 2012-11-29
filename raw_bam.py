from collections import namedtuple
import cStringIO as StringIO

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
            self.find_rec_start()

    def find_rec_start(self):
        rh = self.read_hdr
        s = self.f.read(rh.size)
        hdr = ReadHeader(*rh.unpack(s))
        while not(0 < hdr.block_size < 1000) or not(0 < hdr.l_seq < 200):
            s += self.f.read(1)
            hdr = self.ReadHeader(*rh.unpack(s[-rh.size:]))
        self.f.seek(len(s)-rh.size)

    def _process_header(self):
        s = FileStruct('4sI')
        f = self.f
        rh = StringIO.StringIO()
        magic,l_text = s.unpack(f,rh)
        assert magic=="BAM\1"
        s = FileStruct('%isI' % l_text)
        hdr, n_ref = s.unpack(f,rh)
        refs = []
        for i in range(n_ref):
            l_name, = int32_t.unpack(f,rh)
            s = FileStruct("%iscI" % (l_name - 1))
            name,_,length = s.unpack(f,rh)
            refs.append((name,length))
        self.header = (hdr,n_ref,refs)
        self.rawheader = rh

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

