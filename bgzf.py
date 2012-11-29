from struct import Struct
from collections import namedtuple
import os
import logging

from data_structs import *
from util import FileStruct

BGZFBlock = namedtuple("BGZFBlock", "offset size_u")

class BGZFReader:
    BGZFHdr = namedtuple("BGZFHeader", 
                         "ID1 ID2 CM FLG MTIME XFL OS XLEN SI1 SI2 SLEN BSIZE")
    bgzf_s = FileStruct('4BI2BH2B2H')
    bgzf_ftr_s = FileStruct('2I')

    def __init__(self, bfile):
        self.f = open(bfile,'r')
        self._find_bgzf_blocks()

    def _find_bgzf_blocks(self):
        f = self.f
        seek = 0
        bl = self.blocks = []
        fsz = os.fstat(f.fileno()).st_size
        i = 0
        self.uncompressed_size = 0
        while seek < fsz:
            i += 1
            offset = seek
            f.seek(seek)
            h = self.BGZFHdr(*self.bgzf_s.unpack(f))
            seek += h.BSIZE + 1
            f.seek(seek - uint32_t.size)
            size_u = uint32_t.unpack(f)[0]
            bl = BGZFBlock(offset=offset, size_u=size_u)
            self.blocks.append(bl)
            self.uncompressed_size += bl.size_u

        logging.debug("Uncompressed BAM size: %imb" % 
            (self.uncompressed_size // (1024**2)))

