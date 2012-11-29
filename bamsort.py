#!/usr/bin/env python2.7
import gzip, zlib
import sys, operator, os
from cStringIO import StringIO
from collections import namedtuple
from multiprocessing import Pool, current_process, RawArray, Process
from multiprocessing.queues import Queue
import ctypes
import logging
import os

from bgzf import BGZFReader
from raw_bam import RawBAM
from util import *
from data_structs import *

N_WORKERS = 1
READS_PER_BLOCK = 100
MAX_BLOCK_SIZE = 65536 - BGZFReader.bgzf_s.size - BGZFReader.bgzf_ftr_s.size

CompressedRead = namedtuple("CompressedRead",  "rid pos worker_num ptr bs")

def sort_read_ary(i, bam, blocks, offset, buf, q):
    'Process each read in ary, then store compressed read location and length in bry'
    logging.debug("Sorting %s %i" % (blocks, i))
    with open(bam, 'rb') as f:
        f.seek(blocks[0].offset)
        buf.raw = gzip.GzipFile(fileobj=f).read(len(buf))
    rb = RawBAM(StringIO(buf), i == 0)
    ptr = 0
    mrh = RawBAM.min_read_hdr
    lst = []
    for ptr,read in rb:
        cr = {}
        try:
            cr['bs'],cr['rid'],cr['pos'] = mrh.unpack(read[:mrh.size])
        except:
            raise Error(read)
        cr['ptr'] = ptr
        cr['worker_num'] = i
        lst.append(CompressedRead(**cr))
    # lst.sort()
    logging.debug("Finished sort; queuing list of length %i" % len(lst))
    q.put(lst)
    logging.debug("Finished queue")
    return True

def parallel_sort(bam, out):
    lb = BGZFReader(bam)
    mem = lb.uncompressed_size
    buf = RawArray(ctypes.c_char, mem)
    q = Queue()
    procs = []

    block_allocs = chunk(lb.blocks, N_WORKERS)
    offsets = [0] + list(accumulate(sum(b.offset for b in blocks) 
                                    for blocks in block_allocs))[:-1]
    ary_szs = [sum([b.size_u for b in blocks]) for blocks in block_allocs]
    bufs = [RawArray(ctypes.c_char,mem) for mem in ary_szs]
    z = zip(chunk(lb.blocks, N_WORKERS), offsets, bufs)
    for i,(blocks,off,buf) in enumerate(z):
        args = (i, bam, blocks, off, buf, q)
        # sort_read_ary(*args)
        p = Process(target=sort_read_ary, args=args)
        procs.append(p)
        p.start()

    combined = []
    for _ in procs:
        combined += q.get(True)
    logging.debug("Starting combined sort on %i reads" % len(combined))
    # combined.sort()
    logging.debug("Finished combined sort")

    for p in procs:
        p.join()
        logging.debug("Returned from " + str(p))

    hdr = RawBAM(gzip.GzipFile(bam),header=True).rawheader
    with open(out, 'wb') as f:
        write_bgzf_block(f, hdr)
        for creads in grouper(READS_PER_BLOCK, combined):
            data = ""
            for i,cr in enumerate(creads):
                data += bufs[cr.worker_num][cr.ptr:(cr.ptr+cr.bs+4)] 
            write_bgzf_block(f, data)
        write_bam_eof(f)
 
def deflate(data):
    return zlib.compress(data)[2:-4]

def write_bgzf_block(f, data):
    bgzf_s = BGZFReader.bgzf_s
    bgzf_ftr_s = BGZFReader.bgzf_ftr_s
    for d in ["".join(s) for s in grouper(MAX_BLOCK_SIZE, data, "")]:
        comp = deflate(d)
        bsize = bgzf_s.size + len(comp) + bgzf_ftr_s.size
        hdr = (31,139,8,4,0,0,255,6,66,67,2,bsize-1)
        f.write(bgzf_s.pack(*hdr))
        f.write(comp)
        crc32 = zlib.crc32(d)
        isize = len(d)
        f.write(int32_t.pack(crc32))
        f.write(uint32_t.pack(isize))

def write_bam_eof(f):
    # (undocumented?) EOF block required for samtools to be happy
    f.write("\037\213\010\4\0\0\0\0\0\377\6\0\102\103\2\0\033\0\3\0\0\0\0\0\0\0\0\0")

if __name__=="__main__":
    logging.basicConfig(level=logging.DEBUG, format="%(asctime)s;%(levelname)s;%(message)s")
    parallel_sort('test.bam', 'out.bam')
