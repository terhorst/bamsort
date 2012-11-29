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

from bgzf_reader import BGZFReader
from raw_bam import RawBAM
from util import *

N_WORKERS = 1
READS_PER_BLOCK = 2
logging.basicConfig(level=logging.DEBUG, format="%(asctime)s;%(levelname)s;%(message)s")

CompressedRead = namedtuple("CompressedRead",  "rid pos worker_num ptr bs")

def sort_read_ary(i, bam, blocks, offset, buf, q):
    'Process each read in ary, then store compressed read location and length in bry'
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
            logging.debug(read)
            raise
        cr['ptr'] = ptr
        cr['worker_num'] = i
        lst.append(CompressedRead(**cr))
    logging.debug("Wrote %i bytes" % ptr)
    lst.sort()
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
        sort_read_ary(*args)
        # p = Process(target=sort_read_ary, args=args)
        # procs.append(p)
        # p.start()

    combined = []
    for _ in procs:
        combined += q.get(True)
    logging.debug("Starting combined sort")
    combined.sort()
    logging.debug("Finished combined sort")

    for p in procs:
        p.join()
        logging.debug("Returned from " + str(p))

    hdr = RawBAM(gzip.GzipFile(bam),header=True).rawheader
    with open(out,'w') as f:
        write_bgzf_block(f, deflate(hdr), zlib.crc32(hdr), len(hdr))

        for creads in grouper(READS_PER_BLOCK, combined):
            creads = [cr for cr in creads if cr is not None]
            data = ""
            tlen = 0
            crc32 = 0
            for i,cr in enumerate(creads):
                d = arys[cr.worker_num][cr.ptr:(cr.ptr+cr.cs)] 
                e = byte_t.pack(byte_t.unpack(d[0])[0] | 1) + d[1:]
                data += e if i==len(creads)-1 else d
                crc32 = crc32_combine(crc32,cr.crc32,cr.rlen)
                tlen += cr.rlen
            write_bgzf_block(f,data,crc32,tlen)

def deflate(data):
    return zlib.compress(data)[2:-4]

def write_bgzf_block(f, data, crc32, isize):
    bgzf_s = BGZFReader.bgzf_s
    bgzf_ftr_s = BGZFReader.bgzf_ftr_s
    bsize = bgzf_s.size + len(data) + bgzf_ftr_s.size
    hdr = (31,139,8,4,0,0,255,6,66,67,2,bsize-1)
    f.write(bgzf_s.pack(*hdr))
    dc = zlib.decompress(data)
    assert zlib.crc32(dc)==crc32
    assert len(dc)==isize
    f.write(deflate(dc))
    # f.write(data)
    # FIXME: have to write as int32
    f.write(int32_t.pack(crc32))
    f.write(uint32_t.pack(isize))

if __name__=="__main__":
    parallel_sort('test.bam', 'out.bam')
