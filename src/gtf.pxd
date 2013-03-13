# GTF classes
from cpython cimport bool

cdef extern from "stdint.h":
    ctypedef int uint64_t
    ctypedef int int64_t
    ctypedef int uint32_t

cdef class Feature(object):
    cdef public str seqname
    cdef public str source
    cdef public str feature
    cdef public uint64_t start
    cdef public uint64_t end
    cdef public int64_t score
    cdef public int64_t strand
    cdef public int64_t frame
    cdef public dict attributes

cdef class GTF(object):
    cdef public list features
    cpdef print_gene_table(self, target)