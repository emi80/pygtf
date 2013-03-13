import os

__QUOTES = ['"', "'"]

cdef class Feature(object):
    def __init__(self, str seqname, str source, str feature, uint64_t start, uint64_t end, int64_t score, int64_t strand, int64_t frame):
        """Initialize a new feature"""
        self.seqname = seqname
        self.source = source
        self.feature = feature
        self.start = start
        self.end = end
        self.score = score
        self.strand = strand
        self.frame = frame
        self.attributes = {}


cdef class GTF(object):
    def __init__(self):
        self.features = []

    cpdef print_gene_table(self, target):
        cdef dict genes = {}

        for feature in self.features:
            if feature.feature == "transcript":
                if "gene_id" in feature.attributes:
                    for id in feature.attributes["gene_id"].split(","):
                        id = feature.attributes["gene_id"]
                        reads = feature.attributes["reads"] + genes.get(id, 0)
                        genes[id] = reads
        for i in genes:
            target.write("%s\t%d\n" % (i, genes[i]))


cpdef GTF read(source, features=None, attributes=None):
    of = source
    if isinstance(source, basestring):
        of = open(source, 'r')

    gtf = GTF()
    feature = None
    for line in of:
        feature = parse(line, features=features, attributes=attributes)
        if feature is not None:
            gtf.features.append(feature)

    if isinstance(source, basestring):
        of.close()

    return gtf


cdef class iterate(object):
    cdef source
    cdef list features
    cdef list attributes
    cdef of

    def __init__(self, source, features=None, attributes=None):
        self.source = source
        self.features = features
        self.attributes = attributes
        self.of = source
        if isinstance(source, basestring):
            self.of = open(source, 'r')

    def __iter__(self):
        return self

    def __next__(self):
        feature = None
        for line in self.of:
            feature = parse(line, features=self.features, attributes=self.attributes)
            if feature is not None:
                return feature
        if isinstance(self.source, basestring):
            self.of.close()
        raise StopIteration()


cpdef print_gene_table(source, target, header=True):
    cdef dict genes = {}
    if not isinstance(source, (list, tuple)):
        source = [source]

    for i in range(len(source)):
        for feature in source[i]:
            if feature.feature == "transcript":
                if "gene_id" in feature.attributes:
                    for id in feature.attributes["gene_id"].split(","):
                        reads_list = genes.get(id, [])
                        if len(reads_list) <= i:
                            reads_list.append(0)

                        reads_list[i] = feature.attributes["reads"] + reads_list[i]
                        genes[id] = reads_list
    if header:
        ns = []
        for i,n in enumerate(source):
            if isinstance(source, basestring):
                ns.append(os.path.basename(source))
            else:
                ns.append("sample_%d" % (i + 1))

        target.write("%s\t%s\n" % ("gene_id", "\t".join(ns)))
    for i in genes:
        target.write("%s\t%s\n" % (i, "\t".join( ["%d" % x for x in genes[i]]) )  )

cpdef Feature parse(str line, features=None, attributes=None):
    """Parse GTF line and reaturn feature, or
    return None if the line is excluded

    line -- the source line
    """
    if line is None:
        return None

    # exclude comments
    line = line.strip()
    if line[0] == "#":
        return None

    # split line and check field count
    fields = line.split("\t")
    if len(fields) < 8:
        raise ValueError("Line has insufficient number of fields: %d", len(fields))

    # check for exclusion by feature
    if features is not None and fields[2] not in features:
        return None

    # parse the fields
    #
    # we have to parse score and frame and strand separately as
    # they might be '.'
    cdef int64_t _score = -1
    cdef int64_t _frame = -1
    cdef int64_t _strand = 0
    if fields[5] != '.':
        _score = int(fields[5])

    if fields[6] == '+':
        _strand = 1
    elif fields[6] == '-':
        _strand = -1

    if fields[7] != '.':
        _frame = int(fields[7])

    feature = Feature(fields[0], fields[1], fields[2], int(fields[3]), int(fields[4]), _score, _strand, _frame)

    # parse attributes
    if len(fields) > 8 and (attributes is None or len(attributes) > 0):
        __parse_attributes(feature, fields[8], attributes)

    return feature

cdef __parse_attributes(feature, source, attributes=None):
    attrs = [s.strip() for s in source.split(";")]
    for attr in attrs:
        if len(attr) == 0:
            continue
        # split name and value
        name, value = attr.split(' ', 1)

        # clean the name
        name = __remove_quotes(name.strip())

        # check for attributes inclusion
        if attributes is not None and name not in attributes:
            continue

        # clean the value and
        # try to figure out its type
        # for now, if the value is quoted,
        # we treat it as string, otherwise,
        # we try a float conversion
        value = value.strip()

        if __is_quoted(value):
            value = __remove_quotes(value)
        else:
            try:
                value = float(value)
            except:
                pass
        feature.attributes[name] = value


cdef bool __is_quoted(source):
    """Return true if the source is quoted"""
    return source[0] in __QUOTES and source[-1] in __QUOTES

cdef str __remove_quotes(source):
    """Remove quotes from the source if the source is quoted.
    Return the original if it is not quoted.
    """
    if __is_quoted(source):
        return source[1:-1]
    else:
        return source


