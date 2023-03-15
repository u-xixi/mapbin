class NbRead:
    def __init__(self, header, seq):
        self.header = header
        self.seq = seq


class NbContig():
    # contig is the basic unit of this network
    # use this class as the basis
    # it might help reduce the overhead
    # and allows modifying or adding info to the contig
    __slots__ = ['id', 'length', 'num_id']

    def __init__(self, id:str, length=NotImplemented, num_id=NotImplemented):
        """
        :param id: NODE_1
        :param length: contig length
        :param num_id: numeric id for informap
        """
        self.id = id
        self.length = length
        self.num_id = num_id