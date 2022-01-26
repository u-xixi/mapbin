import gzip


def plain_reader(file_name):
    with open(file_name, "r") as f:
        return f.read().strip().split("\n")
        # return f.readlines()


def plain_writer(file_as_str, file_name):
    with open(file_name, "w") as f:
        f.write(file_as_str)


def file_uncompress(file_path):
    with gzip.open(file_path, "rt") as f:
        return f.read().strip().split("\n")
        # return f.readlines()


def file_compress(file_as_string, file_path):
    with gzip.open(file_path, "wt") as f:
        f.write(file_as_string)


def isgzip_reader(file_path):
    try:
        if file_path.endswith(".gz") or \
                file_path.endswith("gzip"):
            return file_uncompress(file_path)

        else:
            return plain_reader(file_path)
    except ValueError:
        print("File corrupted or file path invalid!")


def format_checker(file_content):
    if file_content[0].startswith("@"):  # improve later!!
        if file_content[2].startswith("+"):
            return "fq"
        else:
            raise SystemExit("ValueError: invalid fastq file")
    elif file_content[0].startswith(">"):
        return "fa"
    else:
        raise SystemExit("ValueError: invalid sequence file")


class SeqFile:
    def __init__(self, file_path):
        # self.gzip = False
        self.seqformat = ""
        self.idf = ""
        self.step = -1
        self.file_content = []
        self.__infile_receptor(file_path)

    def isgzip(self, file_path):
        if file_path.endswith(".gz") or \
                file_path.endswith("gzip"):
            try:
                self.file_content = file_uncompress(file_path)
            except ValueError:
                print("Invalid gz file!")
        else:
            try:
                self.file_content = plain_reader(file_path)
            except ValueError:
                print("File corrupted or file path invalid!")

    def format_checker(self):
        if self.file_content[0].startswith("@"):  # improve later!!
            if self.file_content[2].startswith("+"):
                self.step = 4
                self.idf = "@"
                self.seqformat = "fastq"
            else:
                raise SystemExit("ValueError: invalid fastq file")
        elif self.file_content[0].startswith(">"):
            self.step = 2  # in fact not correct
            self.idf = ">"
            self.seqformat = "fasta"
        else:
            raise SystemExit("ValueError: invalid sequence file")

    def __infile_receptor(self, file_path):
        self.isgzip(file_path)
        self.format_checker()


class SamFile:
    def __init__(self, filepath, header, alignments):
        self.header = header
        self.alignments = alignments

    def issam(self):
        pass

    def __sam_reader(self, filepath):
#         This is temporary!!
        pass


