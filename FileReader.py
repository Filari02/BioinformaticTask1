from ORFFinder import ORFFinder
from Bio.Seq import Seq

class FileReader:
    
    def read_file(self, file_name):
        file = open(file_name)
        file.readline()
        sequence = ""

        for line in file:
            #line = line[:-1]
            sequence += line

        file.close()
        
        seq = Seq(sequence)
        return seq

    