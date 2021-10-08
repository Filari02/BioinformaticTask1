from Bio.Seq import Seq


class ORFFinder:

    start = Seq("ATG")
    stop = [Seq("TAA"), Seq("TAG"), Seq("TGA")]


    def find_all_ORFs(self, seq):
        frames_list = self.find_ORFs(seq)
        frames_list.extend(self.find_ORFs(seq.reverse_complement()))        
        return frames_list


    def find_ORFs(self, seq):
        i = 0
        list = []
        while i < 3:
            codon_started = False
            start_codon = 0
            for j in range(i, len(seq), 3):
                triplet = seq[j:j+3]
                if triplet == self.start and not(codon_started):
                    codon_started = True
                    start_codon = j
                elif any(end in triplet for end in self.stop) and codon_started:
                    codon_started = False
                    list.append(seq[start_codon:j+3])
                            
            i += 1  

        return list

    def filter_sequences(self, seq, number):
        remove_list = []
        for item in seq:
            if len(item) < number:
                remove_list.append(item)

        for remove_item in remove_list:
            seq.remove(remove_item)

        return seq



