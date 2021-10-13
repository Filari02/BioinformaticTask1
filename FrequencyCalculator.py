from Bio.Seq import Seq

class FrequencyCalculator:

    
    codon_list = []
    protein_list = []

    def __init__(self):
        codons = ['T', 'A', 'C', 'G']

        for codon1 in codons:
            for codon2 in codons:
                for codon3 in codons:
                    codon = codon1 + codon2 + codon3
                    self.codon_list.append(codon)

        for codon in self.codon_list:
            protein = Seq(codon).translate()
            if not any(p in protein for p in self.protein_list) and protein != "*":
                self.protein_list.append(protein)       


    def calculate_codon_frequency(self, codon, seq_list):

        counter = 0 
        all_codons = 0 
        for seq in seq_list: 
            for i in range(0, len(seq), 3):             
                if codon == seq[i:i+3]:
                    counter += 1

                all_codons += 1    

        frequency = counter/all_codons * 100
        return frequency


    def calculate_all_codons_frequency(self, seq_list):
        codon_frequences = [0 for i in range(len(self.codon_list))]
        for codon in self.codon_list:
            codon_frequences[self.codon_list.index(codon)] = self.calculate_codon_frequency(codon, seq_list)
        
        return codon_frequences




    def calculate_dicodon_frequency(self, dicodon, seq_list):
        counter = 0
        all_dicodons = 0
        for seq in seq_list:
            for i in range (0, len(seq)-3, 3):
                dicodon1 = Seq(seq[i:i+6]).translate()
                if dicodon == dicodon1:
                    counter += 1
                all_dicodons += 1

        frequency = counter/all_dicodons * 100
        return frequency

    def calculate_all_dicodons_frequency(self, seq_list):
        dicodon_frequences = [[0 for i in range(len(self.protein_list))] for j in range(len(self.protein_list))]
        for protein1 in self.protein_list:
            for protein2 in self.protein_list:
                     
                dicodon = protein1 + protein2


                dicodon_frequences[self.protein_list.index(protein1)][self.protein_list.index(protein2)] = self.calculate_dicodon_frequency(dicodon, seq_list)
        

        return dicodon_frequences


    def compare_codon_frequences(self, codon_frequences1, codon_frequences2):
        distances = [0 for i in range(len(self.codon_list))]

        for i in range(0,len(self.codon_list)):

            
            codon1 = codon_frequences1[i]
            codon2 = codon_frequences2[i]
            
            if codon1 == 0 and codon2 == 0:
                    distances[i] = 1
                    continue


            distances[i] = abs(codon1 - codon2) / max(codon1, codon2)


        return distances    


    def compare_dicodon_frequences(self, protein_frequences1, protein_frequences2):
        distances = [[0 for i in range(len(self.protein_list))] for j in range (len(self.protein_list))]

        for i in range(0,len(self.protein_list)):
            for j in range(0, len(self.protein_list)):

                protein1 = protein_frequences1[i][j]
                protein2 = protein_frequences2[i][j]


                if protein1 == 0 and protein2 == 0:
                    distances[i][j] = 1
                    continue
               
                distances[i][j] = abs(protein1 - protein2) / max(protein1, protein2)


        return distances





        
