from Bio.Seq import Seq
from FileReader import FileReader
from FileWriter import FileWriter
from FrequencyCalculator import FrequencyCalculator
from ORFFinder import ORFFinder
import sys 


if __name__ == "__main__":

    file_names = ["bacterial1.fasta", "bacterial2.fasta", "bacterial3.fasta", "bacterial4.fasta", "mamalian1.fasta", "mamalian2.fasta", "mamalian3.fasta", "mamalian4.fasta"]
    codons_frequences = []
    dicodons_frequences = []
    file_reader = FileReader()
    file_writer = FileWriter()
    
    frequency_calculator = FrequencyCalculator()

    for file_name in file_names:

        s = file_reader.read_file(file_name)

        orf_finder = ORFFinder()

        s_list = orf_finder.find_all_ORFs(s)

        s = orf_finder.filter_sequences(s_list, 100)

        codon_frequency = frequency_calculator.calculate_all_codons_frequency(s)
        dicodon_frequency = frequency_calculator.calculate_all_dicodons_frequency(s)

        codons_frequences.append(codon_frequency)
        dicodons_frequences.append(dicodon_frequency)

    for i in range(0, len(codons_frequences)-1):
        for j in range(i+1, len(codons_frequences)):

            compared_codon_frequences = frequency_calculator.compare_codon_frequences(codons_frequences[i], codons_frequences[j])
            file_writer.write_to_file_codons(file_names[i], file_names[j], compared_codon_frequences, frequency_calculator.codon_list)

            compared_dicodon_frequences = frequency_calculator.compare_dicodon_frequences(dicodons_frequences[i], dicodons_frequences[j])
            file_writer.write_to_file_dicodons(file_names[i], file_names[j], compared_dicodon_frequences, frequency_calculator.protein_list)   
      
    