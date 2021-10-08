from os import write


class FileWriter:

    def write_to_file_codons(self, name1, name2, codon_directions, codon_list):
        file_name = "results/" + name1 + "&" + name2 + ".txt"

        f = open(file_name, "a")


        f.write("codons: \n\n")
        f.write(str(len(codon_list)) + "\n")

        for i in range(0, len(codon_list)):
            f.write(codon_list[i] + "\t" + "{:.5f}".format(codon_directions[i]) + "\n")

        f.close()    


    def write_to_file_dicodons(self, name1, name2, protein_directions, protein_list):
        file_name = "results/" + name1 + "&" + name2 + ".txt"

        f = open(file_name, "a")


        f.write("\ndicodons: \n\n")

        for i in range(0, len(protein_list)):
            f.write("\t" + str(protein_list[i]))

        f.write("\n" + str(len(protein_list)))


        for i in range(0, len(protein_list)):
            f.write("\n" + str(protein_list[i]))

            for j in range(0, len(protein_list)):

                f.write("\t" + "{:.4f}".format(protein_directions[i][j]))

        f.close()          