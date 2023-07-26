# Ishay Eldar

import sys
import csv
import os

class Cell:
    def __init__(self, name, genome):
        """
        Initialize a Cell object.

        :param name: Name of the cell
        :param genome: List of genetic sequences
        """
        self.name = name
        self.genome = genome
        self.num_sequence = len(genome)

    def __str__(self):
        """
        Return a string representation of the Cell object.

        :return: String representation of the Cell object
        """
        return f'<{self.name}, {self.num_sequence}>'

    def find_ssr(self, genome_seq_index):
        """
        Find repeated sequences in the genome.

        The function receives an input of a genetic sequence and returns a dictionary of the number of repetitions
        of each segment up to 6 nucleotides long.

        :param genome_seq_index: Index of the genetic sequence in the genome
        :return: Dictionary containing the repeated sequences and their repetitions
        """
        ssr_dict = {}  # The dictionary that the function will return
        counter = 1  # number of repeats
        string = self.genome[genome_seq_index % self.num_sequence][0]  # the sequence that we check for SSR

        for n in range(1, 7):
            for index in range(len(string) - n):
                counter = 1
                ssr_length = string[index: index + n]  # the SSR that we're checking
                if string[index + n: index + (2 * n)] == ssr_length:
                    for k in range(len(string) // n):
                        if string[index + n + k * n: index + n + (
                                k + 1) * n] == ssr_length:  # We found an SSR with 3 repeats
                            counter += 1
                        else:
                            break
                    if counter >= 3:  # We will add the SSR to the dictionary
                        if ssr_length not in ssr_dict or counter > ssr_dict[ssr_length]:
                            ssr_dict[ssr_length] = counter

        if not any(ssr_dict):  # We didn't find any SSR
            return None

        ssr_list = self.dict_to_list(ssr_dict)

        return ssr_list

    def dict_to_list(self, dict_ssr):
        """
        Convert the dictionary of SSRs and their repetitions to a list of lists.

        The function passes the dictionary to a sorted list that adds ";" between the values.

        :param dict_ssr: Dictionary of SSRs and their repetitions
        :return: List of lists representing the SSRs and their repetitions
        """

        if not any(dict_ssr):  # We didn't find any SSR
            return None
        # lists
        temp = []
        dict_List = []

        # My attempt:
        for key, value in dict_ssr.items():
            aKey = key
            aValue = value
            temp.append(aKey)
            temp.append(aValue)
            # print(temp)
            dict_List.append(temp)
            temp = []
            aKey = ""
            aValue = ""
            # print(dict_List)

        return dict_List

    def ssr_printer(self, genome_seq_index):
        """
        Print the SSRs and their repetitions as a sorted list.

        The function checks whether the dictionary is empty or the input is empty.
        Otherwise, it prints the dictionary values as a sorted list.

        :param genome_seq_index: Index of the genetic sequence in the genome
        :return: String representation of the sorted SSR list
        """
        str_to_print = ""

        if len(self.genome[genome_seq_index % self.num_sequence][0]) == 0:
            return "No simple repeats in DNA sequence"
        if self.find_ssr(genome_seq_index) is None:
            return "No simple repeats in DNA sequence"
        else:
            list_sorted_keys = sorted(self.find_ssr(genome_seq_index))
            # print(list_sorted_keys)
            for i in range(len(list_sorted_keys)):
                # print(key)
                str_to_print += str(list_sorted_keys[i][0])
                str_to_print += ','
                str_to_print += str(list_sorted_keys[i][1])
                str_to_print += ';'

            return str_to_print[:-1]

    def transcribe(self, genome_seq_index):
        """
        Transcribe the DNA sequence to RNA.

        The function receives an input of a DNA sequence and returns the RNA that will be obtained
        in the transcription process.

        :param genome_seq_index: Index of the genetic sequence in the genome
        :return: Transcribed RNA sequence
        """
        dna = self.genome[genome_seq_index % self.num_sequence][0]
        rna_temporary = dna[len(dna) - 1::-1].upper().replace('A', 'U').replace('T', 'A').replace('C', 'g').replace('G', 'C')
        rna = rna_temporary.replace('g', 'G')
        return rna


    def transcribe_printer(self, genome_seq_index):
        """
        Print the RNA sequence obtained from the input DNA sequence.

        The function prints the RNA sequence obtained from the input DNA sequence.

        :param genome_seq_index: Index of the genetic sequence in the genome
        :return: Tuple containing the label "RNA sequence" and the transcribed RNA sequence
        """
        rna_to_print = self.transcribe(genome_seq_index)
        return ("RNA sequence:", rna_to_print)

    def translate(self, genome_seq_index):
        """
        Translates the DNA sequence to amino acids.

        The function receives an input of a DNA sequence, transcribes it into RNA, and translates the RNA sequence
        into the corresponding amino acids based on the genetic code. It identifies the longest possible protein
        sequence and returns the translated amino acids.

        :param genome_seq_index: Index of the genetic sequence in the genome
        :return: String containing the translated amino acids
        """
        str_amino_acids = ""  # A string of the amino acids that the function will return
        reading_frame = self.genome[genome_seq_index % self.num_sequence][1]  # The reading frame of the sequence
        dna_sec = self.genome[genome_seq_index % self.num_sequence][0]  # The DNA sequence
        rna_sec = self.transcribe(genome_seq_index)  # Transcribed RNA sequence

        dict_of_amino_acids = {  # A dictionary of the letters that mark the amino acids
            "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
            "UCU": "S", "UCC": "S", "UCA": "S", "UCG": "S",
            "UAU": "Y", "UAC": "Y", "UAA": "", "UAG": "",
            "UGU": "C", "UGC": "C", "UGA": "*", "UGG": "W",
            "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
            "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
            "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
            "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
            "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
            "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
            "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
            "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
            "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
            "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
            "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
            "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"
        }

        start_codon = "AUG"
        stop_codons = ["UAA", "UAG", "UGA"]

        length_protein = 0  # Saves the length of the current protein
        max_length = 0  # Saves the length of the longest protein
        start_translate = 0  # Saves the location of the start codon
        end_translate = 0  # Saves the location of the stop codon or end of translation

        # Determine the longest possible protein sequence
        for index in range(reading_frame, len(rna_sec) - (3 + reading_frame), 3):
            if rna_sec[index:index + 3] == start_codon:
                length_protein = 1
                for n in range(1, (len(rna_sec) // 3) + 1):
                    if index + (3 * n) > len(rna_sec) - 3:  # Reached the end of the sequence
                        if length_protein > max_length:
                            max_length = length_protein
                            start_translate = index
                            end_translate = index + 3 * (max_length)
                        break
                    elif rna_sec[index + 3 * n: index + 3 * (n + 1)] in stop_codons:  # Reached a stop codon
                        if length_protein > max_length:
                            max_length = length_protein
                            start_translate = index
                            end_translate = index + 3 * (max_length)
                        break
                    else:
                        length_protein += 1

        # Translate the longest sequence into amino acids
        if max_length > 0:
            for j in range(start_translate, end_translate, 3):
                codon = rna_sec[j: j + 3]
                amino_acid = dict_of_amino_acids[codon]
                str_amino_acids += amino_acid

        if str_amino_acids == "":
            return 0

        tag = ";"
        str_amino_acids_with_tags = tag.join(str_amino_acids)
        return str_amino_acids_with_tags

    def translate_printer(self, genome_seq_index):
        """
        Prints the longest possible amino acid sequence from the given RNA sequence.

        This function calls the `translate` method to obtain the translation of the RNA sequence
        at the specified index in the genome. If the translation is determined to be non-coding RNA
        (indicated by the translation being 0), it returns the message "Non-coding RNA". Otherwise,
        it returns the message "Translation: " followed by the obtained translation.

        :param genome_seq_index: Index of the RNA sequence in the genome
        :return: String representing the translation or non-coding RNA message
        """
        translate_to_print = self.translate(genome_seq_index)
        if translate_to_print == 0:
            return "Non-coding RNA"
        else:
            return "Translation: " + translate_to_print

    def repertoire(self):
        """
        Creates a repertoire of tuples containing the output of SSR and translation for each RNA sequence.

        This function iterates over each RNA sequence in the genome and calls the `ssr_printer` and `translate_printer`
        methods to obtain the SSR and translation outputs, respectively. It then creates a tuple with the SSR output
        and translation output and appends it to the `repertoire_list`. Finally, it returns the `repertoire_list`
        containing all the tuples.

        :return: List of tuples representing the repertoire of SSR and translation outputs
        """
        repertoire_list = []
        temp = ()
        for i in range(self.num_sequence):
            temp = (self.ssr_printer(i), self.translate_printer(i))
            repertoire_list.append(temp)
            temp = ()
            # temp_ssr = ""
            # temp_translate = ""
        return repertoire_list


class StemCell(Cell):
        """
        A class representing a stem cell.

        This class inherits from the Cell class and provides additional functionality specific to stem cells.
        """

        class NerveCell(Cell):
            """
            A class representing a nerve cell derived from a stem cell.

            This class inherits from the Cell class and represents a nerve cell that is derived from a stem cell.
            """

            def __init__(self, stem_cell, internal_coefficient):
                """
                Initializes a NerveCell object.

                :param stem_cell: The stem cell object from which the nerve cell is derived.
                :param internal_coefficient: The internal coefficient value for the nerve cell.
                """
                super().__init__("Nerve Cell", stem_cell.genome)
                self.internal_coefficient = internal_coefficient
                self.signal = 0

            def receive(self, signal):
                """
                Receives a signal by the nerve cell.

                :param signal: The signal received by the nerve cell.
                """
                self.signal = signal

            def send(self):
                """
                Sends a signal from the nerve cell.

                :return: The value of the signal multiplied by the internal coefficient.
                """
                return float(self.signal) * float(self.internal_coefficient)

        class MuscleCell(Cell):
            """
            A class representing a muscle cell derived from a stem cell.

            This class inherits from the Cell class and represents a muscle cell that is derived from a stem cell.
            """

            def __init__(self, stem_cell, file_path, threshold):
                """
                Initializes a MuscleCell object.

                :param stem_cell: The stem cell object from which the muscle cell is derived.
                :param file_path: The file path where the muscle cell writes its output.
                :param threshold: The threshold value for the muscle cell.
                """
                super().__init__("Muscle Cell", stem_cell.genome)
                self.threshold = threshold
                self.file_path = file_path

            def receive(self, signal):
                """
                Receives a signal by the muscle cell and writes it to a file if it exceeds the threshold.

                :param signal: The signal received by the muscle cell.
                :return: The received signal.
                """
                if signal >= int(self.threshold):
                    self.write_to_file(signal)
                return signal

            def write_to_file(self, signal):
                """
                Writes the signal to a file.

                :param signal: The signal to be written to the file.
                """
                with open(self.file_path, "a") as file:
                    file.write(f"{signal}, I like to move it\n")

        def __mul__(self, multiplier):
            """
            Performs multiplication of a stem cell object.

            :param multiplier: The number of times the stem cell should be multiplied.
            :return: A list of string representations of the multiplied stem cells.
            """
            cells = [self] + [self.__class__(self.name, self.genome) for _ in range(multiplier - 1)]
            return [str(cell) for cell in cells]

        def mitosis(self):
            """
            Performs mitosis of a stem cell.

            :return: A list containing the original stem cell and a new stem cell.
            """
            cells = [self, self.__class__(self.name, self.genome)]
            return cells

        def differentiate(self, cell_type, parameters):
            """
            Differentiates the stem cell into a specific cell type.

            :param cell_type: The type of cell to differentiate into (either "Nerve Cell" or "Muscle Cell").
            :param parameters: The parameters required for the differentiation.
                               - If the cell type is "Nerve Cell", parameters should be a list containing one element.
                               - If the cell type is "Muscle Cell", parameters should be a string containing comma-separated values.

            :return: The differentiated cell object.
            :raises ValueError: If the provided cell type is unsupported.
            """
            if cell_type == "Nerve Cell":
                parameter = parameters[0]
                new_nerve = self.NerveCell(self, parameter)
                return new_nerve
            elif cell_type == "Muscle Cell":
                parameters_list = parameters.split(",")
                new_muscle = self.MuscleCell(self, *parameters_list)
                return new_muscle
            else:
                raise ValueError("Unsupported cell type")

class NerveNetwork:
            """
            A class representing a nerve network.

            This class encapsulates a group of nerve cells and a muscle cell, and provides methods to send signals through the network.
            """

            def __init__(self, nerve_cells, muscle_cell):
                """
                Initializes a NerveNetwork object.

                :param nerve_cells: The list of nerve cells in the network.
                :param muscle_cell: The muscle cell in the network.
                """
                self.nerve_cells = nerve_cells
                self.muscle_cell = muscle_cell

            def send_signal(self, signal):
                """
                Sends a signal through the nerve network.

                The signal is passed through each nerve cell in the network, and the final signal is received by the muscle cell.

                :param signal: The initial signal to send through the network.
                :return: The final signal received by the muscle cell.
                """
                for nerve_cell in self.nerve_cells:
                    nerve_cell.receive(signal)
                    signal = nerve_cell.send()
                return self.muscle_cell.receive(signal)

            def __str__(self):
                """
                Returns a string representation of the nerve network.

                The string representation includes the string representations of all nerve cells followed by the string representation of the muscle cell.

                :return: A string representation of the nerve network.
                """
                cells_str = '\n'.join(str(cell) for cell in self.nerve_cells)
                return cells_str + '\n' + str(self.muscle_cell)

if __name__ == '__main__':
    name_file = sys.argv[1]
    signals_str = sys.argv[2]
    signals_list = signals_str.split(',')

    with open(name_file, 'r') as input_file:
        reader = csv.DictReader(input_file, delimiter='\t')
        nc_list = []

        for row in reader:
            #print(row['type'])
            type = row['type']
            dna_secs = row['DNA']
            letters_list = dna_secs.split(',')
            reading_frames = row['reading_frames']
            reading_frames_list = reading_frames.split(',')
            parameter = row['parameter']
            parameters_list = parameter.split(',')

            #parameters_list = parameter.split(',')


            # Check Category
            assert type in ['MC', 'NC'], "File illegal"

            # Check Sequence Letters
            for letters in letters_list:
                assert all(letter in 'ACTG' for letter in letters), "File illegal"

            # Check Sequence Digits
            for digits in reading_frames_list:
                assert all(digit in '012' for digit in digits), "File illegal"

            # Check identical number of items in columns
            assert len(letters_list) == len(reading_frames_list), "File illegal"

            genome = [(letters_list[i], int(reading_frames_list[i])) for i in range(len(letters_list))]
            #print(genome)
            source_stemCell = StemCell('Stem Cell', genome)

            if type == 'NC':
                #print (parameter.split(','))
                assert len(parameters_list) == 1, "File illegal"
                #assert isinstance(parameters_list[0], float), "File illegal"
                assert float(parameters_list[0]) > 0, "File illegal"
                #print(source_stemCell.num_sequence)
                temp_list = source_stemCell.mitosis()
                #print(temp_list[0])
                #print(temp_list[1])
                #print(source_stemCell.differentiate("Nerve Cell", parameters_list))
                #print(temp_list[0].differentiate("Nerve Cell", parameters_list))
                #nerve_1 = temp_list[0].differentiate("Nerve Cell", parameters_list)
                #nerve_2 = temp_list[1].differentiate("Nerve Cell", parameters_list)
                nc_list.append(temp_list[0].differentiate("Nerve Cell", parameters_list))
                nc_list.append(temp_list[1].differentiate("Nerve Cell", parameters_list))
                #print(len(nc_list))
                #print(len(nc_list) + 1)

                #nc_list.extend(temp_list[1].differentiate("Nerve Cell", parameters_list[1]))
                temp_list = []
                #print(nc_list)


            elif type == 'MC':
                assert len(parameters_list) == 2, "File illegal"
                assert float(parameters_list[1]) > 0 , "File illegal"
                mc = source_stemCell.differentiate("Muscle Cell", parameter)

    nerve_network_to_print = NerveNetwork(nc_list, mc)
    if len(signals_list) > 0:
        for signal in signals_list:
            nerve_network_to_print.send_signal(signal)


    if len(nc_list) > 0:
        print(nerve_network_to_print)
        repertoire_to_print = mc.repertoire()
        for i in range(len(repertoire_to_print)):
            print('\n'.join(repertoire_to_print[i]))



