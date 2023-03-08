#autor: Sheyla Carmen
from Bio import SeqIO
from random import *
import pandas as pd0

class validation:

    def __init__(self, name_file,kmer, name_database):
        # print(str(name_file))
        if ".fasta" in name_file:
            self.myfile = SeqIO.parse(name_file, "fasta")
        else:
            self.myfile = pd.read_csv(name_file, header=None, sep="\t")
            self.epitopes = list(self.myfile.iloc[:, 1])
        self.kmer = kmer
        self.database = pd.read_csv(name_database)
        self.experimental = list(self.database['Object'])

    def create_epitopes(self):
        count_all=0
        ran_epitopes = {}
        for record in self.myfile:
            sequence=record.seq
            seq_len=len(sequence)
            count_gen=0
            for seq in list(range(seq_len-(self.kmer-1))):
                count_gen += 1
                count_all += 1
                my_kmer = (sequence[seq:seq+self.kmer])
                iden_key = '>' + str(record.id) + '|kmer_gen_' + str(count_gen) + '|kmer_all_' + str(count_all)
                iden_val = str(my_kmer)
                ran_epitopes[iden_key] = iden_val
        return ran_epitopes

    def ran_1000(self, ran_epitopes):
        if isinstance(ran_epitopes, dict):
            ale = sample(list(ran_epitopes.values()), 1000)
        else:
            ale = sample(list(ran_epitopes), 1000)
            long=[]
            for i in ale:
                long.append(len(i))
            print(len(ale))
        return ale

    def find_1000(self, ale):
        com_find = set(ale).intersection(self.experimental)
        cantidad = len(com_find)
        # print(common_ep, len(common_ep))
        return cantidad, com_find

    def inside_1000(self, ale):
        ale = set(ale)
        com_in = []
        for elem in ale:
            for s in set(self.experimental):
                if elem in s:
                    com_in.append(elem)
        cantidad = len(com_in)
        return cantidad, com_in

    def any_1000(self,ale):
        #si al menos hay una coincidencia, asi que si un epitope largo tiene dos encontrados, lo toma como uno
        ale = set(ale)
        com_any = []
        for s in set(self.experimental):
            if any(substring in s for substring in ale):
                com_any.append(s)
        cantidad = len(com_any)
        return cantidad, com_any

    def ran_verification(self):
        total_epitopes =self.create_epitopes()
        samples = self.ran_1000(total_epitopes)
        cam = self.any_1000(samples)
        return cam

    def pip_verification(self):
        samples_pipe = self.ran_1000(self.epitopes)
        cam_pip = self.any_1000(samples_pipe)
        return cam_pip

ran_hla2 = validation("consenso_pertussis.fasta", 15, "epitope_table_export_1677729886.csv" ).ran_verification()
print("ran", ran_hla2)
print("")
pipe_hla2 = validation("epitopes_HLAII_to _validate_pipeline.csv", 15, "epitope_table_export_1677729886.csv").pip_verification()
print("pipe",pipe_hla2)

# ran_hla1 = validation("consenso_pertussis.fasta", 8, "epitope_table_export_1677729916.csv").ran_verification()
# print(ran_hla1)
#



