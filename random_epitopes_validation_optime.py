import random
from Bio import SeqIO
from random import *
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


class EpitopeGenerator:
    def __init__(self, name_file, kmer):
        if ".fasta" in name_file:
            self.myfile = SeqIO.parse(name_file, "fasta")
            self.myfile2 = SeqIO.parse(name_file, "fasta")
        else:
            self.myfile = pd.read_csv(name_file, header=None)
        self.kmer = kmer

    def create_epitopes(self, kmer, myfile):
        count_all = 0
        random_epitopes = {}
        for record in myfile:
            sequence = record.seq
            seq_len = len(sequence)
            count_gen = 0
            for seq in list(range(seq_len - (kmer - 1))):
                count_gen += 1
                count_all += 1
                my_kmer = sequence[seq : seq + kmer]
                iden_key = (
                    ">"
                    + str(kmer)
                    + str(record.id)
                    + "|kmer_gen_"
                    + str(count_gen)
                    + "|kmer_all_"
                    + str(count_all)
                )
                iden_val = str(my_kmer)
                random_epitopes[iden_key] = iden_val
        return random_epitopes

    def kmers_creation(self):
        if isinstance(self.kmer, list):
            group1 = self.create_epitopes(self.kmer[0], self.myfile)
            group2 = self.create_epitopes(self.kmer[1], self.myfile2)
            epitopes = group2 | group1
            new_keys = list(epitopes.keys())
            shuffle(new_keys)
            ran_epitopes = {}
            for clave in new_keys:
                ran_epitopes[clave] = epitopes[clave]
        else:
            ran_epitopes = self.create_epitopes(self.kmer, self.myfile)
        return ran_epitopes


class Validation:
    def __init__(self, name_database):
        self.database = pd.read_csv(name_database)
        self.experimental = list(self.database["Type"])

    def ran_1000(self, ran_epitopes):
        if isinstance(ran_epitopes, dict):
            return sample(list(ran_epitopes.values()), 967)
        else:
            return sample(list(ran_epitopes), 967)

    def inside_1000(self, ale):
        ale = set(
            ale
        )  # Convertimos a set para asegurar que no haya duplicados en los epitopos del pipeline
        com_in = []

        for elem in ale:
            for s in set(self.experimental):
                if elem in s or elem == s or s in elem:
                    com_in.append(elem)
                    break  # Detenemos la búsqueda para evitar contar múltiples coincidencias del mismo epítopo

        com_in = set(com_in)  # Convertimos a set para eliminar duplicados
        cantidad = len(com_in)
        return cantidad, list(com_in)

    def ran_verification(self, ran_epitopes):
        samples = self.ran_1000(ran_epitopes)
        return self.inside_1000(samples)

    def pip_verification(self, epitopes):
        return self.inside_1000(epitopes)


class Analysis:
    @staticmethod
    def save_results_and_calculate_prob(results, upper_limit, filename):
        count = 0
        total = len(results)

        with open(filename, "w") as archivo:
            for elemento in results:
                archivo.write(str(elemento) + "\n")
                if elemento >= upper_limit:
                    count += 1

        probabilidad = (count / total) * 100
        print(
            f"La probabilidad se encuentren la misma cantidad de epitopes al azar como por el pipeline es: {probabilidad:.5f}%"
        )

    @staticmethod
    def plot_results(results, reference_value, output_file):
        df = pd.DataFrame(results, columns=["Numero"])
        conteo = df["Numero"].value_counts().sort_index()
        sns.set(style="darkgrid")
        plt.figure(figsize=(10, 6))
        plt.axvline(reference_value, color="red")
        sns.lineplot(x=conteo.index, y=conteo.values, marker="o")
        plt.fill_between(conteo.index, conteo.values, alpha=0.45)
        plt.savefig(output_file)
        plt.show()


# Configuración
name_file = "files/pertusis.fasta"
kmer = [9, 8]  # O 15 si es hla2
name_database = "files/epitope_table_export_1677729916.csv"
pipeline_file = "files/epitopes_HLAI_to_validate_pipeline.csv"

# Generación de epitopos
generator = EpitopeGenerator(name_file, kmer)
total_epitopes = generator.kmers_creation()

# Validación para el valor de referencia del pipeline
pipeline_validator = Validation(name_database)
pipe_hla2 = pipeline_validator.pip_verification(
    pd.read_csv(pipeline_file, header=None).iloc[:, 1].tolist()
)

# Validación con k-mers generados
validator = Validation(name_database)
n_ran_hla2 = []

for i in range(1000000):
    print(i)
    ran_hla2 = validator.ran_verification(total_epitopes)
    n_ran_hla2.append(ran_hla2[0])

# Guardar y graficar resultados
Analysis.save_results_and_calculate_prob(n_ran_hla2, pipe_hla2[0], "hla1.txt")
Analysis.plot_results(n_ran_hla2, pipe_hla2[0], "hla1.svg")
