from datetime import datetime
class Sequence:
    ids = [] # Enables to check every loaded sequences object
    __alphabets ={
        'dna_alphabet': ['T'], 'rna_alphabet' : ['U'],
        'aa_alphabet': ['H', 'V', 'B', 'R', 'D', 'F', 'E', 'Z', '+', 'X', 'J', 'M', 'I','Q', 'P', 'Y', 'K', 'L', 'S',
                        'W', '-']  # List of AA specific letters and symbols
     }

    def __init__(self, acs_number, sequence, file, file_type, instance_name):
        self.accession_number = acs_number
        self.file_type = file_type
        self.file = file
        self.sequence = sequence
        self.sequence_type = None
        self.profile = dict()
        self.instance_name = instance_name
        self.head_msg = str(f">> zbik.Alin {datetime.now().strftime('%H:%M')} :")
        for base in set(self.sequence):
            if len(set(self.sequence)) > 5:
                self.sequence_type = 'AA'
                break
            self.profile[base] = self.sequence.count(base)
            if base in Sequence.__alphabets['dna_alphabet']:
                self.sequence_type = 'DNA'
            elif base in Sequence.__alphabets['rna_alphabet']:
                self.sequence_type = 'RNA'
            elif base in Sequence.__alphabets['aa_alphabet']:
                self.sequence_type = 'AA'
        print()
        Sequence.ids.append(instance_name)

    def __str__(self):
        if self.sequence_type == 'AA':
            unit = 'AA'
        elif self.sequence_type == 'DNA' or self.sequence_type == 'RNA':
            unit = 'bp'
        else:
            return 'Instance has not been initialized'
        return f'> {self.instance_name}: {self.sequence_type} Sequence instance of {self.accession_number} entry, ' \
               f'length: {len(self.sequence)} {unit}'
