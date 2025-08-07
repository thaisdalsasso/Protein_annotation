import os

def save_effectors_count(table_file):
    output_file = os.path.splitext(table_file)[0] + '_counts.txt'

    effectors_count = 0
    apoplastic_effectors_count = 0
    cytoplasmic_effectors_count = 0
    apoplastic_cytoplasmic_effectors_count = 0
    cytoplasmic_apoplastic_effectors_count = 0

    with open(table_file, 'r') as f:
        next(f)  # Skip the header line
        for line in f:
            fields = line.strip().split('\t')
            prediction = fields[-1]

            if prediction != 'Non-effector':
                effectors_count += 1

                if prediction == 'Apoplastic effector':
                    apoplastic_effectors_count += 1
                elif prediction == 'Cytoplasmic effector':
                    cytoplasmic_effectors_count += 1
                elif prediction == 'Apoplastic/cytoplasmic effector':
                    apoplastic_cytoplasmic_effectors_count += 1
                elif prediction == 'Cytoplasmic/apoplastic effector':
                    cytoplasmic_apoplastic_effectors_count += 1

    with open(output_file, 'w') as output:
        output.write(f"Number of effectors:\t{effectors_count}\n")
        output.write(f"Number of apoplastic effectors:\t{apoplastic_effectors_count}\n")
        output.write(f"Number of cytoplasmic effectors:\t{cytoplasmic_effectors_count}\n")
        output.write(f"Number of apoplastic/cytoplasmic effectors:\t{apoplastic_cytoplasmic_effectors_count}\n")
        output.write(f"Number of cytoplasmic/apoplastic effectors:\t{cytoplasmic_apoplastic_effectors_count}\n")

table_file = 'Zpa796.secretome.effectorp3'
save_effectors_count(table_file)

