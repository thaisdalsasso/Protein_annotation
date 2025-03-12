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
        output.write(f"Number of effectors: {effectors_count}\n")
        output.write(f"Number of apoplastic effectors: {apoplastic_effectors_count}\n")
        output.write(f"Number of cytoplasmic effectors: {cytoplasmic_effectors_count}\n")
        output.write(f"Number of apoplastic/cytoplasmic effectors: {apoplastic_cytoplasmic_effectors_count}\n")
        output.write(f"Number of cytoplasmic/apoplastic effectors: {cytoplasmic_apoplastic_effectors_count}\n")

directory = '/home/dalsasso/annotations/effectorp3'
for filename in os.listdir(directory):
    if filename.endswith("secretome.effectorp3"):
        file_path = os.path.join(directory, filename)
        save_effectors_count(file_path)

