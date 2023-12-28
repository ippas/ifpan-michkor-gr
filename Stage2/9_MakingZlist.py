def extract_phenotypes(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    components = {}
    current_component = None
    phenotypes_started = False

    for line in lines:
        line = line.strip()
        if line.startswith("Top phenotypes for Component"):
            current_component = line.split()[-1][:-1]
            components[current_component] = []
            phenotypes_started = True
        elif phenotypes_started and line.startswith("["):
            phenotypes = list(map(int, line[1:-1].split(', ')))
            components[current_component].extend(phenotypes)
        elif phenotypes_started and not line:
            phenotypes_started = False

    return components

def write_phenotypes_to_file(components, output_file):
    with open(output_file, 'w') as file:
        file.write("Row\tPhenotype_Code\tPhenotype_Description\n")
        for component, phenotypes in components.items():
            file.write(f"{component}:\n")
            for i, phenotype in enumerate(phenotypes[:20], start=1):
                file.write(f"{i}\t{phenotype}\tPhenotype_Description\n")
            file.write("\n")

if __name__ == "__main__":
    input_file_path = "top_phenotypes_outputFreq50.txt"
    output_file_path = "zlistFreq50.txt"

    extracted_components = extract_phenotypes(input_file_path)
    write_phenotypes_to_file(extracted_components, output_file_path)

    print(f"Output saved to {output_file_path}")
