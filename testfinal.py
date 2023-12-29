#PART 1


import os
import sys
import textwrap

# If one of the required files was not found it exits the program
if not os.path.exists("exon_positions.txt"):
     print("exon_positions.txt not found")
     sys.exit()
if not os.path.exists("gene_sequence.txt"):
    print("gene_sequence.txt not found")
    sys.exit()
if not os.path.exists("code.txt"):
    print("code.txt not found")
    sys.exit()


#we read the gene_sequence.txt store the text in var called "content" and remove the \n so it's no longer counted as a char
with open("gene_sequence.txt", "r") as o:
    content = o.read().replace('\n', '')

#Function to read the positions of the Exon    
def CheckExonPos(content):
    exon_positions = []
    with open('exon_positions.txt', 'r') as f:
        #read the file and extract the start and the end
        for line in f:
            start, end = line.strip().split()
            exon_positions.append((int(start), int(end)))

    exon_sequence = ""
    for start, end in exon_positions:
        #extract the data from the string from each start to each end
        exon_sequence += content[start:end]
    return exon_sequence

#Store the exon data in "var"
var = CheckExonPos(content)

# Transcription Function... replace "T" with "U"
def transcription(String):
        newString = String.replace('T', 'U')
        return newString
#Store the data transcripted in "var"
var = transcription(var)


#Function that Finds "AUG" starts to splitting from "AUG" until a group of 3 chars is in stop_codons, and adds all of the groups into a a variable
def find_and_translate(string):
    #If AUG not founds returns empty string
    if "AUG" not in string:
        return ""
    #Search for AUG and assign the start index to it
    start_index = string.index("AUG")
    #Assign the stop codons
    stop_codons = {"UAG", "UGA", "UAA"}
    #Start spliting from AUG
    codons = [string[i:i+3] for i in range(start_index, len(string), 3)]
    translated = ""
    #Read each codon and compare it to the threeletter in "code.txt", and assign the add its corresponding oneletter to the var
    for codon in codons:
        with open("code.txt", "r") as f:
            for line in f :
                args = line.split()
                threeletter = args[0]
                oneletter = args[1]
                if codon == threeletter:
                    translated += oneletter
        if codon in stop_codons:
            #remove the oneletter assigned to the stop codon
            translated = translated[:-1]
            break

    return translated

#Store the data after translation in "allnewgroup"
allnewgroup=find_and_translate(var)


#function to write to the file with 80 chars per line
def write_to_file(string, filename):
    wrapped_lines = textwrap.wrap(string, width=80)
    with open(filename, 'w') as file:
        for line in wrapped_lines:
            file.write(line + '\n')
    print(filename,"has been saved with the desired output !", "\n")

write_to_file(allnewgroup, "output.txt")

print('\n','After conducting a BLASTP search, we have identified a match for this sequence with UniProt ID P09848, which corresponds to Lactase/phlorizin hydrolase.')
##################################################################################################################

#PART 2

while True:

    print("Menu:")
    print("-" * 20)
    print("1. ALL")
    print("2. Per Category")
    print("3. Within Category")
    print("4. Specific AA")
    print("5. Exit")

    choice = input("Select an option (1-5): ")
    # Read the output file and strip the Chars in it
    with open('output.txt', 'r') as fileoutput:
        protein_sequence = fileoutput.read().strip()
        # Define a list of amino acids of interest
        amino_acids = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V"]
    # If choice 1 was made (ALL)
    if choice == "1":
        print("Option 1 selected")


        # Count the frequency of each amino acid in the protein sequence and store the results in a dictionary
        amino_acid_counts = {aa: sum(1 for x in protein_sequence if x == aa) for aa in amino_acids}
        # Calculate the total number of amino acids in the protein sequence
        total_aa = sum(amino_acid_counts.values())
        # Calculate the percentage of each amino acid in the protein sequence and store the results in a dictionary
        amino_acid_percentages = {aa: count / total_aa * 100 for aa, count in amino_acid_counts.items()}
        # Print a table of the amino acid frequencies and percentages
        print("Amino Acid | Frequency | Percentage")
        print("-" * 36)

        # Arranges a dictionary of amino acid percentages in descending order and displays each amino acid's symbol, frequency, and percentage in the protein sequence using formatted string with aligned and padded fields
        frequencies = [(aa, amino_acid_counts[aa], percentage) for aa, percentage in sorted(amino_acid_percentages.items(), key=lambda x: x[1], reverse=True)]
        [print(f"{aa:^10} | {frequency:^9} | {percentage:>9.2f}%") for aa, frequency, percentage in frequencies]


    # If choice 2 was made (Per Category)
    elif choice == "2":
        print("Option 2 selected")
        # Open the file "output.txt" and read its contents as a string
        with open('output.txt', 'r') as fileoutput:
            protein_sequence = fileoutput.read().strip()

        # Define lists of amino acids in each category
        positively_charged = ["R", "H", "K"]
        negatively_charged = ["D", "E"]
        polar = ["N", "C", "Q", "S", "T", "Y"]
        nonpolar = ["A", "I", "L", "M", "F", "P", "W", "V", "G"]

        # Count the frequency of each amino acid in the protein sequence and store the results in a dictionary
        aa_counts = {aa: protein_sequence.count(aa) for aa in amino_acids}

        # Initialize dictionaries to store the counts and percentages for each category
        category_counts = {"Positively charged": 0, "Negatively charged": 0, "Polar": 0, "Non-polar": 0}
        category_percentages = {"Positively charged": 0, "Negatively charged": 0, "Polar": 0, "Non-polar": 0}

        # Iterate through each amino acid and add its count to the appropriate category
        for aa, count in aa_counts.items():
            if aa in positively_charged:
                category_counts["Positively charged"] += count
            elif aa in negatively_charged:
                category_counts["Negatively charged"] += count
            elif aa in polar:
                category_counts["Polar"] += count
            elif aa in nonpolar:
                category_counts["Non-polar"] += count

        # Calculate the total number of amino acids in the protein sequence
        total_aa = sum(category_counts.values())

        # Calculate the percentage of each category in the protein sequence
        category_percentages = {category: count / total_aa * 100 for category, count in category_counts.items()}

        # Print a table of the amino acid category frequencies and percentages
        print("Amino Acid Category | Frequency | Percentage")
        print("-" * 45)
        [print(f"{category:^20} | {category_counts[category]:^9} | {percentage:>9.2f}%") for category, percentage in sorted(category_percentages.items(), key=lambda x: x[1], reverse=True)]


    # If choice 3 was made (Within Category)
    elif choice == "3":
        print("Option 3 selected")
        # Open the file "output.txt" and read its contents as a string
        with open('output.txt', 'r') as fileoutput:
            protein_sequence = fileoutput.read().strip()

        # Define lists of amino acids in each category
        positively_charged = ["R", "H", "K"]
        negatively_charged = ["D", "E"]
        polar = ["N", "C", "Q", "S", "T", "Y"]
        nonpolar = ["A", "I", "L", "M", "F", "P", "W", "V", "G"]

        # Define a dictionary that maps category names to their corresponding lists of amino acids
        categories = {
            "Positively charged": positively_charged,
            "Negatively charged": negatively_charged,
            "Polar": polar,
            "Non-polar": nonpolar
        }
        print("Choose a Category:")
        print("1. Positively charged")
        print("2. Negatively charged")
        print("3. Polar")
        print("4. Non-polar")

        # Prompt the user to select an amino acid category
        choice = input("Please choose an amino acid category: ")
        if choice == "1":
            category="Positively charged"
        elif choice == "2":
            category="Negatively charged"
        elif choice == "3":
            category="Polar"
        elif choice == "4":
            category="Non-polar"
        else:
            print("Category not found")
            continue

        # Check that the user input is valid
        if category not in categories:
            print(f"Invalid category: {category}")
        else:
            # Count the frequency of each amino acid in the selected category and store the results in a dictionary
            aa_counts = {aa: protein_sequence.count(aa) for aa in categories[category]}

            # Calculate the total number of amino acids in the selected category
            total_aa = sum(aa_counts.values())

            # Calculate the percentage of each amino acid in the selected category
            aa_percentages = {aa: count / total_aa * 100 for aa, count in aa_counts.items()}

            # Print a table of the amino acid frequencies and percentages for the selected category
            print(f"Amino Acid Frequencies and Percentages for {category}")
            print("Amino Acid | Frequency | Percentage")
            print("-" * 36)
            [print(f"{aa:^10} | {aa_counts[aa]:^9} | {percentage:>9.2f}%") for aa, percentage in sorted(aa_percentages.items(), key=lambda x: x[1], reverse=True)]


    # If choice 4 was made (Specific AA)
    elif choice == "4":
        print("Option 4 selected")
        # Open the file "output.txt" and read its contents as a string
        with open('output.txt', 'r') as fileoutput:
            protein_sequence = fileoutput.read().strip()

        # Prompt the user to enter a one-letter code for the amino acid of interest
        aa = input("Please enter a one-letter code for the amino acid of interest: ")

        # Check that the user input is valid
        if len(aa) != 1 or aa not in amino_acids:
            print(f"Invalid amino acid code: {aa}")
        else:
            # Count the frequency of the amino acid of interest and store the result in a variable
            aa_count = protein_sequence.count(aa)

            # Calculate the total number of amino acids in the protein sequence
            total_aa = len(protein_sequence)

            # Calculate the percentage of the amino acid of interest in the protein sequence
            aa_percentage = aa_count / total_aa * 100

            # Print the percentage of the amino acid of interest
            print(f"Percentage of {aa}: {aa_percentage:.2f}%")


    # If choice 5 was made (Exit)
    elif choice == "5":
        # Quit the program
        print("Quitting...")
        break
    else:
        print("Invalid choice. Please select an option from 1-5.")

##################################################################################################################

#PART 3
while True:
    CHOICE = input("Do you want to proceed to the mutation part ?  1.Yes    2.No : ")
    if CHOICE == "1":
        #Replace the character[30049] with the mutated nucleotide
        with open('gene_sequence.txt', 'r') as variable:
            chars=variable.read()
            chars = chars.replace('\n', '')
            chars = chars[:30049] + 'A' + chars[30049+1:]
            print("\n","A mutation has occured at position 30049 and replaced T with A... The goal of this program is to check the effect of this mutation.", "\n")


        #Same as part 1

        alltxt3 = CheckExonPos(chars)
        s = transcription(alltxt3)
        groups3= find_and_translate(s) 

        write_to_file(groups3, "part3.txt")

        print("A T nucleotide was swapped out for an A nucleotide at position 30049 of the gene sequence as a consequence of a point mutation.An early stop codon can be produced as a result of point mutations, which are changes to a single nucleotide in the DNA sequence. Early protein synthesis termination brought on by premature stop codons results in a truncated protein product that is shorter than the original protein. Additionally, it may be inferred from the information at hand that the altered gene sequence is shorter than the original gene sequence. The observation is based on the fact that the altered sequence abruptly ended at the A nucleotide at position 30049, but the original sequence went on for a number of further characters after this.")
        break
    elif CHOICE == "2":
        print("Good choice, mutations could be harmful... Have a nice day !")
        break

    else :
        print("Invalid input, Try again !")

