# Plasmid Designer (C)

A simple C program for plasmid/vector design and management.  
It provides a menu-driven interface to create and modify plasmids as doubly linked lists of genes.

## Features
- Initialize a new plasmid
- Add genes at a specific position
- Delete genes by position
- Edit gene name, sequence, and function
- Print plasmid contents or details of a single gene
- Design PCR primers (forward and reverse) for a selected gene
- Save plasmid to a file (FASTA-like format with legend)
- Load plasmid data from a CSV file (`name,sequence,function`)

## Data structures
- **Doubly linked list** to store genes
- **Gene struct** with:
  - `name`
  - `sequence`
  - `function`

## Example workflow
1. Run the program
2. Use the menu to add genes:
   - Name, DNA sequence (5'–3'), function, and position
3. Print plasmid contents
4. Design PCR primers for selected genes
5. Save the plasmid to a `.txt` file or load from a `.csv` file

## Example CSV format
```bash
lacZ,ATGACCATGATTACGCCA,Therapeutic marker
KanR,ATGGATCCTCTAGAGGAT,Antibiotic resistance
ori,GGGCGCCCTTTCGTCTTCA,Replication origin
```

## Compilation
Use `gcc` (or any C compiler) with support for the C99 standard or newer:
```bash
gcc plasmid_designer.c -o plasmid_designer
```

## Running
```bash
./plasmid_designer
```

## Example menu
```bash
--- Plasmid Manager Menu ---
1) Add gene at position
2) Delete gene at position
3) Print plasmid contents
4) Print gene details from position
5) Design PCR primers for gene
6) Save to file
7) Load plasmid from CSV file
8) Edit gene data
0) Exit
```

## Tested on
	•	Linux Ubuntu 
	•	macOS Sequoia 
	•	macOS Tahoe 
	
## License
This project is licensed under the MIT License.
You are free to use, modify, and distribute the code with proper attribution.