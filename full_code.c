#include <stdio.h>
#include <stdlib.h>
#include <string.h> // I use this library to work with strings more easily as char

typedef struct node_t {
    struct node_t *prev;
    struct node_t *next;
    void *data; 
} node_t;

typedef struct {
    node_t *head;
    node_t *tail;
} list_t;

typedef struct {
    char *name;
    char *sequence;
    char *function;
} gene_t;

list_t plasmid;

// FUNCTION FOR COPYING STRING

char *strdup_or_null(const char *s) { // const means s does not change, avoids warnings
    size_t len = strlen(s) + 1; // +1 for terminating character
    char *p = malloc(len); // allocate space for data to be copied
    if (!p) return NULL;
    memcpy(p, s, len); // copy len bytes from s to p
    return p;
}

// DECLARATIONS FOR DOUBLY LINKED LISTS

void listinit(list_t *l) {
    l->head = l->tail = NULL;
}

int listput(list_t *l, void *data) {
    node_t *n = malloc(sizeof(node_t));
    if (!n) return -1;
    n->data = data;
    n->prev = l->tail;
    n->next = NULL;
    if (l->tail) l->tail->next = n;
    else         l->head = n;
    l->tail = n;
    return 0;
}

int listinsert(list_t *l, void *data, unsigned int pos) {
    node_t *new_node = malloc(sizeof(node_t));
    if (!new_node) return -1;
    new_node->data = data;
    if (pos <= 1 || !l->head) {
        new_node->prev = NULL;
        new_node->next = l->head;
        if (l->head) l->head->prev = new_node;
        else         l->tail = new_node;
        l->head = new_node;
        return 0;
    }
    node_t *cur = l->head;
    unsigned int idx = 1;
    while (cur->next && idx < pos - 1) { // move to the node before the desired position (if pos=3, move to node #2)
        cur = cur->next;
        idx++;
    }
    new_node->prev = cur;
    new_node->next = cur->next;
    if (cur->next) cur->next->prev = new_node;
    else           l->tail = new_node;
    cur->next = new_node;
    return 0;
}

node_t *listdel(list_t *l, node_t *n) {
    if (!n) return NULL;
    if (n == l->head) {
        l->head = n->next;
        if (l->head) l->head->prev = NULL;
        else         l->tail = NULL;
    } else if (n == l->tail) {
        l->tail = n->prev;
        if (l->tail) l->tail->next = NULL;
        else         l->head = NULL;
    } else {
        n->prev->next = n->next;
        n->next->prev = n->prev;
    }
    n->prev = n->next = NULL;
    return n;
}

void listfree(list_t *l) { // frees only nodes, not their data
    node_t *n;
    while (l->tail) {
        n = l->tail;
        l->tail = n->prev;
        free(n);
    }
    l->head = NULL;
}

// OPERATIONS ON PLASMID

void initPlasmid(void) {
    listinit(&plasmid);
}

int addGene(const char *name,
            const char *sequence,
            const char *function,
            unsigned int position)
{
    gene_t *g = malloc(sizeof(gene_t));
    if (!g) return -1;
    g->name     = strdup_or_null(name); // copy strings so they can be safely freed later
    g->sequence = strdup_or_null(sequence);
    g->function = strdup_or_null(function);
    if (!g->name || !g->sequence || !g->function) {
        free(g->name); free(g->sequence); free(g->function);
        free(g);
        return -1;
    }
    if (listinsert(&plasmid, g, position) != 0) {
        free(g->name); free(g->sequence); free(g->function);
        free(g);
        return -1;
    }
    return 0;
}

int deleteGene(char *name) {
    node_t *n;
    for (n = plasmid.head; n; n = n->next) {
        gene_t *g = (gene_t *)n->data; // n->data is void, so cast to gene_t
        if (strcmp(g->name, name) == 0) {
            node_t *removed = listdel(&plasmid, n);
            if (!removed) return -1;
            free(g->name); free(g->sequence); free(g->function);
            free(g);
            free(removed);
            return 0;
        }
    }
    return -1;
}

void printPlasmid(void) {
    node_t *n;
    int idx = 1;
    if (!plasmid.head) {
        printf("Plasmid is empty.\n");
        return;
    }
    printf("Plasmid contents:\n");
    for (n = plasmid.head; n; n = n->next, idx++) {
        gene_t *g = (gene_t *)n->data;
        printf(" %2d: %s\n", idx, g->name);
    }
}

void printGeneInfo(unsigned int position) {
    if (!plasmid.head) {
        printf("Plasmid is empty.\n");
        return;
    }
    node_t *n = plasmid.head;
    unsigned int idx = 1;
    while (n && idx < position) {
        n = n->next; idx++;
    }
    if (!n) {
        printf("No gene at position %u.\n", position);
        return;
    }
    gene_t *g = (gene_t *)n->data;
    printf("Gene information at position %u:\n", position);
    printf(" Name: %s\n", g->name);
    printf(" Sequence: %s\n", g->sequence);
    printf(" Function: %s\n", g->function);
}

char complement(char base) {
    switch(base) {
        case 'A': case 'a': return 'T';
        case 'T': case 't': return 'A';
        case 'C': case 'c': return 'G';
        case 'G': case 'g': return 'C';
        default:            return 'N';
    }
}

void designPCR(unsigned int position, unsigned int primer_length) {
    if (!plasmid.head) {
        printf("Plasmid is empty.\n");
        return;
    }
    node_t *n = plasmid.head;
    unsigned int idx = 1;
    while (n && idx < position) { n = n->next; idx++; }
    if (!n) { printf("No gene at position %u.\n", position); return; }
    gene_t *g = (gene_t *)n->data;
    size_t len = strlen(g->sequence);
    if (primer_length == 0 || primer_length > len) {
        printf("Invalid primer length (max %zu).\n", len);
        return;
    }
    char *forward = malloc(primer_length + 1);
    char *reverse = malloc(primer_length + 1);
    if (!forward || !reverse) { free(forward); free(reverse); printf("Memory allocation error.\n"); return; }
    strncpy(forward, g->sequence, primer_length);
    forward[primer_length] = '\0';
    for (unsigned int i = 0; i < primer_length; i++) {
        char base = g->sequence[len - 1 - i];
        reverse[i] = complement(base);
    }
    reverse[primer_length] = '\0';
    printf("PCR design for gene '%s' (position %u, length %u):\n", g->name, position, primer_length);
    printf(" Forward primer: %s\n", forward);
    printf(" Reverse primer: %s\n", reverse);
    free(forward); free(reverse);
}

int savePlasmidToFile(const char *filename) {
    FILE *f = fopen(filename, "w");
    if (!f) return -1;
    node_t *n;
    for (n = plasmid.head; n; n = n->next) {
        gene_t *g = (gene_t*)n->data;
        fputs(g->sequence, f);
    }
    fputc('\n', f);
    fprintf(f, "\n# LEGEND: %4s  %-10s  %6s  %6s\n", "No.", "Name", "Start", "End");
    size_t pos = 1;
    int idx = 1;
    for (n = plasmid.head; n; n = n->next, idx++) {
        gene_t *g = (gene_t*)n->data;
        size_t len = strlen(g->sequence);
        fprintf(f, "%4d  %-10s  %6zu  %6zu\n", idx, g->name, pos, pos + len - 1);
        pos += len;
    }
    fclose(f);
    return 0;
}

int loadPlasmidFromCSV(const char *filename) {
    FILE *f = fopen(filename, "r");
    if (!f) {
        printf("Cannot open CSV file: %s\n", filename);
        return -1;
    }
    char line[2048];
    unsigned int position = 1;
    while (fgets(line, sizeof(line), f)) {
        char *name = strtok(line, ",");
        char *seq  = strtok(NULL, ",");
        char *func = strtok(NULL, "\n");
        if (!name || !seq || !func) {
            printf("CSV format error in line: %s\n", line);
            continue;
        }
        size_t len = strlen(func);
        if (func[len - 1] == '\n') func[len - 1] = '\0';
        if (addGene(name, seq, func, position++) != 0) {
            printf("Failed to add gene: %s\n", name);
        }
    }

    fclose(f);
    return 0;
}

int editGene(unsigned int position, const char *new_name, const char *new_seq, const char *new_func) {
    node_t *n = plasmid.head;
    unsigned int idx = 1;
    while (n && idx < position) { n = n->next; idx++; }
    if (!n) return -1;
    gene_t *g = (gene_t *)n->data;
    char *name_copy = strdup_or_null(new_name);
    char *seq_copy  = strdup_or_null(new_seq);
    char *func_copy = strdup_or_null(new_func);
    if (!name_copy || !seq_copy || !func_copy) {
        free(name_copy); free(seq_copy); free(func_copy);
        return -1;
    }
    free(g->name);
    free(g->sequence);
    free(g->function);
    g->name = name_copy;
    g->sequence = seq_copy;
    g->function = func_copy;
    return 0;
}

int main(void) {
    initPlasmid();
    printf("Creating a new empty plasmid...\n");
    int choice;
    do {
        printf("\n--- Plasmid Manager Menu ---\n");
        printf("1) Add gene at position\n");
        printf("2) Delete gene at position\n");
        printf("3) Print plasmid contents\n");
        printf("4) Print gene details from position\n");
        printf("5) Design PCR primers for gene\n");
        printf("6) Save to file\n");
        printf("7) Load plasmid from CSV file (name,sequence,function)\n");
        printf("8) Edit gene data\n");
        printf("0) Exit\n");
        printf("Choice: ");
        if (scanf("%d", &choice) != 1) break;
        switch (choice) {
            case 1: {
                char name[64], seq[1024], func[256];
                unsigned int pos;
                printf("Gene name: "); scanf("%63s", name);
                printf("Sequence (5'-3'): "); scanf("%1023s", seq);
                printf("Function: "); scanf("%255s", func);
                printf("Position in plasmid: "); scanf("%u", &pos);
                if (addGene(name, seq, func, pos) == 0)
                    printf("Gene added.\n");
                else
                    printf("Error adding gene.\n");
                break;
            }
            case 2: {
                unsigned int pos;
                printf("Gene position to delete: "); scanf("%u", &pos);
                node_t *n = plasmid.head;
                unsigned int idx = 1;
                while (n && idx < pos) { n = n->next; idx++; }
                if (!n) {
                    printf("No gene at position %u.\n", pos);
                } else {
                    gene_t *g = (gene_t*)n->data;
                    char name_buf[64];
                    strncpy(name_buf, g->name, 63); name_buf[63]='\0';
                    if (deleteGene(name_buf)==0)
                        printf("Gene '%s' deleted.\n", name_buf);
                    else
                        printf("Error deleting gene.\n");
                }
                break;
            }
            case 3:
                printPlasmid();
                break;
            case 4: {
                unsigned int pos;
                printf("Gene position: "); scanf("%u", &pos);
                printGeneInfo(pos);
                break;
            }
            case 5: {
                unsigned int pos, length;
                printf("Gene position: "); scanf("%u", &pos);
                printf("Primer length: "); scanf("%u", &length);
                designPCR(pos, length);
                break;
            }
            case 6: {
                char fname[256];
                printf("File name (format .txt): "); scanf("%255s", fname);
                if (savePlasmidToFile(fname) == 0)
                    printf("Plasmid saved to %s\n", fname);
                else
                    printf("Error saving!\n");
                break;
            }
            case 7: {
                char fname[256];
                printf("CSV file name to load: "); scanf("%255s", fname);
                if (loadPlasmidFromCSV(fname) == 0)
                    printf("Plasmid loaded from file %s\n", fname);
                else
                   printf("Error loading file.\n");
                break;
            }
            case 8: {
                unsigned int pos;
                char new_name[64], new_seq[1024], new_func[256];
                printf("Gene position to edit: "); scanf("%u", &pos);
                printf("New name: "); scanf("%63s", new_name);
                printf("New sequence: "); scanf("%1023s", new_seq);
                printf("New function: "); scanf("%255s", new_func);
                if (editGene(pos, new_name, new_seq, new_func) == 0)
                    printf("Gene at position %u updated.\n", pos);
                else
                    printf("Error editing gene.\n");
                break;
            }
            case 0:
                printf("Exiting program.\n");
                break;
            default:
                printf("Invalid choice.\n");
        }
    } while (choice != 0);
    listfree(&plasmid);
    return 0;
}