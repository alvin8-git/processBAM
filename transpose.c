/*
 * transpose.c - Transpose tab-separated data
 * Usage: transpose -t < input.tsv > output.tsv
 *        transpose < input.tsv > output.tsv (pass-through)
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_LINES 50000
#define MAX_LINE_LEN 500000

int main(int argc, char *argv[]) {
    int transpose_mode = 0;
    if (argc > 1 && strcmp(argv[1], "-t") == 0) {
        transpose_mode = 1;
    }

    char *lines[MAX_LINES];
    int nlines = 0;
    char buf[MAX_LINE_LEN];

    while (fgets(buf, sizeof(buf), stdin) && nlines < MAX_LINES) {
        lines[nlines++] = strdup(buf);
    }

    if (!transpose_mode) {
        for (int i = 0; i < nlines; i++) {
            printf("%s", lines[i]);
        }
        for (int i = 0; i < nlines; i++) {
            free(lines[i]);
        }
        return 0;
    }

    /* Count columns from first line */
    int ncols = 0;
    char *p = lines[0];
    while (*p) {
        if (*p == '\t') ncols++;
        p++;
    }
    ncols++;

    /* Transpose */
    for (int c = 0; c < ncols; c++) {
        for (int r = 0; r < nlines; r++) {
            char *s = lines[r];
            char *e;
            int col = 0;

            /* Find start of column c */
            while (col < c && (e = strchr(s, '\t'))) {
                s = e + 1;
                col++;
            }

            /* Find end of this field */
            e = strchr(s, '\t');
            if (!e) e = strchr(s, '\n');
            if (!e) e = s + strlen(s);

            /* Print field */
            printf("%.*s%c", (int)(e - s), s, r < nlines - 1 ? '\t' : '\n');
        }
    }

    for (int i = 0; i < nlines; i++) {
        free(lines[i]);
    }
    return 0;
}
