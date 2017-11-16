#include <stdio.h>
#include <stdlib.h>

#include "s1a.h"

int main(int c, char *v[])
{
	if (c != 3) return fprintf(stderr, "usage:\n\t%s raw.dat n\n", *v);
	//                                             0 1       2
	char *filename_x = v[1];
	int n = atoi(v[2]);

	struct s1a_file x[1];

	s1a_load_whole_datafile_trunc(x, filename_x, n+1);
	printf("%s: %d records\n", filename_x, x ->n);

	s1a_decode_line(NULL, x->t + n);

	return 0;
}
