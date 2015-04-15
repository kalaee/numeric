#include <stdio.h>

int main()
{
	FILE * f = fopen("out.A.txt","w");
	fprintf(f,"Hello, World!\n");
	fclose(f);

	return 0;
}
