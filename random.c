#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void outputMatrix(FILE *fp, int size)
{
    int i, j;

    fprintf(fp, "%d %d\n", size, size);
    for (i = 0; i < size; i++)
    {
        for (j = 0; j < size; j++)
        {
            fprintf(fp, "%d ", rand() % 100);
        }
        fprintf(fp, "\n");
    }
}

int main(int argc, char **argv)
{

    if (argc != 1 && argc != 2)
    {
        printf("usage: %s <size>\n", argv[0]);
        printf("or     %s for default size 64\n", argv[0]);
        return 1;
    }

    srand(time(NULL));
    int size = (argc == 1) ? 64 : atoi(argv[1]); /* default size is 64 */
    FILE *fp = fopen("input.txt", "w");

    outputMatrix(fp, size);
    outputMatrix(fp, size);

    fclose(fp);
    return 0;
}
