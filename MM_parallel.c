#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <memory.h>
#include <omp.h>
#include "setting.h"

double diff_in_second(struct timespec t1, struct timespec t2)
{
    struct timespec diff;
    if (t2.tv_nsec - t1.tv_nsec < 0) {
        diff.tv_sec  = t2.tv_sec - t1.tv_sec - 1;
        diff.tv_nsec = t2.tv_nsec - t1.tv_nsec + 1000000000;
    } else {
        diff.tv_sec  = t2.tv_sec - t1.tv_sec;
        diff.tv_nsec = t2.tv_nsec - t1.tv_nsec;
    }

    return (diff.tv_sec + diff.tv_nsec / 1000000000.0);
}

int main (int argc, char **argv) {

    int num;
    FILE *fp = NULL;
    int i, j, k;
    int a_row, a_col, b_row, b_col;
    int *a, *b, *c;
    char a_buff_row[255], a_buff_col[255], b_buff_row[255], b_buff_col[255], buff[255];
    struct timespec start, end;
    memset(a_buff_row, 0, 255);
    memset(a_buff_col, 0, 255);
    memset(b_buff_row, 0, 255);
    memset(b_buff_col, 0, 255);
    memset(buff, 0, 255);
    
    /* input file : default or cmd arguments */
    fp = (argc == 1) ? fopen("input.txt", "r") : fopen(argv[1], "r");

    /* scan matrix a */
    int status = fscanf(fp,"%s", a_buff_row);
    a_row = atoi(a_buff_row);
    status = fscanf(fp,"%s", a_buff_col);
    a_col = atoi(a_buff_col);
    
    a = (int*)malloc(sizeof(int*) * a_row * a_col);
    
    for(i = 0 ; i < a_row ; i++){
        for(j = 0 ; j < a_col ; j++){
            status = fscanf(fp,"%s", buff);
            num = atoi(buff);
            a[i * a_row + j] = num;
            memset(buff, 0, 255);
        }
    }

    /* scanf matrix b */
    status = fscanf(fp,"%s", b_buff_row);
    b_row = atoi(b_buff_row);
    status = fscanf(fp,"%s", b_buff_col);
    b_col = atoi(b_buff_col);
    
    b = (int*)malloc(sizeof(int*) * b_row * b_col);
    
    for(i = 0 ; i < b_row ; i++){
        for(j = 0 ; j < b_col ; j++){
            status = fscanf(fp,"%s", buff);
            num = atoi(buff);
            b[j * b_row + i] = num;
            memset(buff, 0, 255);
        }
    }
    fclose(fp);

    /* create matrix c to store result */
    c = (int*)malloc(sizeof(int*) * a_row * b_col);

    for(i = 0 ; i < a_row ; i++){
        for(j = 0  ; j < b_col ; j++){
            c[i * a_row + j]=0;
        }
    }
    int cij = 0;
    int ii, jj;
    /* matrix multiply */ 
    clock_gettime(CLOCK_REALTIME, &start);
    
    #pragma omp parallel for private(jj,i,j,k) num_threads(OMP)
    for(ii = 0 ; ii < a_row ; ii+=tile_size)
        for(jj = 0; jj < b_col ; jj+=tile_size)    
            for (i = ii ; i < ii + tile_size ;++i) 
            {
                for (j = jj ; j < jj + tile_size ;++j){
                    cij = c[i * a_row + j];
                    for (k = 0 ; k < a_col ; ++k) 
                        cij +=a[i * a_row + k]*b[j * b_col + k];
                    c[i * a_row + j] = cij; 
                }       
            }
    //printf("size of int: %d", sizeof(int));
    clock_gettime(CLOCK_REALTIME, &end);
    printf("MM_parallel\t  : %f sec\n",diff_in_second(start, end));
    /* export file : parallel.txt */
    FILE *fp_out = NULL;
    fp_out = fopen("MM_parallel.txt", "w");
    for(i = 0 ; i < a_row ; i++){
        for(j = 0 ; j < b_col ; j++){
            fprintf(fp_out,"%d ",c[i * a_row + j]);
        }
        fprintf(fp_out,"\n");
    }
    fclose(fp_out);

    /* free memory */
    free(a);
    free(b);    

}