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

void add (int size, int *a, int *b, int *c) 
{
	for (int i = 0 ; i < size *size ; i++)
		c[i] = a[i] + b[i];		
}

void sub (int size, int *a, int *b, int *c)
{
	for (int i = 0 ; i < size *size ; i++)
		c[i] = a[i] - b[i];
}
void MM (int size, int *a, int *b, int *c)
{
	int i, j, k, ii, jj;
	int cij;
    for(ii = 0 ; ii < size ; ii += tile_size)
        for(jj = 0; jj < size ; jj += tile_size)    
            for (i = ii ; i < ii + tile_size ;++i) 
            {
                for (j = jj ; j < jj + tile_size ;++j){
                    cij = c[i * size + j];
                    for (k = 0 ; k < size ; ++k) 
                        cij +=a[i * size + k]*b[j + size * k];
                    c[i * size + j] = cij; 
                }       
            }
}
int* createMatrix(int size)
{
	int* matrix;
	matrix = (int *)malloc(size * size * sizeof(int));
	return matrix;
}

void freeMatrix(int *matrix, int n)
{
	free(matrix);
}

void strassen (int size, int *a, int *b, int *c)
{
	int *m1, *m2, *m3, *m4, *m5, *m6, *m7, *a11, *a12, *a21, *a22,
	    *b11, *b12, *b21, *b22, *c11, *c12, *c21, *c22, *temp1, *temp2,
		*temp3, *temp4, *temp5, *temp6, *temp7, *temp8, *temp9, *temp10;
	int newSize;

	/* alloc memory */
	newSize = size/2;
	temp1 = createMatrix(newSize);
	temp2 = createMatrix(newSize);
	temp3 = createMatrix(newSize);
	temp4 = createMatrix(newSize);
	temp5 = createMatrix(newSize);
	temp6 = createMatrix(newSize);
	temp7 = createMatrix(newSize);
	temp8 = createMatrix(newSize);
	temp9 = createMatrix(newSize);
	temp10 = createMatrix(newSize);
	m1 = createMatrix(newSize);
	m2 = createMatrix(newSize);
	m3 = createMatrix(newSize);
	m4 = createMatrix(newSize);
	m5 = createMatrix(newSize);
	m6 = createMatrix(newSize);
	m7 = createMatrix(newSize);
	a11 = createMatrix(newSize);
	a12 = createMatrix(newSize);
	a21 = createMatrix(newSize);
	a22 = createMatrix(newSize);
	b11 = createMatrix(newSize);
	b12 = createMatrix(newSize);
	b21 = createMatrix(newSize);
	b22 = createMatrix(newSize);
	c11 = createMatrix(newSize);
	c12 = createMatrix(newSize);
	c21 = createMatrix(newSize);
	c22 = createMatrix(newSize);

	/* seperate matrix */
	for (int i = 0 ; i < newSize; i++) 
		for (int j = 0 ; j < newSize; j++){
			a11[i * newSize + j]=a[i * size + j];
			b11[i * newSize + j]=b[i * size + j];
			a12[i * newSize + j]=a[newSize + i * size + j];
			b12[i * newSize + j]=b[newSize + i * size + j];
			a21[i * newSize + j]=a[newSize * size + i * size + j];
			b21[i * newSize + j]=b[newSize * size + i * size + j];
			a22[i * newSize + j]=a[newSize + newSize * size + i * size + j];
			b22[i * newSize + j]=b[newSize + newSize * size + i * size + j];	
		}

	#pragma omp parallel sections  num_threads(OMP)
	{
		#pragma omp section 
		{
			//printf( "%d", omp_get_thread_num());
			/* M1 = (A11 + A22) (B11 + B22) */
			add(newSize, a11, a22, temp1);
			add(newSize, b11, b22, temp2);
			MM(newSize, temp1, temp2, m1);
	
			/* M2 = (A21 + A22)B11 */
			add(newSize, a21, a22, temp3);
			MM(newSize, temp3, b11, m2);
		
			//printf( "%d", omp_get_thread_num());
			/* M3 = A11(B12 − B22) */
			sub(newSize, b12, b22, temp4);
			MM(newSize, a11, temp4, m3);
			/* M4 = A22(B21 − B11) */
			sub(newSize, b21, b11, temp5);
			MM(newSize, a22, temp5, m4);
		}
		#pragma omp section
		{
			//printf( "%d", omp_get_thread_num());
			/* M5 = (A11 + A12)B22 */
			add(newSize, a11, a12, temp6);
			MM(newSize, temp6, b22, m5);
		
			/* M6 = (A21 − A11)(B11 + B12) */
			sub(newSize, a21, a11, temp7);
			add(newSize, b11, b12, temp8);
			MM(newSize, temp7, temp8, m6);
		
			//printf( "%d", omp_get_thread_num());
			/* M7 = (A12 − A22)(B21 + B22) */
			sub(newSize, a12, a22, temp9);
			add(newSize, b21, b22, temp10);
			MM(newSize, temp9, temp10, m7);
		}
	}
	/* C11 = M1 + M4 − M5 + M7 */
	add(newSize, m1, m4, temp1);
	sub(newSize, m7, m5, temp2);
	add(newSize, temp1, temp2, c11);
	/* C12 = M3 + M5 */
	add(newSize, m3, m5, c12);
	/* C21 = M2 + M4 */
	add(newSize, m2, m4, c21);
	/* C22 = M1 − M2 + M3 + M6 */
	sub(newSize, m1, m2, temp1);
	add(newSize, m3, m6, temp2);
	add(newSize, temp1, temp2, c22);

	/* combine matrix */
	for (int i = 0 ; i < newSize ; i++) {
		for (int j = 0 ; j < newSize ; j++) {
			c[i * size + j] = c11[i * newSize + j];
			c[newSize + i * size + j] = c12[i * newSize + j];
			c[newSize * size + i * size + j] = c21[i * newSize + j];
			c[newSize + newSize * size + i * size + j] = c22[i * newSize + j];
		}
	}

	/* free memory */
	freeMatrix(temp1, newSize);
	freeMatrix(temp2, newSize);
	freeMatrix(m1, newSize);
	freeMatrix(m2, newSize);
	freeMatrix(m3, newSize);
	freeMatrix(m4, newSize);
	freeMatrix(m5, newSize);
	freeMatrix(m6, newSize);
	freeMatrix(m7, newSize);
	freeMatrix(a11, newSize);
	freeMatrix(a12, newSize);
	freeMatrix(a21, newSize);
	freeMatrix(a22, newSize);
	freeMatrix(b11, newSize);
	freeMatrix(b12, newSize);
	freeMatrix(b21, newSize);
	freeMatrix(b22, newSize);
	freeMatrix(c11, newSize);
	freeMatrix(c12, newSize);
	freeMatrix(c21, newSize);
	freeMatrix(c22, newSize);
}


int main (int argc, char **argv)
{
	srand(time(NULL));
	int *a, *b, *c;
	int i, j, k , input_size;
	int num,tmp = 2;
	struct timespec start, end;
	FILE *fp = fopen("input.txt", "r");
	char a_buff_row[7];
	char a_buff_col[7];
	char b_buff_row[7];
	char b_buff_col[7];
	char buff[7];
	memset(a_buff_row, 0, 7);
	memset(a_buff_col, 0, 7);
	memset(b_buff_row, 0, 7);
	memset(b_buff_col, 0, 7);
	memset(buff, 0, 7);

	int status = fscanf(fp,"%s", a_buff_row);
	int a_row = atoi(a_buff_row);
	status = fscanf(fp,"%s", a_buff_col);
	int a_col = atoi(a_buff_col);
	if(a_row != a_col) {
		printf("error a input!!\n");
		return 0;
	}
	
	/* make size of matrix (input_size) = 2^n */
	input_size = a_row;
	while(input_size > tmp)
		tmp *= 2;
	input_size = tmp;
	
	a = createMatrix(input_size);
	for(i = 0 ; i < tmp ; i++) {
		for(j = 0 ; j < tmp ; j++) {
			if((i < a_row) && (j < a_col)){
				status = fscanf(fp,"%s",buff);
				num = atoi(buff);
				a[i * tmp +  j] = num;
				memset(buff, 0, 7);
			}
			else 
			{
				a[i * tmp +  j] = 0;
			}
			
		}
	}

	status = fscanf(fp,"%s", b_buff_row);
	int b_row = atoi(b_buff_row);
	status = fscanf(fp,"%s", b_buff_col);
	int b_col = atoi(b_buff_col);
	if(b_row != b_col || b_row != a_row || b_row != a_col) {
		printf("error b input!!\n");
		return 0;
	}

	b = createMatrix(input_size);
	for(i = 0 ; i < tmp ; i++) {
		for(j = 0 ; j < tmp ; j++) {
			if(i < b_row && j < b_col){
				status = fscanf(fp,"%s",buff);
				num = atoi(buff);
				b[i * tmp +  j] = num;
				memset(buff, 0, 7);
			}
			else
			{
				b[i * tmp +  j] = 0;	
			}
			
		}
	}
	fclose(fp);     /* input file */

	/* create matrix c */
	c = createMatrix(input_size);
	
	/* matrix multiplication */
	clock_gettime(CLOCK_REALTIME, &start);
	strassen(input_size, a, b, c);
	clock_gettime(CLOCK_REALTIME, &end);
	printf("hybrid strassen parallel section   : %f sec\n",diff_in_second(start, end));

	FILE *fp_out = fopen("strassen_hybrid_parallel_section.txt", "w");
	for(i = 0 ; i < tmp ; i++) {
		for(j = 0 ; j < tmp ; j++) {
			if(i < a_row && j < a_col)
				fprintf(fp_out,"%d ",c[i * tmp +  j]);
			else
				break;
		}
		if(i < a_row)
			fprintf(fp_out,"\n");
	}
	fclose(fp_out);
}