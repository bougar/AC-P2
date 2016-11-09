#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "mpi.h"

/*La siguiente función inicializa la matriz asignando valores aleatorios*/
float * getRamdomMatrix(int rows, int columns){
    int i, j;
    srand(time(NULL));
    float * matrix = (float *) malloc(rows * columns * sizeof(int));
    for (i=0; i < rows; i++)
        for (j=0; j < columns; j++)
            matrix[i*columns+j] = (float) rand() + ((float)rand()/(float)RAND_MAX); 
    return matrix;
}

/*La siguiente función inicializa la matriz, siguiendo el ejemplo de matrizVector*/
float * getMatrix(int rows, int columns){
    int i, j;
    float * matrix = (float *) malloc(rows * columns * sizeof(float));
    for (i=0; i < rows; i++)
        for (j=0; j < columns; j++)
            matrix[i*columns+j] = (float) i*columns+j+1;
    return matrix;
}

int main(int argc, char * argv[]){
    int rank, numprocs, n, k, m, alfa;
	float * matrixA = NULL;
	float * matrixB = NULL;

    MPI_Init(&argc, &argv);
    // Determinar el rango del proceso invocado
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Determinar el numero de procesos
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	if (!rank){	
		if ( argc < 5 ){
			fprintf(stderr, "Número de parámetros incorrecto\n");
			MPI_Finalize();
		}

		m = atoi(argv[1]);
		k = atoi(argv[2]);
		n = atoi(argv[3]);
		alfa = atoi(argv[4]);

		if ( n < 1 || k < 1 || m < 1 ){
			fprintf(stderr, "Algún parámetro no tiene un valor aceptable\n");
			MPI_Finalize();
		}

		float * matrixA = getMatrix(m, k);
		float * matrixB = getMatrix(k, n);
		
		if ( matrixA == NULL || matrixB == NULL ){
			fprintf(stderr, "No se ha podido obtener alguna de las matrices.\n");
			MPI_Finalize();
		}
	}

	// Replica matrixB en todos los procesos
    MPI_Bcast(matrixB, k*n, MPI_FLOAT, 0, MPI_COMM_WORLD);

	if (!rank)
		free(matrixA);
	free(matrixB);
	MPI_Finalize();
}
