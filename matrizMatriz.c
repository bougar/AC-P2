#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include "mpi.h"

/*
* Funtion: printMatrix
* --------------------
* Imprime una matriz por pantalla
*
* rows: Número de filas de la matriz
* columns: Número de columns de la matriz
* matrix: La matriz a imprimir
*
*/
void printMatrix(int rows, int columns, float * matrix){
    int i, j;
    for (i=0; i < rows; i++){
        printf("|");
        for (j=0; j < columns; j++)
            printf(" %f |", matrix[i*columns+j]);
        printf("\n");
    }
    printf("\n");
}

/*
* Funtion: getRandomMatrix 
* ------------------------
* Obtiene una matriz inicializada con números aleatorios 
*
* rows: Número de filas de la matriz
* columns: Número de columns de la matriz
*
* returns: Puntero a array de flotantes inicializado con valores
  aleatorios
*/
float * getRamdomMatrix(int rows, int columns){
    int i, j;

    srand(time(NULL));

    float * matrix = (float *) malloc(rows * columns * sizeof(int));

    if ( matrix == NULL ){
		fprintf(stderr, "Error de memoria\n");
		MPI_Finalize();
        exit(1);
    }

    for (i=0; i < rows; i++)
        for (j=0; j < columns; j++)
            matrix[i*columns+j] = (float) rand() + ((float)rand()/(float)RAND_MAX); 

    return matrix;
}

/*
* Funtion: getMatrix 
* ------------------------
* Obtiene una matriz inicializada de forma secuencial
*
* rows: Número de filas de la matriz
* columns: Número de columns de la matriz
*
* returns: Puntero a array de flotantes inicializado 
  de forma secuencial
*/
float * getMatrix(int rows, int columns){
    int i, j;
    float * matrix = (float *) malloc(rows * columns * sizeof(float));

    if ( matrix == NULL ){
		fprintf(stderr, "Error de memoria\n");
		MPI_Finalize();
        exit(1);
    }

    for (i=0; i < rows; i++)
        for (j=0; j < columns; j++)
            matrix[i*columns+j] = (float) i*columns+j+1;

    return matrix;
}

/*
* Funtion: getSendCounts
* ------------------------
* Calcula el reparto de elementos por proceso
*
* rows: Número de filas de la matriz
* columns: Número de columns de la matriz
* numProcs: Número de procesos
*
* returns: Puntero de enteros que guarda en cada
  posición 'k' cuantos elementos le tocan al 'k'
  proceso
*
*/
int * getSendCounts(int rows, int columns, int numProcs){
    int n, i, aux;
    int * sendCounts = (int *) malloc(numProcs * sizeof(int));

    if ( sendCounts == NULL ){
		fprintf(stderr, "Error de memoria\n");
        MPI_Abort(MPI_COMM_WORLD, 1);  
        exit(1);
    }

    aux = numProcs; 
    for (i=0;i < aux; i++){
        n = rows/numProcs;
        if ( (rows % numProcs) != 0 )
            n++;
        sendCounts[i]=n*columns;
        rows -= n;
        numProcs--;
    }
    return sendCounts;
}

/*
* Funtion: getRecvCounts
* ------------------------
* Calcula cuantos elementos le tiene que devolver cada
  proceso al proceso raíz
*
* columns: Número de columns de la matriz
* columnsB: Número de columnas de la segunda matriz
* sendCounts: Reparto de elementos por proceso
* numProcs: Número de procesos
*
* returns: Puntero de enteros que guarda en cada
  posición 'k' cuantos elementos tiene que entregar
  el proceso 'k' al proceso raiz
*
*/
int * getRecvCounts(int columns, int columnsB, int * sendCounts, int numProcs){
    int * recvCounts = (int *) malloc( numProcs * sizeof(int));
    int i;

    if ( recvCounts == NULL ){
		fprintf(stderr, "Error de memoria\n");
        MPI_Abort(MPI_COMM_WORLD, 1);  
        exit(1);
    }

    for (i=0; i < numProcs; i++)
        recvCounts[i] = (sendCounts[i] / columns) * columnsB;
    return recvCounts;
        
}

/*
* Funtion: getDispls
* ------------------------
* Obtiene los desplazamientos relativos de un array. 
*
* rows: Número de filas de la matriz
* columns: Número de columns de la matriz
* numProcs: Número de procesos
* Counts: Elementos por proceso
*
* returns: Puntero de enteros que guarda en cada
  posición 'k' cuantos elementos a partir de que 
  elemento se tiene que enviar en el scatterv o
  a partir de que elemnto se debe guardar en el 
  gatherv, para cada proceso
*
*/
int * getDispls(int rows, int columns, int numProcs, int  * Counts){
    int i, j;
    int * displs = (int *) malloc ( numProcs * sizeof(int));

    if ( displs == NULL ){
		fprintf(stderr, "Error de memoria\n");
        MPI_Abort(MPI_COMM_WORLD, 1);  
        exit(1);
    }

    for (j=0, i=0; i < rows && j < numProcs; j++, i = i + (Counts[j] / columns)){
        displs[j]=i*columns;
    }
    return displs;
}

/*
* Funtion: matrixProduct 
* ------------------------
* Computa la multiplicación de dos matrices
*
* rows: Número de filas de la matriz
* columns: Número de columns de la matriz
* columnsB: Número de columnas de la segunda matriz
* matrixA: La primera matriz
* matrixB: La segunda matriz
* alfa: El valor que multiplica a la matriz a
*
* returns: Un puntero de flotantes que represanta una matriz
  inicializada con el resultado de multiplicar matrixA por
  matrixB
*
*/
float *  matrixProduct(int rows, int columns, int columnsB, float * matrixA, float * matrixB, float alfa){
    int i,j,l;
    float * matrixResult = (float *) malloc ( rows * columnsB * sizeof(float) );
    memset(matrixResult, 0, rows * columnsB * sizeof(float));

    if ( matrixResult == NULL ){
		fprintf(stderr, "Error de memoria\n");
        MPI_Abort(MPI_COMM_WORLD, 1);  
        exit(1);
    }

    for (i=0; i<rows; i++) {
        for(j=0; j<columnsB; j++){
            for (l=0; l<columns; l++) {
                matrixResult[i*columnsB+j] += alfa*matrixA[i*columns+l]*matrixB[l*columnsB+j];
            }
        }
    }
    return matrixResult;
}

int main(int argc, char * argv[]){
    int rank, numprocs, n, k, m;
    float alfa;
    int * displs;
    int * sendCounts;
    int * recvCounts;
	float * matrixA = NULL;
	float * matrixB = NULL;
	float * matrixC = NULL;
    float * matrixResult = NULL;

    MPI_Init(&argc, &argv);
    // Determinar el rango del proceso invocado
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Determinar el numero de procesos

    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	if (!rank){	
		if ( argc < 5 ){
			fprintf(stderr, "Número de parámetros incorrecto\n");
            MPI_Abort(MPI_COMM_WORLD, 1);  
            exit(1);
		}
    
		m = atoi(argv[1]);
		k = atoi(argv[2]);
		n = atoi(argv[3]);
		alfa = atof(argv[4]);

		if ( n < 1 || k < 1 || m < 1 ){
			fprintf(stderr, "Las dimensiones de las matrices no pueden ser menores a 0\n");
            MPI_Abort(MPI_COMM_WORLD, 1);  
            exit(1);
		}
		matrixA = getMatrix(m, k);
		matrixB = getMatrix(k, n);
        matrixC = (float *) malloc(m * n * sizeof(float));

		if ( matrixC == NULL ){
			fprintf(stderr, "Error de memoria\n");
            MPI_Abort(MPI_COMM_WORLD, 1);  
            exit(1);
		}
	}

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&k, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&alfa, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    sendCounts = getSendCounts(m, k, numprocs);
    displs = getDispls(m, k, numprocs, sendCounts);

    if (rank){
        matrixA = (float *) malloc (sendCounts[rank] * k * sizeof(float));
        matrixB = (float *) malloc (k * n * sizeof(float));

        if ( matrixB == NULL || matrixA == NULL ){
			fprintf(stderr, "Error de memoria\n");
            MPI_Abort(MPI_COMM_WORLD, 1);  
            exit(1);
        }
    }

    MPI_Scatterv(matrixA, sendCounts, displs, MPI_INT, matrixA, sendCounts[rank], MPI_INT, 0, MPI_COMM_WORLD);

	// Replica matrixB en todos los procesos
    MPI_Bcast(matrixB, k*n, MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    matrixResult = matrixProduct(sendCounts[rank] / k, k, n, matrixA, matrixB, alfa);

    recvCounts = getRecvCounts(k, n, sendCounts, numprocs);
    free(displs);
    displs = getDispls(m, n, numprocs, recvCounts);
    
    MPI_Gatherv(matrixResult, recvCounts[rank], MPI_FLOAT, matrixC, recvCounts, displs, MPI_FLOAT, 0, MPI_COMM_WORLD);


    if (!rank){
        free(matrixC);
        matrixC = matrixProduct(m,k,n, matrixA, matrixB, alfa);
        free(matrixC);
    }

    free(sendCounts);
	free(matrixA);
	free(matrixB);
    free(recvCounts);
    free(matrixResult);
    free(displs);
	MPI_Finalize();
}
