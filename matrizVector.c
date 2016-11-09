/* Programa que multiplica una matrix por un vector de forma y=alfa*A*x
   La matriz es de tamaño n*n y los vectores de longitud n
   n es múltiplo del número total de procesos
   Todos los datos están al inicio y al final almacenados en el procesador 0
   La matriz se distribuye por filas de forma bloque
   El vector x se replica en todos los procesos
   El vector y se encuentra distribuido entre los procesos y se recoge al final en el 0*/

#include <stdlib.h>
#include <stdio.h>
#include "mpi.h"

/* Función que multiplica una matriz por un vector en la memoria local */
void matrizVectorLocal(float alfa, float *A, int lda, float *x, float *y, int filas, int cols){
    int i, j;
    for (i=0; i<filas; i++) {
        for (j=0; j<cols; j++) {
            y[i] += alfa*A[i*lda+j]*x[j];
        }
    }
    return;
}


int main(int argc, char *argv[]) {

    int numprocs, rank, i, j, n, test;
    float alfa;

    MPI_Init(&argc, &argv);
    // Determinar el rango del proceso invocado
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    // Determinar el numero de procesos
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

    // Proceso 0 lee parámetros de entrada
    // Parámetro 1 -> n
    // Parámetro 2 -> alfa
    // Parámetro 3 -> booleano que nos indica si se desea imprimir matrices y vectores de entrada y salida
    if(!rank){
        if(argc>3){
            n = atoi(argv[1]);
            alfa = atof(argv[2]);
            test = atoi(argv[3]);
        }
        else{
            printf("NUMERO DE PARAMETROS INCORRECTO\n");
            MPI_Finalize();
            return 0;
        }

        if((n<=0) || (n%numprocs)){
            printf("TAMAÑO NO VALIDO\n");
            MPI_Finalize();
            return 0;
        }
    }

    // Proceso 0 envía n a todos los demás con un broadcast
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    // Proceso 0 envía alfa a todos los demás con un broadcast
    MPI_Bcast(&alfa, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);

    float *x = (float *) malloc(n*sizeof(float));
    float *A;
    float *y;
    if(!rank){ // Proceso 0 inicializa la matriz y el vector
        A = (float *) malloc(n*n*sizeof(float));
        y = (float *) malloc(n*sizeof(float));

        for(i=0; i<n; i++){
            for(j=0; j<n; j++){
                A[i*n+j] = 1+i+j;
            }
        }

        for(i=0; i<n; i++){
            x[i] = (1+i);
        }

        if(test){
            printf("\nMatriz A es...\n");
            for(i=0; i<n; i++){
                for(j=0; j<n; j++){
                    printf("%f ", A[i*n+j]);
                }
                printf("\n");
            }

            printf("\nVector x es...\n");
            for(i=0; i<n; i++){
                printf("%f ", x[i]);
            }
            printf("\n");
        }
    }

    

    // Número de filas por procesador. Sólo válido si n múltiplo de numprocs
    int filasBloque = n/numprocs; 
    // Cada proceso reserva espacio para las filas que le corresponden
    float *localA = (float *) malloc(filasBloque*n*sizeof(float));
    float *localy = (float *) malloc(filasBloque*sizeof(float));

    // Todos los procesos inicializan su parte de y a 0
    for(i=0; i<filasBloque; i++){
        localy[i] = 0.0;
    }
 
    double t;
    MPI_Barrier(MPI_COMM_WORLD); // Barrera para garantizar una correcta medida de tiempo
    t = MPI_Wtime();

    // Se envía A desde el proceso 0 a los diferentes procesos dividiendo por filas con un scatter
    MPI_Scatter(A, filasBloque*n, MPI_FLOAT, localA, filasBloque*n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // Replica x en todos los procesos
    MPI_Bcast(x, n, MPI_FLOAT, 0, MPI_COMM_WORLD);

    // Multiplica la submatriz por el vector
    matrizVectorLocal(alfa, localA, n, x, localy, filasBloque, n);

    //Reduce los elementos de y en el proceso 0
    MPI_Gather(localy, filasBloque, MPI_FLOAT, y, filasBloque, MPI_FLOAT, 0, MPI_COMM_WORLD);
    t = MPI_Wtime()-t;

    if(!rank){
        if(test){
            printf("\nAl final vector y es...\n");
            for(i=0; i<n; i++){
                printf("%f ", y[i]);
            }
            printf("\n");

            // Solo el proceso 0 calcula el producto y compara los resultados del programa secuencial con el paralelo
            float *testy = (float *) malloc(n*sizeof(float));
            for(i=0; i<n; i++){
                testy[i] = 0.0;
            }

            matrizVectorLocal(alfa, A, n, x, testy, n, n);
            int errores = 0;
            for(i=0; i<n; i++){
                if(testy[i] != y[i]){
                    errores++;
                    printf("\n Error en la posicion %d porque %f != %f", i, y[i], testy[i]);
                }
            }
            printf("\n%d errores en el producto matriz vector con dimension %d\n", errores, n);
            free(testy);
        }

        free(A);
        free(y);
    }

    free(localA);
    free(x);
    free(localy);

    // Barrera para que no se mezcle la impresión del tiempo con la de los resultados
    fflush(NULL);
    MPI_Barrier(MPI_COMM_WORLD);
    printf("\nProducto matriz vector con dimension %d y %d procesos: Tiempo de ejecucion del proceso %d fue %lf\n", n, numprocs, rank, t);

    MPI_Finalize();
    return 0;
}
