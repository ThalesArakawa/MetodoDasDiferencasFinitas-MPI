#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>
#include <string.h>
#define N 600

struct quad
{
	int cc;
	double a;
	double d;
};

void importCC(struct quad A[N][N])//Importa valores de Contorno
{
	int i,j;
	double aux;
	FILE *im = fopen("uff.dat","r");
  int c;



	for(i=0;i<=N;i++)			
	{
		for(j=0;j<=N;j++)
		{
			c=fgetc(im);
			if(i!=N&&j!=N){
			
				if(!feof(im))
				{
					if(c=='A')//Logo Carregada
					{
						A[i][j].cc=1;
						A[i][j].a=100;
					}
	      				else if(c=='B'){
		 				A[i][j].cc=0;
						A[i][j].a=0;
					}
				}
				//Aterrando laterais
				if(j==0)
				{
				  	A[i][j].a=0;//terra
					A[i][j].cc=1;
				}
				if(j==N-1)
				{
					A[i][j].a=0;//terra
					A[i][j].cc=1;
				}
				if(i==N-1)
				{
					A[i][j].a=0;//terra
					A[i][j].cc=1;
				}
				if(i==0)
				{
					A[i][j].a=0;//terra
					A[i][j].cc=1;
	      			}
			
			}		
		}
	}
	fclose(im);
}
void printMap(struct quad map[N][N], char c, int my_rank)//Imprimir Mapa
{
	int i,j;
	FILE *fp;

	if(c=='c')	//IMPRIME CONDIÇÃO DE CONTORNO
	{	
		fp=fopen("cc.dat","w");

		for(i=0;i<=N;i++)			
		{
			if(i!=N){
				for(j=0;j<=N;j++)
				{
					if(j==N){
						fprintf(fp,"\n");
					}
					else{
						fprintf(fp,"%f\t",map[j][i].a);
					}
				}
			}	
		}
	}
	else	//IMPRIME VALORES DA MALHA
	{	
		fp=fopen("out.dat","w");

		for(i=0;i<=N;i++)			
		{
			if(i!=N){
				for(j=0;j<=N;j++)
				{
					if(j==N){
						fprintf(fp,"\n");
					}
					else{
					fprintf(fp,"%f\t",map[j][i].a);
					}
				
				}
			}
			
		}
	}
	fclose(fp);
}
void opDif(struct quad A[][N],int nLin)//Operador Diferencial
{
	double erro=0.0;
  	int i,j;
  	
  	struct quad aux1[N], aux2[N];
  	
  	int my_rank,size;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Status status;
	
	///// CRIANDO O DADO DO TIPO DE STRUCT PARA O MPI /////////////////////////////////
		MPI_Datatype quad_type;
		int lengths[3] = { 1, 1, 1 };
	     
		// Calculate displacements
		// In C, by default padding can be inserted between fields. MPI_Get_address will allow
		// to get the address of each struct field and calculate the corresponding displacement
		// relative to that struct base address. The displacements thus calculated will therefore
		// include padding if any.
		MPI_Aint displacements[3];
		struct quad dummy_quad;
		MPI_Aint base_address;
		MPI_Get_address(&dummy_quad, &base_address);
		MPI_Get_address(&dummy_quad.cc, &displacements[0]);
		MPI_Get_address(&dummy_quad.a, &displacements[1]);
		MPI_Get_address(&dummy_quad.d, &displacements[2]);
		displacements[0] = MPI_Aint_diff(displacements[0], base_address);
		displacements[1] = MPI_Aint_diff(displacements[1], base_address);
		displacements[2] = MPI_Aint_diff(displacements[2], base_address);
	     
		MPI_Datatype types[3] = { MPI_INT, MPI_DOUBLE, MPI_DOUBLE };
		MPI_Type_create_struct(3, lengths, displacements, types, &quad_type);
		MPI_Type_commit(&quad_type);
		//////////////////////////////////////////////////////////////////////////////////
	
	if(my_rank!=0){
		MPI_Send(&A[0][0],N,quad_type, my_rank-1, 99, MPI_COMM_WORLD);
		MPI_Recv(&aux1[0], N, quad_type, my_rank-1, 99, MPI_COMM_WORLD, &status);
	}
	if(my_rank!=size-1)
	{
		MPI_Send(&A[nLin-1][0],N,quad_type, my_rank+1, 99, MPI_COMM_WORLD);
		MPI_Recv(&aux2[0], N, quad_type, my_rank+1, 99, MPI_COMM_WORLD, &status);
	}
	
	for(i=0;i<nLin;i++)			
	{
		for(j=0;j<N;j++)
		{
			if(A[i][j].cc==0)
			{
				
				if(i==0){
				
					A[i][j].d=0.25*(A[i][j+1].a+A[i][j-1].a+A[i+1][j].a+aux1[j].a);
				
				}
				else if(i==nLin-1){
				
					A[i][j].d=0.25*(A[i][j+1].a+A[i][j-1].a+aux2[j].a+A[i-1][j].a);
				
				}
				else{
					A[i][j].d=0.25*(A[i][j+1].a+A[i][j-1].a+A[i+1][j].a+A[i-1][j].a);
				}
			}
		}
	}
}
double att(struct quad A[][N],int nLin)//Operador Diferencial
{
	int i,j;
	double erro=0.0;

	for(i=0;i<nLin;i++)			
	{
		for(j=0;j<N;j++)
		{
			if(A[i][j].cc==0)
			{
				if(A[i][j].a!=A[i][j].d)
				{
					if(erro<fabs(A[i][j].d-A[i][j].a))
					{
						erro=fabs(A[i][j].d-A[i][j].a);
					}
          A[i][j].a=A[i][j].d;
				} 
			}
		}	
	}

  return erro;

}
int main(int argc, char *argv[])
{

	MPI_Init(&argc, &argv);
	int i,j;
     
        	///// CRIANDO O DADO DO TIPO STRUCT PARA O MPI ////////////////////////////////////////////////////////////
		MPI_Datatype quad_type;
		int lengths[3] = { 1, 1, 1 };
	     
		// Calculate displacements
		// In C, by default padding can be inserted between fields. MPI_Get_address will allow
		// to get the address of each struct field and calculate the corresponding displacement
		// relative to that struct base address. The displacements thus calculated will therefore
		// include padding if any.
		MPI_Aint displacements[3];
		struct quad dummy_quad;
		MPI_Aint base_address;
		MPI_Get_address(&dummy_quad, &base_address);
		MPI_Get_address(&dummy_quad.cc, &displacements[0]);
		MPI_Get_address(&dummy_quad.a, &displacements[1]);
		MPI_Get_address(&dummy_quad.d, &displacements[2]);
		displacements[0] = MPI_Aint_diff(displacements[0], base_address);
		displacements[1] = MPI_Aint_diff(displacements[1], base_address);
		displacements[2] = MPI_Aint_diff(displacements[2], base_address);
	     
		MPI_Datatype types[3] = { MPI_INT, MPI_DOUBLE, MPI_DOUBLE };
		MPI_Type_create_struct(3, lengths, displacements, types, &quad_type);
		MPI_Type_commit(&quad_type);
		///////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	// Get the number of processes
        int size;
        MPI_Comm_size(MPI_COMM_WORLD, &size);
	
	//Obtendo rank do trabalhador
	int my_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	
	//Determinando buffer/submatriz de cada trabalhador
	int trabalho,nLin;
	nLin=(my_rank < N%size) ? N/size + 1 : N/size;
	struct quad own_map[nLin][N];
	trabalho = nLin*N;
	printf("\n P %d nLin %d",my_rank,nLin);
	

	//Trabalho único para o mestre
	if(my_rank == 0){

		struct quad map[N][N];
		//Importa Condição de Contorno de um arquivo em ASCII
		importCC(map); 

		//Vetor com o número de elementos que cada processo possuirá
           	int counts[size];
           	//Vetor com os deslocamentos sobre o buffer a ser aplicado o Scatterv
           	int displacements[size];
           	for(i=0;i<size;i++){
           	
           		counts[i]=(i < N%size) ? N*((int)(N/size) + 1) : N*(int)(N/size);
           		displacements[i]=(i>(N%size)) ? N*(i*(int)(N/size) + N%size) : N*i*((int)(N/size) + 1) ;
           		printf("\n i %d counts %d displa %d",i,counts[i],displacements[i]);
           	
           	}
 
		//Separando matriz por tamanhos diferentes
		MPI_Scatterv(map, counts, displacements, quad_type, own_map, trabalho, quad_type, 0, MPI_COMM_WORLD);

	}
	else
	{
		MPI_Scatterv(NULL, NULL, NULL, quad_type, own_map, trabalho, quad_type, 0, MPI_COMM_WORLD);
	}


	
	double erro=0.0;
	double tol=1e-5;
	double reduction_result;
	
	do{
		opDif(own_map,nLin);
		erro=att(own_map,nLin);
		MPI_Allreduce(&erro, &reduction_result, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	}while(tol<reduction_result);
	
	
	

	//Atribuindo a tarefa de unir a submatriz de cada trabalhador para a thread master
	if(my_rank == 0){
		struct quad B[N][N];
		//Vetor com o número de elementos que cada processo possuirá
           	int counts[size];
           	//Vetor com os deslocamentos sobre o buffer a ser aplicado o Scatterv
           	int displacements[size];
		for(i=0;i<size;i++){
           	
           		counts[i]=(i < N%size) ? N*((int)(N/size) + 1) : N*(int)(N/size);
           		displacements[i]=(i>(N%size)) ? N*(i*(int)(N/size) + N%size) : N*i*((int)(N/size) + 1) ;
           		printf("\n i %d counts %d displa %d",i,counts[i],displacements[i]);
           	
           	}
		MPI_Gatherv(own_map, trabalho, quad_type, B, counts, displacements, quad_type, 0, MPI_COMM_WORLD);
		printMap(B,'m',0);//Imprime a malha num arquivo out.dat
	}
	else
	{
		MPI_Gatherv(own_map, trabalho, quad_type, NULL, NULL, NULL, quad_type, 0, MPI_COMM_WORLD);
	}

	MPI_Finalize();

}

