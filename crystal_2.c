/* Programa inical para geração de estruturas cristalinas - IC 2020 - Mateus Corradini Lopes */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<malloc.h>
#include<string.h>

struct Tatomic{
	int id; //para designar a espécie
	char nome[50]; //nome da espécie atômica
	double coord[3]; //coordenadas do átomo
};

void main()
{
	//Váriaveis
	struct Tatomic *ptr; //ponteiro para alocar dinamicamente o array das espécies na célula unitária
	struct Tatomic *estrutura; //ponteiro pra estrutura já c
	int i,j,k,l,m; //contadores auxiliares
	int naps; //numero de átomos na célula unitária por espécie
 	int scx,scy,scz; //parâmetros para gerar a supercélula
	int natom, nesp, ntot; //numero de atomos na célula unitária, numero de espécies e número total de átomos na estrutura
	//double **runit; //ponteiro para a matriz da célula unitária que será alocado dinâmicamente
	//double **estrutura; //ponteiro pra matriz que vai alocar a estrutura cristalina que será alocada dinâmicamente
	double par[3],pva[3],pvb[3],pvc[3]; //vetores da base e  parâmetro de rede
	double raux[3]; //vetor auxiliar para translação
	
	//Arquivos
	FILE *fpe, *fpc;
	fpe = fopen("estrutura.dat","r"); //Arquivo de entrada com os parâmetros a serem lidos
	fpc = fopen("cristal_GaAs.xyz","w");	//Arquivo de saída com as posições dos átomos na célula unitária

	//leitura do arquivo.
	
	fscanf(fpe,"%d %d %d",&scx,&scy,&scz);
	fscanf(fpe,"%d",&natom);
	fscanf(fpe,"%d",&nesp);

	//Tendo o número de espécies posso alocar o array do tipo Tatom
	ptr = (struct Tatomic *) malloc(nesp * sizeof(struct Tatomic)); // isso criará um vetor de dimensão NESP do tipo Tatom!
	fscanf(fpe,"%lf %lf %lf",&par[0],&par[1],&par[2]);
	fscanf(fpe,"%lf %lf %lf",&pva[0],&pva[1],&pva[2]);
	fscanf(fpe,"%lf %lf %lf",&pvb[0],&pvb[1],&pvb[2]);
	fscanf(fpe,"%lf %lf %lf",&pvc[0],&pvc[1],&pvc[2]);
	
	for(i = 0; i < nesp; i++)
	{
	
		fscanf(fpe,"%s %d",(ptr+i)->nome, &naps);
		for(j = 0 ; j < naps ; j++) //entrando com as posições dos atómos, pressupondo coordenadas reduzidas, nas células unitárias alocando-as no array
		{
			fscanf(fpe,"%d %lf %lf %lf", &(ptr+j)->id, &(ptr+j)->coord[0], &(ptr+j)->coord[1], &(ptr+j)->coord[2]);
		}
	}	
	
	
	ntot = natom*scx*scy*scz;  //número total de átomos que estarão na estrutura final.


	//alocação dinamicamente da matriz das coordenadas da célula unitária
	//runit = (double **) calloc(natom, sizeof(double *)); //aloca linhas correspondente ao numero de atomos pos celula unitária	
	
	//for(k=0;k<natom;k++) //alocando as colunas correspondentes pra célula unitária
	//{
	//	runit[k] = (double *) calloc(3,sizeof(double));
	//}

	//estrutura = (double **) calloc(ntot, sizeof(double *)); //aloca linhas correspondente ao numero total de átomos que estarão na estrutura
	
	//for(k=0;k<ntot;k++) //alocando as colunas correspondentes para a estrutura
	//{
	//	estrutura[k] = (double *) calloc(3,sizeof(double));
	//}

	estrutura = (struct Tatomic *) malloc (ntot * sizeof(struct Tatomic));
	
	fprintf(fpc,"\t%d\n",ntot); //imprimindo o número total de átomos na estrutura no arquivo.
	
	for(i=0;i<natom;i++) //escrevendo no documento de saída a célula unitária em coordenadas ja convertidas para cartesianas e guardando no array estrutura
	{		
		(estrutura + i)->id = (ptr + i)->id;
		strcpy((estrutura + i)->nome , (ptr + i)->nome);
		
		for(j=0;j<3;j++)
		{
			(ptr+i)->coord[j] = (ptr+i)->coord[j]*pva[j]*par[j] + (ptr+i)->coord[j]*pvb[j]*par[j] + (ptr+i)->coord[j]*pvc[j]*par[j];
			(estrutura+i)->coord[j] = (ptr+i)->coord[j];
			
			if(j%3 == 0)
			{
				fprintf(fpc,"\n%s\t%lf ",(estrutura+i)->nome,(ptr+i)->coord[j]);
			}
			else
			{
				fprintf(fpc,"\t%lf ",(ptr+i)->coord[j]);
			}
		}
	}
	m=i; //salva o ultimo átomo alocado na matriz estrutura;


	//Transladando os átomos da célula unitária na primeira direção e escrevendo no arquivo
	for(i=0;i<natom;i++) 
	{
		
		for(j=1 ; j < scx ; j++)
		{
				raux[0] = (ptr+i)->coord[0] + j*(pva[0]*par[0]);
				raux[1] = (ptr+i)->coord[1] + j*(pva[1]*par[0]);
				raux[2] = (ptr+i)->coord[2] + j*(pva[2]*par[0]);
				(estrutura+m)->id = (ptr+i)->id;
				strcpy((estrutura+m)->nome , (ptr+i)->nome);
				(estrutura+m)->coord[0] = raux[0];
				(estrutura+m)->coord[1] = raux[1];
				(estrutura+m)->coord[2] = raux[2];
				m++;
				fprintf(fpc,"\n%s\t%lf\t%lf\t%lf ",(estrutura+i)->nome, raux[0],raux[1],raux[2]);	
		}
	}
	k=m; 
	
	//Transladando os átomos resultantes do último processo segunda direção e escrevendo no arquivo
	for(i=0 ; i < natom*scx ; i++)
	{
		for(j=1 ; j < scy ; j++)
		{
				raux[0] = (estrutura+i)->coord[0] + j*(pvb[0]*par[1]);
				raux[1] = (estrutura+i)->coord[1] + j*(pvb[1]*par[1]);
				raux[2] = (estrutura+i)->coord[2] + j*(pvb[2]*par[1]);
				(estrutura+k)->id = (estrutura+i)->id;
				strcpy((estrutura+k)->nome , (estrutura+i)->nome);
				(estrutura+k)->coord[0] = raux[0];
				(estrutura+k)->coord[1] = raux[1];
				(estrutura+k)->coord[2] = raux[2];
				k++;
				fprintf(fpc,"\n%s\t%lf\t%lf\t%lf ",(estrutura+i)->nome, raux[0],raux[1],raux[2]);	
		}
	}
	l=k;
	
	//idem para a terceira direção
	for(i=0 ; i < natom*scx*scy ; i++)
	{
		for(j=1 ; j < scz ; j++)
		{
				raux[0] = (estrutura+i)->coord[0] + j*(pvc[0]*par[2]);
				raux[1] = (estrutura+i)->coord[1] + j*(pvc[1]*par[2]);
				raux[2] = (estrutura+i)->coord[2] + j*(pvc[2]*par[2]);
				(estrutura+l)->id = (estrutura+i)->id;
				strcpy((estrutura+l)->nome , (estrutura+i)->nome);
				(estrutura+l)->coord[0] = raux[0];
				(estrutura+l)->coord[1] = raux[1];
				(estrutura+l)->coord[2] = raux[2];
				l++;
				fprintf(fpc,"\n%s\t%lf\t%lf\t%lf ",(estrutura+i)->nome, raux[0],raux[1],raux[2]);	
		}
	}

	

	//liberando o espaço alocado
	free(ptr); 
	free(estrutura);

	//fechando os arquivos
	fclose(fpe);
	fclose(fpc);

	return;	

}
