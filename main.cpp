#include <algorithm>
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>
#include <fstream>
#include <string> 
#include "mtwist.h"
#include <sstream>

// Macro para conversão de string para numero
#define SSTR( x ) static_cast< std::ostringstream & >( \
        ( std::ostringstream() << std::dec << x ) ).str()

typedef struct Disk{
	float x, y;
	unsigned int id;
} Disk;

typedef struct Celula{
	std::vector<Disk*> lista_discos;

	struct Celula* right, 
				 *left, 
				 *upper, 
				 *down, 
				 *up_right, 
				 *up_left, 
				 *down_right, 
				 *down_left;

	unsigned int linha, coluna;

}Celula;

#define PI (3.14159265358979323846)

const int 			SAMPLES 	= 10;
const int           Q = 200;
const int         	N = pow(256, 2);
const int         	N_sqrt = sqrt(N) + 0.5;
const float      	eta = 0.72;
const float      	sigma = sqrt(eta / (N * PI));
const float      	D = N_sqrt;
const float      	delxy = 1.0/ (2.0 * N_sqrt);
const float        two_delxy = 2.0 * delxy;


void novoL(mt_state* state, 
	Disk (&disk)[N_sqrt][N_sqrt], 
	Celula (&celula)[N_sqrt][N_sqrt],
	float D, 
	float sigma,
	int N_sqrt);

float event(float b_x, 
	float b_y, 
	float a_x, 
	float a_y, 
	int dirc, 
	float sigma);

/*void encontra_disco(Disk disco, Celula (&celula)[N_sqrt][N_sqrt], int* _i, int* _j, const int N_sqrt);*/

int main(){
	Celula 			celula[N_sqrt][N_sqrt];
	Disk 			disk[N_sqrt][N_sqrt];
	mt_state 		state[50];

	unsigned int id_count = 1;


	//gerando configuração inicial
	for (int i = 0; i < N_sqrt; i++){
		for (int j = 0; j < N_sqrt; j++){
			disk[i][j].x = delxy + i * two_delxy;
			disk[i][j].y = delxy + j * two_delxy;
			disk[i][j].id = id_count++;
		}
	}

	for (int i = 0; i < N_sqrt; ++i)
	{
		for (int j = 0; j < N_sqrt; ++j)
		{
			// Definição dos vizinhos

			// Se for a primeira  célula da grade
			if(i == 0 && j == 0){
				// Armazenando a linha e a coluna da celula
				celula[i][j].linha = i; celula[i][j].coluna = j;
				celula[i][j].left = &celula[i][N_sqrt-1];
				celula[i][j].right = &celula[i][j+1];
				celula[i][j].upper = &celula[i+1][j];
				celula[i][j].down = &celula[N_sqrt-1][j];
				celula[i][j].down_right = &celula[N_sqrt-1][j+1];
				celula[i][j].down_left = &celula[N_sqrt-1][N_sqrt-1];
				celula[i][j].up_right = &celula[i+1][j+1];
				celula[i][j].up_left = &celula[i+1][N_sqrt-1];
			}
			else if(i == 0 && j == N_sqrt-1){
				celula[i][j].linha = i; celula[i][j].coluna = j;
				celula[i][j].left = &celula[i][j-1];
				celula[i][j].right = &celula[0][0];
				celula[i][j].upper = &celula[i+1][j];
				celula[i][j].down = &celula[N_sqrt-1][j];
				celula[i][j].down_right = &celula[N_sqrt-1][0];
				celula[i][j].down_left = &celula[N_sqrt-1][j-1];
				celula[i][j].up_right = &celula[i+1][0];
				celula[i][j].up_left = &celula[i+1][j-1];
			}
			else if(i==N_sqrt-1 && j == 0){
				celula[i][j].linha = i; celula[i][j].coluna = j;
				celula[i][j].left = &celula[i][N_sqrt-1];
				celula[i][j].right = &celula[i][j+1];
				celula[i][j].upper = &celula[0][j];
				celula[i][j].down = &celula[i-1][j];
				celula[i][j].down_right = &celula[i-1][j+1];
				celula[i][j].down_left = &celula[i-1][N_sqrt-1];
				celula[i][j].up_right = &celula[0][j+1];
				celula[i][j].up_left = &celula[N_sqrt-1][N_sqrt-1];
			}
			else if(i == N_sqrt-1 && j == N_sqrt-1){
				celula[i][j].linha = i; celula[i][j].coluna = j;
				celula[i][j].left = &celula[i][j-1];
				celula[i][j].right = &celula[i][0];
				celula[i][j].upper = &celula[0][j];
				celula[i][j].down = &celula[i-1][j];
				celula[i][j].down_right = &celula[i-1][0];
				celula[i][j].down_left = &celula[i-1][j-1];
				celula[i][j].up_right = &celula[0][0];
				celula[i][j].up_left = &celula[0][j-1];
			}
			else if(i > 0 && i <= N_sqrt-2 && j == 0){
				celula[i][j].linha = i; celula[i][j].coluna = j;
				celula[i][j].left = &celula[i][N_sqrt-1];
				celula[i][j].right = &celula[i][j+1];
				celula[i][j].upper = &celula[i+1][j];
				celula[i][j].down = &celula[i-1][j];
				celula[i][j].down_right = &celula[i-1][j+1];
				celula[i][j].down_left = &celula[i-1][N_sqrt-1];
				celula[i][j].up_right = &celula[i+1][j+1];
				celula[i][j].up_left = &celula[i+1][N_sqrt-1];
			}
			else if(j > 0 && j < N_sqrt-1 && i == 0){
				celula[i][j].linha = i; celula[i][j].coluna = j;
				celula[i][j].left = &celula[i][j-1];
				celula[i][j].right = &celula[i][j+1];
				celula[i][j].upper = &celula[i+1][j];
				celula[i][j].down = &celula[N_sqrt-1][j];
				celula[i][j].down_right = &celula[N_sqrt-1][j+1];
				celula[i][j].down_left = &celula[N_sqrt-1][j-1]; 
				celula[i][j].up_right = &celula[i+1][j+1];
				celula[i][j].up_left = &celula[i+1][j-1];
			}
			else if(i > 0 && i <= N_sqrt-2 && j == N_sqrt-1){
				celula[i][j].linha = i; celula[i][j].coluna = j;
				celula[i][j].left = &celula[i][j-1];
				celula[i][j].right = &celula[i][0];
				celula[i][j].upper = &celula[i+1][j];
				celula[i][j].down = &celula[i-1][j];
				celula[i][j].down_right = &celula[i-1][0];
				celula[i][j].down_left = &celula[i-1][j-1];
				celula[i][j].up_right = &celula[i+1][0];
				celula[i][j].up_left = &celula[i+1][j-1];
			}
			else if(j > 0 && j <= N_sqrt-2 && i == N_sqrt-1){
				celula[i][j].linha = i; celula[i][j].coluna = j;
				celula[i][j].left = &celula[i][j-1];
				celula[i][j].right = &celula[i][j+1];
				celula[i][j].upper = &celula[0][j];
				celula[i][j].down = &celula[i-1][j];
				celula[i][j].down_right = &celula[i-1][j+1];
				celula[i][j].down_left = &celula[i-1][j-1];
				celula[i][j].up_right = &celula[0][j+1];
				celula[i][j].up_left = &celula[0][j-1];
			}
			else{
				celula[i][j].linha = i; celula[i][j].coluna = j;
				celula[i][j].left = &celula[i][j-1];
				celula[i][j].right = &celula[i][j+1];
				celula[i][j].upper = &celula[i+1][j];
				celula[i][j].down = &celula[i-1][j];
				celula[i][j].down_right = &celula[i-1][j+1];
				celula[i][j].down_left = &celula[i-1][j-1];
				celula[i][j].up_right = &celula[i+1][j+1];
				celula[i][j].up_left = &celula[i+1][j-1];
			}
		}
	}


	/****************************************************
	 * 			  ADCIONANDO DICOS AS CELULAS  			*
	 * Adciona a cada celula, o seu respectivo disco 	*
	 ****************************************************/
	for (int i = 0; i < N_sqrt; ++i){
		for (int j = 0; j < N_sqrt; ++j){
			celula[i][j].lista_discos.push_back(&disk[j][i]);
		}
	}
	for (int i = 0; i <= Q; i++){
        
            std::cout << "\t\t\tITERAÇÃO" << i << std::endl;
            std::ofstream myfile;
            std::string s = SSTR( i );
            std::string s2 = "resultados/" + s;
            myfile.open (s2.c_str());
            for (int i = 0; i < N_sqrt; ++i)
            {
                for (int j = 0; j < N_sqrt; ++j)
                {
                    myfile << disk[i][j].x << "," << disk[i][j].y << std::endl;
                }
            }
            myfile.close();
		novoL(state, disk, celula, D, sigma, N_sqrt);
  		
	}

}


void novoL(mt_state* state, Disk (&disk)[N_sqrt][N_sqrt], Celula (&celula)[N_sqrt][N_sqrt], float D, float  sigma, int N_sqrt)
{
    float event_min = 0;

    int dirc = (unsigned int) mts_lrand(state) % 2;

    float distance_to_go = D;
    std::cout << "\n" << D << "\n";
    Disk *next_a;

    int rand_i1, rand_i2;

    rand_i1 = (unsigned int)mts_lrand(state) % N_sqrt ;
    rand_i2 = (unsigned int)mts_lrand(state) % N_sqrt ;

    next_a = &disk[rand_i1][rand_i2];
    std::cout << std::endl;
    std::cout << "Informações da amostra: \ndirc="<<dirc<< std::endl;
    std::cout << "==========COMECOU===========" << std::endl;
    while(distance_to_go > 0.0)
    {
    	std::cout << "\n"<<std::endl;
        Disk* a = next_a;
        event_min = distance_to_go;

        unsigned int l = ceil(a->y / two_delxy) - 1; // linha
        unsigned int c = ceil(a->x / two_delxy) - 1; // coluna

        std::cout << "Informações do disco:\n Celula (" << l << ", "<<c<<")"<<"\nid="<<a->id<<std::endl;

    	Celula* aux;
	    for (int i = 0; i < celula[l][c].lista_discos.size(); i++){
	    	if(celula[l][c].lista_discos[i]->id != a->id){
	    		float event_b = event(celula[l][c].lista_discos[i]->x, 
	    			celula[l][c].lista_discos[i]->y, a->x, a->y, dirc, sigma);
	    		if (event_b < event_min){
	    			std::cout << "O disco esta na propria celula" << std::endl;
	    			event_min = event_b;
	    			next_a = celula[l][c].lista_discos[i];
	    		}
	    	}
	    }
	    // Celula da Direita
	    aux =  celula[l][c].right;
	    for (int i = 0; i < 3; ++i)
	    {
	    	for (int i = 0; i < aux->lista_discos.size(); i++){
				float event_b = event(aux->lista_discos[i]->x, 
					aux->lista_discos[i]->y, a->x, a->y, dirc, sigma);
				if (event_b < event_min){
					event_min = event_b;
					std::cout << "Entrou na celula up_right" << std::endl;
					next_a = aux->lista_discos[i];
				}
		    }
		    aux = aux->right;
	    }

	    // Celula upper-right
	    aux =  celula[l][c].up_right;
	    for (int i = 0; i < 3; ++i)
	    {
	    	for (int i = 0; i < aux->lista_discos.size(); i++){
				float event_b = event(aux->lista_discos[i]->x, 
					aux->lista_discos[i]->y, a->x, a->y, dirc, sigma);
				if (event_b < event_min){
					event_min = event_b;
					std::cout << "Entrou na celula up_right" << std::endl;
					next_a = aux->lista_discos[i];
				}
		    }
		    if(dirc == 0){
		    	aux = aux->upper;
			}
			else{
				aux = aux->right;
			}
	    }

	    aux =  celula[l][c].down_right;
	    for (int i = 0; i < 2; ++i)
	    {
	    	for (int i = 0; i < aux->lista_discos.size(); i++){
				float event_b = event(aux->lista_discos[i]->x, 
					aux->lista_discos[i]->y, a->x, a->y, dirc, sigma);
				if (event_b < event_min){
					event_min = event_b;
					std::cout << "Entrou na celula up_right" << std::endl;
					next_a = aux->lista_discos[i];
				}
		    }
		    aux = aux->right;
	    }

	    // Celula acima
	    aux =  celula[l][c].upper;
	    for (int i = 0; i < 2; ++i)
	    {
	    	for (int i = 0; i < aux->lista_discos.size(); i++){
				float event_b = event(aux->lista_discos[i]->x, 
					aux->lista_discos[i]->y, a->x, a->y, dirc, sigma);
				if (event_b < event_min){
					event_min = event_b;
					std::cout << "Entrou na celula up_right" << std::endl;
					next_a = aux->lista_discos[i];
				}
		    }
		    aux = aux->upper;
	    }

	    aux = celula[l][c].down;
	    for (int i = 0; i < 2; ++i)
	    {
	    	for (int i = 0; i < aux->lista_discos.size(); i++){
				float event_b = event(aux->lista_discos[i]->x, 
					aux->lista_discos[i]->y, a->x, a->y, dirc, sigma);
				if (event_b < event_min){
					event_min = event_b;
					std::cout << "Entrou na celula up_right" << std::endl;
					next_a = aux->lista_discos[i];
				}
		    }
		    aux = aux->down;
	    }


	    aux = celula[l][c].up_left;
	    for (int i = 0; i < 2; ++i)
	    {
	    	for (int i = 0; i < aux->lista_discos.size(); i++){
				float event_b = event(aux->lista_discos[i]->x, 
					aux->lista_discos[i]->y, a->x, a->y, dirc, sigma);
				if (event_b < event_min){
					event_min = event_b;
					std::cout << "Entrou na celula up_right" << std::endl;
					next_a = aux->lista_discos[i];
				}
		    }
		    aux = aux->upper;
	    }

	    aux = celula[l][c].left;
	    for (int i = 0; i < 2; ++i)
	    {
	    	for (int i = 0; i < aux->lista_discos.size(); i++){
				float event_b = event(aux->lista_discos[i]->x, 
					aux->lista_discos[i]->y, a->x, a->y, dirc, sigma);
				if (event_b < event_min){
					event_min = event_b;
					std::cout << "Entrou na celula up_right" << std::endl;
					next_a = aux->lista_discos[i];
				}
		    }
		    aux = aux->left;
	    }


	    aux = celula[l][c].down_left;
	    for (int i = 0; i < 2; ++i)
	    {
	    	for (int i = 0; i < aux->lista_discos.size(); i++){
				float event_b = event(aux->lista_discos[i]->x, 
					aux->lista_discos[i]->y, a->x, a->y, dirc, sigma);
				if (event_b < event_min){
					event_min = event_b;
					std::cout << "Entrou na celula up_right" << std::endl;
					next_a = aux->lista_discos[i];
				}
		    }
		    aux = aux->down;
	    }




        if (dirc == 1){
        	unsigned int l_anterior = l;
        	unsigned int c_anterior = c;
        	a->x = std::fmod((a->x + fabs(event_min-0.0001)), 1.0);
        	unsigned int c_atual = ceil(a->x / two_delxy) - 1; // coluna
        	unsigned int l_atual = ceil(a->y / two_delxy) - 1;
        	if(c_anterior != c_atual || l_atual != l_anterior){
        		Disk* d;
        		std::cout << "O Disco atualizou sua posição e mudou de célula:\nCelula Anterior: ("<<l<<", "<<c<<")"<<std::endl;
        		std::cout << "Celula Atual: (" <<l_atual << ", "<<c_atual<<")"<<std::endl; 
        		// Vou na celula anterior e removo o disco
        		for (int i = 0; i < celula[l_anterior][c_anterior].lista_discos.size(); ++i){
        			if(celula[l_anterior][c_anterior].lista_discos[i]->id == a->id){
						d = celula[l_anterior][c_anterior].lista_discos[i];
        				celula[l_anterior][c_anterior].lista_discos.erase(celula[l_anterior][c_anterior].lista_discos.begin() + i);
        				break;
        			}
        		}
        		std::cout << "Verificando se removeu da celula anterior, e atualizou:" << std::endl;
        		std::cout << "	 ==========================================	 " << std::endl;
        		for(int i = 0; i < celula[l_anterior][c_anterior].lista_discos.size(); ++i){
        			std::cout << "ANTERIOR:" <<  celula[l_anterior][c_anterior].lista_discos[i]->id <<  std::endl;
        		}

        		//Vou na célula atual e armazeno o disco:
        		celula[l_atual][c_atual].lista_discos.push_back(d);
        		for(int i = 0; i < celula[l_atual][c_atual].lista_discos.size(); ++i){
        			std::cout << "NOVA:" <<  celula[l_atual][c_atual].lista_discos[i]->id <<  std::endl;
        		}
        		std::cout << "	 ==========================================	 " << std::endl;
        		
        	}
        }
        else{
        	unsigned int l_anterior = l;
        	unsigned int c_anterior = c;
        	a->y = std::fmod((a->y + fabs(event_min-0.0001)), 1.0);
        	unsigned int c_atual = ceil(a->x / two_delxy) - 1; // coluna
        	unsigned int l_atual = ceil(a->y / two_delxy) - 1;
        	if(l_anterior != l_atual || l_atual != l_anterior){
        		Disk* d;
        		std::cout << "O Disco atualizou sua posição e mudou de célula:\nCelula Anterior: ("<<l<<", "<<c<<")"<<std::endl;
        		std::cout << "Celula Atual: (" <<l_atual << ", "<<c_atual<<")"<<std::endl; 
        		// Vou na celula anterior e removo o disco
        		for (int i = 0; i < celula[l_anterior][c_anterior].lista_discos.size(); ++i){
        			if(celula[l_anterior][c_anterior].lista_discos[i]->id == a->id){
						d = celula[l_anterior][c_anterior].lista_discos[i];
        				celula[l_anterior][c_anterior].lista_discos.erase(celula[l_anterior][c_anterior].lista_discos.begin() + i);
        				break;
        			}
        		}
        		std::cout << "Verificando se removeu da celula anterior, e atualizou:" << std::endl;
        		std::cout << "	 ==========================================	 " << std::endl;
        		for(int i = 0; i < celula[l_anterior][c_anterior].lista_discos.size(); ++i){
        			std::cout << "ANTERIOR:" <<  celula[l_anterior][c_anterior].lista_discos[i]->id <<  std::endl;
        			/*if(celula[l_anterior][coluna].lista_discos[i]->id == a->id){
						std::cout << "TA NA ANTERIOR!!!" << std::endl;
         			}*/
        		}

        		//Vou na célula atual e armazeno o disco:
        		celula[l_atual][c_atual].lista_discos.push_back(d);
        		for(int i = 0; i < celula[l_atual][c_atual].lista_discos.size(); ++i){
        			std::cout << "NOVA:" <<  celula[l_atual][c_atual].lista_discos[i]->id <<  std::endl;
        			/*if(celula[l_atual][coluna].lista_discos[i]->id == a->id){
						std::cout << "TA NA NOVA!!!" << std::endl;
         			}*/
        		}
        		std::cout << "	 ==========================================	 " << std::endl;
        	}
        }
		std::cout << "distance to go ANTERIOR = " << distance_to_go << std::endl << "event_min = "<< event_min << std::endl;
        distance_to_go -= event_min;
        std::cout << "dtg = " << distance_to_go << std::endl;
    }

}


float event(float b_x, float b_y, float a_x, float a_y, int dirc, float sigma)
{
    float d_perp = 0.0;
    float d_para = 0.0;

    if (dirc == 1)
        d_perp = std::fmod(fabs(b_y - a_y), 1.0);
    else
        d_perp = std::fmod(fabs(b_x - a_x), 1.0);

    d_perp = fmin(d_perp, 1.0 - d_perp);
    if (d_perp > 2.0 * sigma)
        return (float)INFINITY;
    else{
        d_para = sqrt(fabs(4.0 * sigma*sigma - d_perp*d_perp));

        if (dirc == 1){
            return std::fmod((b_x - a_x - d_para + 1.0), 1.0);
        }
        else{
            return std::fmod((b_y - a_y - d_para + 1.0), 1.0);
        }
    }
}
