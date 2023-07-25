#include <stdio.h>
#include <stdlib.h>
#include <SDL.h> // Para a biblioteca SDL
#include <conio.h>
#include <math.h>

#define pi 3.14159265358979323

typedef struct
{
    double x; 		// define a coordenada x de um ponto no plano cartesiano
    double y;		// define a coordenada y de um ponto no plano cartesiano
} ponto_xy;


typedef struct
{
    ponto_xy A; 	// ponto A de uma reta no plano cartesiano
    ponto_xy B;		// ponto B de uma reta no plano cartesiano
} reta_xy;

typedef struct
{
    ponto_xy C;				// centro da Terra
    double raio; 			//raio da Terra
    double camada_h[10001];	// define a altura de cada camada atmosferica
    double n[10001]; 		// define indice de refracao de cada camada atmosferica
    int atm_num; 			//numero de camadas amosfericas
} Terra;

int teste_2(reta_xy* raio, Terra* T, ponto_xy P, double ang, int luz_cor); //funcao que calcula os raios refratados na atmosfera
double angulo_incidencia(ponto_xy v, ponto_xy P, ponto_xy C); //funcao que calcula o angulo de incidencia em relacao a normal
ponto_xy rotacao_vetorial(ponto_xy v, double ang); // funcao que rotaciona o vetor diretor
char posicao_vetor_raio(ponto_xy v, ponto_xy P, ponto_xy C); // funcao que calcula a posicao relativa do raio incidente em relacao a normal
void viewport_LightRay_render(reta_xy* L_ray, double xpp, double ypp, double zpp, int nray, int luz_cor[3]); // funcao que desenha os raios nas coordenadas da tela
void viewport_Earth_render(Terra T, double xpp, double ypp, double zpp); // funcao que desenha a atmosfera nas coordenadas da tela
void drawLine(SDL_Renderer* renderer, int xi, int yi, int xf, int yf); 
void drawCircle(SDL_Renderer* renderer, int centerX, int centerY, int radius);
void drawCircle2(SDL_Renderer* renderer, int centerX, int centerY, int Raio);


/* Variaveis globais para dados auxiliares. */
ponto_xy Paux, Caux; //Paux salva o ponto em que o raio atinge a primeira camada nas coordenadas da tela. Caux salva o centro da terra nas coordenadas da tela
double Raux, C00; //Raux eh o raio da Terra na escala das coordenadas da tela
int camaux = 0; // determina um ponto dentro do vetor de raios
SDL_Renderer* renderer = NULL;


int main(int argc, char* argv[])
{
    Terra T; // Variavel com os parametros da terra
    int i, nray1, nray2, point; // variaveis de controle
    double altM; //altura maxima da atmosfera
    double esp; 		// espessura de cada camada
    reta_xy L_ray1[7000], L_ray2[7000]; 	// declaração do conjunto de direções que representara o raio de luz na atmosfera
    ponto_xy Po1, Po2, Pcam;			// ponto de teste
    double aux1=0, aux2=0; // variavel de controle para atribuicao dos indices de refracao
    double z;

	// parametros para cor
	int luz_cor1[3], luz_cor2[3];
	luz_cor1[0] = 255;
	luz_cor1[1] = 255;
	luz_cor1[2] = 50;
	luz_cor2[0] = 255;
	luz_cor2[1] = 0;
	luz_cor2[2] = 255;
    
	/********************************************** atribuicao de parametros iniciais reais ********************************************/

    // Unidade de comprimento: Km
    T.C.x = 0;			// coordenada x do centro
    T.C.y = 0;			// coordenada y do centro
	T.raio = 6000;		// raio
    altM = 10;			// altitude máxima da atmosfera considerada
    T.atm_num = 5000;		// numero de camadas
    esp = altM / T.atm_num;	//calculo da espessura de cada camada
    Po1.x = 0;			// atribuindo ponto de Teste x
    Po1.y = T.raio + 0.01;				// atribuindo ponto de Teste y
    Po2.x = 0;
    Po2.y = T.raio + 0.02;


	// Imprime o raio da Terra
    printf("\n\nRaio: %lf \n\n", T.raio);

	// Taxa de variacao do indice de refracao em funcao da altura (unidade m^(-1))
	double dn_dh = 30 * pow(10, -8);

    for (i = 0; i < T.atm_num; i++)
    {
        T.camada_h[i] = T.raio + (i + 1) * esp;			//calculo da altura de cada camada.
		
		// Calculo do indice de refracao de cada camada
		aux1 = 1 + 1000 * (T.atm_num * esp) * dn_dh; 
		aux1 = aux1 - esp * 1000 * dn_dh * i;
		//aux1 = aux1 - esp * 1000 * dn_dh;
		if (aux1 < 1) aux1 = 1;
		T.n[i] = aux1;

        printf("\n  h(%d) = %lf    n(%d) = %lf", i, T.camada_h[i], i, T.n[i]);	//linha de teste - testa se as alturas estão sendo atribuidas corretamente
    }
    T.n[0] = 1; // Atribui o menor indice de refracao para a primeira camada
    /*************************************************************************************************************************************/
  
	// Inicializa o SDL
	SDL_Init(SDL_INIT_VIDEO);

	// Cria uma janela
	SDL_Window* window = SDL_CreateWindow("Desenhar Linha com SDL", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 2500, 1400, SDL_WINDOW_SHOWN);

	// Cria um renderer para a janela
	renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

	// Calcula a trajetória de cada raio na atmosferaa e retorna o numero de raios calculados.
    nray1 = teste_2(L_ray1, &T, Po1, 3.95 * pi / 4 /*4*pi/4*/, 14);
    //nray2 = teste_2(L_ray2, &T, Po2, 3.95001 * pi / 4 /*4*pi/4*/, 14);

	// Imprime na tela alguns parametros recferentes ao raio
	for (i = 0; i <= nray1; i++)
	{
		printf("\n  P%d = (%lf; %lf)    d = %lf", i, L_ray1[i].B.x, L_ray1[i].B.y, sqrt(L_ray1[i].B.x * L_ray1[i].B.x + L_ray1[i].B.y * L_ray1[i].B.y));
	}

	// ajusa as cores e atualiza o renderer
	SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
	SDL_RenderClear(renderer);
	SDL_RenderPresent(renderer);

	Paux.x = 0; // recebe zero como controle de variavel (variavel camaux de isualizacao do raio)
	point = camaux; // determina o ponto em que sera visualizado o raio
	Pcam = Paux;
	//point = 0;
	Pcam.x = -582.8527950206;
	Pcam.y = 5971.6251238118;
	Pcam.x = -200;
	Pcam.y = 5971.6251238118;
	z = 500; // determina a distancia de visualizacao
	//Pcam = L_ray1[point].B; // determinacao do ponto de visualizacao da camera dentro do sistema de coordenadas
	viewport_Earth_render(T, Pcam.x, Pcam.y, z);
	viewport_LightRay_render(L_ray1, Pcam.x, Pcam.y, z, nray1, luz_cor1);

	//viewport_LightRay_render(L_ray2, Pcam.x, Pcam.y, z, nray2, luz_cor2);
	SDL_RenderPresent(renderer);
	SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
	SDL_RenderClear(renderer);
	
	SDL_Event event1;
	while (SDL_WaitEvent(&event1)) {
		if (event1.type == SDL_KEYDOWN) {
			break;
		}
	}
	/*
	for (i = 1; i <= 100; i++)
	{
		if(i <= 5) z = 0.5*i;
		if (5 < i && i <= 10) z = i;
		if (10 < i && i <= 20) z = 2 * i;
		if (20 < i && i <= 30) z = 4 * i;
		if (30 < i && i <= 40) z = 6 * i;
		if (40 < i && i <= 50) z = 8 * i;
		if (50 < i && i <= 60) z = 10 * i;
		if (60 < i && i <= 80) z = 12 * i;
		if (i > 80) z = 14 * i;

		viewport_Earth_render(T, Pcam.x, Pcam.y, z);
		viewport_LightRay_render(L_ray1, Pcam.x, Pcam.y, z, nray1, 14);
		viewport_LightRay_render(L_ray2, Pcam.x, Pcam.y, z, nray2, 12);B

		// Atualiza o Render

		SDL_RenderPresent(renderer);
		SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
		SDL_RenderClear(renderer);
	}*/


	printf("\n\n\n Raio               = %lf \n Camada 0 transform = %lf \n camada 0 calculada = %lf", Raux, C00, sqrt((Paux.x - Caux.x) * (Paux.x - Caux.x) + (Paux.y - Caux.y) * (Paux.y - Caux.y))); /***************************************************************************************************************************************/
	printf("\n\n\n %d", camaux);
	//getch();

	    // Aguarda um evento
    SDL_Event event;
    while (SDL_WaitEvent(&event)) {
        if (event.type == SDL_KEYDOWN) {
            break;
        }
    }
    //getch();

    // Libera recursos
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();

    return 0;
}

int teste_2(reta_xy* L_ray, Terra* T, ponto_xy P, double ang)
{
	int i;
	int camada, aux;		// variaveis de controle para armazenar a camada em que o raio se encontra.
	char pos;				// variavel de controle. Verifica aa posicao do vetor diretor em relacao a a reta normal.
	double a, b, c;			// parametros da equacao gerada para t no calculo das interseccoes.
	double ang_i, ang_r;	// angulos de incidencia e refracao
	double H, h;			// altura da camada suerior e altura da camada inferior ao ponto inicial
	double t;				// parametro para a equação vetorial da reta
	double d;				// distancia do ponto inicial ao centro da Terra (tambem usada como variavel auxiliar nos calculos)
	ponto_xy v;				// vetor diretor do raio de luz




	/************** TRECHO DE CALCULOS PARA O PRIMEIRO RAIO DE LUZ *********************/

	// inicia o vetor diretor a partir do angulo de direcao inicial do raio de luz
	v.x = cos(ang);
	v.y = sin(ang);


	// Atribui o ponto inicial ao primeiro ponto do conjunto de retas do raio
	L_ray[0].A.x = P.x;
	L_ray[0].A.y = P.y;

	// calculo da distancia do ponto inicial ao centro da Terra
	d = sqrt((L_ray[0].A.x - (*T).C.x) * (L_ray[0].A.x - (*T).C.x) + (L_ray[0].A.y - (*T).C.y) * (L_ray[0].A.y - (*T).C.y));


	// calculo do raio (altura) da camada em que o ponto inicial está localizado
	i = 0;
	while (i < (*T).atm_num)		// laco para analisar camada por camada
	{
		if ((*T).camada_h[i] > d)	// se a altura da camada for maior que a distancia do ponto ao centro, o ponto estara nessa camada
		{
			H = (*T).camada_h[i];					// atribuicao da altura da camada superior ao ponto inicial
			if (i == 0) h = (*T).raio;				// se a camada eh a primeira, a camada inferior deve ser a superficie
			else		h = (*T).camada_h[i - 1];	// atribuicao da altura da camada superior ao ponto inicial

			// armazena a camada em que se encontra o raio de luz que está sendo calculado.
			camada = i;

			i = (*T).atm_num + 1;
		}

		else
		{
			i++;
		}
	}


	// calculo do parametro do ponto de interseccao do raio com a camada superior
	a = 1;
	b = 2 * (v.x * (L_ray[0].A.x - (*T).C.x) + v.y * (L_ray[0].A.y - (*T).C.y));
	c = pow((L_ray[0].A.x - (*T).C.x), 2) + pow((L_ray[0].A.y - (*T).C.y), 2) - H * H;

	t = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);		//calculo efetivo do parametro t
	aux = 1;											// variavel controle indicando que a proxima camada esta mais acima

	//calculo do parametro ponto de interseccao do raio com a camada inferior (caso exista)
	c = pow((L_ray[0].A.x - (*T).C.x), 2) + pow((L_ray[0].A.y - (*T).C.y), 2) - h * h;
	d = b * b - 4 * a * c;
	if (d >= 0)	//condicao de exeistencia
	{
		if ((-b - sqrt(b * b - 4 * a * c)) / (2 * a) >= 0) //verifica se a intersecao correu antes do ponto inicial.
		{
			t = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
			aux = -1;									// variavel controle indicando que a proxima camada esta mais abaixo.
		}

	}

	// atribuicao do ponto de interseccao com a camada calculada a partir do paramentro e do vetor diretor
	L_ray[0].B.x = L_ray[0].A.x + v.x * t;
	L_ray[0].B.y = L_ray[0].A.y + v.y * t;

	/*************************** FIM DO TRECHO ***********************************/


	i = 1;

	while (i < 6999)
	{
		/*
		Esse laco calcula as direcoes dos raios de luz atraves das camadas.
		Segue a sequencia basica de passos:
		1) Atribui ao pononto inicial do raio atual o ponto final do raio anterior. Esse ponto ocorre exatamente na transicao entre duas camadas.
		2) calcula o anguo de incidencia em relacao a a reta normal no ponto de incidencia. O ponto de incidencia eh o ponto A (inicial do raio atual).
		3) Verifica se havera refracao ou reflexao total.
			3.1) Caso haja refracao, o angulo de refracao eh calculado e a direcao do raio eh atualizada. Nesse caso, a camada tambem sera atualizada, pois indica que o raio atravessou.
			3.2) Caso haja reflexao total, o angulo de reflexao eh calculado e a direcao do raio de luz eh atualizada. Nesse caso a camada nao sera atualizada, pois o raio refletiu de volta para a mesma camada.
		4) Calcula o ponto de interseccao com a proxima camada.
			4.1) Caso a interseccao ocorra na camada acima, aux recebera 1, atualizando (posteriormente) a camada para a superior.
			4.2) Caso a interseccao ocorra na camada abaixo, aux recebera -1, atualizando (posteriormente) a camada para a inferior.
		5) Atualiza o valor do ponto B (ponto final do raio) conforme o calculo realizado anteriormente.
		6) Repete o processo.
		Obs. A variavel aux indica apenas em qual camada ocorrera a interseccao com o raio de luz. aux = 1 caso ocorra na camada de cima e aux = -1 caso ocorra na camada de baixo.
		Esse controle eh necessario para que as camadas sejam atualizadas corretamente.

		*/

		// o ponto inicial do raio i deve coincidir com o ponto final do raio i-1. pela continuidade dos raios de luz.
		L_ray[i].A = L_ray[i - 1].B;

		// Calculo do ângulo de incidência no ponto pertencente a a superficie dioptrica
		ang_i = angulo_incidencia(v, L_ray[i].A, (*T).C);


		/******** calculo do angulo de refracao pela Lei de Snell-Descartes. *********/
		// indica que haverá refração
		if ((*T).n[camada] * sin(ang_i) / (*T).n[camada + aux] >= 0 && (*T).n[camada] * sin(ang_i) / (*T).n[camada + aux] <= 1)
		{
			//Lei de Snell-Descartes
			ang_r = asin((*T).n[camada] * sin(ang_i) / (*T).n[camada + aux]);

			// analisa em qual posicao o vetor diretor esta em relacao ao vetor normal para determinar  o sentido de rotacao 			
			pos = posicao_vetor_raio(v, L_ray[i].A, (*T).C);

			// Atualizacao do vetor diretor apos a refracao. O calculo eh gerado a partir da rotacao do vetor diretor em termos do angulo de refracao 
			if (pos == 'E')
			{
				ang = ang_r - ang_i;
			}

			else
			{
				ang = ang_i - ang_r;
			}

			v = rotacao_vetorial(v, ang);


			//atualizacao da camada. Armazena a nova camada em que o raio se encontra.
			camada = camada + aux;
		}

		else //indica que haverá reflexão total
		{
			if (camada > (*T).atm_num - 10) printf("\n\n Reflexão total superior na camada : %d \n iteracao: %d \n Coordenadas do raio inicial (%.10lf, %.10lf) \n", camada, i, L_ray[i].A.x, L_ray[i].A.y);
			if (camada < 55) printf("\n\n Reflexão total inferior na camada : %d \n iteracao: %d \n Coordenadas do raio inicial (%.10lf, %.10lf) \n", camada, i, L_ray[i].A.x, L_ray[i].A.y);

			pos = posicao_vetor_raio(v, L_ray[i].A, (*T).C);

			if (pos == 'E')
			{
				ang = pi - 2 * ang_i;
			}
			else
			{
				ang = pi + 2 * ang_i;
			}
			v = rotacao_vetorial(v, ang);
		}



		/****************************  Calculo do ponto de interseccao do raio com as camadas adjacentes ******************************************/
		// Calculo das novas alturas limites na camada em que o raio se encontra
		H = (*T).camada_h[camada];						// atribuicao da altura da camada superior ao ponto inicial
		if (camada == 0) h = (*T).raio;				// se a camada eh a primeira, a camada inferior deve ser a superficie
		else		h = (*T).camada_h[camada - 1];		// atribuicao da altura da camada superior ao ponto inicial




		// calculo do parametro do ponto de interseccao do raio com a camada superior
		a = 1;
		b = 2 * (v.x * (L_ray[i].A.x - (*T).C.x) + v.y * (L_ray[i].A.y - (*T).C.y));
		c = pow((L_ray[i].A.x - (*T).C.x), 2) + pow((L_ray[i].A.y - (*T).C.y), 2) - H * H;

		t = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);		//calculo efetivo do parametro t
		aux = 1;											// variavel controle indicando que a proxima camada esta mais acima

		//calculo do parametro ponto de interseccao do raio com a camada inferior (caso exista)
		c = pow((L_ray[i].A.x - (*T).C.x), 2) + pow((L_ray[i].A.y - (*T).C.y), 2) - h * h;
		d = b * b - 4 * a * c;
		if (d >= 0)
		{
			if ((-b - sqrt(b * b - 4 * a * c)) / (2 * a) > 0)
			{
				t = (-b - sqrt(b * b - 4 * a * c)) / (2 * a);
				aux = -1;									// variavel controle indicando que a proxima camada esta mais abaixo.
			}

		}

		// atribuicao do ponto de interseccao com a camada calculada a partir do paramentro e do vetor diretor
		L_ray[i].B.x = L_ray[i].A.x + v.x * t;
		L_ray[i].B.y = L_ray[i].A.y + v.y * t;

		/*********************************** Fim do calculo ***************************************************************************************/

		i++;

		if (camada + aux < 0 || camada >(*T).atm_num)
		{
			return i;
		}


	}

	return i - 1;
}


double angulo_incidencia(ponto_xy v, ponto_xy P, ponto_xy C)
{
	ponto_xy n;		// declaracao do vetor normal a a superficie no ponto P.
	double ang_i;

	n.x = (P.x - C.x) / sqrt((P.x - C.x) * (P.x - C.x) + (P.y - C.y) * (P.y - C.y));	// calculo da coordenada x do vetor normal.
	n.y = (P.y - C.y) / sqrt((P.x - C.x) * (P.x - C.x) + (P.y - C.y) * (P.y - C.y));	// calculo da cordenada y do vero normtal.

	// calculo do angulo de incidencia pelo produto inerno entre o vetor diretor e o vetor normal a a camada.
	ang_i = acos(v.x * n.x + v.y * n.y);

	// caso o angulo seja maior que 90 graus, deve-se considerar o vetor oposto ao veor normal aneriormente considerado.
	if (ang_i > pi / 2)
	{
		n.x = (-1) * n.x;
		n.y = (-1) * n.y;
		ang_i = acos(v.x * n.x + v.y * n.y);
	}

	// retorna o angulo calculado.p
	return ang_i;
}


ponto_xy rotacao_vetorial(ponto_xy v, double ang)
{
	ponto_xy u;

	u.x = v.x * cos(ang) - v.y * sin(ang);
	u.y = v.x * sin(ang) + v.y * cos(ang);

	return u;
}


char posicao_vetor_raio(ponto_xy v, ponto_xy P, ponto_xy C)
{
	ponto_xy n;		// declaracao do vetor normal a a superficie no ponto P.
	double ang;
	char c;

	n.x = (P.x - C.x) / sqrt((P.x - C.x) * (P.x - C.x) + (P.y - C.y) * (P.y - C.y));	// calculo da coordenada x do vetor normal.
	n.y = (P.y - C.y) / sqrt((P.x - C.x) * (P.x - C.x) + (P.y - C.y) * (P.y - C.y));	// calculo da cordenada y do vero normtal.

	// calculo do angulo de incidencia pelo produto inerno entre o vetor diretor e o vetor normal a a camada.
	ang = acos(v.x * n.x + v.y * n.y);

	// caso o angulo seja maior que 90 graus, deve-se considerar o vetor oposto ao veor normal anteriormente considerado.
	if (ang > pi / 2)
	{
		n.x = (-1) * n.x;
		n.y = (-1) * n.y;
	}

	if (n.y >= 0)
	{
		ang = acos(n.x);
	}
	else
	{
		ang = 2 * pi - acos(n.x);
	}

	v = rotacao_vetorial(v, -ang);

	if (v.y >= 0)
	{
		c = 'E';
	}
	else
	{
		c = 'D';
	}

	return c;
}


void viewport_LightRay_render(reta_xy* L_ray, double xpp, double ypp, double zpp, int nray, int luz_cor[3])
{
	int i;
	double xi, yi, zi, xf, yf, zf;
	double d = 0.05;


	for (i = 0; i <= nray; i++)
	{
		xi = L_ray[i].A.x;
		yi = L_ray[i].A.y;
		xf = L_ray[i].B.x;
		yf = L_ray[i].B.y;

		xi = xi - xpp;
		yi = yi - ypp;
		xf = xf - xpp;
		yf = yf - ypp;

		xi = (d / (d + zpp)) * xi;
		yi = (d / (d + zpp)) * yi;
		xf = (d / (d + zpp)) * xf;
		yf = (d / (d + zpp)) * yf;



		xi = (2000 / d) * xi + 1250;
		yi = -(2000 / d) * yi + 700;
		xf = (2000 / d) * xf + 1250;
		yf = -(2000 / d) * yf + 700;

		//printf("\n  P%d = (%lf; %lf)", i, xf, yf);
		if (i == camaux && Paux.x == 0)
		{
			Paux.x = xf;
			Paux.y = yf; /*****************************************************************************************************************************************************************************************************************************/
		}

		

		// Define a cor da linha do raio
		SDL_SetRenderDrawColor(renderer, luz_cor[0], luz_cor[1], luz_cor[2], 255);

		
		// Desenha uma linha de raio no renderer
		drawLine(renderer, xi, yi, xf, yf);

	}

}


void viewport_Earth_render(Terra T, double xpp, double ypp, double zpp)
{
	int i = 0, int_delta = 0, camada_n = 0, var = 1;
	double d = 0.05;
	double int_set;

	T.C.x = T.C.x - xpp;
	T.C.y = T.C.y - ypp;

	T.C.x = (d / (d + zpp)) * T.C.x;
	T.C.y = (d / (d + zpp)) * T.C.y;


	printf("\n  C_i = (%lf; %lf)", T.C.x, T.C.y);
	T.C.x = (2000 / d) * T.C.x + 1250;
	T.C.y = -(2000 / d) * T.C.y + 700;

	Caux = T.C; /****************************************************************************************************************************************************************************************************************************************/
	printf("\n  C_f = (%lf; %lf)", T.C.x, T.C.y);


	printf("\n  Raio_i = %lf", T.raio);
	T.raio = (d / (d + zpp)) * T.raio;
	T.raio = (2000 / d) * T.raio;
	Raux = T.raio;  /********************************************************************************************************************************************************************************************************************************/
	printf("\n  Raio_f = %lf", T.raio);


	// Define a cor da linha (vermelho)
	SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255);
	drawCircle(renderer, T.C.x, T.C.y, T.raio);
	// Atualiza o renderer
	//SDL_RenderPresent(renderer);

	if (zpp == 0) camada_n = 50;
	if (0 < zpp && zpp <= 0.5) camada_n = 100;
	if (0.5 < zpp && zpp <= 1) camada_n = 1000;
	if (1 < zpp && zpp <= 2) camada_n = 1800;
	if (2 < zpp) camada_n = 2000;
	

	if (7 < zpp && zpp <= 10)
	{
		var = 4;
	}

	if (10 < zpp && zpp <= 13)
	{
		var = 6;
	}

	if (13 < zpp && zpp <= 16)
	{
		var = 10;
	}

	if (16 < zpp)
	{
		var = 15;
	}

	int_set = 0.05;

	for (i = 0; i <= camada_n; i += var  )
	{
		T.camada_h[i] = (d / (d + zpp)) * T.camada_h[i];
		T.camada_h[i] = (2000 / d) * T.camada_h[i];
		
		int_delta = i * int_set;
		
		// Define a cor da linha (vermelho)
		SDL_SetRenderDrawColor(renderer, 100 - int_delta, 100 - int_delta, 100 - int_delta, 0);
		drawCircle2(renderer, T.C.x, T.C.y, T.camada_h[i]);
					
	}

	SDL_SetRenderDrawColor(renderer, 0, 0, 255, 255);
	drawCircle2(renderer, T.C.x, T.C.y, T.raio);

	C00 = T.camada_h[0]; /***********************************************************************************************************************************************************************************************************************/
	//printf("\n  Camada_0 = %lf", T.camada_h[0]);
}

void drawLine(SDL_Renderer* renderer, int xi, int yi, int xf, int yf)
{
	int x1 = xi;
	int x2 = xf;
	int y1 = yi;
	int y2 = yf;
	int dx = abs(x2 - x1);
	int dy = abs(y2 - y1);
	int sx = (x1 < x2) ? 1 : -1;
	int sy = (y1 < y2) ? 1 : -1;
	int err = dx - dy;

	while (1)
	{
		
		if (x1 <= 2500 && x1 >= 0 && y1 <= 1500 && y1 >= 0)
		{
			SDL_RenderDrawPoint(renderer, x1, y1);
		}
	
		if (x1 == x2 && y1 == y2)
			break;

		int e2 = 2 * err;
		if (e2 > -dy)
		{
			err -= dy;
			x1 += sx;
		}
		if (e2 < dx)
		{
			err += dx;
			y1 += sy;
		}
	}
}

void drawCircle(SDL_Renderer* renderer, int centerX, int centerY, int radius) {
	int x = radius - 1;
	int y = 0;
	int dx = 1;
	int dy = 1;
	int error = dx - (radius << 1);


	while (x >= y) {
		// Desenhe os pixels simétricos de oito vias em torno do centro
		if (centerX + x <= 2500 && centerX + x >= 0 && centerY + y <= 1500 && centerY + y >= 0)
		{
			SDL_RenderDrawPoint(renderer, centerX + x, centerY + y);
		}

		if (centerX + x <= 2500 && centerX + x >= 0 && centerY - y <= 1500 && centerY - y >= 0)
		{
			SDL_RenderDrawPoint(renderer, centerX + x, centerY - y);
		}

		if (centerX - x <= 2500 && centerX - x >= 0 && centerY + y <= 1500 && centerY + y >= 0)
		{
			SDL_RenderDrawPoint(renderer, centerX - x, centerY + y);
		}

		if (centerX - x <= 2500 && centerX - x >= 0 && centerY - y <= 1500 && centerY - y >= 0)
		{
			SDL_RenderDrawPoint(renderer, centerX - x, centerY - y);
		}


		if (centerX + y <= 2500 && centerX + y >= 0 && centerY + x <= 1500 && centerY + x >= 0)
		{
			SDL_RenderDrawPoint(renderer, centerX + y, centerY + x);
		}

		if (centerX + y <= 2500 && centerX + y >= 0 && centerY - x <= 1500 && centerY - x >= 0)
		{
			SDL_RenderDrawPoint(renderer, centerX + y, centerY - x);
		}

		if (centerX - y <= 2500 && centerX - y >= 0 && centerY + x <= 1500 && centerY + x >= 0)
		{
			SDL_RenderDrawPoint(renderer, centerX - y, centerY + x);
		}

		if (centerX - y <= 2500 && centerX - y >= 0 && centerY - x <= 1500 && centerY - x >= 0)
		{
			SDL_RenderDrawPoint(renderer, centerX - y, centerY - x);
		}

		// Atualiza a posição e o erro
		if (error <= 0) {
			y++;
			error += dy;
			dy += 2;
		}
		if (error > 0) {
			x--;
			dx += 2;
			error += dx - (radius << 1);
		}

	}
}

void drawCircle2(SDL_Renderer* renderer, int centerX, int centerY, int Raio)
{
	double x, y, xl;
	double cos_ang, sin_ang, ang, R_inv;

	x = -centerX + 2500;
	ang = acos(x / Raio);

	y = Raio * sin(ang);

	R_inv = Raio;
	R_inv = 1.0 / Raio;

	cos_ang = cos(R_inv);
	sin_ang = sin(R_inv);

	while (x >= (-centerX))
	{
		xl = x;
		x = x * cos_ang - y * sin_ang;
		y = xl * sin_ang + y * cos_ang;
		SDL_RenderDrawPoint(renderer, centerX + x, centerY - y);
	}

}
