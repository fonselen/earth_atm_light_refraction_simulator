#include <stdio.h>
#include <stdlib.h>
#include <SDL.h> // Para a biblioteca SDL
#include <conio.h>
#include <math.h>

#define pi 3.14159265358979323
#define EXP 2.71828182845904523536

typedef struct
{
	long double x; 		// define a coordenada x de um ponto no plano cartesiano
	long double y;		// define a coordenada y de um ponto no plano cartesiano
} ponto_xy;


typedef struct
{
	ponto_xy A; 	// ponto A de uma reta no plano cartesiano
	ponto_xy B;		// ponto B de uma reta no plano cartesiano
} reta_xy;

typedef struct
{
	ponto_xy* Data;
	int color[3];
	int TotalReflectionPosition;

} LightRay;

typedef struct
{
	ponto_xy C;				// centro da Terra
	long double raio; 			//raio da Terra
	long double* camada_h;	// define a altura de cada camada atmosferica
	long double* n; 		// define indice de refracao de cada camada atmosferica
	int atm_num; 			//numero de camadas amosfericas
} Terra;

int teste_2(LightRay* L_ray, Terra* T, ponto_xy P, long double ang, long double xM, int* pos_ReflexaoTotal); //funcao que calcula os raios refratados na atmosfera
void LightRay_ResolutionConverter(ponto_xy* L_ray_i, ponto_xy* L_ray_f, int pos_ReflexaoTotal, int nray, int LightRay_Resolution);
long double angulo_incidencia(ponto_xy v, ponto_xy P, ponto_xy C); //funcao que calcula o angulo de incidencia em relacao a normal
ponto_xy rotacao_vetorial(ponto_xy v, long double ang); // funcao que rotaciona o vetor diretor
char posicao_vetor_raio(ponto_xy v, ponto_xy P, ponto_xy C); // funcao que calcula a posicao relativa do raio incidente em relacao a normal
void viewport_LightRay_render(ponto_xy* L_ray, long double xpp, long double ypp, long double zpp, int nray, int luz_cor[3]); // funcao que desenha os raios nas coordenadas da tela
void viewport_Earth_render(Terra* T, long double xpp, long double ypp, long double zpp); // funcao que desenha a atmosfera nas coordenadas da tela
void drawLine(SDL_Renderer* renderer, int xi, int yi, int xf, int yf);
void drawCircle(SDL_Renderer* renderer, int centerX, int centerY, int radius);
void drawCircle2(SDL_Renderer* renderer, int centerX, int centerY, int Raio);
void drawCircle3(SDL_Renderer* renderer, int centerX, int centerY, int Raio);
ponto_xy Earth_tangent_point(ponto_xy PO, Terra* T);
reta_xy Earth_tangent_line(ponto_xy PO, Terra* T);
double FindTargetPointAngle(LightRay* L_ray, Terra* T, ponto_xy Pi, ponto_xy TargetP, long double StartingAng, long double EndingAng);
double FindHorizonAngle(LightRay* L_ray, Terra* T, ponto_xy P, long double angi, long double xM);
void DataShow(ponto_xy Obs, ponto_xy TangentP, ponto_xy TUp, ponto_xy TDown, ponto_xy RefP1, ponto_xy RefP2);


/* Variaveis globais para dados auxiliares. */
SDL_Renderer* renderer = NULL;


int main(int argc, char* argv[])
{
	Terra T; // Variavel com os parametros da terra
	int i, j, k=0, *n_pontos, n_pontos_total, posRT = 0, L_ray_view = 1, L_ray_refracted_view = 0; // variaveis de controle
	long double altM; //altura maxima da atmosfera
	long double esp; 		// espessura de cada camada
	LightRay L_ray1, *L_ray, *L_rayRefracted; 	// declaração do conjunto de direções que representara o raio de luz na atmosfera
	ponto_xy* P_0, Pcam, obs, target, Paux1, Paux2, Tangent_Line[2];			// pontos iniciais
	long double *ang_0;		// angulos de direcao iniciais para os raios de luz
	long double *xM, z;						// profundidade da camera
	long double expoente = 0;			// variavel ara o calculo dos indices de refracao			
	int luz_corT[3]; // parametros para cor
	long double m;

	luz_corT[0] = 0;
	luz_corT[1] = 255;
	luz_corT[2] = 25;

	target.x = 35.0;
	target.y = 6369.9537087717281 + 1.0/1000.0;


	/*************** alocacao de parametros ***********************/
	
	
	L_ray1.Data = (ponto_xy*)malloc(2000000 * sizeof(ponto_xy));
	
	L_ray = (LightRay*)malloc(21 * sizeof(LightRay));
	L_rayRefracted = (LightRay*)malloc(21 * sizeof(LightRay));

	for (i = 0; i < 21; i++)
	{
		L_ray[i].Data = (ponto_xy*)malloc(1013 * sizeof(ponto_xy));
		L_rayRefracted[i].Data = (ponto_xy*)malloc(2 * sizeof(ponto_xy));
	}

	T.camada_h = (long double*)malloc(500000 * sizeof(long double));

	T.n = (long double*)malloc(500000 * sizeof(long double));

	P_0 = (ponto_xy*)malloc(21 * sizeof(ponto_xy));
	ang_0 = (long double*)malloc(21 * sizeof(long double));

	xM = (long double*)malloc(21 * sizeof(long double));
	n_pontos = (int*)malloc(21 * sizeof(int));
	 /****************************************************************/

	/************ Inicializacao de parametros gerais************/

	obs.x = -5;
	obs.y = 6370;

	for (i = 0; i < 21; i++)
	{
		n_pontos[i] = 1013;
		xM[i] = target.x;
		P_0[i] = obs;
		ang_0[i] = 0;
	}
	
	for (i = 0; i < 21; i++)
	{
		L_ray[i].color[0] = 40;
		L_ray[i].color[1] = 230;
		L_ray[i].color[2] = 240;

		L_rayRefracted[i].color[0] = 40;
		L_rayRefracted[i].color[1] = 230;
		L_rayRefracted[i].color[2] = 240;
	}

	/***********************************************************/


	/*********** atribuicao de parametros iniciais reais *************/
		// Unidade de comprimento: Km
	T.C.x = 0;			// coordenada x do centro
	T.C.y = 0;			// coordenada y do centro
	T.raio = 6370;		// raio;
	altM = 0.1;			// altura maxima para as camadas. No caso, as camadas serao consideradas ate essa altura.
	T.atm_num = 100000;		// numero de camadas
	esp = altM / T.atm_num;	//calculo da espessura de cada camada
	/*****************************************************************/

	//Atribuindo as coordenas da camera
	Pcam.x = 35;
	Pcam.y = 6370;
	z = 0;
	
	//Atribuicao dos indices de refracao proximos ao solo
	for (i = 0; i < 200; i++)
	{
		T.n[199 - i] = 1.00035 - i * 0.00025 / 200;
		//if (i < 10 || i > 990) 
		printf("\n n(%d) = %.15lf", 199 - i, T.n[199 - i]);

	}

	// Atribuicao das alturas e indices de refracao das camadas
	for (i = 0; i < T.atm_num; i++)
	{
		T.camada_h[i] = 0.0;
		T.camada_h[i] = T.raio + (i + 1) * esp;			//calculo da altura de cada camada.
		expoente = -0.2 * esp * (i - 2);
		if (i >= 200) T.n[i] = 1 + 0.00035 * pow(EXP, expoente);

	}
	
	ang_0[0] = 0.00148;

	ang_0[20] = FindHorizonAngle(&L_ray1, &T, P_0[20], ang_0[20], xM[20]);
	
		
	L_ray1.color[0] = 40;
	L_ray1.color[1] = 230;
	L_ray1.color[2] = 240;
	//Calculo da variacao angular na varredura


	for (i = 0; i < 10; i++)
	{
		//Raio n refletido
		ang_0[i] = FindTargetPointAngle(&L_ray1, &T, P_0[i], target, 0.0014, ang_0[20]);
		
		// Calcula a trajetoria de cada raio na atmosfera e retorna o numero de raios calculados.
		n_pontos_total = teste_2(&L_ray1, &T, P_0[i], ang_0[i], xM[i], &posRT);
		LightRay_ResolutionConverter(L_ray1.Data, L_ray[i].Data, posRT, n_pontos_total, n_pontos[i]);

		Paux1 = L_ray[i].Data[0];
		Paux2 = L_ray[i].Data[1];
		m = (Paux2.y - Paux1.y) / (Paux2.x - Paux1.x);
		L_rayRefracted[i].Data[0] = Paux1;
		L_rayRefracted[i].Data[1].x = xM[i];
		L_rayRefracted[i].Data[1].y = m * (L_rayRefracted[i].Data[1].x - Paux1.x) + Paux1.y;

		L_ray[i].color[0] = L_ray1.color[0];
		L_ray[i].color[1] = L_ray1.color[1];
		L_ray[i].color[2] = L_ray1.color[2];
		
		L_rayRefracted[i].color[0] = L_ray1.color[0];
		L_rayRefracted[i].color[1] = L_ray1.color[1];
		L_rayRefracted[i].color[2] = L_ray1.color[2];

		L_ray1.color[0] = 40;
		L_ray1.color[1] = 230;
		L_ray1.color[2] = 240;

		posRT = 0;

		//Raio Refletido
		ang_0[i + 10] = FindTargetPointAngle(&L_ray1, &T, P_0[i + 10], target, -0.00025, ang_0[20]);

		// Calcula a trajetoria de cada raio na atmosfera e retorna o numero de raios calculados.
		n_pontos_total = teste_2(&L_ray1, &T, P_0[i + 10], ang_0[i + 10], xM[i + 10], &posRT);
		LightRay_ResolutionConverter(L_ray1.Data, L_ray[i + 10].Data, posRT, n_pontos_total, n_pontos[i]);

		Paux1 = L_ray[i + 10].Data[0];
		Paux2 = L_ray[i + 10].Data[1];
		m = (Paux2.y - Paux1.y) / (Paux2.x - Paux1.x);
		L_rayRefracted[i + 10].Data[0] = Paux1;
		L_rayRefracted[i + 10].Data[1].x = xM[i + 10];
		L_rayRefracted[i + 10].Data[1].y = m * (L_rayRefracted[i + 10].Data[1].x - Paux1.x) + Paux1.y;

		L_ray[i + 10].color[0] = L_ray1.color[0];
		L_ray[i + 10].color[1] = L_ray1.color[1];
		L_ray[i + 10].color[2] = L_ray1.color[2];

		L_rayRefracted[i + 10].color[0] = L_ray1.color[0];
		L_rayRefracted[i + 10].color[1] = L_ray1.color[1];
		L_rayRefracted[i + 10].color[2] = L_ray1.color[2];

		L_ray1.color[0] = 40;
		L_ray1.color[1] = 230;
		L_ray1.color[2] = 240;

		posRT = 0;

		target.y = target.y + 5.0 / 1000.0;
	}

	/*************/
	ang_0[i + 10] = 0.0013;

	ang_0[i + 10] = FindHorizonAngle(&L_ray1, &T, P_0[i + 10], ang_0[i + 10], xM[i + 10]);

	n_pontos_total = teste_2(&L_ray1, &T, P_0[i + 10], ang_0[i + 10], xM[i + 10], &posRT);
	LightRay_ResolutionConverter(L_ray1.Data, L_ray[i + 10].Data, posRT, n_pontos_total, n_pontos[i + 10]);

	Paux1 = L_ray[i + 10].Data[0];
	Paux2 = L_ray[i + 10].Data[1];
	m = (Paux2.y - Paux1.y) / (Paux2.x - Paux1.x);
	L_rayRefracted[i + 10].Data[0] = Paux1;
	L_rayRefracted[i + 10].Data[1].x = xM[i + 10];
	L_rayRefracted[i + 10].Data[1].y = m * (L_rayRefracted[i + 10].Data[1].x - Paux1.x) + Paux1.y;

	L_ray[i + 10].color[0] = L_ray1.color[0];
	L_ray[i + 10].color[1] = L_ray1.color[1];
	L_ray[i + 10].color[2] = L_ray1.color[2];

	L_rayRefracted[i + 10].color[0] = L_ray1.color[0];
	L_rayRefracted[i + 10].color[1] = L_ray1.color[1];
	L_rayRefracted[i + 10].color[2] = L_ray1.color[2];
	/*****************************/

	//Calculo do ponto de tangencia
	Tangent_Line[0] = obs;
	Tangent_Line[1] = Earth_tangent_point(obs, &T);

	// Mostraa as distancias relevantes
	//DataShow(obs, Tangent_Line[1], L_ray, );

	// Calculo do ponto de tangencia a partir do ponto de observacao
	Tangent_Line[1] = Earth_tangent_line(obs, &T).B;

	// Inicializa o SDL
	SDL_Init(SDL_INIT_VIDEO);

	// Cria uma janela
	SDL_Window* window = SDL_CreateWindow("Desenhar Linha com SDL", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 2500, 1400, SDL_WINDOW_SHOWN);
	if (window == NULL) {
		printf("Erro ao criar a janela: %s\n", SDL_GetError());
		return 1;
	}

	// Cria um renderer para a janela
	renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

	// ajusa as cores e atualiza o renderer
	SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
	SDL_RenderClear(renderer);
	SDL_RenderPresent(renderer);

	SDL_Event event;

	while (k == 0)
	{
		while (SDL_PollEvent(&event))
		{


			if (event.type == SDL_MOUSEWHEEL)
			{
				if (event.wheel.y < 0) // scroll down
				{
					z += 0.05;
				}

				if (event.wheel.y > 0) // scroll up
				{
					z -= 0.05;
					if (z < 0) z = 0;
				}

			}

			if (event.type == SDL_KEYDOWN)
			{
				if (event.key.keysym.sym == SDLK_ESCAPE) k = 1;

				if (event.key.keysym.sym == SDLK_w)
				{
					Pcam.y = Pcam.y + 0.002;

				}

				if (event.key.keysym.sym == SDLK_s)
				{
					Pcam.y = Pcam.y - 0.002;

				}

				if (event.key.keysym.sym == SDLK_d)
				{
					if (Pcam.x < target.x) Pcam.x = Pcam.x + 0.02;
					else Pcam.x = target.x;

				}

				if (event.key.keysym.sym == SDLK_a)
				{
					if (Pcam.x > obs.x) Pcam.x = Pcam.x - 0.02;
					else Pcam.x = obs.x;

				}

				if (event.key.keysym.sym == SDLK_r)
				{
					if (L_ray_view == 0) L_ray_view = 1;
					else L_ray_view = 0;

				}

				if (event.key.keysym.sym == SDLK_p)
				{
					if (L_ray_refracted_view == 0) L_ray_refracted_view = 1;
					else L_ray_refracted_view = 0;

				}

			}
		}


		viewport_Earth_render(&T, Pcam.x, Pcam.y, z);

		if (L_ray_view == 1)
		{
			for (i = 0; i < 21; i++)
			{
				viewport_LightRay_render(L_ray[i].Data, Pcam.x, Pcam.y, z, n_pontos[i], L_ray[i].color);
			}
		}

		if (L_ray_refracted_view == 1)
		{
			for (i = 0; i < 21; i++)
			{
				viewport_LightRay_render(L_rayRefracted[i].Data, Pcam.x, Pcam.y, z, 2, L_rayRefracted[i].color);
			}
		}
		
		viewport_LightRay_render(Tangent_Line, Pcam.x, Pcam.y, z, 2, luz_corT); //visualizacao da linha tangente
		SDL_RenderPresent(renderer);
		SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
		SDL_RenderClear(renderer);

		

		/*
		if (L_ray_refracted_view == 1)
		{
			viewport_Earth_render(&T, Pcam.x, Pcam.y, z);

			for (i = 0; i < 21; i++)
			{
				viewport_LightRay_render(L_rayRefracted[i].Data, Pcam.x, Pcam.y, z, 2, L_rayRefracted[i].color);
			}

			viewport_LightRay_render(Tangent_Line, Pcam.x, Pcam.y, z, 2, luz_corT); //visualizacao da linha tangente
			SDL_RenderPresent(renderer);
			SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
			SDL_RenderClear(renderer);
		}	
		*/

	}

	// Libera recursos
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();

	for (i = 0; i < 21; i++)
	{
		free(L_ray[i].Data);
		free(L_rayRefracted[i].Data);
	}

	free(L_ray);

	free(L_rayRefracted);

	free(L_ray1.Data);

	free(T.camada_h);

	free(T.n);

	free(P_0);
	free(ang_0);

	free(xM);
	free(n_pontos);

	return 0;
}

int teste_2(LightRay* L_ray, Terra* T, ponto_xy P, long double ang, long double xM, int* pos_ReflexaoTotal)
{
	int i;
	int camada = 0, aux;		// variaveis de controle para armazenar a camada em que o raio se encontra.
	char pos;				// variavel de controle. Verifica aa posicao do vetor diretor em relacao a a reta normal.
	long double a, b, c;			// parametros da equacao gerada para t no calculo das interseccoes.
	long double ang_i, ang_r;	// angulos de incidencia e refracao
	long double H = 0, h = 0;			// altura da camada suerior e altura da camada inferior ao ponto inicial
	long double t;				// parametro para a equação vetorial da reta
	long double d;				// distancia do ponto inicial ao centro da Terra (tambem usada como variavel auxiliar nos calculos)
	long double ang_reta = 0;
	ponto_xy v;				// vetor diretor do raio de luz


	/************** TRECHO DE CALCULOS PARA O PRIMEIRO RAIO DE LUZ *********************/

	// inicia o vetor diretor a partir do angulo de direcao inicial do raio de luz
	v.x = cosl(ang);
	v.y = sinl(ang);


	// Atribui o ponto inicial ao primeiro ponto do conjunto de retas do raio
	(*L_ray).Data[0] = P;


	// calculo da distancia do ponto inicial ao centro da Terra
	d = sqrtl(((*L_ray).Data[0].x - (*T).C.x) * ((*L_ray).Data[0].x - (*T).C.x) + ((*L_ray).Data[0].y - (*T).C.y) * ((*L_ray).Data[0].y - (*T).C.y));


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
	b = 2 * (v.x * ((*L_ray).Data[0].x - (*T).C.x) + v.y * ((*L_ray).Data[0].y - (*T).C.y));
	//c = pow((L_ray[0].x - (*T).C.x), 2) + pow((L_ray[0].y - (*T).C.y), 2) - H * H;
	c = ((*L_ray).Data[0].x - (*T).C.x) * ((*L_ray).Data[0].x - (*T).C.x) + ((*L_ray).Data[0].y - (*T).C.y) * ((*L_ray).Data[0].y - (*T).C.y) - H * H;

	t = (-b + sqrtl(b * b - 4 * a * c)) / (2 * a);		//calculo efetivo do parametro t
	aux = 1;											// variavel controle indicando que a proxima camada esta mais acima

	//calculo do parametro ponto de interseccao do raio com a camada inferior (caso exista)
	c = ((*L_ray).Data[0].x - (*T).C.x) * ((*L_ray).Data[0].x - (*T).C.x) + ((*L_ray).Data[0].y - (*T).C.y) * ((*L_ray).Data[0].y - (*T).C.y) - h * h;
	d = b * b - 4 * a * c;
	if (d >= 0)	//condicao de exeistencia
	{
		if ((-b - sqrtl(b * b - 4 * a * c)) / (2 * a) >= 0) //verifica se a intersecao correu antes do ponto inicial.
		{
			t = (-b - sqrtl(b * b - 4 * a * c)) / (2 * a);
			aux = -1;									// variavel controle indicando que a proxima camada esta mais abaixo.
		}

	}

	// atribuicao do ponto de interseccao com a camada calculada a partir do paramentro e do vetor diretor
	(*L_ray).Data[1].x = (*L_ray).Data[0].x + v.x * t;
	(*L_ray).Data[1].y = (*L_ray).Data[0].y + v.y * t;

	/*************************** FIM DO TRECHO ***********************************/
	
	i = 1;

	while (i < 3999999)
	{
		/*
		Esse laco calcula as direcoes dos raios de luz atraves das camadas.
		Segue a sequencia basica de passos:
		1) calcula o anguo de incidencia em relacao a a reta normal no ponto de incidencia. O ponto de incidencia eh o ponto A (inicial do raio atual).
		2) Verifica se havera refracao ou reflexao total.
			2.1) Caso haja refracao, o angulo de refracao eh calculado e a direcao do raio eh atualizada. Nesse caso, a camada tambem sera atualizada, pois indica que o raio atravessou.
			2.2) Caso haja reflexao total, o angulo de reflexao eh calculado e a direcao do raio de luz eh atualizada. Nesse caso a camada nao sera atualizada, pois o raio refletiu de volta para a mesma camada.
		3) Calcula o ponto de interseccao com a proxima camada.
			3.1) Caso a interseccao ocorra na camada acima, aux recebera 1, atualizando (posteriormente) a camada para a superior.
			3.2) Caso a interseccao ocorra na camada abaixo, aux recebera -1, atualizando (posteriormente) a camada para a inferior.
		4) Atualiza o valor do ponto B (ponto final do raio) conforme o calculo realizado anteriormente.
		5s) Repete o processo.
		Obs. A variavel aux indica apenas em qual camada ocorrera a interseccao com o raio de luz. aux = 1 caso ocorra na camada de cima e aux = -1 caso ocorra na camada de baixo.
		Esse controle eh necessario para que as camadas sejam atualizadas corretamente.

		*/

		// Calculo do ângulo de incidência no ponto pertencente a a superficie dioptrica
		ang_i = angulo_incidencia(v, (*L_ray).Data[i], (*T).C);


		/******** calculo do angulo de refracao pela Lei de Snell-Descartes. *********/
		// indica que haverá refração
		if ((*T).n[camada] * sinl(ang_i) / (*T).n[camada + aux] >= 0 && (*T).n[camada] * sinl(ang_i) / (*T).n[camada + aux] <= 1)
		{
			//Lei de Snell-Descartes
			ang_r = asinl((*T).n[camada] * sinl(ang_i) / (*T).n[camada + aux]);

			// analisa em qual posicao o vetor diretor esta em relacao ao vetor normal para determinar  o sentido de rotacao 			
			pos = posicao_vetor_raio(v, (*L_ray).Data[i], (*T).C);

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

			//printf("\n\n iteracao = %d \n camda_I = %d | camada_R = %d", i, camada, camada + aux);
			//printf("\n n1 = %.15lf | n2 = %.15lf | ang_I = %.15lf | ang_R = %.15lf", (*T).n[camada], (*T).n[camada + aux], ang_i, ang_r);

			//atualizacao da camada. Armazena a nova camada em que o raio se encontra.
			camada = camada + aux;
		}

		else //indica que haverá reflexão total
		{
			(*L_ray).color[0] = 255;
			(*L_ray).color[1] = 50;
			(*L_ray).color[2] = 50;

			if (camada > (*T).atm_num - 10) printf("\n\n Reflexão total superior na camada : %d \n iteracao: %d \n Coordenadas do ponto de reflexao (%.10lf, %.10lf) \n", camada, i, (*L_ray).Data[i].x, (*L_ray).Data[i].y);
			if (camada < 2000) printf("\n\n Reflexão total inferior na camada : %d \n iteracao: %d \n Coordenadas do ponto de reflexao (%.15lf, %.15lf) \n", camada, i, (*L_ray).Data[i].x, (*L_ray).Data[i].y);

			if ((*pos_ReflexaoTotal) == 0) (*pos_ReflexaoTotal) = i;

			pos = posicao_vetor_raio(v, (*L_ray).Data[i], (*T).C);

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
		b = 2 * (v.x * ((*L_ray).Data[i].x - (*T).C.x) + v.y * ((*L_ray).Data[i].y - (*T).C.y));
		c = ((*L_ray).Data[i].x - (*T).C.x) * ((*L_ray).Data[i].x - (*T).C.x) + ((*L_ray).Data[i].y - (*T).C.y) * ((*L_ray).Data[i].y - (*T).C.y) - H * H;

		t = (-b + sqrtl(b * b - 4 * a * c)) / (2 * a);		//calculo efetivo do parametro t
		aux = 1;											// variavel controle indicando que a proxima camada esta mais acima

		//calculo do parametro ponto de interseccao do raio com a camada inferior (caso exista)
		c = ((*L_ray).Data[i].x - (*T).C.x) * ((*L_ray).Data[i].x - (*T).C.x) + ((*L_ray).Data[i].y - (*T).C.y) * ((*L_ray).Data[i].y - (*T).C.y) - h * h;
		d = b * b - 4 * a * c;
		if (d >= 0)
		{
			if ((-b - sqrtl(b * b - 4 * a * c)) / (2 * a) > 0)
			{
				t = (-b - sqrtl(b * b - 4 * a * c)) / (2 * a);
				aux = -1;									// variavel controle indicando que a proxima camada esta mais abaixo.
			}

		}


		// atribuicao do ponto de interseccao com a camada calculada a partir do paramentro e do vetor diretor
		(*L_ray).Data[i + 1].x = (*L_ray).Data[i].x + v.x * t;
		(*L_ray).Data[i + 1].y = (*L_ray).Data[i].y + v.y * t;

		//printf("\n I: (%.15lf, %.15lf) (%.15lf, %.15lf)", L_ray[i - 1].A.x, L_ray[i - 1].A.y, L_ray[i - 1].B.x, L_ray[i - 1].B.y);
		//printf("\n R: (%.15lf, %.15lf) (%.15lf, %.15lf)", L_ray[i].A.x, L_ray[i].A.y, L_ray[i].B.x, L_ray[i].B.y);

		//vaux.x = (L_ray[i].B.x - L_ray[i].A.x) / sqrt((L_ray[i].B.x - L_ray[i].A.x) * (L_ray[i].B.x - L_ray[i].A.x) + (L_ray[i].B.y - L_ray[i].A.y) * (L_ray[i].B.y - L_ray[i].A.y));
		//vaux.y = (L_ray[i].B.y - L_ray[i].A.y) / sqrt((L_ray[i].B.x - L_ray[i].A.x) * (L_ray[i].B.x - L_ray[i].A.x) + (L_ray[i].B.y - L_ray[i].A.y) * (L_ray[i].B.y - L_ray[i].A.y));
		//printf("\n\n vetor diretor1: (%.15lf; %.15lf) \n vetor diretor2: (%.15lf; %.15lf)", v.x, v.y, vaux.x, vaux.y);

		//Verifica se o raio chegou ao ponto final
		if (((*L_ray).Data[i + 1].x <= xM && xM <= 0) || ((*L_ray).Data[i + 1].x >= xM && xM > 0))
		{
			ang_reta = ((*L_ray).Data[i + 1].y - (*L_ray).Data[i].y) / ((*L_ray).Data[i + 1].x - (*L_ray).Data[i].x);
			(*L_ray).Data[i + 1].x = xM;
			(*L_ray).Data[i + 1].y = ang_reta * (xM - (*L_ray).Data[i].x) + (*L_ray).Data[i].y;
			printf("\nangulo final: %.20lf", atanl(ang_reta) + pi);
			return (i + 2);
		}

		/*********************************** Fim do calculo ***************************************************************************************/

		i++;

		if ( (camada + aux) < 2 || camada > ((*T).atm_num - 5) )
		{
			return (i + 2);
		}


	}

	return (i + 1);
}



void LightRay_ResolutionConverter(ponto_xy* L_ray_i, ponto_xy* L_ray_f, int pos_ReflexaoTotal, int n_pontos_i, int n_pontos_f)
{
	int i = 1, j, * datasort, aux_change;
	long double delta_i;

	datasort = malloc(n_pontos_f * sizeof(int));

	delta_i = (long double)(n_pontos_i) / (long double)(n_pontos_f - 13);

	datasort[0] = 0;
	datasort[1] = 1;
	datasort[2] = 2;
	datasort[3] = 3;

	while (i <= n_pontos_f - 13)
	{
		datasort[i + 3] = (int)round(i * delta_i);
		if (datasort[i + 3] >= n_pontos_i)
		{
			datasort[i + 3] = (int)(datasort[i + 3] + datasort[i + 2]) / 2;
		}

		i++;
	}

	if (pos_ReflexaoTotal >= 2)
	{
		datasort[i + 3] = pos_ReflexaoTotal - 2;
		datasort[i + 4] = pos_ReflexaoTotal - 1;
		datasort[i + 5] = pos_ReflexaoTotal;
		datasort[i + 6] = pos_ReflexaoTotal + 1;
		datasort[i + 7] = pos_ReflexaoTotal + 2;
		datasort[i + 8] = n_pontos_i - 4;
		datasort[i + 9] = n_pontos_i - 3;
		datasort[i + 10] = n_pontos_i - 2;
		datasort[i + 11] = n_pontos_i - 1;
	}
	else
	{
		//caso nã haja reflexão total, o algoritmo preenche com 4, 5, 6, 7 e 8 para ocupar as posicoes que seriam
		//ocupadas pelos pontos da reflexão total. Atribui-se 4 em diante para  não gerar conflito com os pontos iniciais.
		//O array datasort será ordenado de acordo com oss indices.
		datasort[i + 3] = 4;
		datasort[i + 4] = 5;
		datasort[i + 5] = 6;
		datasort[i + 6] = 7;
		datasort[i + 7] = 8;
		datasort[i + 8] = n_pontos_i - 4;
		datasort[i + 9] = n_pontos_i - 3;
		datasort[i + 10] = n_pontos_i - 2;
		datasort[i + 11] = n_pontos_i - 1;
	}

	while (i < n_pontos_f - 3)
	{
		j = i + 3 - 1;

		while (datasort[j + 1] < datasort[j] && j > 0)
		{
			aux_change = datasort[j + 1];
			datasort[j + 1] = datasort[j];
			datasort[j] = aux_change;
			j--;
		}
		i++;
	}


	for (i = 3; i >= 1; i--)
	{
		j = i + 1;

		while (datasort[j] < datasort[j - 1] && j < n_pontos_f - 13)
		{
			aux_change = datasort[j - 1];
			datasort[j - 1] = datasort[j];
			datasort[j] = aux_change;
			j++;
		}
	}


	for (i = 0; i < n_pontos_f; i++)
	{
		L_ray_f[i] = L_ray_i[datasort[i]];
		//printf("\n datasort[%d] = %d", i, datasort[i]);
	}

	free(datasort);

}


long double angulo_incidencia(ponto_xy v, ponto_xy P, ponto_xy C)
{
	ponto_xy n;		// declaracao do vetor normal a a superficie no ponto P.
	long double ang_i;

	n.x = (P.x - C.x) / sqrtl((P.x - C.x) * (P.x - C.x) + (P.y - C.y) * (P.y - C.y));	// calculo da coordenada x do vetor normal.
	n.y = (P.y - C.y) / sqrtl((P.x - C.x) * (P.x - C.x) + (P.y - C.y) * (P.y - C.y));	// calculo da coordenada y do vetor normal.

	// calculo do angulo de incidencia pelo produto inerno entre o vetor diretor e o vetor normal a a camada.
	ang_i = acosl(v.x * n.x + v.y * n.y);

	// caso o angulo seja maior que 90 graus, deve-se considerar o vetor oposto ao veor normal anteriormente considerado.
	if (ang_i > pi / 2)
	{
		n.x = (-1) * n.x;
		n.y = (-1) * n.y;
		ang_i = acosl(v.x * n.x + v.y * n.y);
	}

	// retorna o angulo calculado.
	return ang_i;
}


ponto_xy rotacao_vetorial(ponto_xy v, long double ang)
{
	ponto_xy u;

	u.x = v.x * cosl(ang) - v.y * sinl(ang);
	u.y = v.x * sinl(ang) + v.y * cosl(ang);

	return u;
}


char posicao_vetor_raio(ponto_xy v, ponto_xy P, ponto_xy C)
{
	ponto_xy n;		// declaracao do vetor normal a a superficie no ponto P.
	long double ang;
	char c;

	n.x = (P.x - C.x) / sqrtl((P.x - C.x) * (P.x - C.x) + (P.y - C.y) * (P.y - C.y));	// calculo da coordenada x do vetor normal.
	n.y = (P.y - C.y) / sqrtl((P.x - C.x) * (P.x - C.x) + (P.y - C.y) * (P.y - C.y));	// calculo da cordenada y do vero normtal.

	// calculo do angulo de incidencia pelo produto inerno entre o vetor diretor e o vetor normal a a camada.
	ang = acosl(v.x * n.x + v.y * n.y);

	// caso o angulo seja maior que 90 graus, deve-se considerar o vetor oposto ao veor normal anteriormente considerado.
	if (ang > pi / 2)
	{
		n.x = (-1) * n.x;
		n.y = (-1) * n.y;
	}

	if (n.y >= 0)
	{
		ang = acosl(n.x);
	}
	else
	{
		ang = 2 * pi - acosl(n.x);
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


void viewport_LightRay_render(ponto_xy* L_ray, long double xpp, long double ypp, long double zpp, int nray, int luz_cor[3])
{
	int i;
	long double xi, yi, xf, yf, yauxi, yauxf, m = 0;
	long double d = 0.05;
	
	
	xi = L_ray[0].x;
	yi = L_ray[0].y;

	xi = xi - xpp;
	yi = yi - ypp;

	xi = (d / (d + zpp)) * xi;
	yi = (d / (d + zpp)) * yi;

	xi = (2000 / d) * xi + 1250;
	yi = -(2000 / d) * yi + 700;

	for (i = 1; i < nray; i++)
	{
		xf = L_ray[i].x;
		yf = L_ray[i].y;

		xf = xf - xpp;
		yf = yf - ypp;

		xf = (d / (d + zpp)) * xf;
		yf = (d / (d + zpp)) * yf;

		xf = (2000 / d) * xf + 1250;
		yf = -(2000 / d) * yf + 700;

		if (!(xi < 0 && xf < 0 || xi > 2500 && xf > 2500))
		{
			// Define a cor da linha do raio
			SDL_SetRenderDrawColor(renderer, luz_cor[0], luz_cor[1], luz_cor[2], 255);

			if ((xi < 0 && xf > 2500) || (xf < 0 && xi > 2500))			
			{
				m = (yf - yi) / (xf - xi);

				yauxi = m * (0 - xi) + yi;
				yauxf = m * (2500 - xi) + yi;

				// Desenha uma linha de raio no renderer
				drawLine(renderer, 0, yauxi, 2500, yauxf);

			}
			else
			{
				// Desenha uma linha de raio no renderer
				drawLine(renderer, xi, yi, xf, yf);
			}
		
		}
		xi = xf;
		yi = yf;
	}
		

}


void viewport_Earth_render(Terra* T, long double xpp, long double ypp, long double zpp)
{
	int i = 0, int_delta = 0, camada_n = 0, var = 50000;
	long double d = 0.05;
	long double int_set;
	long double Raio, CamadaH;
	ponto_xy Center;

	Raio = (*T).raio;
	Center = (*T).C;

	Center.x = Center.x - xpp;
	Center.y = Center.y - ypp;

	Center.x = (d / (d + zpp)) * Center.x;
	Center.y = (d / (d + zpp)) * Center.y;


	Center.x = (2000 / d) * Center.x + 1250;
	Center.y = -(2000 / d) * Center.y + 700;


	Raio = (d / (d + zpp)) * Raio;
	Raio = (2000 / d) * Raio;

	SDL_SetRenderDrawColor(renderer, 30, 50, 150, 255);
	drawCircle3(renderer, Center.x, Center.y, Raio);

	int_set = 0.00004;
	camada_n = (*T).atm_num;

	for (i = 0; i <= 10; i++)
	{
		CamadaH = (*T).camada_h[i];
		CamadaH = (d / (d + zpp)) * CamadaH;
		CamadaH = (2000 / d) * CamadaH;

		int_delta = i * int_set;

		SDL_SetRenderDrawColor(renderer, 100 - int_delta, 100 - int_delta, 100 - int_delta, 0);
		drawCircle2(renderer, Center.x, Center.y, CamadaH);

	}

	//var = 1;
	for (i = 0; i < camada_n; i += var)
	{
		CamadaH = (*T).camada_h[i];
		CamadaH = (d / (d + zpp)) * CamadaH;
		CamadaH = (2000 / d) * CamadaH;

		int_delta = i * int_set;

		SDL_SetRenderDrawColor(renderer, 100 - int_delta, 100 - int_delta, 100 - int_delta, 0);
		//SDL_SetRenderDrawColor(renderer, 200 , 200, 200, 0);

		drawCircle2(renderer, Center.x, Center.y, CamadaH);

	}
	CamadaH = (*T).camada_h[camada_n - 1];
	CamadaH = (d / (d + zpp)) * CamadaH;
	CamadaH = (2000 / d) * CamadaH;

	SDL_SetRenderDrawColor(renderer, 100, 100, 100, 0);
	drawCircle2(renderer, Center.x, Center.y, CamadaH);

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
			SDL_RenderDrawPoint(renderer, x1, y1 + 1);
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
	int x = 0;
	int y = radius;
	int d = 1 - radius;
	int deltaE = 3;
	int deltaSE = -2 * radius + 5;

	while (y >= x) {
		if (centerX + x >= 0 && centerX + x < 2500 && centerY + y >= 0 && centerY + y < 1400) {
			SDL_RenderDrawPoint(renderer, centerX + x, centerY + y);
		}

		if (centerX - x >= 0 && centerX - x < 2500 && centerY + y >= 0 && centerY + y < 1400) {
			SDL_RenderDrawPoint(renderer, centerX - x, centerY + y);
		}

		if (centerX + x >= 0 && centerX + x < 2500 && centerY - y >= 0 && centerY - y < 1400) {
			SDL_RenderDrawPoint(renderer, centerX + x, centerY - y);
		}

		if (centerX - x >= 0 && centerX - x < 2500 && centerY - y >= 0 && centerY - y < 1400) {
			SDL_RenderDrawPoint(renderer, centerX - x, centerY - y);
		}

		if (centerX + y >= 0 && centerX + y < 2500 && centerY + x >= 0 && centerY + x < 1400) {
			SDL_RenderDrawPoint(renderer, centerX + y, centerY + x);
		}

		if (centerX - y >= 0 && centerX - y < 2500 && centerY + x >= 0 && centerY + x < 1400) {
			SDL_RenderDrawPoint(renderer, centerX - y, centerY + x);
		}

		if (centerX + y >= 0 && centerX + y < 2500 && centerY - x >= 0 && centerY - x < 1400) {
			SDL_RenderDrawPoint(renderer, centerX + y, centerY - x);
		}

		if (centerX - y >= 0 && centerX - y < 2500 && centerY - x >= 0 && centerY - x < 1400) {
			SDL_RenderDrawPoint(renderer, centerX - y, centerY - x);
		}

		if (d < 0) {
			d += deltaE;
			deltaE += 2;
			deltaSE += 2;
		}
		else {
			d += deltaSE;
			deltaE += 2;
			deltaSE += 4;
			y--;
		}

		x++;
	}
}


void drawCircle2(SDL_Renderer* renderer, int centerX, int centerY, int Raio)
{
	long double x, y, xl;
	long double cos_ang, sin_ang, ang, R_inv;

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


void drawCircle3(SDL_Renderer* renderer, int centerX, int centerY, int Raio)
{
	long double x, y, xl;
	long double cos_ang, sin_ang, ang, R_inv;

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
		SDL_RenderDrawLine(renderer, centerX + x, centerY - y, centerX + x, 1400);
	}

}

ponto_xy Earth_tangent_point(ponto_xy PO, Terra* T)
{
	long double angP = 0, angT = 0, m = 0;
	ponto_xy Tangent_P = { 0.0, 0.0 };


	if (PO.x != 0)
	{
		angP = atan2(PO.y, PO.x); //calcula o angulo do ponto

	}

	else
	{
		angP = pi / 2;
	}

	angT = acos((*T).raio / sqrt(PO.x * PO.x + PO.y * PO.y)); //calcula o angulo do triangulo formado com o ponto, raio e reta tangente

	m = -1 / tan(angP - angT); //calcula o coeficiente angular da reta tangente

	Tangent_P.x = (m * m * PO.x - m * PO.y) / (1 + m * m);
	Tangent_P.y = m * (Tangent_P.x - PO.x) + PO.y;

	return Tangent_P;
}

reta_xy Earth_tangent_line(ponto_xy PO, Terra* T)
{
	reta_xy Tangent_Line;

	Tangent_Line.A = PO;
	Tangent_Line.B = Earth_tangent_point(PO, T);

	printf("\n\n Ponto de tangencia: (%.10lf, %.10lf) \n", Tangent_Line.B.x, Tangent_Line.B.y);


	Tangent_Line.B.x = Tangent_Line.B.x + 15 * (Tangent_Line.B.x - Tangent_Line.A.x);
	Tangent_Line.B.y = Tangent_Line.B.y + 15 * (Tangent_Line.B.y - Tangent_Line.A.y);

	return Tangent_Line;
}

double FindTargetPointAngle(LightRay* L_ray, Terra* T, ponto_xy Pi, ponto_xy TargetP, long double StartingAng, long double EndingAng)
{
	long double Ang_A, Ang_B, ang_aux, delta_ang, factor = 1, Distance, Dist_A, Dist_B;
	int TotalReflection_pos = 0;
	unsigned int n_points;
	ponto_xy Pa, Pb, Paux;

	Ang_A = StartingAng;

	Ang_B = EndingAng;

	n_points = teste_2(L_ray, T, Pi, Ang_A, TargetP.x, &TotalReflection_pos);

	Pa = (*L_ray).Data[n_points - 1];

	n_points = teste_2(L_ray, T, Pi, Ang_B, TargetP.x, &TotalReflection_pos);

	Pb = (*L_ray).Data[n_points - 1];
	
	if ((Pa.y - TargetP.y) * (Pb.y - TargetP.y) > 0)
	{
		return 10;
	}
	
	
	while ( ((Pa.y - Pb.y) > 0.00000001 || (Pb.y - Pa.y) > 0.00000001) && ((Ang_A - Ang_B) > 0.000000000000001 || (Ang_B - Ang_A) > 0.000000000000001))
	{
		ang_aux = (Ang_A + Ang_B) / 2;

		n_points = teste_2(L_ray, T, Pi, ang_aux, TargetP.x, &TotalReflection_pos);

		Paux = (*L_ray).Data[n_points - 1];

		if ((Paux.y - TargetP.y) * (Pb.y - TargetP.y) <= 0)
		{
			Pa = Paux;
			Ang_A = ang_aux;
		}
		else
		{
			Pb = Paux;
			Ang_B = ang_aux;
		}


	}
	
	Dist_A = Pa.y - TargetP.y;
	Dist_B = Pb.y - TargetP.y;

	if (Dist_A < 0) Dist_A = Dist_A * (-1.0);
	if (Dist_B < 0) Dist_B = Dist_B * (-1.0);

	if (Dist_A < Dist_B)
	{
		ang_aux = Ang_A;
	}
	else
	{
		ang_aux = Ang_B;
	}

	return ang_aux;

}

double FindHorizonAngle(LightRay* L_ray, Terra* T, ponto_xy Pi, long double angi, long double xM)
{
	double ang_aux, factor = 1;
	int TotalReflection_pos = 0;
	
	ang_aux = angi; // angulo inicial

	// Verificar se esta acima ou abaixo do angulo critico
	// Caso esteja acima, nao havera reflexao total, caso esteja abaixo, havera reflexao total.
	teste_2(L_ray, T, Pi, ang_aux, xM, &TotalReflection_pos);
	
	//O laco ocorre ate que a precisao seja alcancada 
	while (factor > (1.0 / (1000000000000.0)) )
	{
		//Esse condicional indica que nao houve reflexao total, entao o angulo deve ser diminuido
		if (TotalReflection_pos == 0)
		{
			while (TotalReflection_pos == 0)
			{
				ang_aux = ang_aux - 0.0005 * factor;
				teste_2(L_ray, T, Pi, ang_aux, xM, &TotalReflection_pos);
			}
		}

		factor = factor / 10.0;

		//Esse condicional indica que houve reflexao total, entao o angulo deve ser aumentado
		if (TotalReflection_pos != 0)
		{
			while (TotalReflection_pos != 0)
			{
				(*L_ray).color[0] = 40;
				(*L_ray).color[1] = 40;
				(*L_ray).color[2] = 240;
				TotalReflection_pos = 0; // Recebe zero para calcular o proximo valor.
				ang_aux = ang_aux + 0.0005 * factor;
				teste_2(L_ray, T, Pi, ang_aux, xM, &TotalReflection_pos);
				ang_aux = ang_aux;
			}
		}

		factor = factor / 10.0;
	}
	//Apos atingir a precisao de convergencia, a funcao retorna o angulo calculado.
	return ang_aux;
}

void DataShow(ponto_xy Obs, ponto_xy TangentP, ponto_xy TUp, ponto_xy TDown, ponto_xy RefP1, ponto_xy RefP2)
{
	long double dist, ang;

	//altura do observador
	dist = sqrt(Obs.x * Obs.x + Obs.y * Obs.y) - 6370;
	printf("\n Altura do ponto de observação: %.6lf km", dist);

	//Altura da base do alvo
	dist = sqrt(TDown.x * TDown.x + TDown.y * TDown.y) - 6370;
	printf("\n Altura da base do alvo: %.3lf km", dist);

	//Altura do topo do alvo
	dist = sqrt(TUp.x * TUp.x + TUp.y * TUp.y) - 6370;
	printf("\n Altura do topo do alvo: %.3lf km", dist);

	// Distancia alvo-observador
	ang = acos((Obs.x * TDown.x + Obs.y * TDown.y) / (sqrt(Obs.x * Obs.x + Obs.y * Obs.y) * sqrt(TDown.x * TDown.x + TDown.y * TDown.y)));

	dist = ang * 6370.0;

	printf("\n Distancia observadr-alvo: %.3lf km", dist);


	// distancia ponto de tangencia
	ang = acos((Obs.x * TangentP.x + Obs.y * TangentP.y) / (sqrt(Obs.x * Obs.x + Obs.y * Obs.y) * sqrt(TangentP.x * TangentP.x + TangentP.y * TangentP.y)));

	dist = ang * 6370.0;

	printf("\n Distancia da linha do horizonte: %.3lf km", dist);


	// distancia ponto de tangencia
	ang = acos((Obs.x * RefP2.x + Obs.y * RefP2.y) / (sqrt(Obs.x * Obs.x + Obs.y * Obs.y) * sqrt(RefP2.x * RefP2.x + RefP2.y * RefP2.y)));

	dist = ang * 6370.0;

	printf("\n Distancia do ponto de reflexao do raio inferior: %.3lf km", dist);

	// distancia ponto de tangencia
	ang = acos((Obs.x * RefP1.x + Obs.y * RefP1.y) / (sqrt(Obs.x * Obs.x + Obs.y * Obs.y) * sqrt(RefP1.x * RefP1.x + RefP1.y * RefP1.y)));

	dist = ang * 6370.0;

	printf("\n Distancia do ponto de reflexao do raio inferior: %.3lf km", dist);

}


