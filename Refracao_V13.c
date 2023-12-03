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


/* Definicao das estruturas 3D */
typedef struct
{
	double x; 		// define a coordenada x de um ponto no espaco
	double y;		// define a coordenada y de um ponto no espaco
	double z;		// define a coordenada z de um ponto no espaco
} XYZPoint;

typedef struct
{
	XYZPoint* WorldVertexList; //Lista de vertices no espaco do mundo
	XYZPoint* CamCoordSysVertexList; // Lista de vertices no espaco da camera
	XYZPoint* ProjectionSpaceVertexList; // Lista de vertices no espaco de projecao
	char* VisibleVertices; //Lista de controle de vertices visiveis. Recebe o valor T na posicao i indicanco que o vertice i da lista e visivel.
	int* TriangleList; // Lista de triangulos no espaco do mundo
	int* PSTriangleList; //Lista de triangulos do espaço de projecao
	int* TriangleColor; //Cor da superficie
	int* PSTriangleColor; //Cor da superficie
	int Vertex_N; // numero de vertices no espaco do mundo
	int Triangle_N; // numero de triangulos no espaco do mundo
	int PSVertex_N; // numero de vertices no espaco de projecao
	int PSTriangle_N; // numero de triangulos no espaco de projecao
} XYZObject;

typedef struct
{
	XYZPoint base[3];
	XYZPoint Origin;
} XYZCoordinateSystem;

typedef struct
{
	double Data[3][3];
} XYZTransformatioMatrix;

typedef struct
{
	XYZCoordinateSystem camCS;
	double ProjectionPlane_Z;
	double NearPlane_Z;
	double FarPlane_Z;
	double ProjectioPlaneWidth;
	double ProjectionPlaneHeight;
	double WidthDepthRatio;
	double HeightDepthRatio;
	double BB_X; //Bounding Box x coordinate (-BB_X, 0) and (+BB_X, 0)
	double BB_Y; //Bounding Box y coordinate (0, -BB_y) and (0, +BB_Y)
} XYZCam;

//Funcoes de transformacao geometrica
XYZPoint XYZ_PointRotation(XYZPoint P, XYZPoint O, double ang);
XYZPoint XYZ_PointTranslation(XYZPoint P, XYZPoint Vector, double t);
XYZTransformatioMatrix XYZ_SetRotationMatrix(XYZPoint RotationAxis, double ang);
XYZPoint XYZ_ApplyTransformationMatrix(XYZPoint P, XYZTransformatioMatrix M);

//Funcoes de alocacao
void XYZ_AllocateObjectLists(XYZObject* OBJ, int N_Vertices, int N_Triangles);
void XYZ_DeallocatedObjectLists(XYZObject* OBJ);

//Funcoes de pipeline
void XYZ_SetOBJCamCoordinates(XYZObject* OBJ, XYZCam* Cam);
void XYZ_ApplyTriangleCulling(XYZObject* OBJ, XYZCam* Cam);
void XYZ_SetOBJProjectionSpaceCoordinate(XYZObject* OBJ, XYZCam* Cam);
void XYZ_ViewportOBJRender(XYZObject* OBJ, XYZCam* Cam);
void XYZ_ViewportOBJRender2(XYZObject* OBJ, XYZCam* Cam);
void XYZ_ViewportVisibleGroundRender(XYZObject* OBJ, XYZCam* Cam);
void XYZ_ViewportGeometricHorizon(XYZObject* OBJ, XYZCam* Cam);
void XYZ_drawLine(SDL_Renderer* renderer, int xi, int yi, int xf, int yf);
int LightRay_RenderObjectFrontView(LightRay* L_ray1, LightRay* L_ray2, ponto_xy Obs, ponto_xy GeometricHorizon, ponto_xy ApparentHorizon, int nray);

//Funcoes do Simulador de refracao
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
SDL_Window* window;


int main(int argc, char* argv[])
{
	Terra T; // Variavel com os parametros da terra
	int i, j, k = 0, * n_pontos, n_pontos_total, posRT = 0, L_ray_view = 1, L_ray_refracted_view = 0, FrontView = 0; // variaveis de controle
	long double altM; //altura maxima da atmosfera
	long double esp; 		// espessura de cada camada
	LightRay L_ray1, * L_ray, * L_rayRefracted; 	// declaração do conjunto de direções que representara o raio de luz na atmosfera
	ponto_xy* P_0, Pcam, obs, target, Paux1, Paux2, Tangent_Line[2], TangentPoint;			// pontos iniciais
	long double* ang_0;		// angulos de direcao iniciais para os raios de luz
	long double* xM, z;						// profundidade da camera
	long double expoente = 0;			// variavel ara o calculo dos indices de refracao			
	int luz_corT[3]; // parametros para cor
	long double m;
	long double TargetDistance, TargetSize, AngularDistance;

	luz_corT[0] = 0;
	luz_corT[1] = 255;
	luz_corT[2] = 25;

	//coordenadas do alvo
	target.x = 40;
	target.y = 6369.9537087717281 + 1.0 / 1000.0;

	printf("Target Distance:");
	scanf_s("%lf", &TargetDistance);
	printf("Target Size:");
	scanf_s("%lf", &TargetSize);

	target.x = TargetDistance;
	/*************** alocacao de parametros ***********************/


	L_ray1.Data = (ponto_xy*)malloc(2000000 * sizeof(ponto_xy)); // raio principal que armazenara toda a informacao

	L_ray = (LightRay*)malloc(21 * sizeof(LightRay)); // conjunto de raios resolvidos, de menor resolucao
	L_rayRefracted = (LightRay*)malloc(21 * sizeof(LightRay));  // Projecao do ultimo raio considerado

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

   // Coordenadas do observador
	obs.x = -5;
	obs.y = 6370;


	for (i = 0; i < 21; i++)
	{
		n_pontos[i] = 1013; // numero de pontos de cada raio
		xM[i] = target.x; //posicao x de cada raio
		P_0[i] = obs; // posicao inicial de cada raio
		ang_0[i] = 0; // angulo inicial
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
	Pcam.x = target.x;
	Pcam.y = 6370;
	z = 0;

	//Atribuicao dos indices de refracao proximos ao solo
	for (i = 0; i < 200; i++)
	{
		T.n[199 - i] = 1.00030 - i * 0.00025 / 200;
		//if (i < 10 || i > 990) 
		printf("\n n(%d) = %.15lf", 199 - i, T.n[199 - i]);

	}

	// Atribuicao das alturas e indices de refracao das camadas
	for (i = 0; i < T.atm_num; i++)
	{
		T.camada_h[i] = 0.0;
		T.camada_h[i] = T.raio + (i + 1) * esp;			//calculo da altura de cada camada.
		expoente = -0.14 * esp * (i - 2);
		if (i >= 200) T.n[i] = 1 + 0.00030 * pow(EXP, expoente);

	}

	/****** Neste trecho o raio mais baixo e determinado. Esse raio representara a linha do horizonte refratada *******/
	// A rotina e a mesma para os outros raios
	ang_0[20] = 0.00148;

	ang_0[20] = FindHorizonAngle(&L_ray1, &T, P_0[20], ang_0[20], xM[20]); //encontra o raio de luz do horizonte aparente, ou seja, o raio de luz mais baixo possivel.

	n_pontos_total = teste_2(&L_ray1, &T, P_0[20], ang_0[20], xM[20], &posRT);
	LightRay_ResolutionConverter(L_ray1.Data, L_ray[20].Data, posRT, n_pontos_total, n_pontos[20]);


	Paux1 = L_ray[20].Data[0];
	Paux2 = L_ray[20].Data[1];
	m = (Paux2.y - Paux1.y) / (Paux2.x - Paux1.x);
	L_rayRefracted[20].Data[0] = Paux1;
	L_rayRefracted[20].Data[1].x = xM[20];
	L_rayRefracted[20].Data[1].y = m * (L_rayRefracted[20].Data[1].x - Paux1.x) + Paux1.y;

	L_ray[20].color[0] = L_ray1.color[0];
	L_ray[20].color[1] = L_ray1.color[1];
	L_ray[20].color[2] = L_ray1.color[2];

	L_rayRefracted[20].color[0] = L_ray1.color[0];
	L_rayRefracted[20].color[1] = L_ray1.color[1];
	L_rayRefracted[20].color[2] = L_ray1.color[2];
	/*****************************/

	//Ajustando a coordenada y do alvo
	target.y = L_ray[20].Data[1012].y + 1.0 / 1000.0;


	L_ray1.color[0] = 40;
	L_ray1.color[1] = 230;
	L_ray1.color[2] = 240;

	//Neste laco os raios sao calculados e as suas coordenadas sao armazenadas nos arrays  L_ray[i].Data
	for (i = 0; i < 10; i++)
	{
		//Raio n refletido
		ang_0[i] = FindTargetPointAngle(&L_ray1, &T, P_0[i], target, 0.00958, ang_0[20]); //Encontra o angulo do raio que atingira o alvo no ponto target

		// Calcula a trajetoria de cada raio na atmosfera e retorna o numero de raios calculados.
		n_pontos_total = teste_2(&L_ray1, &T, P_0[i], ang_0[i], xM[i], &posRT);
		LightRay_ResolutionConverter(L_ray1.Data, L_ray[i].Data, posRT, n_pontos_total, n_pontos[i]); //Armazena os pontos importantes nos raios de menor resolucao

		// Neste trecho, a rojecao do raio que chega ao observador e calculada.
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

		// Fim dos calculos para os raios que nao sofrem reflexao total


		//Inicio dos calculos para os raios que sofrem a reflexao total.

		//Identico ao anterior
		posRT = 0;

		//Raio Refletido
		ang_0[i + 10] = FindTargetPointAngle(&L_ray1, &T, P_0[i + 10], target, -0.02, ang_0[20]);

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

	//Calculo do ponto de tangencia
	Tangent_Line[0] = obs;
	Tangent_Line[1] = Earth_tangent_point(obs, &T);
	TangentPoint = Earth_tangent_point(obs, &T);
	// Calculo do ponto de tangencia a partir do ponto de observacao
	Tangent_Line[1] = Earth_tangent_line(obs, &T).B;

	// Inicializa o SDL
	SDL_Init(SDL_INIT_VIDEO);

	// Cria uma janela
	window = SDL_CreateWindow("Desenhar Linha com SDL", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 2500, 1400, SDL_WINDOW_SHOWN);

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

				if (event.key.keysym.sym == SDLK_f)
				{
					LightRay_RenderObjectFrontView(L_ray, L_rayRefracted, obs, TangentPoint, L_rayRefracted[20].Data[1], 21);
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

	while (i < 8999999)
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

		if ((camada + aux) < 2 || camada > ((*T).atm_num - 5))
		{
			return (i + 1);
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


	while (((Pa.y - Pb.y) > 0.00000001 || (Pb.y - Pa.y) > 0.00000001) && ((Ang_A - Ang_B) > 0.000000000000001 || (Ang_B - Ang_A) > 0.000000000000001))
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
	while (factor > (1.0 / (1000000000000.0)))
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






















/************************ Funcoes 3D *******************************/

int LightRay_RenderObjectFrontView(LightRay* L_ray1, LightRay* L_ray2, ponto_xy Obs, ponto_xy GeometricHorizon, ponto_xy ApparentHorizon, int nray) //test
{
	int i, j, k, MouseMH, MouseMV, aux_i, aux_j, RotControl = 0, CubeColorOption;
	//int fat = 10;
	XYZCam Cam1;
	XYZObject RealObject, Object1, Object2, ApparentHorizonOBJ, GeometricHorizonOBJ;
	XYZTransformatioMatrix M, M_Esfera;
	XYZPoint P, RotationAxis, WalkAxis;
	SDL_Event event;
	double RotCamAngAzimutal, RotCamAngPolar, RotationAng;
	XYZCoordinateSystem InertialCoordinateSystem;


	InertialCoordinateSystem.Origin.x = 0;
	InertialCoordinateSystem.Origin.y = 0;
	InertialCoordinateSystem.Origin.z = 0;

	InertialCoordinateSystem.base[0].x = 1;
	InertialCoordinateSystem.base[0].y = 0;
	InertialCoordinateSystem.base[0].z = 0;

	InertialCoordinateSystem.base[1].x = 0;
	InertialCoordinateSystem.base[1].y = 1;
	InertialCoordinateSystem.base[1].z = 0;

	InertialCoordinateSystem.base[2].x = 0;
	InertialCoordinateSystem.base[2].y = 0;
	InertialCoordinateSystem.base[2].z = 1;

	RotationAxis.x = 0;
	RotationAxis.y = 1;
	RotationAxis.z = 0;

	/**************** Alocacao de parametros ****************/

	RealObject.Vertex_N = 20;
	RealObject.Triangle_N = 18;
	XYZ_AllocateObjectLists(&RealObject, RealObject.Vertex_N, RealObject.Triangle_N);

	Object1.Vertex_N = 20;
	Object1.Triangle_N = 18;
	XYZ_AllocateObjectLists(&Object1, Object1.Vertex_N, Object1.Triangle_N);

	Object2.Vertex_N = 20;
	Object2.Triangle_N = 18;
	XYZ_AllocateObjectLists(&Object2, Object2.Vertex_N, Object2.Triangle_N);



	ApparentHorizonOBJ.WorldVertexList = (XYZPoint*)malloc(1 * sizeof(XYZPoint));
	ApparentHorizonOBJ.CamCoordSysVertexList = (XYZPoint*)malloc(1 * sizeof(XYZPoint));
	ApparentHorizonOBJ.ProjectionSpaceVertexList = (XYZPoint*)malloc(1 * sizeof(XYZPoint));
	ApparentHorizonOBJ.VisibleVertices = (char*)malloc(1 * sizeof(char));

	ApparentHorizonOBJ.TriangleList = NULL;
	ApparentHorizonOBJ.PSTriangleList = NULL;
	ApparentHorizonOBJ.TriangleColor = NULL;
	ApparentHorizonOBJ.PSTriangleColor = NULL;

	ApparentHorizonOBJ.Vertex_N = 1;
	ApparentHorizonOBJ.Triangle_N = 0;
	ApparentHorizonOBJ.PSVertex_N = 1;
	ApparentHorizonOBJ.PSTriangle_N = 0;

	ApparentHorizonOBJ.WorldVertexList[0].x = ApparentHorizon.x;
	ApparentHorizonOBJ.WorldVertexList[0].y = ApparentHorizon.y;
	ApparentHorizonOBJ.WorldVertexList[0].z = 0;

	ApparentHorizonOBJ.VisibleVertices[0] = 'T';


	GeometricHorizonOBJ.WorldVertexList = (XYZPoint*)malloc(1 * sizeof(XYZPoint));
	GeometricHorizonOBJ.CamCoordSysVertexList = (XYZPoint*)malloc(1 * sizeof(XYZPoint));
	GeometricHorizonOBJ.ProjectionSpaceVertexList = (XYZPoint*)malloc(1 * sizeof(XYZPoint));
	GeometricHorizonOBJ.VisibleVertices = (char*)malloc(1 * sizeof(char));

	GeometricHorizonOBJ.TriangleList = NULL;
	GeometricHorizonOBJ.PSTriangleList = NULL;
	GeometricHorizonOBJ.TriangleColor = NULL;
	GeometricHorizonOBJ.PSTriangleColor = NULL;

	GeometricHorizonOBJ.Vertex_N = 1;
	GeometricHorizonOBJ.Triangle_N = 0;
	GeometricHorizonOBJ.PSVertex_N = 1;
	GeometricHorizonOBJ.PSTriangle_N = 0;

	GeometricHorizonOBJ.WorldVertexList[0].x = GeometricHorizon.x;
	GeometricHorizonOBJ.WorldVertexList[0].y = GeometricHorizon.y;
	GeometricHorizonOBJ.WorldVertexList[0].z = 0;

	GeometricHorizonOBJ.VisibleVertices[0] = 'T';

	/**************** Fim da alocacao ****************/

	//Definicao manual da camera
	Cam1.camCS.base[0].x = 0;
	Cam1.camCS.base[0].y = 0;
	Cam1.camCS.base[0].z = -1;

	Cam1.camCS.base[1].x = 0;
	Cam1.camCS.base[1].y = 1;
	Cam1.camCS.base[1].z = 0;

	Cam1.camCS.base[2].x = 1;
	Cam1.camCS.base[2].y = 0;
	Cam1.camCS.base[2].z = 0;

	Cam1.camCS.Origin.x = Obs.x;
	Cam1.camCS.Origin.y = Obs.y;
	Cam1.camCS.Origin.z = 0;

	Cam1.ProjectionPlane_Z = 0.0005;
	Cam1.NearPlane_Z = 0.002;
	Cam1.FarPlane_Z = 500;
	Cam1.ProjectioPlaneWidth = 0.0005;
	Cam1.ProjectionPlaneHeight = 0.00028;

	Cam1.WidthDepthRatio = (Cam1.ProjectioPlaneWidth / 2) / Cam1.ProjectionPlane_Z;
	Cam1.HeightDepthRatio = (Cam1.ProjectionPlaneHeight / 2) / Cam1.ProjectionPlane_Z;

	WalkAxis = Cam1.camCS.base[2];





	RealObject.WorldVertexList[0].x = L_ray1[0].Data[1].x;
	RealObject.WorldVertexList[0].y = L_ray1[0].Data[1].y;
	RealObject.WorldVertexList[0].z = -0.02;

	RealObject.WorldVertexList[1].x = L_ray1[0].Data[1].x;
	RealObject.WorldVertexList[1].y = L_ray1[0].Data[1].y;
	RealObject.WorldVertexList[1].z = 0.02;

	RealObject.WorldVertexList[2].x = L_ray1[1].Data[1].x;
	RealObject.WorldVertexList[2].y = L_ray1[1].Data[1].y;
	RealObject.WorldVertexList[2].z = -0.02;

	RealObject.WorldVertexList[3].x = L_ray1[1].Data[1].x;
	RealObject.WorldVertexList[3].y = L_ray1[1].Data[1].y;
	RealObject.WorldVertexList[3].z = 0.02;

	RealObject.WorldVertexList[4].x = L_ray1[2].Data[1].x;
	RealObject.WorldVertexList[4].y = L_ray1[2].Data[1].y;
	RealObject.WorldVertexList[4].z = -0.02;

	RealObject.WorldVertexList[5].x = L_ray1[2].Data[1].x;
	RealObject.WorldVertexList[5].y = L_ray1[2].Data[1].y;
	RealObject.WorldVertexList[5].z = 0.02;

	RealObject.WorldVertexList[6].x = L_ray1[3].Data[1].x;
	RealObject.WorldVertexList[6].y = L_ray1[3].Data[1].y;
	RealObject.WorldVertexList[6].z = -0.02;

	RealObject.WorldVertexList[7].x = L_ray1[3].Data[1].x;
	RealObject.WorldVertexList[7].y = L_ray1[3].Data[1].y;
	RealObject.WorldVertexList[7].z = 0.02;

	RealObject.WorldVertexList[8].x = L_ray1[4].Data[1].x;
	RealObject.WorldVertexList[8].y = L_ray1[4].Data[1].y;
	RealObject.WorldVertexList[8].z = -0.02;

	RealObject.WorldVertexList[9].x = L_ray1[4].Data[1].x;
	RealObject.WorldVertexList[9].y = L_ray1[4].Data[1].y;
	RealObject.WorldVertexList[9].z = 0.02;

	RealObject.WorldVertexList[10].x = L_ray1[5].Data[1].x;
	RealObject.WorldVertexList[10].y = L_ray1[5].Data[1].y;
	RealObject.WorldVertexList[10].z = -0.02;

	RealObject.WorldVertexList[11].x = L_ray1[5].Data[1].x;
	RealObject.WorldVertexList[11].y = L_ray1[5].Data[1].y;
	RealObject.WorldVertexList[11].z = 0.02;

	RealObject.WorldVertexList[12].x = L_ray1[6].Data[1].x;
	RealObject.WorldVertexList[12].y = L_ray1[6].Data[1].y;
	RealObject.WorldVertexList[12].z = -0.02;

	RealObject.WorldVertexList[13].x = L_ray1[6].Data[1].x;
	RealObject.WorldVertexList[13].y = L_ray1[6].Data[1].y;
	RealObject.WorldVertexList[13].z = 0.02;

	RealObject.WorldVertexList[14].x = L_ray1[7].Data[1].x;
	RealObject.WorldVertexList[14].y = L_ray1[7].Data[1].y;
	RealObject.WorldVertexList[14].z = -0.02;

	RealObject.WorldVertexList[15].x = L_ray1[7].Data[1].x;
	RealObject.WorldVertexList[15].y = L_ray1[7].Data[1].y;
	RealObject.WorldVertexList[15].z = 0.02;

	RealObject.WorldVertexList[16].x = L_ray1[8].Data[1].x;
	RealObject.WorldVertexList[16].y = L_ray1[8].Data[1].y;
	RealObject.WorldVertexList[16].z = -0.02;

	RealObject.WorldVertexList[17].x = L_ray1[8].Data[1].x;
	RealObject.WorldVertexList[17].y = L_ray1[8].Data[1].y;
	RealObject.WorldVertexList[17].z = 0.02;

	RealObject.WorldVertexList[18].x = L_ray1[9].Data[1].x;
	RealObject.WorldVertexList[18].y = L_ray1[9].Data[1].y;
	RealObject.WorldVertexList[18].z = -0.02;

	RealObject.WorldVertexList[19].x = L_ray1[9].Data[1].x;
	RealObject.WorldVertexList[19].y = L_ray1[9].Data[1].y;
	RealObject.WorldVertexList[19].z = 0.02;

	j = 0;
	for (i = 0; i < (3 * Object1.Triangle_N); i += 6)
	{
		RealObject.TriangleList[i + 0] = 1 + j;
		RealObject.TriangleList[i + 1] = 0 + j;
		RealObject.TriangleList[i + 2] = 2 + j;

		RealObject.TriangleList[i + 3] = 2 + j;
		RealObject.TriangleList[i + 4] = 3 + j;
		RealObject.TriangleList[i + 5] = 1 + j;
		j += 2;
	}

	for (i = 0; i < (3 * Object1.Triangle_N); i += 3)
	{
		RealObject.TriangleColor[i] = L_ray1[0].color[0];
		RealObject.TriangleColor[i + 1] = L_ray1[0].color[1];
		RealObject.TriangleColor[i + 2] = L_ray1[0].color[2];
	}








	Object1.WorldVertexList[0].x = L_ray2[0].Data[1].x;
	Object1.WorldVertexList[0].y = L_ray2[0].Data[1].y;
	Object1.WorldVertexList[0].z = -0.02;

	Object1.WorldVertexList[1].x = L_ray2[0].Data[1].x;
	Object1.WorldVertexList[1].y = L_ray2[0].Data[1].y;
	Object1.WorldVertexList[1].z = 0.02;

	Object1.WorldVertexList[2].x = L_ray2[1].Data[1].x;
	Object1.WorldVertexList[2].y = L_ray2[1].Data[1].y;
	Object1.WorldVertexList[2].z = -0.02;

	Object1.WorldVertexList[3].x = L_ray2[1].Data[1].x;
	Object1.WorldVertexList[3].y = L_ray2[1].Data[1].y;
	Object1.WorldVertexList[3].z = 0.02;

	Object1.WorldVertexList[4].x = L_ray2[2].Data[1].x;
	Object1.WorldVertexList[4].y = L_ray2[2].Data[1].y;
	Object1.WorldVertexList[4].z = -0.02;

	Object1.WorldVertexList[5].x = L_ray2[2].Data[1].x;
	Object1.WorldVertexList[5].y = L_ray2[2].Data[1].y;
	Object1.WorldVertexList[5].z = 0.02;

	Object1.WorldVertexList[6].x = L_ray2[3].Data[1].x;
	Object1.WorldVertexList[6].y = L_ray2[3].Data[1].y;
	Object1.WorldVertexList[6].z = -0.02;

	Object1.WorldVertexList[7].x = L_ray2[3].Data[1].x;
	Object1.WorldVertexList[7].y = L_ray2[3].Data[1].y;
	Object1.WorldVertexList[7].z = 0.02;

	Object1.WorldVertexList[8].x = L_ray2[4].Data[1].x;
	Object1.WorldVertexList[8].y = L_ray2[4].Data[1].y;
	Object1.WorldVertexList[8].z = -0.02;

	Object1.WorldVertexList[9].x = L_ray2[4].Data[1].x;
	Object1.WorldVertexList[9].y = L_ray2[4].Data[1].y;
	Object1.WorldVertexList[9].z = 0.02;

	Object1.WorldVertexList[10].x = L_ray2[5].Data[1].x;
	Object1.WorldVertexList[10].y = L_ray2[5].Data[1].y;
	Object1.WorldVertexList[10].z = -0.02;

	Object1.WorldVertexList[11].x = L_ray2[5].Data[1].x;
	Object1.WorldVertexList[11].y = L_ray2[5].Data[1].y;
	Object1.WorldVertexList[11].z = 0.02;

	Object1.WorldVertexList[12].x = L_ray2[6].Data[1].x;
	Object1.WorldVertexList[12].y = L_ray2[6].Data[1].y;
	Object1.WorldVertexList[12].z = -0.02;

	Object1.WorldVertexList[13].x = L_ray2[6].Data[1].x;
	Object1.WorldVertexList[13].y = L_ray2[6].Data[1].y;
	Object1.WorldVertexList[13].z = 0.02;

	Object1.WorldVertexList[14].x = L_ray2[7].Data[1].x;
	Object1.WorldVertexList[14].y = L_ray2[7].Data[1].y;
	Object1.WorldVertexList[14].z = -0.02;

	Object1.WorldVertexList[15].x = L_ray2[7].Data[1].x;
	Object1.WorldVertexList[15].y = L_ray2[7].Data[1].y;
	Object1.WorldVertexList[15].z = 0.02;

	Object1.WorldVertexList[16].x = L_ray2[8].Data[1].x;
	Object1.WorldVertexList[16].y = L_ray2[8].Data[1].y;
	Object1.WorldVertexList[16].z = -0.02;

	Object1.WorldVertexList[17].x = L_ray2[8].Data[1].x;
	Object1.WorldVertexList[17].y = L_ray2[8].Data[1].y;
	Object1.WorldVertexList[17].z = 0.02;

	Object1.WorldVertexList[18].x = L_ray2[9].Data[1].x;
	Object1.WorldVertexList[18].y = L_ray2[9].Data[1].y;
	Object1.WorldVertexList[18].z = -0.02;

	Object1.WorldVertexList[19].x = L_ray2[9].Data[1].x;
	Object1.WorldVertexList[19].y = L_ray2[9].Data[1].y;
	Object1.WorldVertexList[19].z = 0.02;

	j = 0;
	for (i = 0; i < (3 * Object1.Triangle_N); i += 6)
	{
		Object1.TriangleList[i + 0] = 1 + j;
		Object1.TriangleList[i + 1] = 0 + j;
		Object1.TriangleList[i + 2] = 2 + j;

		Object1.TriangleList[i + 3] = 2 + j;
		Object1.TriangleList[i + 4] = 3 + j;
		Object1.TriangleList[i + 5] = 1 + j;
		j += 2;
	}

	for (i = 0; i < (3 * Object1.Triangle_N); i += 3)
	{
		Object1.TriangleColor[i] = L_ray2[0].color[0];
		Object1.TriangleColor[i + 1] = L_ray2[0].color[1];
		Object1.TriangleColor[i + 2] = L_ray2[0].color[2];
	}

	Object2.WorldVertexList[0].x = L_ray2[0 + 10].Data[1].x;
	Object2.WorldVertexList[0].y = L_ray2[0 + 10].Data[1].y;
	Object2.WorldVertexList[0].z = -0.02;

	Object2.WorldVertexList[1].x = L_ray2[0 + 10].Data[1].x;
	Object2.WorldVertexList[1].y = L_ray2[0 + 10].Data[1].y;
	Object2.WorldVertexList[1].z = 0.02;

	Object2.WorldVertexList[2].x = L_ray2[1 + 10].Data[1].x;
	Object2.WorldVertexList[2].y = L_ray2[1 + 10].Data[1].y;
	Object2.WorldVertexList[2].z = -0.02;

	Object2.WorldVertexList[3].x = L_ray2[1 + 10].Data[1].x;
	Object2.WorldVertexList[3].y = L_ray2[1 + 10].Data[1].y;
	Object2.WorldVertexList[3].z = 0.02;

	Object2.WorldVertexList[4].x = L_ray2[2 + 10].Data[1].x;
	Object2.WorldVertexList[4].y = L_ray2[2 + 10].Data[1].y;
	Object2.WorldVertexList[4].z = -0.02;

	Object2.WorldVertexList[5].x = L_ray2[2 + 10].Data[1].x;
	Object2.WorldVertexList[5].y = L_ray2[2 + 10].Data[1].y;
	Object2.WorldVertexList[5].z = 0.02;

	Object2.WorldVertexList[6].x = L_ray2[3 + 10].Data[1].x;
	Object2.WorldVertexList[6].y = L_ray2[3 + 10].Data[1].y;
	Object2.WorldVertexList[6].z = -0.02;

	Object2.WorldVertexList[7].x = L_ray2[3 + 10].Data[1].x;
	Object2.WorldVertexList[7].y = L_ray2[3 + 10].Data[1].y;
	Object2.WorldVertexList[7].z = 0.02;

	Object2.WorldVertexList[8].x = L_ray2[4 + 10].Data[1].x;
	Object2.WorldVertexList[8].y = L_ray2[4 + 10].Data[1].y;
	Object2.WorldVertexList[8].z = -0.02;

	Object2.WorldVertexList[9].x = L_ray2[4 + 10].Data[1].x;
	Object2.WorldVertexList[9].y = L_ray2[4 + 10].Data[1].y;
	Object2.WorldVertexList[9].z = 0.02;

	Object2.WorldVertexList[10].x = L_ray2[5 + 10].Data[1].x;
	Object2.WorldVertexList[10].y = L_ray2[5 + 10].Data[1].y;
	Object2.WorldVertexList[10].z = -0.02;

	Object2.WorldVertexList[11].x = L_ray2[5 + 10].Data[1].x;
	Object2.WorldVertexList[11].y = L_ray2[5 + 10].Data[1].y;
	Object2.WorldVertexList[11].z = 0.02;

	Object2.WorldVertexList[12].x = L_ray2[6 + 10].Data[1].x;
	Object2.WorldVertexList[12].y = L_ray2[6 + 10].Data[1].y;
	Object2.WorldVertexList[12].z = -0.02;

	Object2.WorldVertexList[13].x = L_ray2[6 + 10].Data[1].x;
	Object2.WorldVertexList[13].y = L_ray2[6 + 10].Data[1].y;
	Object2.WorldVertexList[13].z = 0.02;

	Object2.WorldVertexList[14].x = L_ray2[7 + 10].Data[1].x;
	Object2.WorldVertexList[14].y = L_ray2[7 + 10].Data[1].y;
	Object2.WorldVertexList[14].z = -0.02;

	Object2.WorldVertexList[15].x = L_ray2[7 + 10].Data[1].x;
	Object2.WorldVertexList[15].y = L_ray2[7 + 10].Data[1].y;
	Object2.WorldVertexList[15].z = 0.02;

	Object2.WorldVertexList[16].x = L_ray2[8 + 10].Data[1].x;
	Object2.WorldVertexList[16].y = L_ray2[8 + 10].Data[1].y;
	Object2.WorldVertexList[16].z = -0.02;

	Object2.WorldVertexList[17].x = L_ray2[8 + 10].Data[1].x;
	Object2.WorldVertexList[17].y = L_ray2[8 + 10].Data[1].y;
	Object2.WorldVertexList[17].z = 0.02;

	Object2.WorldVertexList[18].x = L_ray2[9 + 10].Data[1].x;
	Object2.WorldVertexList[18].y = L_ray2[9 + 10].Data[1].y;
	Object2.WorldVertexList[18].z = -0.02;

	Object2.WorldVertexList[19].x = L_ray2[9 + 10].Data[1].x;
	Object2.WorldVertexList[19].y = L_ray2[9 + 10].Data[1].y;
	Object2.WorldVertexList[19].z = 0.02;

	j = 0;
	for (i = 0; i < (3 * Object1.Triangle_N); i += 6)
	{
		Object2.TriangleList[i + 0] = 2 + j;
		Object2.TriangleList[i + 1] = 0 + j;
		Object2.TriangleList[i + 2] = 1 + j;

		Object2.TriangleList[i + 3] = 1 + j;
		Object2.TriangleList[i + 4] = 3 + j;
		Object2.TriangleList[i + 5] = 2 + j;
		j += 2;
	}



	for (i = 0; i < (3 * Object2.Triangle_N); i += 3)
	{
		Object2.TriangleColor[i] = L_ray2[10].color[0];
		Object2.TriangleColor[i + 1] = L_ray2[10].color[1];
		Object2.TriangleColor[i + 2] = L_ray2[10].color[2];
	}


	SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
	SDL_RenderClear(renderer);
	SDL_RenderPresent(renderer);


	k = 0;
	SDL_SetWindowGrab(window, SDL_TRUE);

	SDL_SetRelativeMouseMode(SDL_TRUE);

	while (k == 0)
	{
		while (SDL_PollEvent(&event))
		{

			if (event.type == SDL_MOUSEMOTION)
			{
				MouseMH = event.motion.xrel;
				MouseMV = event.motion.yrel;

				RotCamAngAzimutal = (MouseMH * pi / 180) / 1000;
				RotCamAngPolar = (MouseMV * pi / 180) / 1000;
				//printf("\n Azimutal = %lf - Polar = %lf", RotCamAngAzimutal, RotCamAngPolar);

				M = XYZ_SetRotationMatrix(InertialCoordinateSystem.base[1], RotCamAngAzimutal);

				Cam1.camCS.base[0] = XYZ_ApplyTransformationMatrix(Cam1.camCS.base[0], M);
				Cam1.camCS.base[1] = XYZ_ApplyTransformationMatrix(Cam1.camCS.base[1], M);
				Cam1.camCS.base[2] = XYZ_ApplyTransformationMatrix(Cam1.camCS.base[2], M);
				WalkAxis = XYZ_ApplyTransformationMatrix(WalkAxis, M);

				M = XYZ_SetRotationMatrix(Cam1.camCS.base[0], RotCamAngPolar);

				Cam1.camCS.base[1] = XYZ_ApplyTransformationMatrix(Cam1.camCS.base[1], M);
				Cam1.camCS.base[2] = XYZ_ApplyTransformationMatrix(Cam1.camCS.base[2], M);

			}


			if (event.type == SDL_KEYDOWN)
			{
				if (event.key.keysym.sym == SDLK_ESCAPE) k = 1;


				if (event.key.keysym.sym == SDLK_r)
				{
					if (RotControl == 0) RotControl = 1;
					else RotControl = 0;
				}



				if (event.key.keysym.sym == SDLK_w)
				{
					Cam1.ProjectionPlane_Z = 0.0005;
					Cam1.NearPlane_Z = 0.002;
					Cam1.FarPlane_Z = 50;
					Cam1.ProjectioPlaneWidth *= 0.95;
					Cam1.ProjectionPlaneHeight *= 0.95;

					Cam1.WidthDepthRatio = (Cam1.ProjectioPlaneWidth / 2) / Cam1.ProjectionPlane_Z;
					Cam1.HeightDepthRatio = (Cam1.ProjectionPlaneHeight / 2) / Cam1.ProjectionPlane_Z;

				}

				if (event.key.keysym.sym == SDLK_s)
				{
					Cam1.ProjectionPlane_Z = 0.0005;
					Cam1.NearPlane_Z = 0.002;
					Cam1.FarPlane_Z = 50;
					Cam1.ProjectioPlaneWidth /= 0.95;
					Cam1.ProjectionPlaneHeight /= 0.95;

					Cam1.WidthDepthRatio = (Cam1.ProjectioPlaneWidth / 2) / Cam1.ProjectionPlane_Z;
					Cam1.HeightDepthRatio = (Cam1.ProjectionPlaneHeight / 2) / Cam1.ProjectionPlane_Z;

				}

				if (event.key.keysym.sym == SDLK_d)
				{
					Cam1.camCS.Origin = XYZ_PointTranslation(Cam1.camCS.Origin, InertialCoordinateSystem.base[0], 0.1);

				}

				if (event.key.keysym.sym == SDLK_a)
				{
					Cam1.camCS.Origin = XYZ_PointTranslation(Cam1.camCS.Origin, InertialCoordinateSystem.base[0], -0.1);

				}



			}
		}


		SDL_SetRenderDrawColor(renderer, 50, 50, 200, 255);
		XYZ_SetOBJCamCoordinates(&ApparentHorizonOBJ, &Cam1);
		XYZ_SetOBJProjectionSpaceCoordinate(&ApparentHorizonOBJ, &Cam1);
		XYZ_ViewportVisibleGroundRender(&ApparentHorizonOBJ, &Cam1);


		SDL_SetRenderDrawColor(renderer, 50, 200, 50, 255);
		XYZ_SetOBJCamCoordinates(&GeometricHorizonOBJ, &Cam1);
		XYZ_SetOBJProjectionSpaceCoordinate(&GeometricHorizonOBJ, &Cam1);
		XYZ_ViewportGeometricHorizon(&GeometricHorizonOBJ, &Cam1);




		SDL_SetRenderDrawColor(renderer, 40, 230, 240, 255);
		XYZ_SetOBJCamCoordinates(&Object1, &Cam1);
		XYZ_ApplyTriangleCulling(&Object1, &Cam1);
		XYZ_SetOBJProjectionSpaceCoordinate(&Object1, &Cam1);
		XYZ_ViewportOBJRender2(&Object1, &Cam1);

		SDL_SetRenderDrawColor(renderer, 40, 230, 240, 255);
		XYZ_SetOBJCamCoordinates(&Object2, &Cam1);
		XYZ_ApplyTriangleCulling(&Object2, &Cam1);
		XYZ_SetOBJProjectionSpaceCoordinate(&Object2, &Cam1);
		XYZ_ViewportOBJRender2(&Object2, &Cam1);

		/*SDL_SetRenderDrawColor(renderer, 40, 230, 240, 255);
		XYZ_SetOBJCamCoordinates(&RealObject, &Cam1);
		XYZ_ApplyTriangleCulling(&RealObject, &Cam1);
		XYZ_SetOBJProjectionSpaceCoordinate(&RealObject, &Cam1);
		XYZ_ViewportOBJRender2(&RealObject, &Cam1);*/


		//delay(40);
		//SDL_Delay(10);
		SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);
		//drawLine(renderer, 1248, 700, 1252, 700);
		//drawLine(renderer, 1250, 698, 1250, 702);
		SDL_RenderPresent(renderer);
		SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
		SDL_RenderClear(renderer);
		//GroundShow(Cam1, 0.5, 0, 0, 0, 0);
		SDL_SetRenderDrawColor(renderer, 255, 255, 255, 255);

	}


	SDL_SetWindowGrab(window, SDL_FALSE);

	XYZ_DeallocatedObjectLists(&RealObject);
	XYZ_DeallocatedObjectLists(&Object1);
	XYZ_DeallocatedObjectLists(&Object2);




	free(ApparentHorizonOBJ.WorldVertexList);
	free(ApparentHorizonOBJ.CamCoordSysVertexList);
	free(ApparentHorizonOBJ.ProjectionSpaceVertexList);
	free(ApparentHorizonOBJ.VisibleVertices);


	free(GeometricHorizonOBJ.WorldVertexList);
	free(GeometricHorizonOBJ.CamCoordSysVertexList);
	free(GeometricHorizonOBJ.ProjectionSpaceVertexList);
	free(GeometricHorizonOBJ.VisibleVertices);




















	return 0;

}

XYZPoint XYZ_PointRotation(XYZPoint P, XYZPoint O, double ang)
{
	XYZPoint Paux;

	Paux = P;

	Paux.x = P.x * cos(ang) - P.z * sin(ang);
	Paux.z = P.x * sin(ang) + P.z * cos(ang);

	return Paux;
}

XYZPoint XYZ_PointTranslation(XYZPoint P, XYZPoint Vector, double t)
{
	XYZPoint Paux;

	Paux = P;

	Paux.x = P.x + t * Vector.x;
	Paux.y = P.y + t * Vector.y;
	Paux.z = P.z + t * Vector.z;

	return Paux;
}


XYZTransformatioMatrix XYZ_SetRotationMatrix(XYZPoint RotationAxis, double ang)
{
	XYZTransformatioMatrix RotMatT;
	double I[3][3], J[3][3], JJ[3][3];
	int i, j;

	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			if (i == j)
			{
				I[i][j] = 1;
				J[i][j] = 0;
			}
			else
			{
				I[i][j] = 0;
			}
		}
	}

	J[0][1] = -RotationAxis.z;
	J[0][2] = RotationAxis.y;
	J[1][0] = RotationAxis.z;
	J[1][2] = -RotationAxis.x;
	J[2][0] = -RotationAxis.y;
	J[2][1] = RotationAxis.x;

	for (i = 0; i < 3; i++)
	{
		for (j = 0; j < 3; j++)
		{
			JJ[i][j] = J[i][0] * J[0][j] + J[i][1] * J[1][j] + J[i][2] * J[2][j];

			RotMatT.Data[i][j] = I[i][j] + ((sin(ang)) * J[i][j]) + ((1 - cos(ang)) * JJ[i][j]);
		}
	}

	return RotMatT;
}

XYZPoint XYZ_ApplyTransformationMatrix(XYZPoint P, XYZTransformatioMatrix M)
{
	XYZPoint Q;

	Q.x = M.Data[0][0] * P.x + M.Data[0][1] * P.y + M.Data[0][2] * P.z;
	Q.y = M.Data[1][0] * P.x + M.Data[1][1] * P.y + M.Data[1][2] * P.z;
	Q.z = M.Data[2][0] * P.x + M.Data[2][1] * P.y + M.Data[2][2] * P.z;

	return Q;
}

void XYZ_AllocateObjectLists(XYZObject* OBJ, int N_Vertices, int N_Triangles)
{
	(*OBJ).WorldVertexList = (XYZPoint*)malloc(N_Vertices * sizeof(XYZPoint));
	(*OBJ).CamCoordSysVertexList = (XYZPoint*)malloc(N_Vertices * sizeof(XYZPoint));
	(*OBJ).ProjectionSpaceVertexList = (XYZPoint*)malloc(N_Vertices * sizeof(XYZPoint));
	(*OBJ).VisibleVertices = (char*)malloc(N_Vertices * sizeof(char));
	(*OBJ).TriangleList = (int*)malloc((3 * N_Triangles) * sizeof(int));
	(*OBJ).PSTriangleList = (int*)malloc((3 * N_Triangles) * sizeof(int));
	(*OBJ).TriangleColor = (int*)malloc((3 * N_Triangles) * sizeof(int));
	(*OBJ).PSTriangleColor = (int*)malloc((3 * N_Triangles) * sizeof(int));

	if ((*OBJ).WorldVertexList == NULL)
	{
		printf("Vertex list allocation failed");
		exit(1);
	}

	if ((*OBJ).CamCoordSysVertexList == NULL)
	{
		printf("Vertex list allocation failed");
		exit(1);
	}

	if ((*OBJ).ProjectionSpaceVertexList == NULL)
	{
		printf("Vertex list allocation failed");
		exit(1);
	}

	if ((*OBJ).VisibleVertices == NULL)
	{
		printf("Visible Vertex list allocation failed");
		exit(1);
	}

	if ((*OBJ).TriangleList == NULL)
	{
		printf("Triangle list allocation failed");
		exit(1);
	}

	if ((*OBJ).PSTriangleList == NULL)
	{
		printf("Triangle list allocation failed");
		exit(1);
	}

	if ((*OBJ).TriangleColor == NULL)
	{
		printf("Triangle Color list allocation failed");
		exit(1);
	}

	if ((*OBJ).PSTriangleColor == NULL)
	{
		printf("Triangle Color list allocation failed");
		exit(1);
	}
}

void XYZ_DeallocatedObjectLists(XYZObject* OBJ)
{
	free((*OBJ).WorldVertexList);
	free((*OBJ).CamCoordSysVertexList);
	free((*OBJ).ProjectionSpaceVertexList);
	free((*OBJ).VisibleVertices);
	free((*OBJ).TriangleList);
	free((*OBJ).PSTriangleList);
	free((*OBJ).TriangleColor);
	free((*OBJ).PSTriangleColor);


	(*OBJ).WorldVertexList = NULL;
	(*OBJ).TriangleList = NULL;
}


/************************   Area atualmente em desenvolvimento  *******************************************/

void XYZ_SetOBJCamCoordinates(XYZObject* OBJ, XYZCam* Cam)
{
	int i;
	XYZPoint Paux;

	for (i = 0; i < (*OBJ).Vertex_N; i++)
	{
		Paux = (*OBJ).WorldVertexList[i];

		Paux.x = Paux.x - (*Cam).camCS.Origin.x;
		Paux.y = Paux.y - (*Cam).camCS.Origin.y;
		Paux.z = Paux.z - (*Cam).camCS.Origin.z;

		(*OBJ).CamCoordSysVertexList[i].x = (*Cam).camCS.base[0].x * Paux.x + (*Cam).camCS.base[0].y * Paux.y + (*Cam).camCS.base[0].z * Paux.z;
		(*OBJ).CamCoordSysVertexList[i].y = (*Cam).camCS.base[1].x * Paux.x + (*Cam).camCS.base[1].y * Paux.y + (*Cam).camCS.base[1].z * Paux.z;
		(*OBJ).CamCoordSysVertexList[i].z = (*Cam).camCS.base[2].x * Paux.x + (*Cam).camCS.base[2].y * Paux.y + (*Cam).camCS.base[2].z * Paux.z;
	}
}

void XYZ_ApplyTriangleCulling(XYZObject* OBJ, XYZCam* Cam)
{
	XYZPoint PA, PB, PC, Paux, Normal, Vector1, Vector2;
	int i, j = 0, TriangleListLength;

	TriangleListLength = (3 * (*OBJ).Triangle_N);

	for (i = 0; i < (*OBJ).Vertex_N; i++)
	{
		(*OBJ).VisibleVertices[i] = 'F';
	}

	(*OBJ).PSVertex_N = 0; // zera o numero de vertices no espaco de projecao
	(*OBJ).PSTriangle_N = 0; // zera o numero de triangulos no espaco de projecao

	for (i = 0; i < TriangleListLength; i += 3)
	{
		PA = (*OBJ).CamCoordSysVertexList[(*OBJ).TriangleList[i]];
		PB = (*OBJ).CamCoordSysVertexList[(*OBJ).TriangleList[i + 1]];
		PC = (*OBJ).CamCoordSysVertexList[(*OBJ).TriangleList[i + 2]];


		/************************************** Culling de Frustum(deve ser paralelizado) **************************************/
		//Verifica se esta entre os planos proximo ou distante
		if ((PA.z > (*Cam).FarPlane_Z && PB.z > (*Cam).FarPlane_Z && PC.z > (*Cam).FarPlane_Z) || (PA.z < (*Cam).NearPlane_Z && PB.z < (*Cam).NearPlane_Z && PC.z < (*Cam).NearPlane_Z))
		{
			continue;
		}


		// WidthDepthRatio define o coeficiente angular da reta definida pela interseccao do plano direito do frustum com o plano xz.
		/*
		==============================================================================================================================

				 ^ (eixo z)
				 |
				 |
				 |
		\        |        /  <-- Interseccao dos Planos laterais do frustum com o plano xz (vista superior).
		 \       |       /       A variavel WidthDepthRatio armazena o coeficiente angular dessa reta.
		  \      |      /        Com isso, fica facil determinar se um ponto esta entre os dois planos (dentro do V)
		   \     |     /         ou se está fora
			\    |    /
			 \   |   /
			  \  |  /
			   \ | /
				\|/
		---------*-----------------------------------> (eixo x)
				 |
				 |

		==============================================================================================================================
		*/

		// Como, nesse caso, a coordenada y e irrelevante, o calculo para determinar de qual lado do plano o ponto esta pode ser reduzido ao plano xz
		// Sera analogo para os outros planos 

		// Verificacao para o plano esquerdo
		if (PA.x < (-(*Cam).WidthDepthRatio * PA.z) && PB.x < (-(*Cam).WidthDepthRatio * PB.z) && PC.x < (-(*Cam).WidthDepthRatio * PC.z)) continue;

		// Verificacao para o plano Direito
		if (PA.x > ((*Cam).WidthDepthRatio * PA.z) && PB.x > ((*Cam).WidthDepthRatio * PB.z) && PC.x > ((*Cam).WidthDepthRatio * PC.z)) continue;

		// Verificacao para o plano superior
		if (PA.y > ((*Cam).HeightDepthRatio * PA.z) && PB.y > ((*Cam).HeightDepthRatio * PB.z) && PC.y > ((*Cam).HeightDepthRatio * PC.z)) continue;

		// Verificacao para o plano inferior
		if (PA.y < (-(*Cam).HeightDepthRatio * PA.z) && PB.y < (-(*Cam).HeightDepthRatio * PB.z) && PC.y < (-(*Cam).HeightDepthRatio * PC.z)) continue;
		/************************************** Fim do Culling de Frustum **************************************/

		/************************************** Back-Face Culling (deve ser paralelizado) **************************************/
		//Vetores definidos pelas arestas do triangulo. O ponto a e a origem dos dois vetores.
		Vector1.x = PB.x - PA.x;
		Vector1.y = PB.y - PA.y;
		Vector1.z = PB.z - PA.z;

		Vector2.x = PC.x - PA.x;
		Vector2.y = PC.y - PA.y;
		Vector2.z = PC.z - PA.z;

		//Calculo do vertor normal a partir do produto vetorial.
		Normal.x = Vector1.y * Vector2.z - Vector1.z * Vector2.y;
		Normal.y = Vector1.z * Vector2.x - Vector1.x * Vector2.z;
		Normal.z = Vector1.x * Vector2.y - Vector1.y * Vector2.x;

		//O produto escalar entre o vetor definido pela origem e pelo ponto PA (vetor PA) e o vetor normal deve ser maior que zero para a face estar voltada para a camera.
		if (PA.x * Normal.x + PA.y * Normal.y + PA.z * Normal.z <= 0) continue;


		/************************************** Fim do Back-Face Culling **************************************/


		// Se o triangulo passou por todos os testes, entao o triangulo esta dentro do frustum e voltado para tela.
		//Sendo assim, ele deve ser colocado na lista de triangulos do espaco de projecao e seus vertices deverao ser transformados

		//Aqui a lista de triangulos que estao dentro do frustum e gerada. Esses serao os triangulos do espaco de projecao.
		// o indice i percorre a lista de triangulos original e o indice j percorre a lista de triangulos do espaco de projecao
		(*OBJ).PSTriangleList[j] = (*OBJ).TriangleList[i]; //O i-esimo triangulo da lista original devera entrar como o j-esimo triangulo da lista de triangulos do espaco de projecao
		(*OBJ).PSTriangleList[j + 1] = (*OBJ).TriangleList[i + 1];
		(*OBJ).PSTriangleList[j + 2] = (*OBJ).TriangleList[i + 2];

		(*OBJ).PSTriangleColor[j] = (*OBJ).TriangleColor[i]; //O i-esimo triangulo da lista original devera entrar como o j-esimo triangulo da lista de triangulos do espaco de projecao
		(*OBJ).PSTriangleColor[j + 1] = (*OBJ).TriangleColor[i + 1];
		(*OBJ).PSTriangleColor[j + 2] = (*OBJ).TriangleColor[i + 2];

		j += 3;
		//printf("\nO triangulo %d apareceu", i/3);

		/*
		==================================== Descricao do controle de vertices visiveis ====================================
		Os indices dos vertices visiveis serao salvos no array unidimensional (*OBJ).VisibleVertices. A cada tres posicoes
		consecutivas tem-se um triangulo.
		O array com os vertices visiveis mapeia as posicoes dos vertices originais para transforma-los nas coordenadas da tela.
		Teremos entao:
		Lista de vertices no espaco da camera: (*OBJ).CamCoordSysVertexList
		Lista de triangulo no espaco da camera: (*OBJ).TriangleList
		Lista de Vertices no espaco de projecao: (*OBJ).ProjectionSpaceVertexList
		Lista de triangulos no espaco de projecao: (*OBJ).PSTriangleList
		Lista de vertices visiveis: (*OBJ).VisibleVertices

		Se o controle de vertices na posicao i e true, entao o vertice da posicao i e visivel, caso contrario, ele e invisivel.

		Segue um exemplo de execucao:

		Array de vertices no espaco da camera (CamCoordSysVertexList)
		| Vertex1 | Vertex2 | Vertex3 | Vertex4 | Vertex5 | Vertex6 | Vertex7 | Vertex8 | Vertex9 | Vertex10 | Vertex11 | Vertex12 | Vertex13 | Vertex14 | Vertex15 |
			 0         1         2         3         4         5         6         7         8         9          10         11          12         13        14

		Array de triangulos no espaco da camera ((*OBJ).TriangleList) (eh a mesma lista para o espaco do mundo)
		| 4 | 1 | 5 | 10 | 12 | 7 |  4 | 3 | 5 | 0 | 2 | 7 |  1 | 3 | 5 | 10 | 12 | 9 | 3 | 2 | 5 | 10 | 14 | 6 |  8 | 9 | 11 | 13 | 7 | 14 |  10 | 11 | 12 | 1 | 5 | 7 |


		Array de controle de vertices visiveis ((*OBJ).VisibleVertices)
		Esse e o estado inicial, indica que nenhum vertice esta visivel.
		| false | false | false | false | false | false | false | false | false | false | false | false | false | false | false |
			 0      1       2       3       4       5       6       7       8       9       10      11      12      13      14

		Sao 12 triangulos. Suponha que sejam visieis apenas os seguintes triangulos
		(*OBJ).TriangleList[3]
		(*OBJ).TriangleList[4]
		(*OBJ).TriangleList[5]

		(*OBJ).TriangleList[12]
		(*OBJ).TriangleList[13]
		(*OBJ).TriangleList[14]

		(*OBJ).TriangleList[30]
		(*OBJ).TriangleList[31]
		(*OBJ).TriangleList[32]

		Sendo assim, os vertices visiveis serao:
		(*OBJ).VisibleVertices[10] = true;
		(*OBJ).VisibleVertices[12] = true;
		(*OBJ).VisibleVertices[7] = true;

		(*OBJ).VisibleVertices[1] = true;
		(*OBJ).VisibleVertices[3] = true;
		(*OBJ).VisibleVertices[5] = true;

		(*OBJ).VisibleVertices[10] = true;
		(*OBJ).VisibleVertices[11] = true;
		(*OBJ).VisibleVertices[12] = true;

		O array de controle de vertices visiveis ((*OBJ).VisibleVertices) vai ter a seguinte configuracao:
		| false | true | false | true | false | true | false | true | false | false | true | true | true | false | false |
			 0      1      2       3      4       5      6       7      8       9      10      11    12      13      14

		==============================================================================================================================

		*/

		// triangulo visivel! O contador deve ser incrementado.
		(*OBJ).PSTriangle_N += 1;

		//Marca o Vertice da posicao (*OBJ).TriangleList[i] (e as duas seguintes) como visivel, pois o trianglo foi considerado visivel.
		(*OBJ).VisibleVertices[(*OBJ).TriangleList[i]] = 'T';
		(*OBJ).VisibleVertices[(*OBJ).TriangleList[i + 1]] = 'T';
		(*OBJ).VisibleVertices[(*OBJ).TriangleList[i + 2]] = 'T';

	}
}

void XYZ_SetOBJProjectionSpaceCoordinate(XYZObject* OBJ, XYZCam* Cam)
{
	int i;


	for (i = 0; i < (*OBJ).Vertex_N; i++)
	{
		if ((*OBJ).VisibleVertices[i] == 'T')
		{
			//Transformacao do ponto do espaco da camera para o espaco de projecao
			(*OBJ).ProjectionSpaceVertexList[i].x = ((*Cam).ProjectionPlane_Z / ((*OBJ).CamCoordSysVertexList[i].z)) * (*OBJ).CamCoordSysVertexList[i].x;
			(*OBJ).ProjectionSpaceVertexList[i].y = ((*Cam).ProjectionPlane_Z / ((*OBJ).CamCoordSysVertexList[i].z)) * (*OBJ).CamCoordSysVertexList[i].y;
			(*OBJ).ProjectionSpaceVertexList[i].z = (*Cam).ProjectionPlane_Z;
		}

	}

}

void XYZ_ViewportOBJRender(XYZObject* OBJ, XYZCam* Cam)
{
	XYZPoint PA, PB, PC;
	int i;
	int TriangleListLength;
	double Pixel_Width_Ratio, Pixel_Height_Ratio;

	Pixel_Width_Ratio = 2500 / (*Cam).ProjectioPlaneWidth;
	Pixel_Height_Ratio = 1400 / (*Cam).ProjectionPlaneHeight;

	// variavel criada para otimizacao
	TriangleListLength = 3 * (*OBJ).PSTriangle_N;

	for (i = 0; i < TriangleListLength; i += 3)
	{
		PA = (*OBJ).ProjectionSpaceVertexList[(*OBJ).PSTriangleList[i]];
		PB = (*OBJ).ProjectionSpaceVertexList[(*OBJ).PSTriangleList[i + 1]];
		PC = (*OBJ).ProjectionSpaceVertexList[(*OBJ).PSTriangleList[i + 2]];

		PA.x = Pixel_Width_Ratio * PA.x + 1250;
		PB.x = Pixel_Width_Ratio * PB.x + 1250;
		PC.x = Pixel_Width_Ratio * PC.x + 1250;

		PA.y = (-Pixel_Height_Ratio) * PA.y + 700;
		PB.y = (-Pixel_Height_Ratio) * PB.y + 700;
		PC.y = (-Pixel_Height_Ratio) * PC.y + 700;

		XYZ_drawLine(renderer, PA.x, PA.y, PB.x, PB.y);
		//drawLine(renderer, PA.x, PA.y, PC.x, PC.y);
		XYZ_drawLine(renderer, PB.x, PB.y, PC.x, PC.y);


	}
}


void XYZ_ViewportOBJRender2(XYZObject* OBJ, XYZCam* Cam)
{
	XYZPoint PA, PB, PC;
	int i, x, y;
	int TriangleListLength;
	int PixelA_x, PixelA_y, PixelB_x, PixelB_y, PixelC_x, PixelC_y;
	double Pixel_Width_Ratio, Pixel_Height_Ratio, aux, m_AB = 0, m_AC = 0, m_BC = 0, ymin, yMAX;

	Pixel_Width_Ratio = 2500 / (*Cam).ProjectioPlaneWidth;
	Pixel_Height_Ratio = 1400 / (*Cam).ProjectionPlaneHeight;

	// variavel criada para otimizacao
	TriangleListLength = 3 * (*OBJ).PSTriangle_N;

	for (i = 0; i < TriangleListLength; i += 3)
	{
		SDL_SetRenderDrawColor(renderer, (*OBJ).PSTriangleColor[i], (*OBJ).PSTriangleColor[i + 1], (*OBJ).PSTriangleColor[i + 2], 255);

		PA = (*OBJ).ProjectionSpaceVertexList[(*OBJ).PSTriangleList[i]];
		PB = (*OBJ).ProjectionSpaceVertexList[(*OBJ).PSTriangleList[i + 1]];
		PC = (*OBJ).ProjectionSpaceVertexList[(*OBJ).PSTriangleList[i + 2]];

		PA.x = Pixel_Width_Ratio * PA.x + 1250;
		PB.x = Pixel_Width_Ratio * PB.x + 1250;
		PC.x = Pixel_Width_Ratio * PC.x + 1250;

		PA.y = (-Pixel_Height_Ratio) * PA.y + 700;
		PB.y = (-Pixel_Height_Ratio) * PB.y + 700;
		PC.y = (-Pixel_Height_Ratio) * PC.y + 700;


		if (PA.x > PB.x)
		{
			aux = PA.x;
			PA.x = PB.x;
			PB.x = aux;

			aux = PA.y;
			PA.y = PB.y;
			PB.y = aux;
		}

		if (PA.x > PC.x)
		{
			aux = PA.x;
			PA.x = PC.x;
			PC.x = aux;

			aux = PA.y;
			PA.y = PC.y;
			PC.y = aux;
		}

		if (PC.x < PB.x)
		{
			aux = PC.x;
			PC.x = PB.x;
			PB.x = aux;

			aux = PC.y;
			PC.y = PB.y;
			PB.y = aux;
		}

		PixelA_x = (int)round(PA.x);
		PixelA_y = (int)round(PA.y);
		PixelB_x = (int)round(PB.x);
		PixelB_y = (int)round(PB.y);
		PixelC_x = (int)round(PC.x);
		PixelC_y = (int)round(PC.y);

		if (PixelA_x != PixelC_x)
		{
			m_AC = (double)(PixelC_y - PixelA_y) / (double)(PixelC_x - PixelA_x);
		}


		if (PixelB_x != PixelA_x)
		{
			m_AB = (double)(PixelB_y - PixelA_y) / (double)(PixelB_x - PixelA_x);
		}

		if (PixelB_x != PixelC_x)
		{
			m_BC = (double)(PixelC_y - PixelB_y) / (double)(PixelC_x - PixelB_x);
		}


		for (x = PixelA_x; x <= PixelB_x; x++)
		{
			ymin = m_AB * (x - PB.x) + PB.y;
			yMAX = m_AC * (x - PA.x) + PA.y;

			if (yMAX < ymin)
			{
				aux = yMAX;
				yMAX = ymin;
				ymin = aux;


			}

			for (y = (int)ymin; y <= (int)yMAX; y++)
			{
				SDL_RenderDrawPoint(renderer, x, y);

			}
			y = y;
		}

		while (x <= PixelC_x)
		{
			ymin = m_BC * (x - PC.x) + PC.y;
			yMAX = m_AC * (x - PA.x) + PA.y;

			if (yMAX < ymin)
			{
				aux = yMAX;
				yMAX = ymin;
				ymin = aux;


			}

			for (y = (int)ymin; y <= (int)yMAX; y++)
			{
				SDL_RenderDrawPoint(renderer, x, y);

			}

			x++;
		}
	}
}


void XYZ_ViewportVisibleGroundRender(XYZObject* OBJ, XYZCam* Cam)
{
	XYZPoint PA;
	int i;
	int TriangleListLength;
	double Pixel_Width_Ratio, Pixel_Height_Ratio;

	Pixel_Width_Ratio = 2500 / (*Cam).ProjectioPlaneWidth;
	Pixel_Height_Ratio = 1400 / (*Cam).ProjectionPlaneHeight;

	PA = (*OBJ).ProjectionSpaceVertexList[0];


	PA.x = Pixel_Width_Ratio * PA.x + 1250;

	PA.y = (-Pixel_Height_Ratio) * PA.y + 700;


	/*XYZ_drawLine(renderer, 0, PA.y, 2500, PA.y);
	XYZ_drawLine(renderer, 0, PA.y + 1, 2500, PA.y + 1);*/

	for (i = 0; i < 2500; i++)
	{
		XYZ_drawLine(renderer, i, PA.y, i, 1400);
	}
}

void XYZ_ViewportGeometricHorizon(XYZObject* OBJ, XYZCam* Cam)
{
	XYZPoint PA;
	int i;
	int TriangleListLength;
	double Pixel_Width_Ratio, Pixel_Height_Ratio;

	Pixel_Width_Ratio = 2500 / (*Cam).ProjectioPlaneWidth;
	Pixel_Height_Ratio = 1400 / (*Cam).ProjectionPlaneHeight;

	PA = (*OBJ).ProjectionSpaceVertexList[0];


	PA.x = Pixel_Width_Ratio * PA.x + 1250;

	PA.y = (-Pixel_Height_Ratio) * PA.y + 700;


	XYZ_drawLine(renderer, 0, PA.y, 2500, PA.y);
	XYZ_drawLine(renderer, 0, PA.y + 1, 2500, PA.y + 1);

}



/************************ Fim da area atualmente em desenvolvimento  ************************************/




void XYZ_drawLine(SDL_Renderer* renderer, int xi, int yi, int xf, int yf)
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

		if (x1 <= 2500 && x1 >= 0 && y1 <= 1400 && y1 >= 0)
		{
			SDL_RenderDrawPoint(renderer, x1, y1);
			//SDL_RenderDrawPoint(renderer, x1, y1 + 1);
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

