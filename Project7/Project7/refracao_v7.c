#include <stdio.h>
#include <stdlib.h>
#include <SDL.h> // Para a biblioteca SDL
#include <conio.h>
#include <math.h>

#define pi 3.14159265358979323
#define EXP 2.71828182845904523536

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

int teste_2(reta_xy* L_ray, Terra* T, ponto_xy P, double ang, double xM); //funcao que calcula os raios refratados na atmosfera
double angulo_incidencia(ponto_xy v, ponto_xy P, ponto_xy C); //funcao que calcula o angulo de incidencia em relacao a normal
ponto_xy rotacao_vetorial(ponto_xy v, double ang); // funcao que rotaciona o vetor diretor
char posicao_vetor_raio(ponto_xy v, ponto_xy P, ponto_xy C); // funcao que calcula a posicao relativa do raio incidente em relacao a normal
void viewport_LightRay_render(reta_xy* L_ray, double xpp, double ypp, double zpp, int nray, int luz_cor[3]); // funcao que desenha os raios nas coordenadas da tela
void viewport_Earth_render(Terra* T, double xpp, double ypp, double zpp); // funcao que desenha a atmosfera nas coordenadas da tela
void drawLine(SDL_Renderer* renderer, int xi, int yi, int xf, int yf);
void drawCircle(SDL_Renderer* renderer, int centerX, int centerY, int radius);
void drawCircle2(SDL_Renderer* renderer, int centerX, int centerY, int Raio);
void drawCircle3(SDL_Renderer* renderer, int centerX, int centerY, int Raio);
ponto_xy Earth_tangent_point(ponto_xy PO, Terra* T);
reta_xy Earth_tangent_line(ponto_xy PO, Terra* T);
void DataShow(ponto_xy Obs, ponto_xy TangentP, ponto_xy TUp, ponto_xy TDown, ponto_xy RefP1, ponto_xy RefP2);


/* Variaveis globais para dados auxiliares. */
ponto_xy Paux, Paux2;
SDL_Renderer* renderer = NULL;


int main(int argc, char* argv[])
{
	Terra T; // Variavel com os parametros da terra
	int i, nray1, nray2, nray3, nray4; // variaveis de controle
	double altM; //altura maxima da atmosfera
	double esp; 		// espessura de cada camada
	reta_xy L_ray1[6000], L_ray2[6000], L_ray3[6000], L_ray4[6000], Tangent_Line[1]; 	// declaração do conjunto de direções que representara o raio de luz na atmosfera
	ponto_xy Po1, Po2, Po3, Po4, Pcam, obs;			// pontos iniciais
	double ang1 = 0, ang2 = 0, ang3 = 0, ang4 = 0;		// angulos de direcao iniciais para os raios de luz
	double xM, z;						// profundidade da camera
	double expoente = 0;			// variavel ara o calculo dos indices de refracao			
	int luz_cor1[3], luz_cor2[3], luz_cor3[3], luz_cor4[3], luz_corT[3]; // parametros para cor

	luz_cor1[0] = 255;
	luz_cor1[1] = 255;
	luz_cor1[2] = 50;

	luz_cor2[0] = 255;
	luz_cor2[1] = 0;
	luz_cor2[2] = 255;

	luz_cor3[0] = 0;
	luz_cor3[1] = 55;
	luz_cor3[2] = 255;

	luz_cor4[0] = 0;
	luz_cor4[1] = 255;
	luz_cor4[2] = 25;

	luz_corT[0] = 0;
	luz_corT[1] = 255;
	luz_corT[2] = 25;

	/********************************************** atribuicao de parametros iniciais reais ********************************************/

	// Unidade de comprimento: Km
	T.C.x = 0;			// coordenada x do centro
	T.C.y = 0;			// coordenada y do centro
	T.raio = 6370;		// raio;
	altM = 0.1;			// altura maxima para as camadas. No caso, as camadas serao consideradas ate essa altura.
	T.atm_num = 5000;		// numero de camadas
	esp = altM / T.atm_num;	//calculo da espessura de cada camadax

	/*   Parametros 1  convergencia em x = -6 dn/dh = - 70.10^(-6) km^-1 (expoente = -0.2 * esp * (i - 1);) */
	/* Pontos de tangencia L_ray1[6447].A, L_ray2[4748].A */
	/*xM = -6;
	Po1.x = 29.998955343212295;		// atribuindo ponto de Teste x
	Po1.y = 6369.993836600876421; 	// atribuindo ponto de Teste y
	Po2.x = 30.000546808318862;
	Po2.y = 6369.976838917313216;

	ang1 = 3.140263745777878;
	ang2 = 3.139763422836058;*/



	/* Parametros 2 - convergencia em (-5, 6370)  dn/dh= - 70.10^(-6) km^-1 (expoente = -0.2 * esp * (i - 1);) */
	// Pontos de tangencia L_ray1[6765].A, L_ray2[6155].A
	/*xM = -5;
	Po1.x = 29.998718633117562;		// atribuindo ponto de Teste x
	Po1.y = 6369.997017750913074; 	// atribuindo ponto de Teste y
	Po2.x = 29.999178255219043;
	Po2.y = 6369.990915518699694;

	ang1 = 3.140387436073874;
	ang2 = 3.140205801214904;*/


	/*   Parametros 3 - convergencia em X = -10 dn/dh = - 70.10^(-6) km^-1 (expoente = -0.2 * esp * (i - 1);) */
	// Pontos de tangencia L_ray1[7028].A, L_ray2[5229].A
	/*xM = -10;
	Po1.x = 30;						// atribuindo ponto de Teste x
	Po1.y = 6369.999646487176506; 	// atribuindo ponto de Teste y
	Po2.x = 30;
	Po2.y = 6369.981651808505376;

	ang1 = 3.140406628350369;
	ang2 = 3.139906161607758;*/


	/*   Parametros 4 - convergencia em (-5, 6370)  dn/dh = - 50.10^(-6) km^-1 (expoente = -0.143 * esp * (i - 1);) */
	// Pontos de tangencia L_ray1[9590].A, L_ray2[8370].A
	/*xM = -5;
	Po1.x = 34.999486426587538;						// atribuindo ponto de Teste x
	Po1.y = 6369.999754697519165; 	// atribuindo ponto de Teste y
	Po2.x = 35.000154753828681;
	Po2.y = 6369.987550841238772;

	ang1 = 3.140670863122173;
	ang2 = 3.140349678914915;*/


	/*   Parametros 5 - convergencia em X = -5  dn/dh = - 70.10^(-6) km^-1 (expoente = -0.2 * esp * (i - 1);) */
	// Pontos de tangencia L_ray1[9527].A, L_ray2[7927].A
	/*xM = -5;
	Po1.x = 34.999783820349080;						// atribuindo ponto de Teste x
	Po1.y = 6369.999123053975381; 	// atribuindo ponto de Teste y
	Po2.x = 35.001276623758166;
	Po2.y = 6369.983114610075063;

	ang1 = 3.140273834419551;
	ang2 = 3.139860657934644;*/

	/*   Parametros 6 - convergencia em X = -10  dn/dh = - 105.10^(-6) km^-1 (expoente = -0.3 * esp * (i - 1);) */
	// Pontos de tangencia L_ray1[9601].A, L_ray2[7724].A
	/*xM = -10;
	Po1.x = 33.001040579043206;						// atribuindo ponto de Teste x
	Po1.y = 6370.010531517002164; 	// atribuindo ponto de Teste y
	Po2.x = 33.000465783526906;
	Po2.y = 6369.991764242880890;

	ang1 = 3.139898756854048;
	ang2 = 3.139444052893975;*/

	/*   Parametros 7 - convergencia em X = -20  dn/dh = - 140.10^(-6) km^-1 (expoente = -0.4 * esp * (i - 1);) */
	// Pontos de tangencia L_ray1[9939].A, L_ray2[6880].A
	/*xM = -20;
	Po1.x = 24.998248465233139;						// atribuindo ponto de Teste x
	Po1.y = 6370.050344365408819; 	// atribuindo ponto de Teste y
	Po2.x = 25.000889136201327;
	Po2.y = 6370.019743766170905;

	ang1 = 3.140473891548510;
	ang2 = 3.139764450772945;*/

	/*   Parametros 8 - convergencia em (-10, 6370)  dn/dh = - 140.10^(-6) km^-1 (expoente = -0.4 * esp * (i - 1);) */
	// Pontos de tangencia L_ray1[9860].A, L_ray2[7578].A
	/*xM = -10;
	Po1.x = 36.000911675480737;						// atribuindo ponto de Teste x
	Po1.y = 6369.996873765185228; 	// atribuindo ponto de Teste y
	Po2.x = 36.000778098658301;
	Po2.y = 6369.974054155538397;

	ang1 = 3.138716654746058;
	ang2 = 3.138197687193504;*/

	/*   Parametros 9 - Refracao critica  dn/dh = - 157.10^(-6) km^-1 (expoente = -0.448697 * esp * (i - 1);) */
	/*xM = -33;
	Po1.x = 0;						// atribuindo ponto de Teste x
	Po1.y = 6370 + 0.00005;; 	// atribuindo ponto de Teste y
	Po2.x = 0;
	Po2.y = 6370 + 0.00005;;

	ang1 = 0;
	ang2 = 0;*/

	/*   Parametros 8 - Convergencia em (-8, 6370)  Condicao de duto (armadilha)  dn/dh = - 175.0^(-6) km^-1 (expoente = -0.5 * esp * (i - 1);) */
	// Pontos de tangencia L_ray1[9364].A, L_ray2[3910].A
	/*xM = -8;
	Po1.x = 35.001331470185363;						// atribuindo ponto de Teste x
	Po1.y = 6369.997484525296386; 	// atribuindo ponto de Teste y
	Po2.x = 34.998158729139845;
	Po2.y = 6369.942961131062475;

	ang1 = 3.138123667182083;
	ang2 = 3.136794539305753;*/

	/*   Parametros 9 - Convergencia em (-5, 6370)  Condicao de duto (armadilha)  dn/dh = - 70.10^(-6) km^-1 (expoente = -0.2 * esp * (i - 1);) */
	// Pontos de tangencia L_ray1[9527].A, L_ray2[7927].A
	xM = -5;
	obs.x = -5;
	obs.y = 6370;
	Po1.x = 34.999628410493600;						// atribuindo ponto de Teste x
	Po1.y = 6369.996843873450416; 	// atribuindo ponto de Teste y
	Po2.x = 34.998889843632981;
	Po2.y = 6369.983877735638089;
	Po3.x = 34.999628410493600;						
	Po3.y = 6369.996843873450416; 	
	Po4.x = 34.998889843632981;
	Po4.y = 6369.983877735638089;
	
	ang1 = 3.1402150;
	ang2 = 3.139880;
	ang3 = 3.1401284;
	ang4 = 3.1398018;

	/*  Definicao de parametros  */
	/*xM = -15;
	Po1.x = -5;						// atribuindo ponto de Teste x
	Po1.y = 6370; 	// atribuindo ponto de Teste y
	Po2.x = -5;
	Po2.y = 6370;

	ang1 = -0.00021;
	ang2 = -0.00001;*/


	Paux.x = 35;
	Paux.y = 6370;
	z = 0.1;
	//z =  -0.04406;
	Pcam = Paux; // determinacao do ponto de visualizacao da camera dentro do sistema de coordenadas
	Paux2.x = 35;
	Paux2.y = 6370;

	T.n[0] = 1.0;

	for (i = 0; i < T.atm_num; i++)
	{
		T.camada_h[i] = T.raio + (i + 1) * esp;			//calculo da altura de cada camada.
		expoente = -0.2 * esp * (i - 1);
		if (i >= 1) T.n[i] = 1 + 0.00035 * pow(EXP, expoente);

		//if (i < 5000) printf("\n  h(%d) = %.10lf    n(%d) = %.10lf", i, T.camada_h[i], i, T.n[i]);	//linha de teste - testa se as alturas estão sendo atribuidas corretamente
		//if (i >= 1) printf("   camada %d  dn/dt = %.15lf", i, (T.n[i] - T.n[i - 1]) / (1000*(T.camada_h[i] - T.camada_h[i - 1])));
	}
	//printf("\n  h(10) = %.10lf    n(10) = %.10lf \n  dn/dh = %.10lf \n espReal = %.10lf  \n espCalc = %.10lf", T.camada_h[3], T.n[3], (T.n[3] - T.n[2]) / (T.camada_h[3] - T.camada_h[2]), esp, T.camada_h[3] - T.camada_h[2]);

	/*************************************************************************************************************************************/


	// Inicializa o SDL
	SDL_Init(SDL_INIT_VIDEO);

	// Cria uma janela
	SDL_Window* window = SDL_CreateWindow("Desenhar Linha com SDL", SDL_WINDOWPOS_CENTERED, SDL_WINDOWPOS_CENTERED, 2500, 1400, SDL_WINDOW_SHOWN);

	// Cria um renderer para a janela
	renderer = SDL_CreateRenderer(window, -1, SDL_RENDERER_ACCELERATED);

	// Calcula a trajetoria de cada raio na atmosfera e retorna o numero de raios calculados.
	nray1 = teste_2(L_ray1, &T, Po1, ang1, xM);
	nray2 = teste_2(L_ray2, &T, Po2, ang2, xM);
	nray3 = teste_2(L_ray3, &T, Po3, ang3, xM);
	nray4 = teste_2(L_ray4, &T, Po4, ang4, xM);

	//Calculo do ponto de tangencia
	Tangent_Line[0].B = Earth_tangent_point(obs, &T);

	// Mostraa as distancias relevantes
	DataShow(obs, Tangent_Line[0].B, L_ray1[0].A, L_ray2[0].A, L_ray1[4649].A, L_ray2[4001].A);

	// Calculo do ponto de tangencia a partir do ponto de observacao
	Tangent_Line[0] = Earth_tangent_line(obs, &T);
	//Tangent_Line[0] = Earth_tangent_line(Paux2, &T);

		// Imprime na tela alguns parametros recferentes ao raio
	for (i = 0; i <= nray1; i++)
	{
		//if (i < 1000 || i >= 6000) printf("\n  P%d = (%.15lf; %.15lf)    d = %.10lf", i, L_ray1[i].B.x, L_ray1[i].B.y, sqrt(L_ray1[i].B.x * L_ray1[i].B.x + L_ray1[i].B.y * L_ray1[i].B.y));
	}

	for (i = 0; i <= nray2; i++)
	{
		//if (i < 1000 || i >= 7000) printf("\n  P%d = (%.15lf; %.15lf)    d = %.10lf", i, L_ray2[i].B.x, L_ray2[i].B.y, sqrt(L_ray2[i].B.x * L_ray2[i].B.x + L_ray2[i].B.y * L_ray2[i].B.y));
	}

	/*printf("\n  Altura do ponto inicial    d = %.10lf", sqrt(Paux2.x * Paux2.x + Paux2.y * Paux2.y) - 6370);
	printf("\n\nangulo de visada = %.15lf", atan2((Tangent_Line[0].B.y - Tangent_Line[0].A.y), (Tangent_Line[0].B.x - Tangent_Line[0].A.x)));


	printf("\n\nM1 = (%.15lf, %.15lf)", (L_ray1[9493].B.x + L_ray1[9494].B.x) / 2, (L_ray1[9493].B.y + L_ray1[9494].B.y) / 2);
	printf("\n\nM2 = (%.15lf, %.15lf)", (L_ray2[8196].B.x + L_ray2[8197].B.x) / 2, (L_ray2[8196].B.y + L_ray2[8197].B.y) / 2);

	printf("\n\nang1 = %.15lf", atan2((L_ray1[9493].B.y - L_ray1[9494].B.y), (L_ray1[9493].B.x - L_ray1[9494].B.x)));
	printf("\n\nang2 = %.15lf", atan2((L_ray2[8196].B.y - L_ray2[8197].B.y), (L_ray2[8196].B.x - L_ray2[8197].B.x)));*/

	// ajusa as cores e atualiza o renderer
	SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
	SDL_RenderClear(renderer);
	SDL_RenderPresent(renderer);

	viewport_Earth_render(&T, Pcam.x, Pcam.y, z);
	SDL_RenderPresent(renderer);

	viewport_LightRay_render(L_ray1, Pcam.x, Pcam.y, z, nray1, luz_cor1);
	viewport_LightRay_render(L_ray2, Pcam.x, Pcam.y, z, nray2, luz_cor2);
	viewport_LightRay_render(L_ray3, Pcam.x, Pcam.y, z, nray3, luz_cor3);
	viewport_LightRay_render(L_ray4, Pcam.x, Pcam.y, z, nray4, luz_cor4);
	viewport_LightRay_render(Tangent_Line, Pcam.x, Pcam.y, z, 0, luz_corT); //visualizacao da linha tangente


	SDL_RenderPresent(renderer);
	SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
	SDL_RenderClear(renderer);

	i = 0;

	SDL_Event event1;

	while (1) {
		SDL_PollEvent(&event1);
		if (event1.type == SDL_KEYDOWN) {
			break;
		}


	}

	while (1) {
		SDL_PollEvent(&event1);
		if (event1.type == SDL_KEYDOWN) {
			break;
		}

		if (Pcam.x > xM) Pcam.x = Pcam.x - 0.02;
		//if (Pcam.x < 35) Pcam.x = Pcam.x + 0.02;
		else Pcam.x = xM;
		//else Pcam.x = 35;

		viewport_Earth_render(&T, Pcam.x, Pcam.y, z);
		viewport_LightRay_render(L_ray1, Pcam.x, Pcam.y, z, nray1, luz_cor1);
		viewport_LightRay_render(L_ray2, Pcam.x, Pcam.y, z, nray2, luz_cor2);
		viewport_LightRay_render(L_ray3, Pcam.x, Pcam.y, z, nray3, luz_cor3);
		viewport_LightRay_render(L_ray4, Pcam.x, Pcam.y, z, nray4, luz_cor4);
		//viewport_LightRay_render(Tangent_Line, Pcam.x, Pcam.y, z, 0, luz_corT); //visualizacao da linha tangente
		SDL_RenderPresent(renderer);
		SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
		SDL_RenderClear(renderer);


	}

	while (1) {
		SDL_PollEvent(&event1);
		if (event1.type == SDL_KEYDOWN) {
			break;
		}

		if (z < 0) z = 0;
		viewport_Earth_render(&T, Pcam.x, Pcam.y, z);
		viewport_LightRay_render(L_ray1, Pcam.x, Pcam.y, z, nray1, luz_cor1);
		viewport_LightRay_render(L_ray2, Pcam.x, Pcam.y, z, nray2, luz_cor2);
		viewport_LightRay_render(L_ray3, Pcam.x, Pcam.y, z, nray3, luz_cor3);
		viewport_LightRay_render(L_ray4, Pcam.x, Pcam.y, z, nray4, luz_cor4);
		viewport_LightRay_render(Tangent_Line, Pcam.x, Pcam.y, z, 0, luz_corT); //visualizacao da linha tangente
		SDL_RenderPresent(renderer);
		SDL_SetRenderDrawColor(renderer, 0, 0, 0, 255);
		SDL_RenderClear(renderer);
		z = z - 0.05;

	}

	// Libera recursos
	SDL_DestroyRenderer(renderer);
	SDL_DestroyWindow(window);
	SDL_Quit();


	return 0;
}

int teste_2(reta_xy* L_ray, Terra* T, ponto_xy P, double ang, double xM)
{
	int i;
	int camada = 0, aux;		// variaveis de controle para armazenar a camada em que o raio se encontra.
	char pos;				// variavel de controle. Verifica aa posicao do vetor diretor em relacao a a reta normal.
	double a, b, c;			// parametros da equacao gerada para t no calculo das interseccoes.
	double ang_i, ang_r;	// angulos de incidencia e refracao
	double H = 0, h = 0;			// altura da camada suerior e altura da camada inferior ao ponto inicial
	double t;				// parametro para a equação vetorial da reta
	double d;				// distancia do ponto inicial ao centro da Terra (tambem usada como variavel auxiliar nos calculos)
	double ang_reta = 0;
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

	while (i < 11999)
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
			if (camada > (*T).atm_num - 10) printf("\n\n Reflexão total superior na camada : %d \n iteracao: %d \n Coordenadas do ponto de reflexao (%.10lf, %.10lf) \n", camada, i, L_ray[i].A.x, L_ray[i].A.y);
			if (camada < 55) printf("\n\n Reflexão total inferior na camada : %d \n iteracao: %d \n Coordenadas do ponto de reflexao (%.10lf, %.10lf) \n", camada, i, L_ray[i].A.x, L_ray[i].A.y);

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

		//Verifica se o raio chegou ao ponto final
		if (L_ray[i].B.x <= xM)
		{
			ang_reta = (L_ray[i].B.y - L_ray[i].A.y) / (L_ray[i].B.x - L_ray[i].A.x);
			L_ray[i].B.x = xM;
			L_ray[i].B.y = ang_reta * (xM - L_ray[i].A.x) + L_ray[i].A.y;
			return i;
		}

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

		// Define a cor da linha do raio
		SDL_SetRenderDrawColor(renderer, luz_cor[0], luz_cor[1], luz_cor[2], 255);


		// Desenha uma linha de raio no renderer
		drawLine(renderer, xi, yi, xf, yf);


	}

}


void viewport_Earth_render(Terra* T, double xpp, double ypp, double zpp)
{
	int i = 0, int_delta = 0, camada_n = 0, var = 250;
	double d = 0.05;
	double int_set;
	double Raio, CamadaH;
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

	int_set = 0.01;
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


void drawCircle3(SDL_Renderer* renderer, int centerX, int centerY, int Raio)
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
		SDL_RenderDrawLine(renderer, centerX + x, centerY - y, centerX + x, 1400);
	}

}

ponto_xy Earth_tangent_point(ponto_xy PO, Terra* T)
{
	double angP, angT, m;
	ponto_xy TangentP;

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

	TangentP.x = (m * m * PO.x - m * PO.y) / (1 + m * m);
	TangentP.y = m * (TangentP.x - PO.x) + PO.y;

	return TangentP;
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

void DataShow(ponto_xy Obs, ponto_xy TangentP, ponto_xy TUp, ponto_xy TDown, ponto_xy RefP1, ponto_xy RefP2)
{
	double dist, ang;

	//altura do observador
	dist = sqrt(Obs.x * Obs.x + Obs.y * Obs.y) - 6370;
	printf("\n Altura do ponto de observação: %.3lf km", dist);

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