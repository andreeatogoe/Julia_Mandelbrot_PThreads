#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <math.h>


char *in_filename_julia;
char *in_filename_mandelbrot;
char *out_filename_julia;
char *out_filename_mandelbrot;
int threads_nr;
pthread_barrier_t barrier;

// structura pentru un numar complex
typedef struct _complex {
	double a;
	double b;
} complex;

// structura pentru parametrii unei rulari
typedef struct _params {
	int is_julia, iterations;
	double x_min, x_max, y_min, y_max, resolution;
	complex c_julia;
} params;

typedef struct _threads_arg {
	int thread_id;
	params par_j;
	params par_m;
	int width_j, width_m, height_j, height_m;
	int **result_j, **result_m;
} threads_arg, *arg_t;

params par_j;
params par_m;
int width_j, width_m, height_j, height_m;
int **result_j, **result_m;

// citeste argumentele programului
void get_args(int argc, char **argv)
{
	if (argc < 6) {
		printf("Numar insuficient de parametri:\n\t"
				"./tema1 fisier_intrare_julia fisier_iesire_julia "
				"fisier_intrare_mandelbrot fisier_iesire_mandelbrot "
				"numar de thread-uri\n");
		exit(1);
	}

	in_filename_julia = argv[1];
	out_filename_julia = argv[2];
	in_filename_mandelbrot = argv[3];
	out_filename_mandelbrot = argv[4];
	threads_nr = *argv[5] - '0';
}

// citeste fisierul de intrare
void read_input_file(char *in_filename, params* par)
{
	FILE *file = fopen(in_filename, "r");
	if (file == NULL) {
		printf("Eroare la deschiderea fisierului de intrare!\n");
		exit(1);
	}

	fscanf(file, "%d", &par->is_julia);
	fscanf(file, "%lf %lf %lf %lf",
			&par->x_min, &par->x_max, &par->y_min, &par->y_max);
	fscanf(file, "%lf", &par->resolution);
	fscanf(file, "%d", &par->iterations);

	if (par->is_julia) {
		fscanf(file, "%lf %lf", &par->c_julia.a, &par->c_julia.b);
	}

	fclose(file);
}

// scrie rezultatul in fisierul de iesire
void write_output_file(char *out_filename, int **result, int width, int height)
{
	int i, j;

	FILE *file = fopen(out_filename, "w");
	if (file == NULL) {
		printf("Eroare la deschiderea fisierului de iesire!\n");
		return;
	}

	fprintf(file, "P2\n%d %d\n255\n", width, height);
	for (i = 0; i < height; i++) {
		for (j = 0; j < width; j++) {
			fprintf(file, "%d ", result[i][j]);
		}
		fprintf(file, "\n");
	}

	fclose(file);
}

// aloca memorie pentru rezultat
int **allocate_memory(int width, int height)
{
	int **result;
	int i;

	result = malloc(height * sizeof(int*));
	if (result == NULL) {
		printf("Eroare la malloc!\n");
		exit(1);
	}

	for (i = 0; i < height; i++) {
		result[i] = malloc(width * sizeof(int));
		if (result[i] == NULL) {
			printf("Eroare la malloc!\n");
			exit(1);
		}
	}

	return result;
}

// elibereaza memoria alocata
void free_memory(int **result, int height)
{
	int i;

	for (i = 0; i < height; i++) {
		free(result[i]);
	}
	free(result);
}

//impl alg Julia paralelizat
void run_julia_par(params *par, int **result, int width, int height, int id) {

	int w, h, i;
	int start, end;

	start = id * (double) width / threads_nr;
	end = fmin((id + 1) * (double) width / threads_nr, width);
	for (w = start; w < end; w++) {
		for (h = 0; h < height; h++) {
			int step = 0;
			complex z = { .a = w * par->resolution + par->x_min,
							.b = h * par->resolution + par->y_min };

			while (sqrt(pow(z.a, 2.0) + pow(z.b, 2.0)) < 2.0 && step < par->iterations) {
				complex z_aux = { .a = z.a, .b = z.b };

				z.a = pow(z_aux.a, 2) - pow(z_aux.b, 2) + par->c_julia.a;
				z.b = 2 * z_aux.a * z_aux.b + par->c_julia.b;

				step++;
			}

			result[h][w] = step % 256;
		}
	}
    //bariera
    pthread_barrier_wait(&barrier);

	// transforma rezultatul din coordonate matematice in coordonate ecran
	start = id * (double) (height / 2) / threads_nr;
	end = fmin((id + 1) * (double) (height / 2) / threads_nr, (height / 2));
	for (i = start; i < end; i++) {
		int *aux = result[i];
		result[i] = result[height - i - 1];
		result[height - i - 1] = aux;
	}
}

//impl alg Mandelbrot paralelizat
void run_mandelbrot_par(params *par, int **result, int width, int height, int id) {

	int w, h, i;
	int start, end;
	start = id * (double) width / threads_nr;
	end = fmin((id + 1) * (double) width / threads_nr, width);

	for (w = start; w < end; w++) {
		for (h = 0; h < height; h++) {
			complex c = { .a = w * par->resolution + par->x_min,
							.b = h * par->resolution + par->y_min };
			complex z = { .a = 0, .b = 0 };
			int step = 0;

			while (sqrt(pow(z.a, 2.0) + pow(z.b, 2.0)) < 2.0 && step < par->iterations) {
				complex z_aux = { .a = z.a, .b = z.b };

				z.a = pow(z_aux.a, 2.0) - pow(z_aux.b, 2.0) + c.a;
				z.b = 2.0 * z_aux.a * z_aux.b + c.b;

				step++;
			}

			result[h][w] = step % 256;
		}
	}
	//bariera
    pthread_barrier_wait(&barrier);

	// transforma rezultatul din coordonate matematice in coordonate ecran
	start = id * (double) (height / 2) / threads_nr;
	end = fmin((id + 1) * (double) (height / 2) / threads_nr, (height / 2));
	for (i = start; i < end; i++) {
		int *aux = result[i];
		result[i] = result[height - i - 1];
		result[height - i - 1] = aux;
	}
}

void *f(void *arg)
{
	int *id = (int*)arg;
    run_julia_par(&par_j, result_j, width_j, height_j, *id);

	run_mandelbrot_par(&par_m, result_m, width_m, height_m, *id);
    pthread_exit(NULL);
}


int main(int argc, char *argv[]) {
	
	int  id, r;

	// se citesc argumentele programului
	get_args(argc, argv);

	read_input_file(in_filename_julia, &par_j);

	width_j = (par_j.x_max - par_j.x_min) / par_j.resolution;
	height_j = (par_j.y_max - par_j.y_min) / par_j.resolution;
	result_j = allocate_memory(width_j, height_j);

	read_input_file(in_filename_mandelbrot, &par_m);

	width_m = (par_m.x_max - par_m.x_min) / par_m.resolution;
	height_m = (par_m.y_max - par_m.y_min) / par_m.resolution;
	result_m = allocate_memory(width_m, height_m);

	pthread_t threads[threads_nr];
	pthread_barrier_init(&barrier, NULL, threads_nr);
	int arg[threads_nr];

	for (id = 0; id < threads_nr; id++) {
        
        arg[id] = id;
        r = pthread_create(&threads[id], NULL, f, &arg[id]);
 
        if (r) {
            printf("Eroare la crearea thread-ului %d\n", id);
            exit(-1);
        }
    }

    for (id = 0; id < threads_nr; id++) {
        r = pthread_join(threads[id], NULL);
        if (r) {
            printf("Eroare la asteptarea thread-ului %d\n", id);
            exit(-1);
        }
    }
    write_output_file(out_filename_julia, result_j, width_j, height_j);
    write_output_file(out_filename_mandelbrot, result_m, width_m, height_m);
    free_memory(result_j, height_j);
    free_memory(result_m, height_m);
    pthread_barrier_destroy(&barrier);

 
    pthread_exit(NULL);
	
	return 0;
}