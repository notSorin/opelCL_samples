//Copyright (C) Sorin Draghici <dsorin95@gmail.com>
#define MAX_WINDOW_SIZE 5*5

void buble_sort(float array[], int size)
{
	int i, j;
	float tmp;

	for (i = 1; i < size; i++)
		for (j = 0 ; j < size - i; j++)
			if (array[j] > array[j + 1])
			{
				tmp = array[j];
				array[j] = array[j + 1];
				array[j + 1] = tmp;
			}
}

__kernel void remove_noise(__global float *im, __global float *image_out, float thredshold, int window_size, int height, int width)
{
	int filaGlobal = get_global_id(0), columnaGlobal = get_global_id(1), ii, jj, ws2 = (window_size - 1) >> 1; //ws2 = (3 - 1) / 2 = 1
	float window[MAX_WINDOW_SIZE], median;

	if(filaGlobal >= ws2 && filaGlobal < height - ws2 && columnaGlobal >= ws2 && columnaGlobal < width - ws2)
	{
		for (ii = -ws2; ii <= ws2; ii++)
			for (jj = -ws2; jj <= ws2; jj++)
				window[(ii + ws2) * window_size + jj + ws2] = im[(filaGlobal + ii) * width + columnaGlobal + jj];

		buble_sort(window, window_size * window_size);
		median = window[(window_size * window_size - 1) >> 1];

		if (fabs(((float)(median - im[filaGlobal * width + columnaGlobal])) / median) <= thredshold)
			image_out[filaGlobal * width + columnaGlobal] = im[filaGlobal * width + columnaGlobal];
		else
			image_out[filaGlobal * width + columnaGlobal] = median;
	}
	else
		if(filaGlobal == 0 || filaGlobal == height - 1 || columnaGlobal == 0 || columnaGlobal == width - 1)
			image_out[filaGlobal * width + columnaGlobal] = im[filaGlobal * width + columnaGlobal];
}
