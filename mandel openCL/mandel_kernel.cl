//Copyright (C) Sorin Draghici <dsorin95@gmail.com>
#pragma OPENCL EXTENSION cl_khr_fp64 : enable

typedef struct
{
	unsigned char r, g, b;
}rgb_t;

#define VAL 255

void hsv_to_rgb(int hue, int min, int max, rgb_t *p, int invert, int saturation, int color_rotate)
{
	if (min == max)
		max = min + 1;

	if (invert)
		hue = max - (hue - min);

	if (!saturation)
	{
		p->r = p->g = p->b = 255 * (max - hue) / (max - min);

		return;
	}

	double h = fmod((double)(color_rotate + 0.0001 + 4.0 * (hue - min) / (max - min)), (double)6),
		c = VAL * saturation, X = c * (1 - fabs(fmod((double)h, (double)2) - 1));

	p->r = p->g = p->b = 0;

	switch ((int)h)
	{
	case 0:
	{
		p->r = c;
		p->g = X;

		return;
	}
	case 1:
	{
		p->r = X;
		p->g = c;

		return;
	}
	case 2:
	{
		p->g = c;
		p->b = X;

		return;
	}
	case 3:
	{
		p->g = X;
		p->b = c;

		return;
	}
	case 4:
	{
		p->r = X;
		p->b = c;

		return;
	}
	default:
	{
		p->r = c;
		p->b = X;
	}
	}
}

__kernel void mandelKernel(__global rgb_t *texIn, int width, int height, double scale, 
	double cy, double cx, int max_iter, int invert, int saturation, int color_rotate, int tex_w)
{
	int i = get_global_id(0), j = get_global_id(1), iter = 0;
	double y, x, zx, zy, zx2, zy2;

	y = (i - height / 2) * scale + cy;
	x = (j - width / 2) * scale + cx;
	zx = zy = zx2 = zy2 = 0;
	
	for(iter = 0; iter < max_iter && zx2 + zy2 <= max_iter; iter++)
	{
		zy = 2 * zx * zy + y;
		zx = zx2 - zy2 + x;
		zx2 = zx * zx;
		zy2 = zy * zy;
	}
	
	rgb_t pixel = texIn[i * tex_w + j];
	
	hsv_to_rgb(iter, 0, max_iter, &pixel, invert, saturation, color_rotate);
	
	texIn[i * tex_w + j] = pixel;
}
