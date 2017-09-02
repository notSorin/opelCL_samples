//Copyright (C) Sorin Draghici <dsorin95@gmail.com>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <GL/glut.h>
//#include <GL/GL.h>
//#include <GL/GLU.h>
#include "CL/cl.h"

//#ifndef DEVICE
//#define DEVICE CL_DEVICE_TYPE_GPU
//#endif

#define RUN_SERIAL     0
#define RUN_OPENCL_CPU 1
#define RUN_OPENCL_GPU 2
#define VAL 255

void set_texture();

typedef struct
{
	unsigned char r, g, b;
}rgb_t;

int shots = 1, run_mode, gwin, width, height, tex_w, tex_h,
color_rotate = 0, saturation = 1, invert = 0, max_iter = 256;
double scale = 1. / 256, cx = -.6, cy = 0;
GLuint texture;
rgb_t **tex = 0;

/* Time */
#include <sys/time.h>
#include <sys/resource.h>

static struct timeval tv0;

double getMicroSeconds()
{
double t;
gettimeofday(&tv0, (struct timezone*)0);
t = ((tv0.tv_usec) + (tv0.tv_sec)*1000000);

return (t);
}

void render(void)
{
	double	x = (double)width / tex_w, y = (double)height / tex_h;

	glClear(GL_COLOR_BUFFER_BIT);
	glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
	glBindTexture(GL_TEXTURE_2D, texture);
	glBegin(GL_QUADS);

	glTexCoord2f(0, 0);
	glVertex2i(0, 0);

	glTexCoord2f(x, 0);
	glVertex2i(width, 0);

	glTexCoord2f(x, y);
	glVertex2i(width, height);

	glTexCoord2f(0, y);
	glVertex2i(0, height);

	glEnd();
	glFlush();
	glFinish();
}

void screen_shot()
{
	char fn[100];
	int i;
	FILE *fp;

	//sprintf_s(fn, 100, "screen%03d.ppm", shots++);
	sprintf(fn, "screen%03d.ppm", shots++);
	//fopen_s(&fp, fn, "w");
	fp = fopen(fn, "w");
	fprintf(fp, "P6\n%d %d\n255\n", width, height);

	for (i = height - 1; i >= 0; i--)
		fwrite(tex[i], 1, width * 3, fp);

	fclose(fp);
	printf("%s written\n", fn);
}

void keypress(unsigned char key, int x, int y)
{
	switch (key)
	{
	case 'q':
		glFinish();
		glutDestroyWindow(gwin);

		return;
	case 27:
		scale = 1. / 256;
		cx = -.6;
		cy = 0;

		break;
	case 'r':
		color_rotate = (color_rotate + 1) % 6;

		break;
	case '>':
	case '.':
		max_iter += 64;

		if (max_iter > 1 << 15)
			max_iter = 1 << 15;

		printf("max iter: %d\n", max_iter);

		break;
	case '<':
	case ',':
		max_iter -= 64;

		if (max_iter < 64)
			max_iter = 64;

		printf("max iter: %d\n", max_iter);

		break;
	case 'm':
		saturation = 1 - saturation;

		break;
	case 'i':
		screen_shot();

		return;
	case 'z':
		max_iter = 4096;

		break;
	case 'x':
		max_iter = 128;

		break;
	case 's':
		run_mode = RUN_SERIAL;

		break;
	case 'c':
		run_mode = RUN_OPENCL_CPU;

		break;
	case 'g':
		run_mode = RUN_OPENCL_GPU;

		break;
	case ' ':
		invert = !invert;
	}

	set_texture();
}

void hsv_to_rgb(int hue, int min, int max, rgb_t *p)
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

	double h = fmod(color_rotate + 1e-4 + 4.0 * (hue - min) / (max - min), 6), c = VAL * saturation, X = c * (1 - fabs(fmod(h, 2) - 1));

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

int initCL = 0;
cl_device_id device_id = NULL;
cl_int err;
cl_context context;
cl_command_queue command_queue;
cl_program program;
cl_kernel kernel;

void initOCL()
{
	//variables used to read kernel source file
	FILE *fp;
	size_t filelen, readlen;
	char *kernel_src; //char string to hold kernel source

	//read the kernel
	//fopen_s(&fp, "mandel_kernel.cl", "r");
	fp = fopen("mandel_kernel.cl", "r");
	fseek(fp, 0L, SEEK_END);

	filelen = ftell(fp);

	rewind(fp);

	kernel_src = malloc(sizeof(char) * (filelen + 1));
	readlen = fread(kernel_src, 1, filelen, fp);
	fclose(fp);

	if (readlen != filelen)
	{
		printf("error reading file\n");

		exit(EXIT_FAILURE);
	}

	//ensure the string is NULL terminated
	kernel_src[filelen] = '\0';

	//Set up platform and GPU device
	cl_uint numPlatforms;

	//Find number of platforms
	err = clGetPlatformIDs(0, NULL, &numPlatforms);

	if (err != CL_SUCCESS || numPlatforms <= 0)
	{
		printf("Error: Failed to find a platform!\n\n");

		exit(EXIT_FAILURE);
	}

	//Get all platforms
	cl_platform_id *platform = (cl_platform_id *)malloc(numPlatforms * sizeof(cl_platform_id));

	err = clGetPlatformIDs(numPlatforms, platform, NULL);

	if (err != CL_SUCCESS || numPlatforms <= 0)
	{
		printf("Error: Failed to get the platform!\n\n");

		exit(EXIT_FAILURE);
	}

	//Secure a GPU
	int i;

	for (i = 0; i < (int)numPlatforms; i++)
	{
		err = clGetDeviceIDs(platform[i], run_mode, 1, &device_id, NULL);

		if (err == CL_SUCCESS)
			break;
	}

	if (device_id == NULL)
	{
		printf("Error: Failed to create a device group!\n\n");

		exit(EXIT_FAILURE);
	}

	//Create a compute context 
	context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);

	if (!context)
	{
		printf("Error: Failed to create a compute context!\n\n");

		exit(EXIT_FAILURE);
	}

	//create command queue 
	command_queue = clCreateCommandQueue(context, device_id, 0, &err);

	if (err != CL_SUCCESS)
	{
		printf("Unable to create command queue. Error Code=%d\n", err);

		exit(EXIT_FAILURE);
	}

	//create program object from source.
	// kernel_src contains source read from file earlier
	program = clCreateProgramWithSource(context, 1, (const char **)&kernel_src, NULL, &err);

	free(kernel_src);

	if (err != CL_SUCCESS)
	{
		printf("Unable to create program object. Error Code=%d\n", err);

		exit(EXIT_FAILURE);
	}

	err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);

	char buffer[2048];
	size_t len;
	
	if (err != CL_SUCCESS)
	{
		printf("Build failed. Error Code=%d\n", err);


		// get the build log
		clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
		printf("--- Build Log -- \n %s\n", buffer);
		//system("pause");
		exit(EXIT_FAILURE);
	}
	
	//clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
	//printf("--- Build Log -- \n %s\n", buffer);

	kernel = clCreateKernel(program, "mandelKernel", &err);

	if (err != CL_SUCCESS)
	{
		printf("Unable to create kernel object. Error Code=%d\n", err);

		exit(EXIT_FAILURE);
	}

	initCL = 1;
}

double calc_mandel_opencl()
{
	if (!initCL)
		initOCL();	

	cl_mem deviceTexIn;// , deviceTexOut;
	// create buffer objects to input and output args of kernel function
	deviceTexIn = clCreateBuffer(context, CL_MEM_READ_WRITE, tex_w * tex_h * sizeof(rgb_t), NULL, NULL);
	err = clEnqueueWriteBuffer(command_queue, deviceTexIn, CL_TRUE, 0, tex_w * tex_h * sizeof(rgb_t), &tex[tex_h], 0, NULL, NULL);
	
	if (err != CL_SUCCESS)
	{
		printf("Error enqueuing write buffer command. Error Code=%d\n", err);

		exit(1);
	}

	//deviceTexOut = clCreateBuffer(context, CL_MEM_WRITE_ONLY, tex_w * tex_h * sizeof(rgb_t), NULL, NULL);

		/*clSetKernelArg(kernel, 1, sizeof(cl_mem), &deviceTexOut) ||*/
	if (clSetKernelArg(kernel, 0, sizeof(cl_mem), &deviceTexIn) ||
		clSetKernelArg(kernel, 1, sizeof(cl_int), &width) ||
		clSetKernelArg(kernel, 2, sizeof(cl_int), &height) ||
		clSetKernelArg(kernel, 3, sizeof(cl_double), &scale) ||
		clSetKernelArg(kernel, 4, sizeof(cl_double), &cy) ||
		clSetKernelArg(kernel, 5, sizeof(cl_double), &cx) ||
		clSetKernelArg(kernel, 6, sizeof(cl_int), &max_iter) ||
		clSetKernelArg(kernel, 7, sizeof(cl_int), &invert) ||
		clSetKernelArg(kernel, 8, sizeof(cl_int), &saturation) ||
		clSetKernelArg(kernel, 9, sizeof(cl_int), &color_rotate) ||
		clSetKernelArg(kernel, 10, sizeof(cl_int), &tex_w) != CL_SUCCESS)
	{
		printf("Unable to set kernel arguments. Error Code=%d\n", err);

		exit(EXIT_FAILURE);
	}

	/*printf("width: %d\n", width);
	printf("height: %d\n", height);
	printf("scale: %f\n", scale);
	printf("cy: %f\n", cy);
	printf("cx: %f\n", cx);
	printf("max_iter: %d\n", max_iter);
	printf("invert: %d\n", invert);
	printf("saturatio: %d\n", saturation);
	printf("color_rotate: %d\n\n", color_rotate);*/

	size_t global[2] = {height, width};
	
	double t0 = getMicroSeconds(), t1;
	
	t0 = getMicroSeconds();
	err = clEnqueueNDRangeKernel(command_queue, kernel, 2, NULL, global, NULL, 0, NULL, NULL);
	// wait for the command to finish
	clFinish(command_queue);
	
	t1 = getMicroSeconds();
	
	if (err != CL_SUCCESS)
	{
		printf("Unable to enqueue kernel command. Error Code=%d\n", err);
		exit(1);
	}

	// read the output back to host memory
	err = clEnqueueReadBuffer(command_queue, deviceTexIn, CL_TRUE, 0, tex_w * tex_h * sizeof(rgb_t), &tex[tex_h], 0, NULL, NULL);
		
	if (err != CL_SUCCESS)
	{
		printf("Error enqueuing read buffer command. Error Code=%d\n", err);
		exit(1);
	}

	return t1 - t0;//getMicroSeconds() - t0;
	//return(1.0);
}

double calc_mandel()
{
	int i, j, iter;
	double x, y, zx, zy, zx2, zy2, t0 = getMicroSeconds();

	for (i = 0; i < height; i++)
		for (j = 0; j < width; j++)
		{
			y = (i - height / 2) * scale + cy;
			x = (j - width / 2) * scale + cx;
			zx = zy = zx2 = zy2 = 0;

			for (iter = 0; iter < max_iter && zx2 + zy2 <= max_iter; iter++)
			{
				zy = 2 * zx * zy + y;
				zx = zx2 - zy2 + x;
				zx2 = zx * zx;
				zy2 = zy * zy;
			}

			hsv_to_rgb(iter, 0, max_iter, &tex[i][j]);
		}

	return getMicroSeconds() - t0;
	//return 1.0;
}

void alloc_tex()
{
	int i, ow = tex_w, oh = tex_h;

	for (tex_w = 1; tex_w < width; tex_w *= 2/*tex_w <<= 1*/);
	for (tex_h = 1; tex_h < height; tex_h *= 2/*tex_h <<= 1*/);

	if (tex_h != oh || tex_w != ow)
		tex = realloc(tex, tex_h * tex_w * 3 + tex_h * sizeof(rgb_t*));

	tex[0] = (rgb_t *)(tex + tex_h); //tex[0] apunta al primer rgb_t válido en tex

	for (i = 1; i < tex_h; i++)
		tex[i] = tex[i - 1] + tex_w;
}

void set_texture()
{
	double t;
	char title[128];

	alloc_tex();

	switch (run_mode)
	{
	case RUN_SERIAL:
	{
		t = calc_mandel();

		break;
	}
	case RUN_OPENCL_CPU:
	{
		t = calc_mandel_opencl();

		break;
	}
	case RUN_OPENCL_GPU:
	{
		t = calc_mandel_opencl();
	}
	}

	glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D, texture);
	glTexImage2D(GL_TEXTURE_2D, 0, 3, tex_w, tex_h, 0, GL_RGB, GL_UNSIGNED_BYTE, tex[0]);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	render();
	//sprintf_s(title, 128, "Mandelbrot: %5.2f fps (%ix%i)", 1000000 / t, width, height);
	sprintf(title, "Mandelbrot: %5.2f fps (%ix%i)", 1000000 / t, width, height);
	glutSetWindowTitle(title);
}

void mouseclick(int button, int state, int x, int y)
{
	if (state != GLUT_UP)
		return;

	cx += (x - width / 2) * scale;
	cy -= (y - height / 2) * scale;

	switch (button)
	{
	case GLUT_LEFT_BUTTON: /* zoom in */
		if (scale > fabs(x) * 1e-16 && scale > fabs(y) * 1e-16)
			scale /= 2;
		break;
	case GLUT_RIGHT_BUTTON: /* zoom out */
		scale *= 2;
		break;
		/* any other button recenters */
	}

	set_texture();
}


void resize(int w, int h)
{
	//printf("resize %d %d\n", w, h);
	width = w;
	height = h;

	glViewport(0, 0, w, h);
	glOrtho(0, w, 0, h, -1, 1);

	set_texture();
}

void init_gfx(int *c, char **v)
{
	glutInit(c, v);
	glutInitDisplayMode(GLUT_RGB);
	glutInitWindowSize(640, 480);

	gwin = glutCreateWindow("Mandelbrot");

	glutReshapeFunc(resize);
	glutKeyboardFunc(keypress);
	glutDisplayFunc(render);
	glutMouseFunc(mouseclick);
	glGenTextures(1, &texture);
	set_texture();
}

int main(int c, char **v)
{
	init_gfx(&c, v);
	printf("keys:\n\tr: color rotation\n\tm: monochrome\n\ti: screen shot\n\t"
		"s: serial code\n\tc: OpenCL CPU\n\tg: OpenCL GPU\n\t"
		"<, >: decrease/increase max iteration\n\tq: quit\n\tmouse buttons to zoom\n");
	glutMainLoop();

	return 0;
}
