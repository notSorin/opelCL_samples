//Copyright (C) Sorin Draghici <dsorin95@gmail.com>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>
#include <limits.h>
#include <CL/cl.h>

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

typedef unsigned char uchar;
typedef unsigned short ushort;
typedef long long int lli;

#define CENSUS_HEIGHT 9
#define CENSUS_WIDTH  7
#define SMIN(a, b) ((a < b) ? a : b)

uchar* read_image_bin(char *name, int height, int width)
{
	uchar *im;
	FILE *fd;

	im = (uchar *)malloc(height * width * sizeof(uchar));
	fd = fopen(name, "r+");
	fread(im, height * width, sizeof(uchar), fd);
	fclose(fd);

	return(im);
}

void write_image_bin(char *name, uchar *im, int height, int width)
{
	FILE *fd;

	fd = fopen(name, "w+");
	fwrite(im, height * width, sizeof(uchar), fd);
	fclose(fd);
}

uchar *read_image_BMP(char *file_name, int *width, int *height)
{
    FILE *file;
    unsigned long size;                 // size of the image in bytes.
    unsigned long i;                    // standard counter.
    unsigned short int planes;          // number of planes in image (must be 1) 
    unsigned short int bpp;             // number of bits per pixel (must be 24)
    char temp;                          // temporary color storage for bgr-rgb conversion.

    // make sure the file is there.
    if ((file = fopen(file_name, "rb"))==NULL)
    {
		printf("File Not Found : %s\n",file_name);
		return 0;
    }

    // seek through the bmp header, up to the width/height:
    fseek(file, 18, SEEK_CUR);

    // read the width
    fread(width, 4, 1, file);
    
    // read the height 
    fread(height, 4, 1, file);
    
    // calculate the size (assuming 24 bits or 3 bytes per pixel).
    size = *width * *height * 3;

    // read the planes
    fread(&planes, 2, 1, file);
    
	if (planes != 1)
	{
		printf("Planes from %s is not 1: %u\n", file_name, planes);
		
		return 0;
    }

    // read the bits per pixel
    fread(&bpp, 2, 1, file);
    
	if (bpp != 24)
	{
		printf("Bpp from %s is not 24: %u\n", file_name, bpp);
		
		return 0;
    }
	
    // seek past the rest of the bitmap header.
    fseek(file, 24, SEEK_CUR);

    // read the data.
    unsigned char *image = (char *)malloc(size);
    
	if (image == NULL)
	{
		printf("Error allocating memory for color-corrected image data");
		
		return 0;	
    }

    if ((i = (unsigned long)fread(image, size, 1, file)) != 1)
	{
		printf("Error reading image data from %s.\n", file_name);
		
		return 0;
    }

    for (i = 0; i < size; i += 3)
	{ // reverse all of the colors. (bgr -> rgb)
		temp = image[i];
		image[i] = image[i + 2];
		image[i + 2] = temp;
    }

    // we're done.
    return image;
}

void write_image_BMP(uchar *data, char *filename, int width, int height)
{
	FILE *fd = fopen(filename, "wb");

	static unsigned char header[54] = {66,77,0,0,0,0,0,0,0,0,54,0,0,0,40,0,0,0,0,0,0,0,0,0,0,0,1,0,24}; //rest is zeroes
	unsigned int pixelBytesPerRow = width * 3 * sizeof(uchar);
	unsigned int *sizeOfFileEntry = (unsigned int *)&header[2];
	
	*sizeOfFileEntry = 54 + pixelBytesPerRow * height;  
	
	unsigned int *widthEntry = (unsigned int *)&header[18];    
	
	*widthEntry = width;
	
	unsigned int *heightEntry = (unsigned int *)&header[22];
	
	*heightEntry = height;
	
	fwrite(header, 1, 54, fd);
	
	static unsigned char zeroes[3] = {0,0,0}; //for padding    

	uchar *image = (uchar *)malloc(3 * sizeof(uchar) * width * height);
	int i;

	for (i = 0; i < width * height; i++)
	{
		image[3 * i]     = data[i]; //R 
		image[3 * i + 1] = data[i]; //G
		image[3 * i + 2] = data[i]; //B
	}

	fwrite(image, 1, pixelBytesPerRow * height, fd);
	fclose(fd);
	free(image);
}

uchar *image_RGB2BW(uchar *image_in, int height, int width)
{
	int i, j;
	uchar *imageBW = (uchar *)malloc(sizeof(uchar) * width * height);
	float R, B, G;

	for (i = 0; i < height; i++)
		for (j = 0; j < width; j++)
		{
			R = (float)(image_in[3 * (i * width + j)]);
			G = (float)(image_in[3 * (i * width + j) + 1]);
			B = (float)(image_in[3 * (i * width + j) + 2]);

			imageBW[i * width + j] = (uchar)(0.2989 * R + 0.5870 * G + 0.1140 * B);
		}

	return imageBW;
}

uchar *padding_image(uchar *im, int height, int width, int vert_pad, int hor_pad)
{
	uchar *im_pad;
	int i, j;
	int widthPAD = width + hor_pad;
	int vp_div_2 = (vert_pad - 1) / 2;
	int hp_div_2 = (hor_pad - 1) / 2;

	im_pad = (uchar *)malloc((height + vert_pad) * widthPAD * sizeof(uchar));

	/* upper */ 
	for (i = 0; i < vp_div_2; i++)
		for (j = 0; j < width + hor_pad; j++)
			im_pad[i * widthPAD + j] = 0;

	/* central */ 
	for (i = vp_div_2; i < height + vp_div_2; i++)
	{
		for (j = 0; j < hp_div_2; j++)
			im_pad[i * widthPAD + j] = 0;
		
		for (j = hp_div_2; j < width + hp_div_2; j++)
			im_pad[i * widthPAD + j] = im[(i - vp_div_2) * width + (j - hp_div_2)];

		for (j = width + hp_div_2; j < width + hor_pad; j++)
			im_pad[i * widthPAD + j] = 0;
	}

	/* lower */ 
	for (i = height + vp_div_2; i < height + vert_pad; i++)
		for (j = 0; j < width + hor_pad; j++)
			im_pad[i * widthPAD + j] = 0;

	return im_pad;
}

void censusTransform(uchar *maching_cost, lli *census_im1, lli *census_im2, uchar *im1, uchar *im2, int height, int width, int disp_min, int disp_max)
{
	int i, ii, j, jj, d_index;
	int widthPAD = width + CENSUS_WIDTH;
	int disparity_len  = disp_max - disp_min + 1;

	// 9x7 Census transform
	int top = (CENSUS_HEIGHT - 1) / 2;
	int left = (CENSUS_WIDTH - 1) / 2;
	uchar center1, center2, curr1, curr2;

	/* Compute with a window, make sure it has valid values */
	for(i = top; i < top + height; i++)
		for(j = left; j < left + width; j++)
		{
			// We need 63 bits -> unsigned long long int
			unsigned long long int census1 = 0;
			unsigned long long int census2 = 0;
			
			center1 = im1[i * widthPAD + j];
			center2 = im2[i * widthPAD + j];

			for(ii = -top; ii <= top; ii++)
				for(jj = -left; jj <= left; jj++)
				{
					curr1 = im1[(i + ii) * widthPAD + (j + jj)];

					// Compare to the center
					unsigned long long int tmp = curr1 >= center1;
					
					// Shift to the desired position
					tmp = tmp << ((ii + top) * CENSUS_WIDTH + (jj + left));
					
					// Add it to its place
					census1 = census1 | tmp;

					// Same proccess with the second image
					curr2 = im2[(i + ii) * widthPAD + (j + jj)];
					tmp = curr2 >= center2;
					tmp = tmp << ((ii + top) * CENSUS_WIDTH + (jj + left));
					census2 = census2 | tmp;
				}
			
			census_im1[(i - top) * width + (j - left)] = census1;
			census_im2[(i - top) * width + (j - left)] = census2;
		}

	for(i = 0; i < height; i++)
		for(j = 0; j < width; j++)
			for(d_index = 0; d_index < disparity_len; d_index++)
			{
				short d = d_index + disp_min;
			
				if (j - d >= 0 && j - d < width) // Compute hamming distance
				{
					maching_cost[disparity_len * (i * width + j) + d_index] = __builtin_popcount(census_im1[i * width + j] ^ census_im2[i * width + j - d]);
				}
				else
					maching_cost[disparity_len * (i * width + j) + d_index] = UCHAR_MAX;
			}
}

void SGM(uchar *maching_cost, ushort *S, ushort *L, ushort *min_L, int height, int width, int disp_min, int disp_max)
{
	int i, j, d_index;
	int disparity_len  = disp_max - disp_min + 1;
	ushort *L_lr = &(L[0 * (disparity_len * height * width)]); /* Left->Right */
	ushort *L_rl = &(L[1 * (disparity_len * height * width)]); /* Right->Left */
	ushort *L_ud = &(L[2 * (disparity_len * height * width)]); /* Up->Down */
	ushort *L_du = &(L[3 * (disparity_len * height * width)]); /* Down->Up */
	ushort prev, curr, next;
	ushort p1 = 30;
	ushort p2 = 30;

	/***************/
	/* Left->Right */
	/***************/
	for (i = 0; i < height; i++) //height x disparity_len
		for(d_index = 0; d_index < disparity_len; d_index++)
		{
			L_lr[disparity_len * i * width + d_index] = maching_cost[disparity_len * i * width + d_index];
			min_L[i * width] = USHRT_MAX;//SMIN(maching_cost[disparity_len * i * width + d_index], USHRT_MAX);
		}

	for (i = 0; i < height; i++) //height x (width - 1)
		for(j = 1; j < width; j++)
		{
			min_L[i * width + j] = USHRT_MAX;
			prev = USHRT_MAX;
			curr = L_lr[disparity_len * (i * width + j - 1)];

			for(d_index = 0; d_index < disparity_len - 1; d_index++)
			{
				next = L_lr[disparity_len * (i * width + j - 1) + d_index + 1]; // OJO cambiar aqui
				
				ushort cost_value = maching_cost[disparity_len * (i * width + j) + d_index]
					+ SMIN(SMIN(next + p1, prev + p1), SMIN(curr, min_L[i * width + j - 1] + p2))
					- min_L[i * width + j - 1]; // NOTE: min_L subtractiong to avoid overflow by adding-accumulation

				L_lr[disparity_len * (i * width + j) + d_index] = cost_value;

				min_L[i * width + j] = SMIN(cost_value, min_L[i * width + j]);
				prev = curr;
				curr = next;
			}
			
			ushort cost_value = maching_cost[disparity_len * (i * width + j) + disparity_len - 1]
				+ SMIN(prev + p1, SMIN(curr, min_L[i * width + j - 1] + p2))
				- min_L[i * width + j - 1];

			L_lr[disparity_len * (i * width + j) + disparity_len - 1] = cost_value;
			min_L[i * width + j] = SMIN(cost_value, min_L[i * width + j]);
		}

	/***************/
	/* Right->Left */
	/***************/
	for (i = 0; i < height; i++) //height x disparity_len
		for(d_index = 0; d_index < disparity_len; d_index++)
		{
			L_rl[disparity_len * (i * width + width - 1) + d_index]  = maching_cost[disparity_len * (i * width + width - 1) + d_index];
			min_L[i * width + width - 1] = USHRT_MAX; //SMIN(maching_cost[disparity_len * (i * width + width - 1) + d_index], USHRT_MAX);
		}

	for (i = height - 1; i >= 0; i--) //height x (width - 1)
		for(j = width - 2; j >= 0; j--)
		{
			min_L[i * width + j] = USHRT_MAX;
			prev = USHRT_MAX;
			curr = L_rl[disparity_len * (i * width + j + 1) + width];

			for(d_index = 0; d_index < disparity_len - 1; d_index++)
			{
				next = L_rl[disparity_len * (i * width + j + 1) + d_index + 1];
				
				ushort cost_value = maching_cost[disparity_len * (i * width + j) + d_index]
					+ SMIN(SMIN(next + p1, prev + p1), SMIN(curr, min_L[i * width + j + 1] + p2))
					- min_L[i * width + j + 1]; // NOTE: min_L subtractiong to avoid overflow by adding-accumulation

				L_rl[disparity_len * (i * width + j) + d_index] = cost_value;

				min_L[i * width + j] = SMIN(cost_value, min_L[i * width + j]);
				prev = curr;
				curr = next;			
			}
			
			ushort cost_value = maching_cost[disparity_len * (i * width + j) + disparity_len - 1]
				+ SMIN(prev + p1, SMIN(curr, min_L[i * width + j + 1] + p2))
				- min_L[i * width + j + 1];
						
			L_rl[disparity_len * (i * width + j) + disparity_len - 1] = cost_value;
			min_L[i * width + j] = SMIN(cost_value, min_L[i * width + j]);
		}

	/***************/
	/* Up->Down    */
	/***************/
	for(j = 0; j < width; j++) //width x disparity_len
		for(d_index = 0; d_index < disparity_len; d_index++)
		{
			L_ud[disparity_len * j + d_index]  = maching_cost[disparity_len * j + d_index];
			min_L[j] = USHRT_MAX; //SMIN(maching_cost[disparity_len * j + d_index], USHRT_MAX);
		}

	for (i = 1; i < height; i++) //(height - 1) x width
		for(j = 0; j < width; j++)
		{
			min_L[i * width + j] = USHRT_MAX;
			prev = USHRT_MAX;
			curr = L_ud[disparity_len * ((i - 1) * width + j)];

			for(d_index = 0; d_index < disparity_len - 1; d_index++)
			{
				next = L_ud[disparity_len * ((i - 1) * width + j) + d_index + 1];

				ushort cost_value = maching_cost[disparity_len * (i * width + j) + d_index]
					+ SMIN(SMIN(next + p1, prev + p1), SMIN(curr, min_L[(i - 1) * width + j] + p2))
					- min_L[(i - 1) * width + j]; // NOTE: min_L subtractiong to avoid overflow by adding-accumulation

				L_ud[disparity_len * (i * width + j) + d_index] = cost_value;

				min_L[i * width + j] = SMIN(cost_value, min_L[i * width + j]);
				prev = curr;
				curr = next;			
			}

			ushort cost_value = maching_cost[disparity_len * (i * width + j) + disparity_len - 1]
				+ SMIN(prev + p1, SMIN(curr, min_L[(i - 1) * width + j] + p2))
				- min_L[(i - 1) * width + j];

			L_ud[disparity_len * (i * width + j) + disparity_len - 1] = cost_value;
			min_L[i * width + j] = SMIN(cost_value, min_L[i * width + j]);
		}

	/***************/
	/* Down->Up    */
	/***************/
	for(j = 0; j < width; j++) //width x disparity_len
		for(d_index = 0; d_index < disparity_len; d_index++)
		{
			L_du[disparity_len * ((height - 1) * width + j) + d_index] = maching_cost[disparity_len * ((height - 1) * width + j) + d_index];
			min_L[(height - 1) * width + j] = USHRT_MAX; //SMIN(maching_cost[disparity_len * ((height - 1) * width + j) + d_index], USHRT_MAX);
		}

	for (i = height - 2; i >= 0; i--) //(height - 1) x width
		for(j = width - 1; j >= 0; j--)
		{
			min_L[i * width + j] = USHRT_MAX;
			prev = USHRT_MAX;
			curr = L_du[disparity_len * ((i + 1) * width + j)];

			for(d_index = 0; d_index < disparity_len - 1; d_index++)
			{
				next = L_du[disparity_len * ((i + 1) * width + j) + d_index + 1];

				ushort cost_value = maching_cost[disparity_len * (i * width + j) + d_index]
					+ SMIN(SMIN(next + p1, prev + p1), SMIN(curr, min_L[(i + 1) * width + j] + p2) )
					- min_L[(i+1)*width+j]; // NOTE: min_L subtractiong to avoid overflow by adding-accumulation

				L_du[disparity_len * (i * width + j) + d_index] = cost_value;

				min_L[i * width + j] = SMIN(cost_value, min_L[i * width + j]);
				prev = curr;
				curr = next;			
			}

			ushort cost_value = maching_cost[disparity_len * (i * width + j) + d_index]
				+ SMIN(prev + p1, SMIN(curr, min_L[(i + 1) * width + j] + p2))
				- min_L[(i + 1) * width + j];

			L_du[disparity_len * (i * width + j) + d_index] = cost_value;
			min_L[i * width + j] = SMIN(cost_value, min_L[i * width + j]);
		}

	/**************/
	/* Compute S  */
	/**************/
	for (i = 0; i < height; i++) //height x width
		for(j = 0; j < width; j++)
			for(d_index = 0; d_index < disparity_len; d_index++)
				S[disparity_len * (i * width + j) + d_index] = L_lr[disparity_len * (i * width + j) + d_index] + L_rl[disparity_len * (i * width + j) + d_index]
					+ L_ud[disparity_len * (i * width + j) + d_index] + L_du[disparity_len * (i * width + j) + d_index];
}

void writeCost(uchar *imOUT, uchar *maching_cost, int height, int width, int disp_min, int disp_max)
{
	int i, j, d_index;
	int disparity_len  = disp_max - disp_min + 1;

	// 9x7 Census transform
	int top = (CENSUS_HEIGHT - 1) / 2;
	int left = (CENSUS_WIDTH - 1) / 2;

	for(d_index = 0; d_index < disparity_len; d_index++)
	{
		for(i = 0; i < height; i++)
			for(j = 0; j < width; j++) 
				imOUT[i * width + j] = maching_cost[disparity_len * (i * width + j) + d_index];

		char file[100];
		sprintf(file,"costs/image%i.bin", d_index);
		write_image_bin(file, imOUT, height, width);
	}
}

void minDisparity(uchar *min, ushort *cost_aggretation, int height, int width, int disp_min, int disp_max)
{
	int i, j, d_index;
	int disparity_len  = disp_max - disp_min + 1;

	for(i = 0; i < height; i++)
		for(j = 0; j < width; j++)
		{
			ushort min_MC = USHRT_MAX;
			uchar d_min;
		
			for(d_index = 0; d_index < disparity_len; d_index++)
			{
				if (j - d_index >= 0 && j + d_index < width)
					if (min_MC > cost_aggretation[disparity_len * (i * width + j) + d_index])
					{
						min_MC = cost_aggretation[disparity_len * (i * width + j) + d_index];
						d_min = d_index;
					}
			}

			min[i * width + j] = d_min;
		}
}

void medianFilter(uchar *out, uchar *in, int height, int width)
{
	int i, ii, j, jj;

	for(i = 0; i < height; i++)
		for(j = 0; j < width; j++)
		{
			ushort median = 0;
		
			for (ii = -1; ii <= 1; ii++)
				for (jj = -1; jj <= 1; jj++)
					if (i + ii >= 0 && i + ii < height && j + jj >= 0 && j + jj < width)
						median += in[(i + ii) * width + (j + jj)];

			out[i * width + j] = 4 * median / 9;
		}
}

void executeOCL(uchar *im_pad1, uchar *im_pad2, uchar *Dm, int height, int width,
	int vert_pad, int hor_pad, int disp_min, int disp_max)
{
	//Variables used to read kernel source file
	FILE *fp = fopen("kernels.cl", "r");
	size_t filelen, readlen;
	char *kernel_src;

	//Read the kernel
	fseek(fp, 0L, SEEK_END);

	filelen = ftell(fp);

	rewind(fp);

	kernel_src = (char *)malloc(sizeof(char) * (filelen + 1));
	readlen = fread(kernel_src, 1, filelen, fp);
	kernel_src[filelen] = '\0';

	if (filelen != readlen)
	{
		printf("Failed to read kernel file\n");
		exit(EXIT_FAILURE);
	}

	//Set up platform and GPU device
	cl_uint numPlatforms;

	//Find number of platforms
	clGetPlatformIDs(0, NULL, &numPlatforms);

	//Get all platforms
	cl_platform_id *platform = (cl_platform_id *)malloc(numPlatforms * sizeof(cl_platform_id));

	clGetPlatformIDs(numPlatforms, platform, NULL);

	//Secure a GPU
	cl_device_id device_id = NULL;
	cl_int err = 1;

	int i = 0;

	while (i < (int)numPlatforms && err != CL_SUCCESS)
	{
		err = clGetDeviceIDs(platform[i], CL_DEVICE_TYPE_GPU, 1, &device_id, NULL);

		i++;
	}

	cl_context context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
	cl_command_queue command_queue = clCreateCommandQueue(context, device_id, 0, &err);
	cl_program program = clCreateProgramWithSource(context, 1, (const char **)&kernel_src, NULL, &err);

	if (err != CL_SUCCESS)
	{
		printf("Unable to create program object. Error Code=%d\n", err);
		exit(EXIT_FAILURE);
	}

	clBuildProgram(program, 0, NULL, NULL, NULL, NULL);

	char buffer[2048];
	size_t len;

	clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
	printf("--- Build Log -- \n %s\n", buffer);

	cl_kernel censusTransformImagesKernel = clCreateKernel(program, "censusTransformImages", &err),
		censusTransformMachingCostKernel = clCreateKernel(program, "censusTransformMachingCost", &err),
		sgmHorizontalInitKernel = clCreateKernel(program, "initSGMHorizontal", &err),
		sgmVerticalInitKernel = clCreateKernel(program, "initSGMVertical", &err),
		sgmHorizontalKernel = clCreateKernel(program, "sgmHorizontal", &err),
		sgmVerticalKernel = clCreateKernel(program, "sgmVertical", &err),
		computeSKernel = clCreateKernel(program, "computeS", &err),
		minDisparityKernel = clCreateKernel(program, "minDisparity", &err),
		medianFilterKernel = clCreateKernel(program, "medianFilter", &err);

	if (err != CL_SUCCESS)
	{
		printf("Unable to create kernel object. Error Code=%d\n", err);
		exit(EXIT_FAILURE);
	}

	int widthPAD = width + hor_pad;
	int disp_len = disp_max - disp_min + 1;

	cl_mem census_im1 = clCreateBuffer(context, CL_MEM_READ_WRITE, height * width * sizeof(lli), NULL, &err),
		census_im2 = clCreateBuffer(context, CL_MEM_READ_WRITE, height * width * sizeof(lli), NULL, &err),
		im1pad = clCreateBuffer(context, CL_MEM_READ_ONLY, (height + vert_pad) * widthPAD * sizeof(uchar), NULL, &err),
		im2pad = clCreateBuffer(context, CL_MEM_READ_ONLY, (height + vert_pad) * widthPAD * sizeof(uchar), NULL, &err),
		devMachingCost = clCreateBuffer(context, CL_MEM_READ_WRITE, height * width * disp_len * sizeof(uchar), NULL, &err),
		devL = clCreateBuffer(context, CL_MEM_READ_WRITE, 4 * height * width * disp_len * sizeof(ushort), NULL, &err),
		devMinL = clCreateBuffer(context, CL_MEM_READ_WRITE, height * width * sizeof(ushort), NULL, &err),
		devS = clCreateBuffer(context, CL_MEM_READ_WRITE, height * width * disp_len * sizeof(ushort), NULL, &err),
		devD = clCreateBuffer(context, CL_MEM_READ_WRITE, height * width * sizeof(uchar), NULL, &err),
		devDm = clCreateBuffer(context, CL_MEM_READ_WRITE, height * width * sizeof(uchar), NULL, &err);

	err = clEnqueueWriteBuffer(command_queue, im1pad, CL_TRUE, 0, (height + vert_pad) * widthPAD * sizeof(uchar), im_pad1, 0, NULL, NULL);
	err = clEnqueueWriteBuffer(command_queue, im2pad, CL_TRUE, 0, (height + vert_pad) * widthPAD * sizeof(uchar), im_pad2, 0, NULL, NULL);

	if (clSetKernelArg(censusTransformImagesKernel, 0, sizeof(cl_mem), &census_im1) ||
		clSetKernelArg(censusTransformImagesKernel, 1, sizeof(cl_mem), &census_im2) ||
		clSetKernelArg(censusTransformImagesKernel, 2, sizeof(cl_mem), &im1pad) ||
		clSetKernelArg(censusTransformImagesKernel, 3, sizeof(cl_mem), &im2pad) || 
		clSetKernelArg(censusTransformImagesKernel, 4, sizeof(cl_int), &width) ||
		clSetKernelArg(censusTransformImagesKernel, 5, sizeof(cl_int), &hor_pad) ||
		clSetKernelArg(censusTransformImagesKernel, 6, sizeof(cl_int), &vert_pad) != CL_SUCCESS)
	{
		printf("Unable to set kernel arguments. Error Code=%d\n", err);
		exit(EXIT_FAILURE);
	}

	if (clSetKernelArg(censusTransformMachingCostKernel, 0, sizeof(cl_mem), &devMachingCost) ||
		clSetKernelArg(censusTransformMachingCostKernel, 1, sizeof(cl_mem), &census_im1) ||
		clSetKernelArg(censusTransformMachingCostKernel, 2, sizeof(cl_mem), &census_im2) ||
		clSetKernelArg(censusTransformMachingCostKernel, 3, sizeof(cl_int), &disp_min) ||
		clSetKernelArg(censusTransformMachingCostKernel, 4, sizeof(cl_int), &disp_max) ||
		clSetKernelArg(censusTransformMachingCostKernel, 5, sizeof(cl_int), &width) != CL_SUCCESS)
	{
		printf("Unable to set kernel arguments. Error Code=%d\n", err);
		exit(EXIT_FAILURE);
	}

	if (clSetKernelArg(sgmHorizontalInitKernel, 0, sizeof(cl_mem), &devL) ||
		clSetKernelArg(sgmHorizontalInitKernel, 1, sizeof(cl_mem), &devMachingCost) ||
		clSetKernelArg(sgmHorizontalInitKernel, 2, sizeof(cl_mem), &devMinL) ||
		clSetKernelArg(sgmHorizontalInitKernel, 3, sizeof(cl_int), &height) ||
		clSetKernelArg(sgmHorizontalInitKernel, 4, sizeof(cl_int), &width) ||
		clSetKernelArg(sgmHorizontalInitKernel, 5, sizeof(cl_int), &disp_min) ||
		clSetKernelArg(sgmHorizontalInitKernel, 6, sizeof(cl_int), &disp_max) != CL_SUCCESS)
	{
		printf("Unable to set kernel arguments. Error Code=%d\n", err);
		exit(EXIT_FAILURE);
	}

	if (clSetKernelArg(sgmVerticalInitKernel, 0, sizeof(cl_mem), &devL) ||
		clSetKernelArg(sgmVerticalInitKernel, 1, sizeof(cl_mem), &devMachingCost) ||
		clSetKernelArg(sgmVerticalInitKernel, 2, sizeof(cl_mem), &devMinL) ||
		clSetKernelArg(sgmVerticalInitKernel, 3, sizeof(cl_int), &height) ||
		clSetKernelArg(sgmVerticalInitKernel, 4, sizeof(cl_int), &width) ||
		clSetKernelArg(sgmVerticalInitKernel, 5, sizeof(cl_int), &disp_min) ||
		clSetKernelArg(sgmVerticalInitKernel, 6, sizeof(cl_int), &disp_max) != CL_SUCCESS)
	{
		printf("Unable to set kernel arguments. Error Code=%d\n", err);
		exit(EXIT_FAILURE);
	}

	if (clSetKernelArg(sgmHorizontalKernel, 0, sizeof(cl_mem), &devMachingCost) ||
		clSetKernelArg(sgmHorizontalKernel, 1, sizeof(cl_mem), &devMinL) ||
		clSetKernelArg(sgmHorizontalKernel, 2, sizeof(cl_mem), &devL) ||
		clSetKernelArg(sgmHorizontalKernel, 3, sizeof(cl_int), &height) ||
		clSetKernelArg(sgmHorizontalKernel, 4, sizeof(cl_int), &width) ||
		clSetKernelArg(sgmHorizontalKernel, 5, sizeof(cl_int), &disp_min) ||
		clSetKernelArg(sgmHorizontalKernel, 6, sizeof(cl_int), &disp_max) != CL_SUCCESS)
	{
		printf("Unable to set kernel arguments. Error Code=%d\n", err);
		exit(EXIT_FAILURE);
	}

	if (clSetKernelArg(sgmVerticalKernel, 0, sizeof(cl_mem), &devMachingCost) ||
		clSetKernelArg(sgmVerticalKernel, 1, sizeof(cl_mem), &devMinL) ||
		clSetKernelArg(sgmVerticalKernel, 2, sizeof(cl_mem), &devL) ||
		clSetKernelArg(sgmVerticalKernel, 3, sizeof(cl_int), &height) ||
		clSetKernelArg(sgmVerticalKernel, 4, sizeof(cl_int), &width) ||
		clSetKernelArg(sgmVerticalKernel, 5, sizeof(cl_int), &disp_min) ||
		clSetKernelArg(sgmVerticalKernel, 6, sizeof(cl_int), &disp_max) != CL_SUCCESS)
	{
		printf("Unable to set kernel arguments. Error Code=%d\n", err);
		exit(EXIT_FAILURE);
	}

	if (clSetKernelArg(computeSKernel, 0, sizeof(cl_mem), &devS) ||
		clSetKernelArg(computeSKernel, 1, sizeof(cl_mem), &devL) ||
		clSetKernelArg(computeSKernel, 2, sizeof(cl_int), &height) ||
		clSetKernelArg(computeSKernel, 3, sizeof(cl_int), &width) ||
		clSetKernelArg(computeSKernel, 4, sizeof(cl_int), &disp_min) ||
		clSetKernelArg(computeSKernel, 5, sizeof(cl_int), &disp_max) != CL_SUCCESS)
	{
		printf("Unable to set kernel arguments. Error Code=%d\n", err);
		exit(EXIT_FAILURE);
	}

	if (clSetKernelArg(minDisparityKernel, 0, sizeof(cl_mem), &devD) ||
		clSetKernelArg(minDisparityKernel, 1, sizeof(cl_mem), &devS) ||
		clSetKernelArg(minDisparityKernel, 2, sizeof(cl_int), &width) ||
		clSetKernelArg(minDisparityKernel, 3, sizeof(cl_int), &disp_min) ||
		clSetKernelArg(minDisparityKernel, 4, sizeof(cl_int), &disp_max) != CL_SUCCESS)
	{
		printf("Unable to set kernel arguments. Error Code=%d\n", err);
		exit(EXIT_FAILURE);
	}

	if (clSetKernelArg(medianFilterKernel, 0, sizeof(cl_mem), &devDm) ||
		clSetKernelArg(medianFilterKernel, 1, sizeof(cl_mem), &devD) ||
		clSetKernelArg(medianFilterKernel, 2, sizeof(cl_int), &height) ||
		clSetKernelArg(medianFilterKernel, 3, sizeof(cl_int), &width) != CL_SUCCESS)
	{
		printf("Unable to set kernel arguments. Error Code=%d\n", err);
		exit(EXIT_FAILURE);
	}

	int disparity_len = disp_max - disp_min + 1;
	size_t global[2] = { height, width };
	double t0 = getMicroSeconds();

	//Census Transform
	err = clEnqueueNDRangeKernel(command_queue, censusTransformImagesKernel, 2, NULL, global, NULL, 0, NULL, NULL);
	err = clEnqueueNDRangeKernel(command_queue, censusTransformMachingCostKernel, 2, NULL, global, NULL, 0, NULL, NULL);

	//SGM
	size_t dim = height;

	err = clEnqueueNDRangeKernel(command_queue, sgmHorizontalKernel, 1, NULL, &dim, NULL, 0, NULL, NULL);
	
	dim = width;
	err = clEnqueueNDRangeKernel(command_queue, sgmVerticalKernel, 1, NULL, &dim, NULL, 0, NULL, NULL);

	global[0] = height;
	global[1] = width;

	err = clEnqueueNDRangeKernel(command_queue, computeSKernel, 2, NULL, global, NULL, 0, NULL, NULL);

	//Min Disparity
	err = clEnqueueNDRangeKernel(command_queue, minDisparityKernel, 2, NULL, global, NULL, 0, NULL, NULL);

	//Median Filter
	err = clEnqueueNDRangeKernel(command_queue, medianFilterKernel, 2, NULL, global, NULL, 0, NULL, NULL);

	//Copy result back to host
	err = clEnqueueReadBuffer(command_queue, devDm, CL_TRUE, 0, height *width * sizeof(uchar), Dm, 0, NULL, NULL);
	
	double t1 = getMicroSeconds();
    printf("\nThe kernels ran in %lf seconds\n",(t1-t0)/1000000);

	clReleaseProgram(program);
	clReleaseKernel(censusTransformImagesKernel);
	clReleaseKernel(censusTransformMachingCostKernel);
	clReleaseKernel(sgmHorizontalKernel);
	clReleaseKernel(sgmVerticalKernel);
	clReleaseKernel(computeSKernel);
	clReleaseKernel(minDisparityKernel);
	clReleaseKernel(medianFilterKernel);
	clReleaseCommandQueue(command_queue);
	clReleaseContext(context);
	free(kernel_src);
}

void main(int argc, char **argv)
{
	uchar *im1, *im2, *im_pad1, *im_pad2;
	uchar *imtmp1, *imtmp2;

	/* Only accept a concrete number of arguments */
	if(argc != 6)
	{
		printf("./exec disp_min disp_max image1.bmp image2.bmp [c, g]\n");
		exit(-1);
	}

	int disp_min = atoi(argv[1]);
	int disp_max = atoi(argv[2]);
	
	if (disp_max < disp_min || disp_min < 0 || disp_max >= UCHAR_MAX)
	{
		printf("Disparity min and max must be positive and less than %i\n", UCHAR_MAX); 
		exit(-1);
	}

	char device = argv[5][0];
	int disp_len = disp_max - disp_min + 1;
	int height;
	int width;

	/* Read images */
	imtmp1 = read_image_BMP(argv[3], &width, &height);
	im1 = image_RGB2BW(imtmp1, height, width);

	imtmp2 = read_image_BMP(argv[4], &width, &height);
	im2 = image_RGB2BW(imtmp2, height, width);

	im_pad1 = padding_image(im1, height, width, CENSUS_HEIGHT, CENSUS_WIDTH);
	im_pad2 = padding_image(im2, height, width, CENSUS_HEIGHT, CENSUS_WIDTH);

	uchar *maching_cost = (uchar *)malloc(height * width * disp_len * sizeof(uchar));
	lli *census_im1     = (lli *)malloc(height * width * sizeof(lli));
	lli *census_im2     = (lli *)malloc(height * width * sizeof(lli));
	ushort *L			= (ushort *)malloc(4 * height * width * disp_len * sizeof(ushort));
	ushort *S			= (ushort *)malloc(height * width * disp_len * sizeof(ushort));
	ushort *min_L		= (ushort *)malloc(height * width * sizeof(ushort));
	uchar *D			= (uchar *)malloc(height *width * sizeof(uchar));
	uchar *Dm			= (uchar *)malloc(height *width * sizeof(uchar));

	switch (device)
	{
	case 'c':
	{

		/****************************/
		/* Compute Census Transform */
		/****************************/
		double t0 = getMicroSeconds();

		censusTransform(maching_cost, census_im1, census_im2, im_pad1, im_pad2, height, width, disp_min, disp_max);

		/****************************/
		/* SGM: cost aggregation    */
		/****************************/
		SGM(maching_cost, S, L, min_L, height, width, disp_min, disp_max);

		/****************************/
		/* min disparity &          */
		/* median Blur 3x3          */
		/****************************/
		minDisparity(D, S, height, width, disp_min, disp_max);
		medianFilter(Dm, D, height, width);
		
		double t1 = getMicroSeconds();
		printf("\nThe kernel ran in %lf seconds\n",(t1-t0)/1000000);

		/****************************/
		/* salida                   */
		/****************************/
		write_image_BMP(Dm, "imageOUT.bmp", width, height);

		break;
	}
	case 'g':
	{
		executeOCL(im_pad1, im_pad2, Dm, height, width, CENSUS_HEIGHT, CENSUS_WIDTH, disp_min, disp_max);
		write_image_BMP(Dm, "imageOUT.bmp", width, height);

		break;
	}
	default:
		printf("Incorrect device, use c or g\n");
	}

	/* Free */
	free(im1); 
	free(imtmp1);
	free(im2);
	free(imtmp2);
	free(im_pad1);
	free(im_pad2);
	free(maching_cost);
	free(census_im1);
	free(census_im2);
	free(L);
	free(S);
	free(min_L);
	free(D);
	free(Dm);
}
