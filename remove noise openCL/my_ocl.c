//Copyright (C) Sorin Draghici <dsorin95@gmail.com>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "my_ocl.h"
#include "CL/cl.h"

#ifndef DEVICE
#define DEVICE CL_DEVICE_TYPE_GPU
#endif

#include <sys/time.h>
#include <sys/resource.h>

double get_time_(){
	static struct timeval 	tv0;
	double time_, time;

	gettimeofday(&tv0,(struct timezone*)0);
	time_=(double)((tv0.tv_usec + (tv0.tv_sec)*1000000));
	time = time_/1000000;
	return(time);
}

int remove_noiseOCL(float *im, float *image_out, float thredshold, int window_size, int height, int width)
{
	// OpenCL host variables
	cl_device_id device_id = NULL;
	cl_int err;
	cl_context context;
	cl_command_queue command_queue;
	cl_program program;
	cl_kernel kernel;
	size_t global[2];

	// variables used to read kernel source file
	FILE *fp;
	size_t filelen, readlen;
	char *kernel_src;  // char string to hold kernel source
	
	// read the kernel
	fp = fopen("remove_noise_kernel.cl", "r");
	fseek(fp, 0L, SEEK_END);
	filelen = ftell(fp);
	rewind(fp);

	kernel_src = malloc(sizeof(char) * (filelen + 1));
	readlen = fread(kernel_src, 1, filelen, fp);
	
	if(readlen!= filelen)
	{
		printf("error reading file\n");
		exit(1);
	}
	
	// ensure the string is NULL terminated
	kernel_src[filelen]='\0';
	
	// Set up platform and GPU device

	cl_uint numPlatforms;

	// Find number of platforms
	err = clGetPlatformIDs(0, NULL, &numPlatforms);
	if (err != CL_SUCCESS || numPlatforms <= 0)
	{
		printf("Error: Failed to find a platform!\n\n");
		return EXIT_FAILURE;
	}

	// Get all platforms
	cl_platform_id *Platform = (cl_platform_id *)malloc(numPlatforms * sizeof(cl_platform_id));

	err = clGetPlatformIDs(numPlatforms, Platform, NULL);
	if (err != CL_SUCCESS || numPlatforms <= 0)
	{
		printf("Error: Failed to get the platform!\n\n");
		return EXIT_FAILURE;
	}
	
	// Secure a GPU
	int i;
	for (i = 0; i < (int)numPlatforms; i++)
	{
		err = clGetDeviceIDs(Platform[i], DEVICE, 1, &device_id, NULL);
		if (err == CL_SUCCESS)
		{
			break;
		}
	}

	if (device_id == NULL)
	{
		printf("Error: Failed to create a device group!\n\n");
		return EXIT_FAILURE;
	}

	// Create a compute context 
	context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
	if (!context)
	{
		printf("Error: Failed to create a compute context!\n\n");
		return EXIT_FAILURE;
	}

	// Create a command queue
	cl_command_queue commands = clCreateCommandQueue(context, device_id, 0, &err);
	if (!commands)
	{
		printf("Error: Failed to create a command commands!\n\n");
		return EXIT_FAILURE;
	}

	// create command queue 
	command_queue = clCreateCommandQueue(context,device_id, 0, &err);
	if (err != CL_SUCCESS)
	{	
		printf("Unable to create command queue. Error Code=%d\n",err);
		exit(1);
	}
	
	// create program object from source. 
	// kernel_src contains source read from file earlier
	program = clCreateProgramWithSource(context, 1 ,(const char **)&kernel_src, NULL, &err);
	if (err != CL_SUCCESS)
	{	
		printf("Unable to create program object. Error Code=%d\n",err);
		exit(1);
	}
	
	err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
	if (err != CL_SUCCESS)
	{
        	printf("Build failed. Error Code=%d\n", err);

		size_t len;
		char buffer[2048];
		// get the build log
		clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buffer), buffer, &len);
		printf("--- Build Log -- \n %s\n",buffer);
		exit(1);
	}
	
	kernel = clCreateKernel(program, "remove_noise", &err);
	if (err != CL_SUCCESS)
	{	
		printf("Unable to create kernel object. Error Code=%d\n",err);
		exit(1);
	}
	
	cl_mem dIm, dImOut;
	// create buffer objects to input and output args of kernel function
	dIm = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float) * width * height, NULL, NULL);
	dImOut = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float) * width * height, NULL, NULL);

	err = clEnqueueWriteBuffer(commands, dIm, CL_TRUE, 0, sizeof(float) * width * height, im, 0, NULL, NULL);
    if (err != CL_SUCCESS)
    {
	printf("Error: Failed to write im to dIm array!\n\n");
	exit(1);
    }
	
	
	// set the kernel arguments
	if ( clSetKernelArg(kernel, 0, sizeof(cl_mem), &dIm) ||
         clSetKernelArg(kernel, 1, sizeof(cl_mem), &dImOut) ||
         clSetKernelArg(kernel, 2, sizeof(cl_float), &thredshold) ||
	 	 clSetKernelArg(kernel, 3, sizeof(cl_int), &window_size) ||
	 	 clSetKernelArg(kernel, 4, sizeof(cl_int), &height) ||
	 	 clSetKernelArg(kernel, 5, sizeof(cl_int), &width) != CL_SUCCESS)
	{
		printf("Unable to set kernel arguments. Error Code=%d\n",err);
		exit(1);
	}

	// set the global work dimension size
	global[0]= height;
	global[1]= width;

	//Variables para calcular el tiempo
	double t0, t1;
	double cpu_time_used = 0.0;
	
	t0 = get_time_();
	err = clEnqueueNDRangeKernel(command_queue, kernel, 2, NULL, global, NULL, 0, NULL, NULL);
	
	// wait for the command to finish
	clFinish(command_queue);
	
	t1 = get_time_();
	
	printf("OCL Exection time %f s.\n", t1-t0);

	if (err != CL_SUCCESS)
	{	
		printf("Unable to enqueue kernel command. Error Code=%d\n",err);
		exit(1);
	}

	// read the output back to host memory
	err = clEnqueueReadBuffer(commands, dImOut, CL_TRUE, 0, sizeof(float) * width * height, image_out, 0, NULL, NULL);
	if (err != CL_SUCCESS)
	{
		printf("Error enqueuing read buffer command. Error Code=%d\n",err);
		exit(1);
	}
	
	// clean up
	clReleaseProgram(program);
	clReleaseKernel(kernel);
	clReleaseCommandQueue(command_queue);
	clReleaseContext(context);
	free(kernel_src);

	return 0;
}
