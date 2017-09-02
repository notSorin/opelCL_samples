typedef long long int lli;

#define SMIN(a, b) ((a < b) ? a : b)

__kernel void medianFilter(__global uchar *out, __global uchar *in, int height, int width)
{
	int ii = 0, jj = 0, i = get_global_id(0), j = get_global_id(1);
	ushort median = 0;
		
	for (ii = -1; ii <= 1; ii++)
		for (jj = -1; jj <= 1; jj++)
			if (i + ii >= 0 && i + ii < height && j + jj >= 0 && j + jj < width)
				median += in[(i + ii) * width + (j + jj)];

	out[i * width + j] = 4 * median / 9;
}

__kernel void minDisparity(__global uchar *min, __global ushort *cost_aggretation, int width, int disp_min, int disp_max)
{
	int d_index = 0, disparity_len  = disp_max - disp_min + 1, i = get_global_id(0), j = get_global_id(1);
	ushort min_MC = USHRT_MAX;
	uchar d_min;

	for(d_index = 0; d_index < disparity_len; d_index++)
		if (j - d_index >= 0 && j + d_index < width)
			if (min_MC > cost_aggretation[disparity_len * (i * width + j) + d_index])
			{
				min_MC = cost_aggretation[disparity_len * (i * width + j) + d_index];
				d_min = d_index;
			}

	min[i * width + j] = d_min;
}

__kernel void initSGMHorizontal(__global ushort *L, __global uchar *maching_cost, __global ushort *min_L, int height, int width, int disp_min, int disp_max)
{
	int disparity_len  = disp_max - disp_min + 1, i = get_global_id(0), d_index = get_global_id(1);
	__global ushort *L_lr = &(L[0 * (disparity_len * height * width)]); /* Left->Right */
	__global ushort *L_rl = &(L[1 * (disparity_len * height * width)]); /* Right->Left */

	L_lr[disparity_len * i * width + d_index] = maching_cost[disparity_len * i * width + d_index];
	min_L[i * width] = maching_cost[disparity_len * i * width + d_index];
	
	L_rl[disparity_len * (i * width + width - 1) + d_index]  = maching_cost[disparity_len * (i * width + width - 1) + d_index];
	min_L[i * width + width - 1] = maching_cost[disparity_len * (i * width + width - 1) + d_index];
}

__kernel void initSGMVertical(__global ushort *L, __global uchar *maching_cost, __global ushort *min_L, int height, int width, int disp_min, int disp_max)
{
	int disparity_len  = disp_max - disp_min + 1, j = get_global_id(0), d_index = get_global_id(1);	
	__global ushort *L_ud = &(L[2 * (disparity_len * height * width)]); /* Up->Down */
	__global ushort *L_du = &(L[3 * (disparity_len * height * width)]); /* Down->Up */
	
	L_ud[disparity_len * j + d_index]  = maching_cost[disparity_len * j + d_index];
	min_L[j] = maching_cost[disparity_len * j + d_index];

	L_du[disparity_len * ((height - 1) * width + j) + d_index] = maching_cost[disparity_len * ((height - 1) * width + j) + d_index];
	min_L[(height - 1) * width + j] = maching_cost[disparity_len * ((height - 1) * width + j) + d_index];
}

__kernel void sgmHorizontal(__global uchar *maching_cost, __global ushort *min_L, __global ushort *L, int height, int width, int disp_min, int disp_max)
{
	int i = get_global_id(0), j = 0, d_index = 0, disparity_len  = disp_max - disp_min + 1;
	__global ushort *L_lr = &(L[0 * (disparity_len * height * width)]); /* Left->Right */
	__global ushort *L_rl = &(L[1 * (disparity_len * height * width)]); /* Right->Left */
	ushort prev, curr, next, p1 = 30, p2 = 30;
	
	//Left->Right
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
	
	//Right->Left
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
}

__kernel void sgmVertical(__global uchar *maching_cost, __global ushort *min_L, __global ushort *L, int height, int width, int disp_min, int disp_max)
{
	int i = 0, j = get_global_id(0), disparity_len  = disp_max - disp_min + 1, d_index = 0;
	__global ushort *L_ud = &(L[2 * (disparity_len * height * width)]); /* Up->Down */
	__global ushort *L_du = &(L[3 * (disparity_len * height * width)]); /* Down->Up */
	ushort prev, curr, next, p1 = 30, p2 = 30;
	
	//Up->Down
	for(i = 1; i < height; i++)
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
	
	//Down-Up
	for(i = height - 2; i >= 0; i--)
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
}

__kernel void computeS(__global ushort *S, __global ushort *L, int height, int width, int disp_min, int disp_max)
{
	int d_index = 0, disparity_len  = disp_max - disp_min + 1, i = get_global_id(0), j = get_global_id(1);
	__global ushort *L_lr = &(L[0 * (disparity_len * height * width)]); /* Left->Right */
	__global ushort *L_rl = &(L[1 * (disparity_len * height * width)]); /* Right->Left */
	__global ushort *L_ud = &(L[2 * (disparity_len * height * width)]); /* Up->Down */
	__global ushort *L_du = &(L[3 * (disparity_len * height * width)]); /* Down->Up */
	
	for(d_index = 0; d_index < disparity_len; d_index++)
		S[disparity_len * (i * width + j) + d_index] = L_lr[disparity_len * (i * width + j) + d_index] + L_rl[disparity_len * (i * width + j) + d_index]
			+ L_ud[disparity_len * (i * width + j) + d_index] + L_du[disparity_len * (i * width + j) + d_index];
}

__kernel void censusTransformImages(__global lli *census_im1, __global lli *census_im2, 
__global uchar *im1, __global uchar *im2, int width, const int CENSUS_WIDTH, const int CENSUS_HEIGHT)
{
	int top = (CENSUS_HEIGHT - 1) / 2, left = (CENSUS_WIDTH - 1) / 2,
		i = get_global_id(0) + top, j = get_global_id(1) + left, ii, jj, widthPAD = width + CENSUS_WIDTH;
	uchar center1, center2, curr1, curr2;
	unsigned long long int census1 = 0, census2 = 0;
	
	center1 = im1[i * widthPAD + j];
	center2 = im2[i * widthPAD + j];
	
	for(ii = -top; ii <= top; ii++)
		for(jj = -left; jj <= left; jj++)
		{
			curr1 = im1[(i + ii) * widthPAD + (j + jj)];
			unsigned long long int tmp = curr1 >= center1;
			tmp = tmp << ((ii + top) * CENSUS_WIDTH + (jj + left));
			census1 = census1 | tmp;

			curr2 = im2[(i + ii) * widthPAD + (j + jj)];
			tmp = curr2 >= center2;
			tmp = tmp << ((ii + top) * CENSUS_WIDTH + (jj + left));
			census2 = census2 | tmp;
		}
	
	census_im1[(i - top) * width + (j - left)] = census1;
	census_im2[(i - top) * width + (j - left)] = census2;
}

__kernel void censusTransformMachingCost(__global uchar *maching_cost,
__global lli *census_im1, __global lli *census_im2, int disp_min, int disp_max, int width)
{
	int i = get_global_id(0), j = get_global_id(1), d_index = 0, disparity_len = disp_max - disp_min + 1;
	
	for(d_index = 0; d_index < disparity_len; d_index++)
	{
		short d = d_index + disp_min;
	
		if (j - d >= 0 && j - d < width)
		{
			lli x = census_im1[i * width + j] ^ census_im2[i * width + j - d];
			
			//Count set bits in x
			x = x - ((x >> 1) & 0x55555555);
			x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
			x = (x + (x >> 4)) & 0x0F0F0F0F;
			x = x + (x >> 8);
			x = x + (x >> 16);
			
			maching_cost[disparity_len * (i * width + j) + d_index] = x & 0x0000003F;
		}
		else
			maching_cost[disparity_len * (i * width + j) + d_index] = UCHAR_MAX;
	}
}