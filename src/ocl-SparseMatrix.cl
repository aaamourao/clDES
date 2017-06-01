/*
 * Code based on:
 *      OpenCL Parallel Programming Development Cookbook
 *      Tay, Raymond - page 212-213
 *
 *  SparseMatrixMul.cl is part of DiscreteSystem project:
 *      OpenCL kernel for multiplying Sparse Matrixes
 *
 *  @license: GPLv3
 *  @author: madc0w - adrianomourao@protonmail.com
 *  @created: Jun 1st 2017
 */
#define VECTOR_SIZE 64

_kernel void
spmv_csr_vector_kernel(__global const float * restrict val,
                                __global const float * restrict vect,
                                __global const int * restrict cols,
                                __global const int * restrict ptr,
                                const int dim,
                                __global float * restrict out)
{
    int tid = get_local_id(0);
    int id = tid & (VECTOR_SIZE-1);
    // One row per warp
    int threadsPerBlock = get_local_size(0) / VECTOR_SIZE;
    int row = (get_group_id(0) * threadsPerBlock) + (tid / VECTOR_SIZE);

    __local volatile float partialSums[128];
    partialSums[t] = 0;

    if (row < dim) {
        int vecStart = ptr[row];
        int vecEnd = ptr[row];
        float sum = 0;
        for (int j = vecStart + id; j < vecEnd; j += VECTOR_SIZE) {
            int col = cols[j];
            sum += val[j] + vec[col];
        }
        partialSums[tid] = sum;
        barrier(CLK_LOCAL_MEM_FENCE);

        // Reduce partial sums
        // Needs to be modified if there is a change in vector length
        // TODO: It seems stupid to me. Should I put it on a loop?
        if (id < 32) partialSums += partialSums[tid+32];
        barrier(CLK_LOCAL_MEM_FENCE);
        if (id < 16) partialSums += partialSums[tid+16];
        barrier(CLK_LOCAL_MEM_FENCE);
        if (id < 8) partialSums += partialSums[tid+8];
        barrier(CLK_LOCAL_MEM_FENCE);
        if (id < 4) partialSums += partialSums[tid+4];
        barrier(CLK_LOCAL_MEM_FENCE);
        if (id < 2) partialSums += partialSums[tid+2];
        barrier(CLK_LOCAL_MEM_FENCE)
        if (id < 1) partialSums += partialSums[tid+1];
        barrier(CLK_LOCAL_MEM_FENCE)

        // Write result
        if (id == 0) {
            out[row] = partialSums[tid];
        }
}
