/*
 =========================================================================
 This file is part of clDES

 clDES: an OpenCL library for Discrete Event Systems computing.

 clDES is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 clDES is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with clDES.  If not, see <http://www.gnu.org/licenses/>.

 Copyright (c) 2018 - Adriano Mourao <adrianomourao@protonmail.com>
                      madc0ww @ [https://github.com/madc0ww]

 LacSED - Laboratório de Sistemas a Eventos Discretos
 Universidade Federal de Minas Gerais

 File: include/backend/kernels.cl
 Description: OpenCL kernels for the custom clDES operations.
 =========================================================================
*/

typedef struct StatesTuple {
    unsigned int x0;
    unsigned int x1;
} StatesTuple;

float gcd(float a, float b) {
    if (a == 0 || b == 0) {
        return 1.0f;
    }
    float mod;
    while (b != 0) {
        mod = fmod(a, b);
        a = b;
        b = mod;
    }
    return a;
}

__kernel void Synchronize_Stage1(__global StatesTuple *aSyncTuples,
                                 unsigned int aG0NumberStates) {
    // Workitems gets its index within index space
    int ix0 = get_global_id(0);
    int ix1 = get_global_id(1);

    unsigned int index = ix1 * aG0NumberStates + ix0;

    aSyncTuples[index].x0 = ix0;
    aSyncTuples[index].x1 = ix1;
};

__kernel void Synchronize_Stage2(__global StatesTuple *aTable,
                                 __global const unsigned int *aG0RowIndices,
                                 __global const unsigned int *aG0ColIndices,
                                 __global const float *aG0Elements,
                                 unsigned int aG0Size, float aG0Private,
                                 __global const unsigned int *aG1RowIndices,
                                 __global const unsigned int *aG1ColIndices,
                                 __global const float *aG1Elements,
                                 float aG1Private, __global float *aSync,
                                 unsigned int aPad) {
    int row = get_global_id(0);
    StatesTuple state = aTable[row];

    unsigned int g0_row_start = aG0RowIndices[state.x0];
    unsigned int g0_row_stop = aG0RowIndices[state.x0 + 1];

    if (get_global_id(1) < g0_row_stop - g0_row_start){
        unsigned int i = g0_row_start + get_global_id(1);

        unsigned int g1_row_start = aG1RowIndices[state.x1];
        unsigned int g1_row_stop = aG1RowIndices[state.x1 + 1];

        float g0_elem = aG0Elements[i];
        if (aG0Private > 1.0f) {
            float g0_gcd_priv = gcd(aG0Private, g0_elem);
            if (g0_gcd_priv > 1.0f) {
                if (aSync[(state.x1 * aG0Size + aG0ColIndices[i]) * aPad +
                          row] > 1.0f) {
                    aSync[(state.x1 * aG0Size + aG0ColIndices[i]) * aPad +
                          row] *= g0_gcd_priv;
                } else {
                    aSync[(state.x1 * aG0Size + aG0ColIndices[i]) * aPad +
                          row] = g0_gcd_priv;
                }
                aG0Private = aG0Private / g0_gcd_priv;
                g0_elem = g0_elem / g0_gcd_priv;
            }
        }
        if (g0_elem > 1.0f) {
            for (unsigned int j = g1_row_start; j < g1_row_stop; ++j) {
                float g1_elem = aG1Elements[j];
                if (aG1Private > 1.0f) {
                    float g1_gcd_priv = gcd(aG1Private, g1_elem);
                    if (g1_gcd_priv > 1.0f) {
                        if (aSync[(aG1ColIndices[j] * aG0Size + state.x0) *
                                      aPad +
                                  row] > 1.0f) {
                            aSync[(aG1ColIndices[j] * aG0Size + state.x0) *
                                      aPad +
                                  row] *= g1_gcd_priv;
                        } else {
                            aSync[(aG1ColIndices[j] * aG0Size + state.x0) *
                                      aPad +
                                  row] = g1_gcd_priv;
                        }
                        aG1Private = aG1Private / g1_gcd_priv;
                        g1_elem = g1_elem / g1_gcd_priv;
                    }
                }
                float sync_gcd = gcd(g0_elem, g1_elem);
                if (sync_gcd > 1.0f) {
                    if (aSync[(aG1ColIndices[j] * aG0Size + aG0ColIndices[i]) *
                                  aPad +
                              row] > 1.0f) {
                        aSync[(aG1ColIndices[j] * aG0Size + aG0ColIndices[i]) *
                                  aPad +
                              row] *= sync_gcd;
                    } else {
                        aSync[(aG1ColIndices[j] * aG0Size + aG0ColIndices[i]) *
                                  aPad +
                              row] = sync_gcd;
                    }
                }
            }
        }
    }
/*
    for (unsigned int i = g0_row_start; i < g0_row_stop; ++i) {
        float g0_elem = aG0Elements[i];
        if (aG0Private > 1.0f) {
            float g0_gcd_priv = gcd(aG0Private, g0_elem);
            if (g0_gcd_priv > 1.0f) {
                if (aSync[(state.x1 * aG0Size + aG0ColIndices[i]) * aPad +
                          row] > 1.0f) {
                    aSync[(state.x1 * aG0Size + aG0ColIndices[i]) * aPad +
                          row] *= g0_gcd_priv;
                } else {
                    aSync[(state.x1 * aG0Size + aG0ColIndices[i]) * aPad +
                          row] = g0_gcd_priv;
                }
                aG0Private = aG0Private / g0_gcd_priv;
                g0_elem = g0_elem / g0_gcd_priv;
            }
        }
        if (g0_elem > 1.0f) {
            for (unsigned int j = g1_row_start; j < g1_row_stop; ++j) {
                float g1_elem = aG1Elements[j];
                if (aG1Private > 1.0f) {
                    float g1_gcd_priv = gcd(aG1Private, g1_elem);
                    if (g1_gcd_priv > 1.0f) {
                        if (aSync[(aG1ColIndices[j] * aG0Size + state.x0) *
                                      aPad +
                                  row] > 1.0f) {
                            aSync[(aG1ColIndices[j] * aG0Size + state.x0) *
                                      aPad +
                                  row] *= g1_gcd_priv;
                        } else {
                            aSync[(aG1ColIndices[j] * aG0Size + state.x0) *
                                      aPad +
                                  row] = g1_gcd_priv;
                        }
                        aG1Private = aG1Private / g1_gcd_priv;
                        g1_elem = g1_elem / g1_gcd_priv;
                    }
                }
                float sync_gcd = gcd(g0_elem, g1_elem);
                if (sync_gcd > 1.0f) {
                    if (aSync[(aG1ColIndices[j] * aG0Size + aG0ColIndices[i]) *
                                  aPad +
                              row] > 1.0f) {
                        aSync[(aG1ColIndices[j] * aG0Size + aG0ColIndices[i]) *
                                  aPad +
                              row] *= sync_gcd;
                    } else {
                        aSync[(aG1ColIndices[j] * aG0Size + aG0ColIndices[i]) *
                                  aPad +
                              row] = sync_gcd;
                    }
                }
            }
        }
    }
    */
}

/*
__kernel void AddAccessedStates(__global const unsigned int *aRowIndices,
                                __global const unsigned int *aColIndices,
                                volatile __global unsigned int *aVector,
                                volatile __global unsigned int *aVSize,
                                volatile __global unsigned int *aNAccessed,
                                unsigned int aNSearches) {
    unsigned int search_id = aColIndices[get_global_id(0)];

    unsigned int row = 0;
    while (aRowIndices[row] < get_global_id(0) + 1) {
        ++row;
    }
    unsigned int accessed_state = row - 1;

    if (aVector[search_id * aNSearches + accessed_state] == 0) {
        aVector[search_id * aNSearches + accessed_state] = 1;
        atomic_inc(aNAccessed);
    }
}

__kernel void FilterAccessedStates() {
    unsigned int instant_vec_size = aVSize[search_id];
    bool already_inserted = false;
    for (unsigned int i = 0; i < instant_vec_size; ++i) {
        if (accessed_state == aVector[search_id * aNSearches + i]) {
            already_inserted = true;
            break;
        }
    }

    if (already_inserted == false) {
        atomic_inc(aNAccessed);
        aVector[search_id * aNSearches + atomic_inc(&aVSize[search_id])] =
            accessed_state;
    }
}
*/
/*
__kernel void SetIntersection(__global const unsigned int *aG0RowIndices,
                              __global const unsigned int *aG0ColIndices,
                              __global const float *aG0Elements,
                              __global unsigned int *vector,
                              __global unsigned int *vsize) {
    if (*vsize == 0) {
        // TODO: Insert matrix elements in vector
    } else {
        // TODO: Check if element is in vector, insert if it is not
    }
}
*/
