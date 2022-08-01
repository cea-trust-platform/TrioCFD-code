#include <fftw3.h>
#pragma GCC diagnostic push
#if __GNUC__ > 5
#pragma GCC diagnostic ignored "-Wsuggest-override"
#endif
#include <fftw3-mpi.h>
#pragma GCC diagnostic pop

int main(int argc, char ** argv)
{
void* test = fftw_malloc(sizeof(fftw_complex)*3);
fftw_free(test);

const ptrdiff_t N0 = 7, N1 = 8;
fftw_plan plan;
fftw_complex *data;
ptrdiff_t alloc_local, local_n0, local_0_start;

MPI_Init(&argc, &argv);
fftw_mpi_init();

/* get local data size and allocate */
alloc_local = fftw_mpi_local_size_2d(N0, N1, MPI_COMM_WORLD,
                                     &local_n0, &local_0_start);
data = fftw_alloc_complex(alloc_local);

/* create plan for in-place forward DFT */
plan = fftw_mpi_plan_dft_2d(N0, N1, data, data, MPI_COMM_WORLD,
                            FFTW_FORWARD, FFTW_ESTIMATE);

/* compute transforms, in-place, as many times as desired */
fftw_execute(plan);

fftw_destroy_plan(plan);

MPI_Finalize();


return 0;
}
