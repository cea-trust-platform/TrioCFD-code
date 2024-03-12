#include <fftw3.h>

int main(int argc, char ** argv)
{
void* test = fftw_malloc(sizeof(fftw_complex)*3);
fftw_free(test);
return 0;
}
