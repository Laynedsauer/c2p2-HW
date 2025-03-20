#ifndef FIT_H
#define FIT_H
#include <hls_stream.h>

extern "C" {
	void fit(int* tempxs,int* tempys,int* tempsigmas,int* templasts,int tempN);
}
#endif
