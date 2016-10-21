#pragma once

enum { ERR_OK, ERR_BAD_HEADER };

#define INPUT_FRAME_SIZE 690 
#define L_FRAME_SIZE 240
#define OUTPUT_FRAME_SIZE (L_FRAME_SIZE*2)

#ifdef LIBTETRADEC_EXPORTS 
#define DLL __declspec(dllexport)
#else
#define DLL
#endif

//DLL int tetra_decode(int first_pass, short input[], short output[]);
//DLL void tetra_decode_init(void);