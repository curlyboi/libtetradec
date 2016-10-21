#include <stdio.h>
#include <stdlib.h>
#include <memory.h>
#include "source.h"
#include "channel.h"

#include "libtetradec.h"

#define PRM_SIZE    24
#define SERIAL_SIZE 138
#define SERIAL_SIZE_NOHEADER (138-1)


static void read_tetra_blocks(short block[], short array[]) {
	short          *ptr_block;
	short i;
	short *start_of_array;

	/* Return value: -1  EOF
	*                0  array filled with TETRA frame
	*/
	start_of_array = array;

	/* Copy first valid block */
	ptr_block = block + 1;

	for (i = 0; i < 114; i++)
		*array++ = *ptr_block++;


	/* Copy second valid block */

	ptr_block = block + 161 - 45;

	for (i = 0; i < 114; i++)
		*array++ = *ptr_block++;

	/* Copy third valid block */

	ptr_block = block + 321 - 45 - 45;

	for (i = 0; i < 114; i++)
		*array++ = *ptr_block++;

	/* Copy fourth valid block */

	ptr_block = block + 481 - 45 - 45 - 45;

	for (i = 0; i < 90; i++)
		*array++ = *ptr_block++;

	array = start_of_array;
	for (i = 0; i < 432; i++)
	{

		if ((array[i] & 0x0080) == 0x0080)  array[i]
			= array[i] | 0xFF00;
		if ((array[i] > 127) || (array[i] < -127))
			printf("Input soft bit out of range\n");

	}
}

int frame_stealing = 0;

void sdec(short sinput[], short l_frame[]) {
	short parm[PRM_SIZE];		    /* Synthesis parameters   */

	Bits2prm_Tetra(sinput, parm);	/* serial to parameters */
	Decod_Tetra(parm, l_frame);		/* decoder */
	Post_Process(l_frame, L_FRAME_SIZE);	/* Post processing of synthesis  */
}

DLL void tetra_decode_init(void) {
	Init_Decod_Tetra();
}

DLL int tetra_cdec(int first_pass, short input[], short output[]) {
	short interleaved_coded[432];
	short reordered[286];
	short coded[432];

	if (input[0] != 0x6b21)
		return ERR_BAD_HEADER;
	
	read_tetra_blocks(input, interleaved_coded);

	if (frame_stealing)
	{
		Desinterleaving_Signalling(interleaved_coded + 216, coded + 216);
		/* When Frame Stealing occurs, recopy first half slot : */
		for (int i = 0; i < 216; i++)
			coded[i] = interleaved_coded[i];
	}
	else
	{
		Desinterleaving_Speech(interleaved_coded, coded);
	}

	int bfi1 = frame_stealing;

	/* Channel Decoding */
	int bfi2 = Channel_Decoding(first_pass, frame_stealing, coded, reordered);

	if ((frame_stealing == 0) && (bfi2 == 1))
		bfi1 = 1;

	/* writing  Reordered_array to output */
	output[0] = bfi1;
	memcpy(output + 1, reordered, SERIAL_SIZE_NOHEADER * sizeof(short));
	output[SERIAL_SIZE] = bfi2;
	memcpy(output + SERIAL_SIZE + 1, reordered + SERIAL_SIZE_NOHEADER, SERIAL_SIZE_NOHEADER * sizeof(short));
	
	return (ERR_OK);
}
DLL void tetra_sdec(short input[], short output[]) {
	Word16  serial[SERIAL_SIZE];

	memcpy(serial, input, SERIAL_SIZE * sizeof(short));
	sdec(serial, output);

	memcpy(serial, input + SERIAL_SIZE, SERIAL_SIZE * sizeof(short));
	sdec(serial, output + L_FRAME_SIZE);
}