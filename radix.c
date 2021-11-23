// from http://codercorner.com/RadixSortRevisited.htm


#include <stdio.h>

int main() {

	int NbItems = 5;
	double InputValues[] = { -0.001, 1000, 0.01, -100, 0};
	int Indexes[NbItems];
	int IndexBuffer[NbItems];
	double OutputValues[NbItems];
	for(int i = 0; i < NbItems ; i++){
		printf("%3.2f, ", InputValues[i]);
		Indexes[i] = i;
		IndexBuffer[i] = i;
	}
	printf("\n");
		
	int Counters[256];
	int Offsets[256];
	unsigned char Radix;
	int Npass = 8;
	
	for(int Pass = 0; Pass < Npass; Pass++){
		
		printf("Pass(%i) \n", Pass);
		for(int i = 0; i < 256; i++ ) {
			Counters[i] = 0;
		}	
		
		for(int i = 0 ; i < NbItems ; i++){       
			printf("binary: %li\n", *((unsigned long*) &InputValues[IndexBuffer[i]]));
			printf("right shift: %i\n", (Pass<<3));
			printf("after shift: %li\n", (*((unsigned long*) &InputValues[IndexBuffer[i]])>>(Pass<<3)));
			Radix = (*((unsigned long*) &InputValues[IndexBuffer[i]])>>(Pass<<3)) & 0xFF;	// Get current byte…
			printf("radix: %i\n", Radix);
			Counters[Radix]++;                    		// …and update counter
		}
		printf("\n");
		if(Pass == (Npass - 1)){
			printf("Last Iteration\n");
			int NbNegativeValues = 0;
			for(int i = 128; i < 256; i++){
				NbNegativeValues += Counters[i];
			}
			printf("NbNegative: %i\n", NbNegativeValues);
			Offsets[0] = NbNegativeValues;
			Offsets[255] = 0;
			for(int i = 1; i < 128; i++){
				Offsets[i] = Offsets[i-1] + Counters[i-1];
			}
			for(int i = 0; i < 127; i++){
				Offsets[254-i] = Offsets[255-i] + Counters[255-i];
			}
			for(int i = 1; i < 256; i++){
				printf("%i, ", Offsets[i]);
			}
			printf("\n");
		} else {
			Offsets[0] = 0;
			for(int i = 1; i < 256; i++){
				Offsets[i] = Offsets[i-1] + Counters[i-1];
			}
		}

		for(int i = 0; i < NbItems ; i++){
			Radix = (*((unsigned long*) &InputValues[IndexBuffer[i]])>>(Pass<<3)) & 0xFF;	// Get current byte…
			Indexes[Offsets[Radix]++] = IndexBuffer[i];
		}
		for(int i = 0; i < NbItems ; i++){
			IndexBuffer[i] = Indexes[i];
		}
	}

	for(int i = 0; i < NbItems ; i++){
		printf("%3.3f,\t ", InputValues[i]);
	}
	printf("\n");
	for(int i = 0; i < NbItems ; i++){
		printf("%3.3f,\t ", InputValues[Indexes[i]]);
	}
	printf("\n");
	for(int i = 0; i < NbItems ; i++){
		printf("%i,\t", Indexes[i]);
	}
	printf("\n");
}
