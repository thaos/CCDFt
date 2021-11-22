// from http://codercorner.com/RadixSortRevisited.htm


#include <stdio.h>

int main() {

	int NbItems = 5;
	float InputValues[] = { 0.001, 1000, 0.01, 100, 0};
	float DestinationBuffer[NbItems];
	float OutputValues[NbItems];
	for(int i = 0; i < NbItems ; i++){
		printf("%3.2f, ", InputValues[i]);
		DestinationBuffer[i] = InputValues[i];
	}
	printf("\n");
		
	int Counters[256];
	int Offsets[256];
	unsigned char Radix;
	
	for(unsigned int Pass = 0; Pass < 4; Pass++){
		
		printf("Pass(%i) \n", Pass);
		for(int i = 0; i < 256; i++ ) {
			Counters[i] = 0;
		}	
		
		for(int i = 0 ; i < NbItems ; i++){       
			printf("binary: %i\n", *((unsigned int*) &DestinationBuffer[i]));
			printf("right shift: %i\n", (Pass<<3));
			printf("after shift: %i\n", (*((unsigned int*) &DestinationBuffer[i])>>(Pass<<3)));
			Radix = (*((unsigned int*) &DestinationBuffer[i])>>(Pass<<3)) & 0xFF;	// Get current byte…
			printf("radix: %i\n", Radix);
			Counters[Radix]++;                    		// …and update counter
		}
		printf("\n");
		
		Offsets[0] = 0;
		for(int i = 1; i < 256; i++){
			Offsets[i] = Offsets[i-1] + Counters[i-1];
		}


		for(int i = 0; i < NbItems ; i++){
			Radix = (*((unsigned int*) &DestinationBuffer[i])>>(Pass<<3)) & 0xFF;	// Get current byte…
			OutputValues[Offsets[Radix]++] = DestinationBuffer[i];
		}
		for(int i = 0; i < NbItems ; i++){
			DestinationBuffer[i] = OutputValues[i];
		}
	}

	for(int i = 0; i < NbItems ; i++){
		printf("%3.3f, ", OutputValues[i]);
	}
	printf("\n");
}
