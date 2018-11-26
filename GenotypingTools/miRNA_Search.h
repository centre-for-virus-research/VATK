/*Functions for unlimited long long array*/
typedef struct { long long *array; size_t used; size_t size; } IntArray;
void initIntArray(IntArray *i, size_t initialSize) { i->array = (long long *)malloc(initialSize * sizeof(long long)); i->used = 0; i->size = initialSize; }
void insertIntArray(IntArray *i, long long element) { if (i->used == i->size) { i->size *= 2; i->array = (long long *)realloc(i->array, i->size * sizeof(long long)); } i->array[i->used++] = element; }
void freeIntArray(IntArray *i) { free(i->array); i->array = NULL; i->used = i->size = 0;}

/*Functions for unlimited char array*/
typedef struct { char *array; size_t used; size_t size;} CharArray;
void initCharArray(CharArray *c, size_t initialSize) { c->array = (char *)malloc(initialSize * sizeof(char)); c->used = 0; c->size = initialSize; }
void insertCharArray(CharArray *c, char element) { if (c->used == c->size) { c->size *= 2; c->array = (char *)realloc(c->array, c->size * sizeof(char)); } c->array[c->used++] = element; }
void freeCharArray(CharArray *c) { free(c->array); c->array = NULL; c->used = c->size = 0; }

/*Karp-Rabin algorithm*/
#define REHASH(a, b, h) ((((h) - (a)*d) << 1) + (b))
long long KR2(char *x, long long m, char *y, long long n) { 
	long long d, hx, hy, i, j, found=0;
	/* Preprocessing */ 
	/* computes d = 2^(m-1) with 
	the left-shift operator */ 
	for (d = i = 1; i < m; ++i) d = (d<<1); 
	for (hy = hx = i = 0; i < m; ++i) { 
		hx = ((hx<<1) + x[i]); 
		hy = ((hy<<1) + y[i]); 
	} 
	/* Searching */ 
	j = 0; 
	while (j <= n-m) { 
		if (hx == hy && memcmp(x, y + j, m) == 0){
			//prlong longf("%d\t",j+1); 
			found++;
		}
		hy = REHASH(y[j], y[j + m], hy); 
		++j; 
	} 
	return found;
}
/*End of Karp-Rabin algorithm*/

// PRINTING TEMPLATE SEQUENCE
// USAGE TempLate (sequence_variable, rev_sequence_variable);
void TempLate(char *seq, char *rev_seq){ 
	long long m=0, n=0; 
	for(m=strlen(seq)-1; m>=0; m--){ 
		if (seq[m]=='A') rev_seq[n++]='T'; 
		else if (seq[m]=='T') rev_seq[n++]='A'; 
		else if (seq[m]=='G') rev_seq[n++]='C'; 
		else if (seq[m]=='C') rev_seq[n++]='G'; 
		else rev_seq[n++]=seq[m]; 
	} 
	rev_seq[n]='\0'; 
}
